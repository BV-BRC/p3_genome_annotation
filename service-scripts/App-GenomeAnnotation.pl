#
# The Genome Annotation application.
# temp copy that rips out the checkm/eval code
#

use Bio::KBase::AppService::AppScript;
use Bio::P3::GenomeAnnotationApp::GenomeAnnotationCore;
use Bio::KBase::AppService::AppConfig qw(data_api_url db_host db_user db_pass db_name seedtk application_backend_dir);
use Bio::KBase::AppService::FastaParser 'parse_fasta';
use Bio::KBase::AppService::LongestCommonSubstring qw(BuildString BuildTree LongestCommonSubstring);
use IPC::Run;
use SolrAPI;
use DBI;

use strict;
use Data::Dumper;
use gjoseqlib;
use File::Basename;
use File::Temp;
use LWP::UserAgent;
use JSON::XS;
use IPC::Run;
use IO::File;
use Module::Metadata;
use GenomeTypeObject;

#
# skani organism prediction configuration.
# Database and taxon map are loaded from the application backend directory
# at /vol/bvbrc/production/application-backend/genome_annotation/skani/current/.
#

my $skani_data_dir = application_backend_dir . "/genome_annotation/skani/current";
my $SKANI_DB       = "$skani_data_dir/bvbrc_ref_sketches";
my $SKANI_TAXON_MAP = "$skani_data_dir/genome_taxon_map.tsv";
my $SKANI_MIN_ANI  = 80.0;
my $SKANI_MIN_AF   = 30.0;
my $SKANI_THREADS  = 8;

my $skani_enabled = -d $SKANI_DB;
if ($skani_enabled)
{
    print STDERR "skani organism prediction enabled; database at $SKANI_DB\n";
}
else
{
    print STDERR "skani database not found at $SKANI_DB; organism prediction disabled\n";
}

my $script = Bio::KBase::AppService::AppScript->new(\&process_genome, \&preflight);

my $rc = $script->run(\@ARGV);

exit $rc;

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    #
    # Ensure the contigs are valid, and look up their size.
    #

    my $ctg = $params->{contigs};
    $ctg or die "Contigs must be specified\n";

    my $res = $app->workspace->stat($ctg);
    $res->size > 0 or die "Contigs not found\n";

    #
    # Size estimate based on conservative 500 bytes/second aggregate
    # compute rate for contig size, with a minimum allocated
    # time of 60 minutes (to account for non-annotation portions).
    #
    my $time = $res->size / 500;
    $time = 7200 if $time < 7200;

    my $ram = "16G";
    if ($res->size > 10_000_000)
    {
	$ram = "128G";
    }

    #
    # Request 8 cpus for some of the fatter bits of the compute.
    #
    return {
	cpu => 8,
	memory => $ram,
	runtime => int($time),
	storage => 10 * $res->size,
    };
}

sub process_genome
{
    my($app, $app_def, $raw_params, $params) = @_;

    print "Proc genome ", Dumper($app_def, $raw_params, $params);

    my $json = JSON::XS->new->pretty(1)->canonical(1);

    #
    # Do some sanity checking on params.
    #
    # Both recipe and workflow may not be specified.
    #
    if ($params->{workflow} && $params->{recipe})
    {
	die "Both a workflow document and a recipe may not be supplied to an annotation request";
    }

    my $core = Bio::P3::GenomeAnnotationApp::GenomeAnnotationCore->new(app => $app,
								       app_def => $app_def,
								       params => $params);

    if (exists($raw_params->{tax_id}) && !exists($params->{taxonomy_id}))
    {
	print STDERR "Fixup incorrect taxid in parameters\n";
	$params->{taxonomy_id} = $raw_params->{tax_id};
    }

    if ($params->{taxonomy_id} && $params->{taxonomy_id} !~ /^\s*(\d+)\s*$/)
    {
	die "Invalid taxonomy id (must be an integer)\n";
    }
    $params->{taxonomy_id} = $1;

    #
    # If taxonomy_id was not provided, run skani organism prediction.
    # We need the contigs downloaded to a local temp file first.
    #
    if (!$params->{taxonomy_id} && $skani_enabled)
    {
	print STDERR "No taxonomy_id supplied; running skani organism prediction\n";

	my $skani_temp = File::Temp->new(SUFFIX => '.fasta');
	my $ws = $app->workspace();
	$ws->copy_files_to_handles(1, $core->token,
				   [[$params->{contigs}, $skani_temp]]);
	close($skani_temp);

	#
	# Handle gzipped contigs.
	#
	my $skani_input = "$skani_temp";
	open(my $fh_check, "<", $skani_input) or die "Cannot open $skani_input: $!";
	my $magic;
	read($fh_check, $magic, 2);
	close($fh_check);
	if ($magic eq "\037\213")
	{
	    my $gunzip_temp = File::Temp->new(SUFFIX => '.fasta');
	    IPC::Run::run(["gzip", "-d", "-c", $skani_input],
			  ">", "$gunzip_temp")
		or die "Failed to decompress contigs for skani: $!\n";
	    $skani_input = "$gunzip_temp";
	}

	my $prediction = run_skani_prediction($skani_input);

	print STDERR sprintf(
	    "skani prediction: taxon_id=%s name='%s' (ANI=%.2f%%, AF=%.2f%%)\n",
	    $prediction->{taxonomy_id},
	    $prediction->{scientific_name},
	    $prediction->{ani},
	    $prediction->{af},
	);

	$params->{taxonomy_id}     = $prediction->{taxonomy_id};
	$params->{scientific_name} //= $prediction->{scientific_name};
	$params->{_skani_prediction} = $prediction;

	#
	# Prepend the predicted scientific name to the output_file
	# so it matches the naming convention used when the user
	# selects the organism on the frontend.
	#
	if ($params->{output_file})
	{
	    $params->{output_file} = $params->{scientific_name} . " " . $params->{output_file};
	}
	else
	{
	    $params->{output_file} = $params->{scientific_name};
	}
    }
    elsif (!$params->{taxonomy_id})
    {
	die "No taxonomy_id supplied and skani organism prediction is not configured. " .
	    "Please supply taxonomy_id.\n";
    }

    my $user_id = $core->user_id;

    #
    # If we are missing domain and/or genetic code, look up from taxonomy.
    #
    my $api = P3DataAPI->new();
    my $def = $api->get_taxon_metadata($params->{taxonomy_id});

    print Dumper($params);
    if (!$params->{domain} ||
	$params->{domain} eq 'auto' ||
	!$params->{code})
    {
	print Dumper($def);
	$params->{domain} = $def->{domain} if $def->{domain};
	$params->{code} = $def->{genetic_code} if $def->{genetic_code};
    }
	
    print Dumper($params);
    #
    # Construct genome object metadata and create a new genome object.
    #

    #
    # Ensure we have a scientific_name at this point (either user-supplied
    # or from skani prediction).
    #
    if (!$params->{scientific_name})
    {
	die "No scientific_name provided and organism prediction did not produce one. " .
	    "Please supply scientific_name or omit taxonomy_id to enable prediction.\n";
    }

    my $meta = {
	scientific_name => $params->{scientific_name},
	genetic_code => $params->{code},
	domain => $params->{domain},
	($def->{taxon_lineage} ? (ncbi_lineage => $def->{taxon_lineage}) : ()),
	($params->{taxonomy_id} ? (ncbi_taxonomy_id => $params->{taxonomy_id}) : ()),
	($user_id ? (owner => $user_id) : ()),
	($params->{_skani_prediction} ? (
	    organism_prediction_method => "skani",
	    organism_prediction_ani    => $params->{_skani_prediction}{ani},
	    organism_prediction_af     => $params->{_skani_prediction}{af},
	    organism_prediction_ref    => $params->{_skani_prediction}{genome_id},
	) : ()),
    };
    my $genome = $core->impl->create_genome($meta);

    #
    # Determine workspace paths for our input and output
    #

    my $ws = $app->workspace();

    my($input_path) = $params->{contigs};

    my $output_folder = $app->result_folder();

    my $output_base = $params->{output_file};

    if (!$output_base)
    {
	$output_base = basename($input_path);
    }

    #
    # Read contig data
    #

    my $temp = File::Temp->new();

    $ws->copy_files_to_handles(1, $core->token, [[$input_path, $temp]]);
    
    close($temp);
    my $contig_data_fh = open_contigs("$temp");


    #
    # Use the state-machine based parser ported from RAST.
    #
    # We do an initial read of the first 64 ids of the file. If any of these are too long,
    # we compute the longest common substring to use to replace clip out and replace with
    # "contigs_" to shorten the ids. We keep a mapping from new name to original name.
    #
    # This should hit the common cases where we have IDs like
    # SA-B-8-4-Neg_un-mapped_reads_[SA-B-8-4-Neg_S4_L001_R1_001]_(paired)_contig_22
    # or
    # Ahmed6_GGACTCC_L005_R1_001_(paired)_trimmed_(paired)_merged_contig_1
    #

    my $n = 0;
    my @ids;
    my $too_long = 0;
    my $max_id_len = 60;
    parse_fasta($contig_data_fh, undef, sub {
	my($id, $seq) = @_;
	push(@ids, $id);
	$too_long++ if length($id) > $max_id_len;
	$n++;
	return ($n < 64);
    });
    close($contig_data_fh);

    if ($n == 0)
    {
	die "No contigs loaded from $temp $input_path\n";
    }

    my %name_map;
    my %name_rev_map;
    my $lcs;
    my $qlcs;
    if ($too_long)
    {
	if ($n == 1)
	{
	    $name_map{'contig'} = $ids[0];
	    $name_rev_map{$ids[0]} = 'contig';
	}
	else
	{
	    BuildString(@ids);
	    my $tree = BuildTree();
	    $lcs = LongestCommonSubstring($tree);
	    $qlcs = quotemeta($lcs);
	    print STDERR "Shortening contig names using longest substring '$lcs'\n";
	}
    }

    #
    # Reopen the file and load data, remapping names if necessary.
    #

    $contig_data_fh = open_contigs("$temp");

    my $n = 0;
    my @contigs;
    parse_fasta($contig_data_fh, undef, sub {
	my($id, $seq) = @_;

	my $orig_id = $id;
	my @orig;
	if ($name_rev_map{$id})
	{
	    $id = $name_rev_map{$id};
	    @orig = (original_id => $orig_id);
	}
	elsif ($lcs)
	{
	    $id =~ s/$qlcs/contig_/;
	    @orig = (original_id => $orig_id);
	}
	if (length($id) > $max_id_len)
	{
	    die "Contig id $orig_id too long even after shortening to $id via longest substring $lcs\n";
	}

	push(@contigs, { id => $id, dna => $seq, @orig });
	$n++;
	return 1;
    });
    close($contig_data_fh);
    $core->impl->add_contigs($genome, \@contigs);

    local $Bio::KBase::GenomeAnnotation::Service::CallContext = $core->ctx;

    my $result;

    #
    # Set overrides. This was the method used before we stashed
    # params in the context; unclear if this is better or worse from
    # an engineering perspective.
    #
    my $override;
    if (my $ref = $params->{reference_genome_id})
    {
	$override->{evaluate_genome} =  {
	    evaluate_genome_parameters => { reference_genome_id => $ref },
	};
    }
    if (my $ref = $params->{reference_virus_name})
    {
	$override->{call_features_vigor4} =  {
	    vigor4_parameters => { reference_name => $ref },
	};
    }
    
    $result = $core->run_pipeline($genome, $params->{workflow}, $params->{recipe}, $override);

    my($gto_path, $index_queue_id) = $core->write_output($genome, $result, {}, undef,
							 $params->{public} ? 1 : 0,
							 $params->{queue_nowait} ? 1 : 0,
							 $params->{skip_indexing} ? 1 : 0);

    #
    # Write the genome quality data as a standalone JSON file to aid in downstream
    # quality summarization.
    #

    if (ref($result->{quality}) && !$app->donot_create_result_folder())
    {
	$ws->save_data_to_file($json->encode($result->{quality}),
		           {}, "$output_folder/quality.json", "json", 1, 1, $core->token);
    }

    #
    # Determine if we are one of a peer group of jobs that was started
    # on behalf of a parent job. If we are, and if we are the last job running,
    # invoke post-parent-job processing.
    #

    my $parent_output_folder;
    my $run_last;
    
    if (my $parent = $params->{_parent_job})
    {
	my $dsn = "DBI:mysql:database=" . db_name . ";host=" . db_host;
	my $dbh = DBI->connect($dsn, db_user, db_pass, { RaiseError => 1, AutoCommit => 0 });
	
	#
	# Save information about this genome.
	#
	$dbh->do(qq(INSERT INTO GenomeAnnotation_JobDetails (job_id, parent_job, genome_id, genome_name, gto_path)
		    VALUES (?, ?, ?, ?, ?)), undef,
		 $app->task_id, $parent, $genome->{id}, $genome->{scientific_name}, $gto_path);
	$dbh->commit();

	my $sth = $dbh->prepare(qq(SELECT children_created, children_completed, parent_app, app_spec, app_params
				   FROM JobGroup
				   WHERE parent_job = ?
				   FOR UPDATE));
	my $res = $sth->execute($parent);
	my $last_job = 0;
	my($created, $completed, $app, $spec, $params);
	if ($res != 1)
	{
	    warn "Missing parent job $parent in database\n";
	}
	else
	{
	    ($created, $completed, $app, $spec, $params) = $sth->fetchrow_array();
	    print "Created=$created completed=$completed\n";

	    eval {
		my $params_dat = decode_json($params);
		$parent_output_folder = $params_dat->{output_path} . "/." . $params_dat->{output_file};
		print STDERR "Found parent output folder $parent_output_folder\n";
	    };
	    if ($@)
	    {
		warn "Error parsing parent job params data : $@\n$params\n";
	    }

	    if ($completed == $created - 1)
	    {
		print "We are the last one out!\n";
		$last_job = 1;
	    }
	    elsif ($completed < $created - 1)
	    {
		print "Not so many gone\n";
	    }
	    else
	    {
		warn "completed=$completed created=$created - should not happen here\n";
	    }
	    my $n = $dbh->do(qq(UPDATE JobGroup
				SET children_completed = children_completed + 1
				WHERE parent_job = ? AND children_completed = ?), undef,
			     $parent, $completed);
	    if ($n == 0)
	    {
		print "Failed on update - not us!\n";
	    }
	    else
	    {
		print "Completed with n=$n\n";
	    }
	}

	$dbh->commit();
	if ($last_job)
	{
	    #
	    # We defer this until we are outside this block so we can ensure
	    # the genome report is written first.
	    #
	    $run_last = sub { $core->run_last_job_processing($parent, $app, $spec, $params); };
	}
    }

    
    #
    # Do last-job processing if needed.
    #
    &$run_last if $run_last;

    $core->ctx->stderr(undef);

    return {
	gto_path => $gto_path,
	index_queue_id => $index_queue_id,
	genome_id => $genome->{id},
    };
}

sub open_contigs
{
    my($file) = @_;
    my $contig_data_fh;
    open($contig_data_fh, "<", $file) or die "Cannot open contig file $file: $!";

    #
    # Read first block to see if this is a gzipped file.
    #
    my $block;
    $contig_data_fh->read($block, 256);
    if ($block =~ /^\037\213/)
    {
	#
	# Gzipped. Close and reopen from gunzip.
	#
	
	close($contig_data_fh);
	undef $contig_data_fh;
	open($contig_data_fh, "-|", "gzip", "-d", "-c", $file) or die "Cannot open gzip from $file: $!";
    }
    else
    {
	$contig_data_fh->seek(0, 0);
    }

    return $contig_data_fh;
}

sub run_seedtk_cmd
{
    my(@cmd) = @_;
    local $ENV{PATH} = seedtk . "/bin:$ENV{PATH}";
    my $ok = IPC::Run::run(@cmd);
    $ok or die "Failure $? running seedtk cmd: " . Dumper(\@cmd);
}

#
# skani organism prediction subroutines.
#

sub run_skani_prediction
{
    my($contig_file) = @_;

    print STDERR "Running skani organism prediction on $contig_file\n";
    print STDERR "  Database: $SKANI_DB\n";
    print STDERR "  Thresholds: ANI >= $SKANI_MIN_ANI%, AF >= $SKANI_MIN_AF%\n";

    my @cmd = ("skani", "search",
	       "--qi", $contig_file,
	       "-d", $SKANI_DB,
	       "-n", "1",
	       "-t", $SKANI_THREADS);

    my ($out, $err);
    my $ok = IPC::Run::run(\@cmd, ">", \$out, "2>", \$err);

    if (!$ok)
    {
	die "skani search failed: $err\n";
    }

    print STDERR "skani stderr: $err\n" if $err;

    #
    # Parse skani output. Tab-separated columns:
    # Ref_file  Query_file  ANI  Align_fraction_ref  Align_fraction_query  Ref_name  Query_name
    #
    my @lines = grep { !/^Ref_file/ && /\S/ } split(/\n/, $out);

    if (!@lines)
    {
	die "skani prediction failed: no hits returned. " .
	    "Please supply taxonomy_id manually.\n";
    }

    my $top = $lines[0];
    my @fields = split(/\t/, $top);

    my $ref_file = $fields[0];
    my $ani      = $fields[2];
    my $af_query = $fields[4];
    my $ref_name = $fields[5] // "";

    #
    # Extract genome_id from the reference file path.
    # Expected pattern: .../fasta/<genome_id>.fna
    #
    my $genome_id;
    if ($ref_file =~ m{/([^/]+)\.fn[a-z]*$})
    {
	$genome_id = $1;
    }
    else
    {
	die "Cannot parse genome_id from skani reference path: $ref_file\n";
    }

    print STDERR sprintf("skani top hit: genome=%s ANI=%.2f%% AF=%.2f%% name=%s\n",
			 $genome_id, $ani, $af_query, $ref_name);

    #
    # Apply thresholds.
    #
    if ($ani < $SKANI_MIN_ANI || $af_query < $SKANI_MIN_AF)
    {
	die sprintf(
	    "Cannot determine organism from contigs. " .
	    "Closest reference: %s (%s) with ANI=%.2f%%, AF=%.2f%%. " .
	    "Thresholds: ANI>=%.1f%%, AF>=%.1f%%. " .
	    "Please supply taxonomy_id and scientific_name manually.\n",
	    $genome_id, $ref_name, $ani, $af_query,
	    $SKANI_MIN_ANI, $SKANI_MIN_AF
	);
    }

    #
    # Look up taxonomy from the mapping file.
    #
    my $taxon_info = lookup_genome_taxonomy($genome_id, $SKANI_TAXON_MAP);

    return {
	genome_id       => $genome_id,
	ani             => $ani,
	af              => $af_query,
	ref_name        => $ref_name,
	taxonomy_id     => $taxon_info->{taxon_id},
	scientific_name => $taxon_info->{scientific_name},
    };
}

sub lookup_genome_taxonomy
{
    my($genome_id, $taxon_map_file) = @_;

    open(my $fh, "<", $taxon_map_file)
	or die "Cannot open taxonomy map $taxon_map_file: $!\n";

    my $header = <$fh>;
    while (<$fh>)
    {
	chomp;
	my @f = split(/\t/);
	if ($f[0] eq $genome_id)
	{
	    close($fh);
	    return {
		taxon_id        => $f[1],
		scientific_name => $f[2],
	    };
	}
    }
    close($fh);

    die "Genome $genome_id not found in taxonomy map $taxon_map_file\n";
}
