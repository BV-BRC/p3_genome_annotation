#
# The Genome Annotation application. Genbank input variation.
#

use Bio::KBase::AppService::AppScript;
use Bio::P3::GenomeAnnotationApp::GenomeAnnotationCore;
use Bio::KBase::AppService::AppConfig 'data_api_url';
use SolrAPI;

use strict;
use Data::Dumper;
use gjoseqlib;
use File::Basename;
use File::Slurp;
use File::Temp;
use LWP::UserAgent;
use JSON::XS;
use IPC::Run 'run';
use IO::File;

my $script = Bio::KBase::AppService::AppScript->new(\&process_genome, \&preflight_cb);

my $rc = $script->run(\@ARGV);

exit $rc;

sub preflight_cb
{
    my($app, $app_def, $raw_params, $params) = @_;

    #
    # Ensure the contigs are valid, and look up their size.
    #

    my $gb = $params->{genbank_file};
    $gb or die "Genbank file must be specified\n";

    my $res = $app->workspace->stat($gb);
    $res->size > 0 or die "Genbank file $gb not found\n";

    #
    # Size estimate based on conservative 500 bytes/second aggregate
    # compute rate for contig size, with a minimum allocated
    # time of 5 minutes.
    #
    my $time = $res->size / 500;
    $time = 3600 if $time < 3600;

    my $ram = "16G";
    if ($res->size > 10_000_000)
    {
	$ram = "128G";
    }

    #
    # Request 8 cpus for some of the fatter bits of the compute.
    #
    return {
	cpu => 2,
	memory => $ram,
	runtime => int($time),
	storage => 10 * $res->size,
    };
}


sub process_genome
{
    my($app, $app_def, $raw_params, $params) = @_;

    print "Proc genome ", Dumper($app_def, $raw_params, $params);

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

    my $user_id = $core->user_id;

    #
    # Determine workspace paths for our input and output
    #

    my $ws = $app->workspace();

    my($input_path) = $params->{genbank_file};

    my $output_folder = $app->result_folder();

    my $output_base = $params->{output_file};

    if (!$output_base)
    {
	$output_base = basename($input_path);
    }

    #
    # Read genbank file data
    #
    # If the genbank file is compressed, uncompress and use that. Downstream
    # code in rast2solr needs the uncompressed data.
    #

    my $gb_temp = File::Temp->new();

    $ws->copy_files_to_handles(1, $core->token, [[$input_path, $gb_temp]]);
    
    my $genbank_data_fh;
    close($gb_temp);
    open($genbank_data_fh, "<", $gb_temp) or die "Cannot open contig temp $gb_temp: $!";

    #
    # Read first block to see if this is a gzipped file.
    #
    my $block;
    $genbank_data_fh->read($block, 256);

    my $gb_file;
    if ($block =~ /^\037\213/)
    {
	#
	# Gzipped. Uncompress into temp.
	#
	$gb_file = File::Temp->new();
	my $ok = run(["gunzip", "-d", "-c",  $gb_temp],
	    ">", $gb_file);
	$ok or die "Could not gunzip $gb_temp: $!";
	close($gb_file);
		
	close($genbank_data_fh);
	undef $genbank_data_fh;
	open($genbank_data_fh, "<", $gb_file) or die "Cannot open $gb_file: $!";
    }
    else
    {
	$genbank_data_fh->seek(0, 0);
	$gb_file = $gb_temp;
    }
    
#     my $gb_data = read_file($genbank_data_fh);
#     close($genbank_data_fh);

    #
    # Use the hidden impl method so that we can pass the filehandle.
    # Published interface expects a string.
    #
    my $genome = $core->impl->_create_genome_from_genbank_impl($genbank_data_fh, $params->{skip_contigs} ? 1 : 0);
    $genome = $genome->prepare_for_return();
    close($genbank_data_fh);

    #
    # Overrides from optional parameters.
    #
    for my $override (['code', 'genetic_code'],
		      ['scientific_name', 'scientific_name'],
		      ['taxonomy_id', 'ncbi_taxonomy_id'],
		      ['domain', 'domain'])
    {
	my($param_name, $gto_name) = @$override;
	
	if (exists $params->{$param_name})
	{
	    $genome->{$gto_name} = $params->{$param_name};
	}
    }

    #
    # Add owner field from token
    #
    if ($core->user_id)
    {
	$genome->{owner} = $core->user_id;
    }

    #
    # See if we are missing a genetic code, and if so, look it up in the
    # data api.
    #
    if (!$genome->{genetic_code})
    {
	my $gc = '';
	my $tax = $genome->{ncbi_taxonomy_id};
	if ($tax =~ /^\d+$/)
	{
	    my $solr = SolrAPI->new(data_api_url);
	    my $res = $solr->query_solr('taxonomy', '/?q=taxon_id:' . $tax, '&fl=genetic_code');
	    if (ref($res) eq 'ARRAY' && @$res)
	    {
		$gc = $res->[0]->{genetic_code};
		warn "Setting GC=$gc based on tax id $tax\n";
	    }
	}
	if ($gc !~ /^\d+$/)
	{
	    my @list = ('Acholeplasma',
			'Candidatus Hepatoplasma',
			'Candidatus Hodgkinia',
			'Entomoplasma',
			'Mesoplasma',
			'Mycoplasma',
			'Spiroplasma',
			'Ureaplasma');
	    if (grep { $genome->{scientific_name} =~ /^$_\s/ } @list)
	    {
		warn "Falling back to GC=4 based on name $genome->{scientific_name}\n";
		$gc = 4;
	    }
	    else
	    {
		warn "Setting default GC=11 based on failed lookup for '$tax' and no match for $genome->{scientific_name}\n";
		$gc = 11;
	    }
	}
	$genome->{genetic_code} = $gc;
    }

    #
    # If our domain is unknown but we have a valid bacterial genetic code, set that.
    #
    if (!$genome->{domain} || $genome->{domain} eq 'Unknown')
    {
	if ($genome->{genetic_code} == 4 || $genome->{genetic_code} == 11)
	{
	    $genome->{domain} = 'Bacteria';
	}
    }

    my $workflow = $params->{workflow};
    if ($params->{import_only})
    {
	if ($workflow)
	{
	    die "Invalid input: workflow and import_only may not both be specified";
	}
	$workflow = JSON::XS->new->pretty(1)->encode($core->import_workflow());
    }
    elsif ($params->{raw_import_only})
    {
	if ($workflow)
	{
	    die "Invalid input: workflow and import_only may not both be specified";
	}
	$workflow = JSON::XS->new->pretty(1)->encode($core->raw_import_workflow());
    }

    my $pipeline_override;
    if (my $ref = $params->{reference_virus_name})
    {
	$pipeline_override->{call_features_vigor4} =  {
	    vigor4_parameters => { reference_name => $ref },
	};
    }

    my $result = $core->run_pipeline($genome, $workflow, $params->{recipe}, $pipeline_override);
    
    #
    # TODO fill in metadata?
    my($gto_path, $index_queue_id) = $core->write_output($genome,
							 $result, {},
							 "$gb_file",
							 $params->{public} ? 1 : 0,
							 $params->{queue_nowait} ? 1 : 0,
							 $params->{skip_indexing} ? 1 : 0);
							 

    $core->ctx->stderr(undef);
    return {
	gto_path => $gto_path,
	index_queue_id => $index_queue_id,
    };

}
