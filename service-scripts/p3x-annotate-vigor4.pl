#
# Use vigor4 to annotate the given genome.
#


use strict;
use Data::Dumper;
use Time::HiRes 'gettimeofday';
use gjogenbank;
use gjoseqlib;
use GenomeTypeObject;
use Getopt::Long::Descriptive;
use File::Copy;
use IPC::Run qw(run);
use File::SearchPath qw(searchpath);
use Bio::KBase::GenomeAnnotation::Config qw(vigor_reference_db_directory);
use Bio::P3::GenomeAnnotationApp::VigorTaxonMap;
use P3DataAPI;

use Cwd;

my($opt, $usage) = describe_options("%c %o",
				    ["reference=s" => "Vigor4 reference name"],
				    ["taxon=i" => "Taxon identifier"],
				    ["remove-existing" => "Remove existing CDS and mat_peptide features if vigor4 run is successful"],
				    ["input|i=s" => "Input file"],
				    ["output|o=s" => "Output file"],
				    ["debug|d" => "Enable debugging"],
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 0;

chomp(my $hostname = `hostname`);

my $tempdir = File::Temp->newdir(CLEANUP => ($opt->debug ? 0 : 1));

print STDERR "Tempdir=$tempdir\n" if $opt->debug;

my $here = getcwd;

my $genome_in = GenomeTypeObject->create_from_file($opt->input);
$genome_in or die "Error reading and parsing input";

#
# Determine our reference database.
#
# If --reference passed, use that.
#
# Otherwise get the taxon id from --taxon parameter or from the GTO.
#

my $reference_name = $opt->reference;

if (!$reference_name)
{
    my $taxon = $opt->taxon // $genome_in->taxonomy_id;
    if ($taxon)
    {
	my $api = P3DataAPI->new;
	my @res = $api->query("taxonomy", ['eq', 'taxon_id', $taxon], ['select', 'taxon_name', 'lineage_ids', 'lineage_names']);
	print Dumper(\@res);
	my $res = $res[0];
	my $ids = $res->{lineage_ids};
	my $names = $res->{lineage_names};
	
	for (my $i = $#$names; $i >= 0; $i--)
	{
	    my $taxon = $ids->[$i];
	    my $name = $names->[$i];
	    my $db = $Bio::P3::GenomeAnnotationApp::VigorTaxonMap::map->{$taxon};
	    if ($db)
	    {
		$reference_name = $db->{db};
		last;
	    }
	}
	if (!$reference_name)
	{
	    warn "No reference found for taxon $taxon\n";
	}
    }
}

if (!$reference_name)
{
    warn "No reference found\n";
    $genome_in->destroy_to_file($opt->output);
    exit 0;
}

#
# Invoke vigor4 to annotate viral genome.
#
# Write contigs as fasta file.
# Invoke vigor4, using the reference that was passed in as the reference_name parameter;
# We parse the .pep file that is generated. It contains two feature types.
# This is a CDS:
# >NC_045512.1 location=266..13468,13471..21555 codon_start=1 gene="orf1ab" ref_db="covid19" ref_id="YP_009724389.1"
# This is a mature peptide:
# >NC_045512.1.1 mat_peptide location=266..805 gene="orf1ab" product="leader protein" ref_db="covid19_orf1ab_mp" ref_id="YP_009725297.1"
#
# We create features of type CDS and mat_peptide.
#

my $sequences_file = $genome_in->extract_contig_sequences_to_temp_file();

my $ref_dir = vigor_reference_db_directory;
$ref_dir ne '' or die "Vigor reference directory not configured";
-d $ref_dir or die "Vigor reference directory '$ref_dir' not found";

my @vigor_params = ("-i", $sequences_file,
		    "--reference-database-path", $ref_dir,
		    "-d", $reference_name,
		    "-o", "$here/vigor_out");

print STDERR Dumper(\@vigor_params);
my $ok = run(["echo", "vigor4", @vigor_params],
	     init => sub { chdir $tempdir; system("pwd"); },
	     ">", "$here/vigor4.stdout.txt",
	     "2>", "$here/vigor4.stderr.txt");
if (!$ok)
{
    print STDERR "Vigor run failed with rc=$?. Stdout:\n";
    copy("$here/vigor4.stdout.txt", \*STDERR);
    print STDERR "Stderr:\n";
    copy("$here/vigor4.stderr.txt", \*STDERR);
}
    
my $event = {
    tool_name => "vigor4",
    execution_time => scalar gettimeofday,
    parameters => \@vigor_params,
    hostname => $hostname,
};

my $event_id = $genome_in->add_analysis_event($event);

#
# Parse the generated peptide file. We collect the CDS and mature_peptides, then
# add features so that we can register the counts.
#
if (open(my $pep_fh, "<", "$here/vigor_out.pep"))
{
    my %features;
    while (my($id, $def, $seq) = read_next_fasta_seq($pep_fh))
    {
	my $fq = { truncated_begin => 0, truncated_end => 0 };
	
	my $type;
	my $ctg;
	if ($def =~ s/^mat_peptide\s+//)
	{
	    ($ctg) = $id =~ /^(.*)\.[^.]+\.[^.]+$/;
	    $type = 'mat_peptide';
	}
	elsif ($def =~ s/^pseudogene\s+//)
	{
	    ($ctg) = $id =~ /^(.*)\.[^.]+\.[^.]+$/;
	    $type = 'pseudogene';
	}
	else
	{
	    ($ctg) = $id =~ /^(.*)\.[^.]+$/;
	    $type = 'CDS';
	}
	
	if (!$ctg)
	{
	    print STDERR "Falling back to prefix of id for contig name from $id\n";
	    ($ctg) = $id =~ /^(.*?)\./;
	}
	
	my $feature = {
	    quality => $fq,
	    type => $type,
	    contig => $ctg,
	    aa_sequence => $seq,
	};
	push(@{$features{$type}}, $feature);
	
	while ($def =~ /([^=]+)=((\"([^\"]+)\")|([^\"\s]+))\s*/mg) 
	{
	    my $key = $1;
	    my $val = $4 ? $4 : $5;
	    
	    my @loc;
	    # print "key=$key val=$val\n";

	    if ($key eq 'location')
	    {
		$feature->{genbank_feature} = { genbank_type => $type, genbank_location  => $val, values => {}};
		
		# location=266..13468,13471..21555
		for my $ent (split(/,/, $val))
		{
		    if (my($s_frag, $s, $e_frag, $e) = $ent =~ /^(<?)(\d+)\.\.(>?)(\d+)$/)
		    {
			$fq->{truncated_begin} = 1 if $s_frag;
			$fq->{truncated_end} = 1 if $e_frag;
			
			my $len = abs($s - $e) + 1;
			my $strand = $s < $e ? '+' : '-';
			push(@loc, [$ctg, $s, $strand, $len]);
		    }
		    else
		    {
			die "error parsing location '$ent'\n";
		    }
		}
		$feature->{location} = \@loc;
	    }
	    else
	    {
		$feature->{$key} = $val;
	    }
	}
	$feature->{product} //= $feature->{gene};
	if (!$feature->{location})
	{
	    warn "No location for feature $def " . Dumper($feature);
	}
    }
    #print Dumper(\%features);

    #
    # If we have features, remove any existing CDS or mat_peptide features.
    #

    if (%features && $opt->remove_existing)
    {
	my @to_del = $genome_in->fids_of_type('CDS', 'mat_peptide', 'pseudogene');
	print "Delete @to_del\n";
	$genome_in->delete_feature($_) foreach @to_del;
    }
    
    for my $type (keys %features)
    {
	my $feats = $features{$type};
	my $n = @$feats;
	my $id_type = $type;
	
	for my $feature (@$feats)
	{
	    my $p = {
		-id	     => $genome_in->new_feature_id($id_type),
		-type 	     => $type,
		-location 	     => $feature->{location},
		-analysis_event_id 	     => $event_id,
		-annotator => 'vigor4',
		-protein_translation => $feature->{aa_sequence},
		-alias_pairs => [[gene => $feature->{gene}]],
		-function => $feature->{product},
		-quality_measure => $feature->{quality},
		-genbank_feature => $feature->{genbank_feature},
	    };
	    #die Dumper($p);
	    
	    $genome_in->add_feature($p);
	}
    }
    
}
else
{
    warn "Could not read $here/vigor_out.pep\n";
}

$genome_in->destroy_to_file($opt->output);
