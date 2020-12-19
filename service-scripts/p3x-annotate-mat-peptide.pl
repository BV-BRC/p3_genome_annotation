#
# Use vipr_mat_peptide to annotate the given CDS.
#
# We write a genbank file with rast_export_genome, invoke vipr_mat_peptide,
# and add the generated mat_peptide features.
#

use strict;
use Data::Dumper;
use Time::HiRes 'gettimeofday';
use gjogenbank;
use gjoseqlib;
use GenomeTypeObject;
use IDclient;
use Getopt::Long::Descriptive;
use File::Copy;
use IPC::Run qw(run);
use File::SearchPath qw(searchpath);

my($opt, $usage) = describe_options("%c %o",
				    ["remove-existing" => "Remove existing CDS and mat_peptide features if vipr_mat_peptide run is successful"],
				    ["input|i=s" => "Input file"],
				    ["output|o=s" => "Output file"],
				    ["debug|d" => "Enable debugging"],
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 0;

#my $tempdir = "/tmp/rbtZ1smxfr";
#goto x;

my $tempdir = File::Temp->newdir(CLEANUP => ($opt->debug ? 0 : 1));

print STDERR "Tempdir=$tempdir\n" if $opt->debug;

#
# Ugh. Need to set up a bindir with clustalw pointing at clustalw2.
#
mkdir("$tempdir/bin") or die "mkdir $tempdir/bin failed: $!";
my @p = searchpath("clustalw2");
if (@p)
{
    symlink($p[0], "$tempdir/bin/clustalw");
}

$ENV{PATH} = "$tempdir/bin:$ENV{PATH}";

#
# We copy our input gto to disk if it arrived from
# stdin since we need to both write a genbank from it
# as well as load it to add features later.
#

my $input_gto;
if ($opt->input)
{
    $input_gto = $opt->input;
}
else
{
    $input_gto = "$tempdir/genome.gto";
    copy(\*STDIN, $input_gto);
}
				    
my $gto = GenomeTypeObject->create_from_file($input_gto);
$gto or die "Could not parse input gto\n";

#
# Read again, delete mat_peptides, and write out in order to create genbank
# without these features to force the caller to find them afresh.
#

{
    my $gto_for_gb = GenomeTypeObject->create_from_file($input_gto);
    $gto_for_gb or die "Could not parse input gto\n";

    #
    # Check for vigor4 features.
    #
    my $vigor_ae;
    for my $ae ($gto_for_gb->analysis_events)
    {
	if ($ae->{tool_name} eq 'vigor4')
	{
	    $vigor_ae = $ae;
	    last;
	}
    }

    my @to_del;

    for my $f ($gto_for_gb->features)
    {
	if ($f->{type} eq 'mat_peptide')
	{
	    if ($f->{feature_creation_event} eq $vigor_ae->{id})
	    {
		warn "Already annotated by vigor4. Skipping annotation\n";
		$gto->destroy_to_file($opt->output);
		exit 0;
	    }
	    push(@to_del, $f->{id});
	}
    }
    print "Delete @to_del from original file\n";
    
    $gto_for_gb->delete_feature($_) foreach @to_del;

    $gto_for_gb->destroy_to_file("$tempdir/genome_no_mat_peptide.gto");
}

my $ok = run(["rast_export_genome",
	      "--genbank-roundtrip",
	      "-o", "$tempdir/genome.gb",
	      "-i", "$tempdir/genome_no_mat_peptide.gto",
	      "genbank"]);
$ok or die "Failed to export: $?\n";

my $cmd = "vipr_mat_peptide";
my @params = ("-d", "$tempdir",
	      "-i", "genome.gb");

$ok = run([$cmd, @params],
	   ">", "$tempdir/out.txt",
	  "2>", "$tempdir/err.txt");
$ok or die "vipr_mat_peptide command '$cmd @params' failed: $?\n";

#
# Determine output file; named based on accession but always
# ending with _matpept_msagbk.faa
#

x:

my(@out) = glob("$tempdir/*_matpept_*.faa");

if (@out == 0)
{
    warn "Failed to write an output file\n";
    $gto->destroy_to_file($opt->output);
    exit 0;
}
elsif (@out > 1)
{
    warn "Multiple output files written; this should not happen. Tmp dir reused?\n";
    $gto->destroy_to_file($opt->output);
    exit 0;
}
my $out = $out[0];

if (!open(O, "<", $out))
{
    warn "Cannot open output file $out: $!";
    $gto->destroy_to_file($opt->output);
    exit 0;
}

my $hostname = `hostname`;
chomp $hostname;
my $event = {
    tool_name => $cmd,
    execute_time => scalar gettimeofday,
    parameters => \@params,
    hostname => $hostname,
};
my $event_id = $gto->add_analysis_event($event);
my $id_client = IDclient->new($gto);
my $id_prefix = "fig|$gto->{id}";

if (-s O > 0 && $opt->remove_existing)
{
    my @to_del = $gto->fids_of_type('mat_peptide');
    print "Delete @to_del\n";
    $gto->delete_feature($_) foreach @to_del;
}

while (my($id, $def, $seq) = read_next_fasta_seq(\*O))
{
    # this isn't a strict fasta header; product names may have spaces in which case
    # the data is split over id and def, so just join them with a space before parsing.
    my $dat = join(" " , $id, $def);
    #
    # CDS ids might have vertical bars which munge parsing. Split them out.
    #
    $dat =~ s/\|CDS=(.*?)\|ref/|ref/;
    my $cds_id = $1;
    my %vals = map { /^([A-za-z]+)=(.*)/ } split(/\|/, $dat);
    print Dumper(\%vals);

    my $loc = genbank_loc_2_cbdl($vals{Loc}, $vals{ACC});

    my $f = $gto->add_feature({ -id_client => $id_client,
				    -id_prefix => $id_prefix,
				    -type => "mat_peptide",
				    -location => $loc,
				    -function => $vals{product},
				    -annotator => $cmd,
				    -annotation => "Add mature peptide from reference $vals{ref}",
				    -analyis_event_id => $event_id,
				    -protein_translation => $seq,
				    -genbank_type => 'mat_peptide',
				    -alias_pairs => [[gene => $vals{symbol}]],
				    -analysis_event_id => $event_id,
				    
				});
    print Dumper($f);
}

$gto->destroy_to_file($opt->output);
