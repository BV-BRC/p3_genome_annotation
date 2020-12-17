#
# Use vipr_mat_peptide to annotate the given CDS.
#
# We write a genbank file with rast_export_genome, invoke vipr_mat_peptide,
# and add the generated mat_peptide features.
#

use strict;
use GenomeTypeObject;
use Getopt::Long::Descriptive;
use IPC::Run qw(run);
use File::SearchPath qw(searchpath);

my($opt, $usage) = describe_options("%c %o",
				    ["input|i=s" => "Input file"],
				    ["output|o=s" => "Output file"],
				    ["debug|d" => "Enable debugging"],
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 0;

my $tempdir = File::Temp->newdir(CLEANUP => ($opt->debug ? 0 : 1));

print STDERR "Tempdir=$tempdir\n" if $opt->debug;

#
# Ugh. Need to set up a bindir with clustalw poitning at clustalw2.
#
mkdir("$tempdir/bin") or die "mkdir $tempdir/bin failed: $!";
my @p = searchpath("clustalw2");
if (@p)
{
    symlink($p[0], "$tempdir/bin/clustalw");
}

$ENV{PATH} = "$tempdir/bin:$ENV{PATH}";
				    
my $ok = run(["rast_export_genome",
	      "-o", "$tempdir/genome.gb",
	      ($opt->input ? ("-i", $opt->input) : ()),
	      "genbank"]);
$ok or die "Failed to export: $?\n";

$ok = run(["vipr_mat_peptide",
	   "-d", $tempdir,
	   "-i", "genome.gb"],
	   ">", "$tempdir/out.txt",
	  "2>", "$tempdir/err.txt");
$ok or die "vipr_mat_peptide failed: $?\n";
