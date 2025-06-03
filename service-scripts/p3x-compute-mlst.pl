#
# Compute genome MLST using https://github.com/tseemann/mlst
#

use strict;
use Getopt::Long::Descriptive;  
use File::Temp;
use File::Copy;
use GenomeTypeObject;
use Cwd;
use Time::HiRes qw(gettimeofday);
use File::Slurp;
use JSON::XS;
use Data::Dumper;
use Text::CSV_XS qw(csv);
use IPC::Run qw(run);

my ($opt, $usage) = describe_options("%c %o [< in] [> out]",
				     ["in|i=s", "Input GTO"],
				     ["out|o=s", "Output GTO"],
				     ["parallel=i", "Number of threads to use", { default => 1 }],
				     ["help|h", "Print this help message"]);
				     
print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if (@ARGV != 0);

chomp(my $hostname = `hostname -f`);

# mlst /vol/patric3/downloads/genomes/83332.12/83332.12.fna

my $gto_in;

if ($opt->in)
{
    $gto_in = GenomeTypeObject->new({file => $opt->in});
    $gto_in or die "Error reading GTO from " . $opt->in . "\n";
}
else
{
    $gto_in = GenomeTypeObject->new({file => \*STDIN});
    $gto_in or die "Error reading GTO from standard input\n";
}

my $in_file = File::Temp->new();

$gto_in->write_contigs_to_file($in_file);

my @cmd = ("mlst", $in_file);

print STDERR "Invoke mlst: @cmd\n";

my($stdout, $stderr);
my $ok = run(\@cmd,
	     ">", \$stdout,
	     "2>", \$stderr);

$ok or die "Error $? running @cmd\n";

chomp $stdout;

my($fn, $scheme, $st, @alleles) = split("\t", $stdout);

my($sw_version) = $stderr =~ /This is mlst\s+(\S+)/m;
print STDERR $stderr;

my $event = {
    tool_name => "p3x-compute-mlst sw_version=$sw_version",
    parameters => [map { "$_" } @cmd],
    execution_time => scalar gettimeofday,
    hostname => $hostname,
};
my $event_id = $gto_in->add_analysis_event($event);

my $typing = {
    typing_method => "MLST",
    database => $scheme,
    tag => $st,
    event_id => $event_id,
};
push(@{$gto_in->{typing}}, $typing);

if ($opt->out)
{
    $gto_in->destroy_to_file($opt->out);
}
else
{
    $gto_in->destroy_to_file(\*STDOUT);
}
