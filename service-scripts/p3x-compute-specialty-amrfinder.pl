
# Compute specialty gene assignments using AMRFinder
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

#
# offsets into similarity_association tuple
#
use constant {
    a_source => 0,
    a_source_id => 1,
    a_query_coverage => 2,
    a_subject_coverage => 3,
    a_identity => 4,
    a_e_value => 5,
    a_notes => 6,
};

my ($opt, $usage) = describe_options("%c %o [< in] [> out]",
				     ["in|i=s", "Input GTO"],
				     ["out|o=s", "Output GTO"],
				     ["text=s", "Output text file"],
				     ["parallel=i", "Number of threads to use", { default => 1 }],
				     ["help|h", "Print this help message"]);
				     
print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if (@ARGV != 0);

chomp(my $hostname = `hostname -f`);

# amrfinder -p /vol/patric3/downloads/genomes/83332.12/83332.12.PATRIC.faa -O "Mycobacterium_tuberculosis" --threads 24

my $orig_dir = getcwd();

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

#
# AMRFinder leaves crud in /tmp. Set up our own tmpdir and use that.
#

my $work_dir = File::Temp->newdir();

my $in_file = "$work_dir/amr.in";
my $out_file;
if ($opt->text)
{
    $out_file = $opt->text;
}
else
{
    $out_file = "$work_dir/amr.out";
}

$gto_in->write_protein_translations_to_file($in_file);

#
# Determine org list for AMRFinder and check to see if this organism is one of them
# Also pluck the database version.
#

my $db_version = "Undetermined";
my $sw_version = "Undetermined";
my %amr_orgs;
{
    my($stderr, $stdout);
    run(["amrfinder", "--list_organisms"],
	">", \$stdout,
	"2>", \$stderr);

    if ($stdout =~ /^Available --organism options:\s+(.*)$/m)
    {
	%amr_orgs = map { $_ => 1 } split(/,\s+/, $1);
    }

    if ($stderr =~ /Software version:\s+(\S+)/m)
    {
	$sw_version = $1;
    }
    if ($stderr =~ /Database version:\s+(\S+)/m)
    {
	$db_version = $1;
    }
}

my @organism_flag;

my($genus, $species) = $gto_in->{scientific_name} =~ /^(\S+)\s+(\S+)/;
my $org = "${genus}_$species";

if ($amr_orgs{$org})
{
    @organism_flag = ("-O" => $org);
}
elsif ($amr_orgs{$genus})
{
    @organism_flag = ("-O" => $genus);
}

my @cmd = ("amrfinder",
	   "--threads", $opt->parallel,
	   "--protein", $in_file,
	   "--output", $out_file,
	   @organism_flag);

print STDERR "Invoke AMRFinder: @cmd\n";

my $ok = run(\@cmd,
	     init => sub { $ENV{TMPDIR} = "$work_dir"; });

$ok or die "Error $? running @cmd\n";


my $event = {
    tool_name => "p3x-compute-specialty-amrfinder db_version=$db_version sw_version=$sw_version",
    parameters => [map { "$_" } @cmd],
    execution_time => scalar gettimeofday,
    hostname => $hostname,
};

my $event_id = $gto_in->add_analysis_event($event);

my $data = csv(in => $out_file, sep_char => "\t", headers => 'auto');
for my $ent (@$data)
{
    my $fid = $ent->{'Protein id'};
    my $feat = $gto_in->find_feature($fid);
    if (ref($feat))
    {
	my $assoc = [];
	$assoc->[a_source] = "NDARO:($ent->{Method})";
	$assoc->[a_source_id] = $ent->{'Closest reference accession'};
	my $ali_len = $ent->{'Alignment length'};
	my $target_len = $ent->{'Target length'};
	$assoc->[a_query_coverage] = 0 + sprintf("%.2f", 100 * $ali_len / $target_len);
	$assoc->[a_subject_coverage] = 0 + $ent->{'% Coverage of reference'};
	$assoc->[a_identity] = 0 + $ent->{'% Identity to reference'};
	my @note_keys = ('Class', 'Closest reference name', 'Element symbol', 'Subclass', 'Subtype', 'Type');
	$assoc->[a_notes] = {  (map { $_ => $ent->{$_} } @note_keys),
				   software_version => $sw_version,
				   db_version => $db_version,
				   event_id => $event_id,
			       };
	push(@{$feat->{similarity_associations}}, $assoc);
    }
}

if ($opt->out)
{
    $gto_in->destroy_to_file($opt->out);
}
else
{
    $gto_in->destroy_to_file(\*STDOUT);
}
