
# Compute CARD specialty gene assignments using RGI
#
# RGI wants the database directory to be called localDB in the working directory
# when the program is invoked. We will create a new tempdir to work in
# and symlink the database directory's localdB into the temp space.
#
# We always write output to our temp dir. Default output is the txt file,
# we also allow saving the JSON output file using --json.
#

use strict;
use Getopt::Long::Descriptive;  
use File::Temp;
use File::Copy;
use GenomeTypeObject;
use Cwd;
use IPC::Run qw(run);
use Time::HiRes qw(gettimeofday);
use File::Slurp;
use JSON::XS;

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
				     ["json=s", "Output JSON file"],
				     ["text=s", "Output text file"],
				     ["db-dir=s", "Database directory"],
				     ["parallel=i", "Number of threads to use", { default => 1 }],
				     ["help|h", "Print this help message"]);
				     
print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if (@ARGV != 0);

chomp(my $hostname = `hostname`);

# rgi main -n 24 --local --include_nudge -t protein --clean -a BLAST -i /vol/patric3/downloads/genomes/83332.12/83332.12.PATRIC.faa  -o ~/BV-BRC/card-test/h37rv-appbackend.rgi

my $orig_dir = getcwd();
my $work_dir = File::Temp->newdir();

my $gto_in;

if (!$opt->db_dir)
{
    die "--db-dir must be specified\n";
}

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

my $in_file = "$work_dir/rgi_input";
my $out_file = "$work_dir/rgi_output";

symlink($opt->db_dir . "/localDB", "$work_dir/localDB") or
    die "Cannot symlink " . $opt->db_dir . "/localDB to $work_dir/localDB: $!";

$gto_in->write_protein_translations_to_file($in_file);

my @cmd = ("rgi", "main",
	   "-n", $opt->parallel,
	   "--local",
	   "--include_nudge",
	   "-t", "protein",
	   "--clean",
	   "-a", "BLAST",
	   "-i", $in_file,
	   "-o", $out_file);
print STDERR "Invoke RGI: @cmd\n";
my $ok = run(\@cmd,
	     init => sub { chdir("$work_dir"); });
$ok or die "Error $? running @cmd\n";

my $out_text = "$out_file.txt";
my $out_json = "$out_file.json";

-f $out_text or die "Output text file $out_text missing\n";
-f $out_json or die "Output json file $out_json missing\n";

my $rgi_data = decode_json(scalar read_file($out_json));

while (my($fid, $res) = each %$rgi_data)
{
    #
    # Each fid may have multiple hits. We sort based on the bitscore and evalue. 
    #
    my @hits = sort { $b->{bitscore} <=> $a->{bitscore} or $a->{evalue} <=> $b->{evalue} } values %$res;
    next unless @hits;
    my $best_hit = $hits[0];
    my $assoc = compute_assoc($best_hit);

    my $feat = $gto_in->find_feature($fid);
    if (ref($feat))
    {
	push(@{$feat->{similarity_associations}}, $assoc);
    }
}

my $event = {
    tool_name => "p3x-compute-specialty-rgi",
    parameters => \@cmd,
    execution_time => scalar gettimeofday,
    hostname => $hostname,
};
$gto_in->add_analysis_event($event);

if ($opt->out)
{
    $gto_in->destroy_to_file($opt->out);
}
else
{
    $gto_in->destroy_to_file(\*STDOUT);
}

#
# Copy RGI outputs if desired
#
if ($opt->json)
{
    # copy($out_json, $opt->json) or die "Error copying $out_json to " . $opt->json . ": $!";
    write_file($opt->json, JSON::XS->new->pretty->canonical->encode($rgi_data));
}

if ($opt->text)
{
    copy($out_text, $opt->text) or die "Error copying $out_text to " . $opt->text . ": $!";
}

sub compute_assoc
{
    my($hit_data) = @_;
    
    my $assoc = [];
    $assoc->[a_source] = 'CARD';
    $assoc->[a_source_id] = $hit_data->{model_id};
    $assoc->[a_e_value] = $hit_data->{evalue};
    $assoc->[a_identity] = $hit_data->{perc_identity};
    
    my $match = $hit_data->{match};
    my $match_len = length($match);
    my $query_len = length($hit_data->{orf_prot_sequence});
    my $subject_len = length($hit_data->{sequence_from_broadstreet});
    # print "Len $match_len $query_len $subject_len\n";
    $assoc->[a_query_coverage]  = 0 + sprintf("%.2f", 100 * $match_len / $query_len);
    $assoc->[a_subject_coverage] = 0 + sprintf("%.2f", 100 * $match_len / $subject_len);
    my @note_keys = qw(model_name model_type type_match partial);
    $assoc->[a_notes] = {  map { $_ => $hit_data->{$_} } @note_keys };

    return $assoc;
}

