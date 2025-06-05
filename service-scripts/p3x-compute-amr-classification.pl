#
# Compute AMR classifications using ANL classifiers.
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
use IPC::Run qw(run start);
use File::Basename;
use File::Path qw(make_path);
use POSIX;

my ($opt, $usage) = describe_options("%c %o [< in] [> out]",
				     ["data-dir=s", "Classifier data directory"],
				     ["in|i=s", "Input GTO"],
				     ["out|o=s", "Output GTO"],
				     ["mic-text=s", "Save MIC output text table here"],
				     ["sir-text=s", "Save SIR output text table here"],
				     ["parallel=i", "Number of threads to use", { default => 1 }],
				     ["help|h", "Print this help message"]);
				     
print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if (@ARGV != 0);

chomp(my $hostname = `hostname -f`);

# /opt/patric-common/runtime/bin/python3 predict.py -t /disks/tmp/tmp1 -n 40 -s Staphylococcus_aureus -f ~/rast_comparison_genomes/staph.fa -o /disks/tmp/out1

#
# Determine database version
#
$opt->data_dir or die "--data-dir must be specified\n";

my $data_dir = $opt->data_dir;
-d $data_dir or die "Data directory $data_dir does not exist\n";

my @model_dirs = grep { -d $_ } <$data_dir/*>;
my $model_version = basename($data_dir);
$model_version =~ s/^OLD_//;
$model_version =~ s/_//g;

#
# Map the phenotype from the model to the full name
#
my %phenotype_map = (S => 'Susceptible',
		     R => 'Resistant');
#
# Read antibiotic name lookup table
#
my %antibiotic_name_mapping;
open(NAMES, "<", "$data_dir/antibiotic-names.txt") or die "Cannot open $data_dir/antibiotic-names.txt: $!";
while (<NAMES>)
{
    chomp;
    my($brc_name, $amr_name) = split(/\t/);
    $antibiotic_name_mapping{$amr_name} = $brc_name;
}
close(NAMES);

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

my($genus, $species) = $gto_in->{scientific_name} =~ /^(\S+)\s+(\S+)/;
our $genome_name = "${genus} $species";

my $debug = 0;
my $work_dir;
if ($debug)
{
    $work_dir = "/tmp/3CxubhVI2C";
}
else
{
    $work_dir = File::Temp->newdir(CLEANUP => 1);
}
print "Using $work_dir\n";

my $in_file = "$work_dir/contigs.fa";

$gto_in->write_contigs_to_file($in_file);

my $n_models = @model_dirs;
my $threads_per_model = ceil($opt->parallel / $n_models);

my @models;
for my $model_dir (@model_dirs)
{
    my $model = new Model($model_dir, $threads_per_model);
    push(@models, $model);
    $model->start();
}
for my $model (@models)
{
    $model->wait();

    my $event = {
	tool_name => "p3x-compute-amr-classification model=$model->{model} model_version=$model->{model_version}",
	parameters => [map { "$_" } @{$model->{cmd}}],
	execution_time => scalar gettimeofday,
	hostname => $hostname,
    };
    my $event_id = $gto_in->add_analysis_event($event);

    $model->process_output($gto_in, $event_id, $opt)
}

if ($opt->out)
{
    $gto_in->destroy_to_file($opt->out);
}
else
{
    $gto_in->destroy_to_file(\*STDOUT);
}

package Model;

use Data::Dumper;

use strict;
use File::Basename;
use IPC::Run;
use File::Path qw(make_path);
use File::Copy;
use Text::CSV_XS qw(csv);
#use Carp::Always;

sub new
{
    my($class, $model_dir, $threads) = @_;

    my $model = basename($model_dir);
    my $model_work_dir = "$work_dir/$model";
    my $tmp_dir = "$model_work_dir/tmp";
    make_path($tmp_dir);
    my $out_file = "$model_work_dir/amr.out";
    
    my @cmd = ("python3", "predict.py",
	       "-t", $tmp_dir,
	       "-n", $threads,
	       "-s", $genome_name,
	       "-f", $in_file,
	       "-o", $out_file);

    unshift(@cmd, "echo") if $debug;

    my $self = {
	model_dir => $model_dir,
	model => $model,
	model_work_dir => $model_work_dir,
	tmp_dir => $tmp_dir,
	out_file => $out_file,
	cmd => \@cmd,
    };
    return bless $self, $class;
}

sub start
{
    my($self) = @_;

    print STDERR "Start model run for $self->{model}\n";
    my $handle = IPC::Run::start($self->{cmd},
		       init => sub {
			   chdir $self->{model_dir};
			   $ENV{PATH} = ".:$ENV{PATH}";
		       });
    $handle or die "Cannot start @{$self->{cmd}}: $!";
    $self->{handle} = $handle;
}

sub wait
{
    my($self) = @_;
    print STDERR "Awaiting result from $self->{model}\n";
    my $ok = $self->{handle}->finish();
    if (!$ok)
    {
	die "Model $self->{model} failed with $?\n";
    }
}

sub process_output
{
    my($self, $gto_in, $event_id, $opt) = @_;
    
    my $mic_out = "$self->{out_file}/contigsamr.mic.pred.tab";
    my $sir_out = "$self->{out_file}/contigs.amr.sir.pred.tab";
    if (-s $mic_out)
    {
	if ($opt->mic_text)
	{
	    copy($mic_out, $opt->mic_text);
	}
	$self->process_MIC($mic_out, $gto_in, $event_id);
    }
    elsif (-s $sir_out)
    {
	if ($opt->sir_text)
	{
	    copy($sir_out, $opt->sir_text);
	}
	$self->process_SIR($sir_out, $gto_in, $event_id);
    }
    else
    {
	print STDERR "No model $self->{model} found for $genome_name in $sir_out or $mic_out\n";
    }
}

#
# MIC predictor
#
# Antibiotic      Median MIC      Pred1   Pred2   Pred3   Pred4   Pred5   W1 Avg  W1 CI Low       W1 CI High
# amikacin        0.25    -2.0    -2.0    -2.0    -2.0    -2.0    0.92020202020202        0.9089840601810188      0.9314199802230212
# bedaquiline     0.03125 -5.0    -5.0    -5.0    -5.0    -5.0    0.9131313131313131      0.8886823360892013      0.9375802901734249
    
sub process_MIC
{
    my($self, $data_file, $gto, $event_id) = @_;
    my $data = csv(in => $data_file, sep_char => "\t", headers => "auto");

    #            {
    #              'Pred3' => '7.0',
    #              'Pred2' => '7.0',
    #              'Pred4' => '7.0',
    #              'W1 CI Low' => '1.0',
    #              'Pred5' => '7.0',
    #              'Antibiotic' => 'pyrazinamide',
    #              'Median MIC' => '128.0',
    #              'Pred1' => '7.0',
    #              'W1 Avg' => '1.0',
    #              'W1 CI High' => '1.0'
    #            },

    for my $ent (@$data)
    {
	my $aname = map_antibiotic_name($ent->{Antibiotic});
	my $w1 = chop_float($ent->{'W1 Avg'});
	my $w1_low = chop_float($ent->{'W1 CI Low'});
	my $w1_high = chop_float($ent->{'W1 CI High'});
	my $assertion = {
	    antibiotic_name => $aname,
	    model_antibiotic_name => $ent->{Antibiotic},
	    evidence => "Computational Method",
	    computational_method => "MIC XGBoost Model",
	    computational_method_performance => "W1 score: $w1, CI[$w1_low, $w1_high]",
	    computational_method_version => $model_version,
	    measurement_unit => 'mg/L',
	    measurement_value => $ent->{'Median MIC'},
	    event_id => $event_id,
	};
	push(@{$gto->{amr_assertions}}, $assertion);
    }
}


# Antibiotic      Prediction      SR      F1 Avg  F1 CI Low       F1 CI High
# amikacin        0.0	     S       0.8775605494467064      0.8401388917347086      0.9149822071587042
# capreomycin     0.0     S       0.8882634433906793      0.8588606088840278      0.9176662778973308
# ciprofloxacin   0.0     S       0.9224356279918939      0.8633464678791597      0.9815247881046282

    
sub process_SIR
{
    my($self, $data_file, $gto, $event_id) = @_;
    my $data = csv(in => $data_file, sep_char => "\t", headers => "auto");

    for my $ent (@$data)
    {
	my $aname = map_antibiotic_name($ent->{Antibiotic});
	my $f1 = chop_float($ent->{'F1 Avg'});
	my $f1_low = chop_float($ent->{'F1 CI Low'});
	my $f1_high = chop_float($ent->{'F1 CI High'});
	my $assertion = {
	    antibiotic_name => $aname,
	    model_antibiotic_name => $ent->{Antibiotic},
	    evidence => "Computational Method",
	    computational_method => "SIR XGBoost Model",
	    computational_method_performance => "F1 score: $f1, CI[$f1_low, $f1_high]",
	    computational_method_version => $model_version,
	    resistant_phenotype => map_phenotype($ent->{SR}),
	    event_id => $event_id,
	};
	push(@{$gto->{amr_assertions}}, $assertion);
    }
    
    
# {
#     'F1 CI Low' => '0.8401388917347086',
#     'SR' => 'S',
#     'Prediction' => '0.0',
#     'F1 Avg' => '0.8775605494467064',
#     'F1 CI High' => '0.9149822071587042',
#     'Antibiotic' => 'amikacin'
#     },

}

sub chop_float
{
    my($val) = @_;
    return 0 + sprintf("%.2f", $val);
}

sub map_antibiotic_name
{
    my($model_name) = @_;
    my $aname = $antibiotic_name_mapping{$model_name};
    if (!$aname)
    {
	warn "Cannot map name '$model_name'\n";
	$aname = $model_name;
    }
    return $aname;
}

sub map_phenotype
{
    my($val) = @_;
    my $mapped = $phenotype_map{$val};
    if ($mapped)
    {
	return $mapped;
    }
    else
    {
	warn "Unknown phenotype value '$val'\n";
	return $val;
    }
}

