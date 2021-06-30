
use strict;
use Getopt::Long::Descriptive;
use File::Basename;
use Bio::P3::Workspace::WorkspaceClientExt;
use JSON::XS;
use Fcntl ':mode';
use File::Temp;
use File::Slurp;
use Cwd qw(abs_path getcwd);
use Data::Dumper;
use Proc::ParallelLoop;
use IPC::Run qw(run);
use Time::HiRes 'gettimeofday';

my($opt, $usage) = describe_options("%c %o output-path genbank-file [genbank-file ...]",
				    ["parallel|j=i", "Run this many parallel threads", { default => 1 }],
				    ["workflow-file=s", "Workflow definition"],
				    ["log-dir=s", "Log directory", { default => "." }],
				    ["public", "Mark genomes public"],
				    ["overwrite", "Overwrite existing workspace files"],
				    ["indexing-url=s", "Indexing url"],
				    ["help|h", "Show this help message"]);
print($usage->text) if $opt->help;
die($usage->text) if @ARGV < 2;

my $output_path = shift;

$ENV{P3_ALLOCATED_CPU} = 2;

my $ws = Bio::P3::Workspace::WorkspaceClientExt->new;
my $json = JSON::XS->new->pretty->canonical;
my $app = "GenomeAnnotationGenbank";
my $app_spec = find_app_spec($app);

my $stat = $ws->stat($output_path);

if (!$stat || !S_ISDIR($stat->mode))
{
    die "Output path " . $output_path . " does not exist\n";
}

#
# Handle custom workflow. This may be either local or in
# the workspace.
#
my $workflow;
my $workflow_txt;
if ($opt->workflow_file)
{
    my $wf = $opt->workflow_file;
    my($wf_ws) = $wf =~ /^ws:(.*)/;

    if ($wf_ws)
    {
	$wf_ws = expand_workspace_path($wf_ws);
	$workflow_txt = $ws->download_file_to_string($wf_ws);
    }
    else
    {
	open(F, "<", $opt->workflow_file) or die "Cannot open workflow file " . $opt->workflow_file . ": $!\n";
	local $/;
	undef $/;
	$workflow_txt = <F>;
	close(F);
    }
    
    eval {
	$workflow = decode_json($workflow_txt);
    };
    if (!$workflow)
    {
	die "Error parsing workflow file " . $opt->workflow_file . ": $@\n";
    }

    if (ref($workflow) ne 'HASH' ||
	!exists($workflow->{stages}) ||
	ref($workflow->{stages}) ne 'ARRAY')
    {
	die "Invalid workflow document (must be a object containing a list of stage definitions)\n";
    }
}

my $base_params = {
    output_path => $output_path,
    queue_nowait => 1,
    (defined($workflow) ? (workflow => $workflow_txt) : ()),
    ($opt->indexing_url ? (indexing_url => $opt->indexing_url) : ()),
    ($opt->public ? (public => 1) : ()),
};

my $errs;
my @work;

for my $gb (@ARGV)
{
    if (! -s $gb)
    {
	warn "Genbank file $gb is missing or zero length\n";
	$errs++;
	next;
    }
    my $params = { %$base_params };
    my $base = basename($gb);
    my $out = $base;
    $out =~ s/\.[^.]+$//;
    $params->{genbank_file} = "$params->{output_path}/$base";
    $params->{output_file} = $out;
    push(@work, [abs_path($gb), $params]);
    
}

my $n = @work;
print STDERR "Running $n genomes on " . $opt->parallel . " threads\n";

pareach \@work, sub {
    my($work) = @_;
    my($file, $params) = @$work;

    my $work_dir = File::Temp->newdir(CLEANUP => 0);
    my $here = getcwd();
    print STDERR "Working in $work_dir\n";
    chdir($work_dir) or die "Cannot chdir $work_dir: $!";

    eval {
	$ws->save_file_to_file($file, { genbank_batch => time },
			       $params->{genbank_file},
			       "contigs", $opt->overwrite ? 1 : 0, 1);
    };
    if ($@)
    {
	die "Failed to upload $file to $params->{genbank_file}: $@";
    }

    my $ptmp = File::Temp->new;
    print $ptmp $json->encode($params);
    close($ptmp);
    run(["cat", "$ptmp"]);

    my $log_base = $opt->log_dir . "/$params->{output_file}";
    my $cmd = ["App-$app", "xx", $app_spec, "$ptmp"];

    print STDERR "@$cmd\n";
    my $start = gettimeofday;
    my $ok = run($cmd,
		 ">", "$log_base.out",
		 "2>", "$log_base.err",
		);
    my $end = gettimeofday;
    my $elap = $end - $start;
    chdir($here);
    $ok or die "Failed running $file: $? $!\n";
    print STDERR "Elapsed for $file: $elap\n";
    write_file("$log_base.elapsed", "$elap\n");
    
}, { Num_Workers => $opt->parallel };


sub find_app_spec
{
    my($app) = @_;
    
    my $top = $ENV{KB_TOP};
    my $specs_deploy = "$top/services/app_service/app_specs";

    my $app_spec;
    
    if (-d $specs_deploy)
    {
	$app_spec = "$specs_deploy/$app.json";
    }
    else
    {
	my @specs_dev = glob("$top/modules/*/app_specs");
	for my $path (@specs_dev)
	{
	    my $s = "$path/$app.json";
	    if (-s $s)
	    {
		$app_spec = $s;
		last;
	    }
	}
	-s $app_spec or 
	    die "cannot find specs file for $app in $specs_deploy or @specs_dev\n";
    }
    -s $app_spec or die "Spec file $app_spec does not exist\n";
    return $app_spec;
}
