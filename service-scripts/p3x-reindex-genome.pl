#
# Resubmit the load files from the given genome to Solr.
#
# We verify that the genome ID does not already exist.
#

use strict;
use Data::Dumper;
use File::Slurp;
use IPC::Run qw(run);
use JSON::XS;
use SolrAPI;
use P3DataAPI;
use P3AuthToken;
use Bio::KBase::AppService::AppConfig 'data_api_url';
use Getopt::Long::Descriptive;
use File::Temp;
use Bio::P3::Workspace::WorkspaceClientExt;

my($opt, $usage) = describe_options("%c %o ws-path",
				    ["execute" => "Perform the reload"],
				    ['skip-file=s@', "Skip loading the given file", { default => [] }],
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 1 if $opt->help;
die($usage->text) if @ARGV != 1;

my %skip = map { $_ => 1 } @{$opt->skip_file};

my $ws_path = shift;

my $api = P3DataAPI->new();
my $ws = Bio::P3::Workspace::WorkspaceClientExt->new();

my $path = "$ws_path/load_files";
my $res = $ws->ls({paths => [$path]});
my $files = $res->{$path};

my $tmp = File::Temp->newdir(CLEANUP => 1);
print "$tmp\n";

for my $file (@$files)
{
    my $wsp = "$path/$file->[0]";
    $ws->download_file($wsp, "$tmp/$file->[0]", 1);
}

my $genome = decode_json(read_file("$tmp/genome.json"));
if (!$genome)
{
    die "Could not find genome info\n";
}

my $genome_id = $genome->[0]->{genome_id};

my @res = $api->query("genome", ["eq", "genome_id", $genome_id], ["select", "genome_id,genome_name"]);
if (@res)
{
    die "Genome already loaded: " . Dumper(\@res);
}


if (!$opt->execute)
{
    print "Would reindex $genome_id if --execute was specified\n";
    exit;
}


print "Reloading $genome_id\n";

my $genome_url = data_api_url . "/indexer/genome";
my $token = P3AuthToken->new();
my $token_str = $token->token;
my @opts;
push(@opts, "-H", "Authorization: $token_str");
push(@opts, "-H", "Content-Type: multipart/form-data");

my @files = ([genome => "genome.json"],
	     [genome_feature => "genome_feature.json"],
	     [genome_amr => "genome_amr.json"],
	     [genome_sequence => "genome_sequence.json"],
	     [pathway => "pathway.json"],
	     [subsystem => "subsystem.json"],
	     [feature_sequence => "feature_sequence.json"],
	     [sp_gene => "sp_gene.json"],
	     [taxonomy => "taxonomy.json"]);

for my $tup (@files)
{
    my($key, $file) = @$tup;
    my $path = "$tmp/$file";
    if (-f $path)
    {
	if ($skip{$key})
	{
	    print STDERR "Skipping $key - $path\n";
	    next;
	}
	push(@opts, "-F", "$key=\@$path");
    }
}

push(@opts, "-D", "-", "-v");

push(@opts, $genome_url);

#curl -H "Authorization: AUTHORIZATION_TOKEN_HERE" -H "Content-Type: multipart/form-data" -F "genome=@genome.json" -F "genome_feature=@genome_feature_patric.json" -F "genome_feature=@genome_feature_refseq.json" -F "genome_feature=@genome_feature_brc1.json" -F "genome_sequence=@genome_sequence.json" -F "pathway=@pathway.json" -F "sp_gene=@sp_gene.json"  

my($stdout, $stderr);

print Dumper(\@opts);
my $ok = run(["curl", @opts], '>', \$stdout);
if (!$ok)
{
    warn "Error $? invoking curl @opts\n";
}

my $json = JSON->new->allow_nonref;

my $data;
eval {
    $data = $json->decode($stdout);
};
if ($@)
{
    die "Bad submission - output did not parse: $@\n$stdout\n";
}

my $queue_id = $data->{id};

print "Submitted indexing job $queue_id\n";

my $solr = SolrAPI->new(data_api_url);

#
# For now, wait up to an hour for the indexing to complete.
#
my $wait_until = time + 3600;

while (time < $wait_until)
{
    my $status = $solr->query_rest("/indexer/$queue_id");
    if (ref($status) ne 'HASH')
    {
	warn "Parse failed for indexer query for $queue_id: " . Dumper($status);
    }
    else
    {
	my $state = $status->{state};
	print STDERR "status for $queue_id (state=$state): " . Dumper($status);
	if ($state ne 'queued')
	{
	    print STDERR "Finishing with state $state\n";
	    last;
	}
    }
    sleep 60;
}

