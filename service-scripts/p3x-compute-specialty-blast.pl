#
# Compute specialty gene assignments using blast or blat.
#


use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;  
use File::Temp qw(tempfile);
use File::Copy;
use GenomeTypeObject;
use Cwd;
use IPC::Run qw(run);
use Time::HiRes qw(gettimeofday);
use File::Slurp;
use JSON::XS;
use File::Basename;
use Bio::SearchIO;

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
				     ["tabular=s", "Tabular output file"],
				     ["report|r=s", "Raw BLAST output file"],
				     ["program=s", "Program to use (blast or blat)", { default => 'blat' }],
				     ["description=s", "Show descriptions", { default => 'N' }],
				     ["db=s@", "Database name"],
				     ["db-dir=s", "Database directory"],
				     ["filter=s", "Filter outputs (Y or N)", { default => 'Y' }],
				     ["evalue=s", "Base evalue to use for analysis", { default => 0.0001 }],
				     ["parallel=i", "Number of threads to use", { default => 1 }],
				     ["help|h", "Print this help message"],
				    { show_defaults => 1 });

print($usage->text), exit 0 if $opt->help;
die($usage->text) if (@ARGV != 0);

#
# Since we don't use human homologs any more, the required parameters are as follows.
#
my $min_qcov = 80;
my $min_scov = 80;
my $min_iden = 80;

if ($opt->program ne 'blast' && $opt->program ne 'blat' && $opt->program ne 'diamond')
{
    die "Invalid value for --program option (must be either 'blast', 'blat', or 'diamond')\n";
}

my $db_dir = $opt->db_dir;

$db_dir or die "Database directory must be specified.\n";
-d $db_dir or die "Database directory $db_dir not found.\n";

chomp(my $hostname = `hostname -f`);

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

my $in_file = "$work_dir/blast_input";

$gto_in->write_protein_translations_to_file($in_file);

my @db = @{$opt->db || []};

#
# Enumerate available databases;
#

my @avail_dbs = <$db_dir/*.faa>;
my @avail_names = map { basename($_, '.faa') } @avail_dbs;
my %avail_names = map { $_ => 1 } @avail_names;

if (@db == 0 || grep { $_ eq 'all' } @db)
{
    @db = @avail_names;
}
else
{
    my @err;
    for my $d (@db)
    {
	if (!$avail_names{$d})
	{
	    push @err, "Database $d is not available\n";
	    next;
	}
    }
    die join("", @err) if @err;
}

my $raw_blast;
my $raw_blast_fh;
if ($opt->report)
{
    $raw_blast = $opt->report;
    unlink($raw_blast);
    open($raw_blast_fh, ">". $raw_blast) or die "Cannot write $raw_blast: $!";
}
else
{
    ($raw_blast_fh, $raw_blast) = tempfile();
}

my($tabular_fh, $tabular_file);
if ($opt->tabular)
{
    $tabular_file = $opt->tabular;
    open($tabular_fh, ">", $tabular_file) or die "Cannot write $tabular_file: $!";

    my $header;
    if ($opt->description eq 'Y') {
	$header = "QID\tQAnnotation\tQOrganism\tDatabase\tSubID\tSubAnnotation\tSubOrganism\tQueryCoverage\tSubCoverage\tIdentity\tP-Value\n";
    } else {
	$header = "QID\tDatabase\tSubID\tQueryCoverage\tSubCoverage\tIdentity\tP-Value\n" ;
    }
    print $tabular_fh $header;
}

my %db_event;

for my $db (@db)
{
    my $db_file = "$db_dir/$db";
    -f "$db_file.faa" or die "DB file $db_file.faa missing\n";

    print STDERR "Compute $db\n";

    my @cmd;
    if ($opt->program eq 'blast')
    {
	@cmd = ("blastp", "-query", $in_file, "-db", "$db_file.faa",
		   "-num_threads", $opt->parallel,
		   "-num_alignments", 1,
		   "-num_descriptions", 1,
		   "-evalue", $opt->evalue,
		   "-outfmt", 0);
    }
    elsif ($opt->program eq 'diamond')
    {
	@cmd = ("diamond", "blastp",
		"--query", $in_file,
		"--db", "$db_file.faa.dmnd",
		"--outfmt", "0",
		"--evalue", $opt->evalue,
		"--id", $min_iden,
		"--sensitive",
		"--threads", $opt->parallel);
    }
    else
    {
	@cmd = ("blat", "-prot", "$db_file.faa", $in_file, "-out=blast", "/dev/fd/1");
    }

    my $event = {
	tool_name => "p3x-compute-specialty-blast",
	parameters => \@cmd,
	execution_time => scalar gettimeofday,
	hostname => $hostname,
    };
    $db_event{$db} = $gto_in->add_analysis_event($event);

    print STDERR "@cmd\n";
    my $ok = run(\@cmd, ">>", $raw_blast_fh);
    
    if (!$ok)
    {
	die "Blat failed with $?: $@cmd\n";
    }
}
close($raw_blast_fh);

my $blast_data = new Bio::SearchIO(-format => 'blast', -file=> $raw_blast);

my %val;
while( my $result = $blast_data->next_result )
{
    # $result is a Bio::Search::Result::ResultI compliant object

    my $hit = $result->next_hit;
    if ($hit)
    {
	# $hit is a Bio::Search::Hit::HitI compliant object
	
	my $hsp = $hit->next_hsp;
	if ($hsp)
	{
	    # $hsp is a Bio::Search::HSP::HSPI compliant object
	    
	    my $qid = $result->query_name;
	    my ($qAnnot, $qOrg) = ($result->query_description, "");
	    ($qAnnot, $qOrg) = $result->query_description=~/(.*)\s*\[(.*)\]/ if $result->query_description=~/(.*)\s*\[(.*)\]/;
	    
	    my $database=$result->database_name;

	    my $sid = $hit->name;
	    my $dbname;
	    #
	    # Diamond doesn't emit database name
	    #
	    if (!$database)
	    {
		($dbname) = $sid =~ /^([^|]+)/;
	    }
	    else
	    {
		$database =~ s/\s*$//;
		$database =~ s/^\s*//;
		
		$dbname = basename($database, '.faa');
	    }


	    my ($sAnnot, $sOrg) = ($hit->description, "");
	    ($sAnnot, $sOrg) = $hit->description=~/(.*)\s*\[(.*)\]/ if $hit->description=~/(.*)\s*\[(.*)\]/;
	    
	    my $iden=int(abs($hsp->percent_identity));	
	    my $qcov=int(abs( $hsp->length('query') * 100 / $result->query_length ));
	    my $scov=int(abs( $hsp->length('hit') * 100 / $hit->length ));
	    my $pvalue = $hsp->significance;

	    # print "INIT $dbname $qid $sid $qcov $scov $iden $pvalue " . $hsp->length('query') . " " . $hsp->length('hit') . " " . $result->query_length ." " . $hit->length . "\n";
	    $val{$dbname}->{$qid} = [$qid, $dbname, $sid, $qcov, $scov, $iden, $pvalue, $qAnnot, $qOrg, $sAnnot, $sOrg];
	}
    }
}
undef $blast_data;
if (!$opt->report)
{
    unlink($raw_blast);
}
# print STDERR Dumper(\%val);

for my $db (@db)
{
    my $v = $val{$db};
    my $event_id = $db_event{$db};

    my $lookup = eval { decode_json(scalar read_file("$db_dir/$db.json")); };
    $lookup or die "Failure reading json: $@";
    my %lookup = map { $_->{source_id} => $_ } @$lookup;
    
    while (my($id, $ent) = each %$v)
    {
	my(undef, undef, $sid, $qcov, $scov, $iden, $pvalue, $qAnnot, $qOrg, $sAnnot, $sOrg) = @$ent;
	#  print join("\t", $id, $sid, $qcov, $scov, $iden), "\n";
	if ($sid)
	{
	    if ($opt->filter eq 'Y')
	    {
		next unless ($qcov >= $min_qcov || $scov >= $min_scov) && ($iden >= $min_iden);
		    
		# if ($db =~ /Human/i)
		# {
		#     next unless ($qcov>=25 || $scov>=25) && $iden>=40;
		# }
		# else
		# {
		#     next unless ($qcov>=80 || $scov>=80) && $iden>=80 ;
		# }
	    }
	    if ($tabular_fh)
	    {
		my $sim = "";
		if ($opt->description eq 'Y') {
		    $sim = "$id\t$qAnnot\t$qOrg\t$db\t$sid\t$sAnnot\t$sOrg\t$qcov\t$scov\t$iden\t$pvalue\n";
		} else {
		    $sim = "$id\t$db\t$sid\t$qcov\t$scov\t$iden\t$pvalue\n";
		}
		print $tabular_fh $sim;
	    }
	}

	
	my $feat = $gto_in->find_feature($id);
	if (ref($feat))
	{
	    my($ident) = $sid =~ /\|(.*)$/;
	    my $info = $lookup{$ident};
	    my @keys = qw(source locus_tag property organism source_id product classification);
	    my $assoc = [$db, $sid, $qcov, $scov, $iden, $pvalue,
		     { (map { $_ => $info->{$_} } grep { $info->{$_} }  @keys),
		       event_id => $event_id, tool => $opt->program }];
	    push(@{$feat->{similarity_associations}}, $assoc);
	}
    }
}
	
close($tabular_fh) if $tabular_fh;


if ($opt->out)
{
    $gto_in->destroy_to_file($opt->out);
}
else
{
    $gto_in->destroy_to_file(\*STDOUT);
}
