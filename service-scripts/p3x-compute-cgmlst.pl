#
# Compute genome MLST using https://github.com/tseemann/mlst
#

use Bio::KBase::AppService::AppConfig qw(application_backend_dir);
use P3DataAPI;
use strict;
use Getopt::Long::Descriptive;  
use File::Temp;
use File::Copy;
use File::Slurp;
use GenomeTypeObject;
use Cwd;
use Data::Dumper;
use Time::HiRes qw(gettimeofday);
use JSON::XS;
use Text::CSV_XS qw(csv);
use IPC::Run qw(run);

my ($opt, $usage) = describe_options("%c %o [< in] [> out]",
                     ["in|i=s", "Input GTO"],
                     ["out|o=s", "Output GTO"],
                     ["parallel=i", "Number of threads to use", { default => 1 }],
                     ["dry_run|n", "Dry run - print commands without executing"],
                     ["help|h", "Print this help message"]);
                     
print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if (@ARGV != 0);

chomp(my $hostname = `hostname -f`);

# helper sub to run or dry-run a command
sub run_cmd {
    my ($cmd_ref, $label) = @_;
    print STDERR "Invoke $label command: @$cmd_ref\n";
    if ($opt->dry_run) {
        print STDERR "DRY RUN - would execute: @$cmd_ref\n";
        return 1;
    }
    my($stdout, $stderr);
    my $ok = run($cmd_ref,
                 ">", \$stdout,
                 "2>", \$stderr);
    $ok or die "Error $? running $label: @$cmd_ref\n";
    return $ok;
}

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

my %taxon_to_schema = (
    470 => "acinetobacter_baumannii",
    1392 => "bacillus_anthracis",
    520 => "bordetella_pertussis",
    235 => "brucella_spp",
    36855 => "brucella_spp",
    120577 => "brucella_spp",
    1218315 => "brucella_spp",
    29459 => "brucella_melitensis",
    444163 => "brucella_spp",
    29460 => "brucella_spp",
    236 => "brucella_spp",
    120576 => "brucella_spp",
    29461 => "brucella_spp",
    13373 => "burkholderia_mallei_fli",
    # burkholderia_mallei_rki => 13373,
    111527 => "burkholderia_pseudomallei",
    197 => "campylobacter_jejuni_coli",
    195 => "campylobacter_jejuni_coli",
    1496 => "clostridioides_difficile",
    1502 => "clostridium_perfringens",
    1717 => "corynebacterium_diphtheriae",
    1719 => "corynebacterium_pseudotuberculosis",
    413503 => "cronobacter_sakazakii_malonaticus",
    28141 => "cronobacter_sakazakii_malonaticus",
    1351 => "enterococcus_faecalis",
    1352 => "enterococcus_faecium",
    562 => "escherichia_coli",
    263 => "francisella_tularensis",
    2058152 => "Klebsiella_oxytoca_grimontii_michiganensis_pasteurii",
    1134687 => "Klebsiella_oxytoca_grimontii_michiganensis_pasteurii",
    571 => "Klebsiella_oxytoca_grimontii_michiganensis_pasteurii",
    2587529 => "Klebsiella_oxytoca_grimontii_michiganensis_pasteurii",
    573 => "Klebsiella_pneumoniae_variicola_quasipneumoniae",
    1463165 => "Klebsiella_pneumoniae_variicola_quasipneumoniae",
    244366 => "Klebsiella_pneumoniae_variicola_quasipneumoniae",
    446 => "legionella_pneumophila",
    1639 => "listeria_monocytogenes",
    33894 => "mycobacterium_tuberculosis_bovis_africanum_canettii",
    1765 => "mycobacterium_tuberculosis_bovis_africanum_canettii",
    78331 => "mycobacterium_tuberculosis_bovis_africanum_canettii",
    77643 => "mycobacterium_tuberculosis_bovis_africanum_canettii",
    36809 => "mycobacteroides_abscessus",
    2096 => "mycoplasma_gallisepticum",
    1464 => "paenibacillus_larvae",
    287 => "pseudomonas_aeruginosa",
    28901 => "salmonella_enterica",
    615 => "serratia_marcescens",
    #staphylococcus_argenteus => 985002,
    1280 => "staphylococcus_aureus",
    29388 => "staphylococcus_capitis",
    1314 => "spyogenes",
    630 => "yersinia_enterocolitica"
);


# grab all lineage 
my $api = P3DataAPI->new();
my $taxon_id = $gto_in->{ncbi_taxonomy_id};
my @res = $api->query("taxonomy", ['eq', 'taxon_id', $taxon_id], ['select', 'taxon_name', 'lineage_ids', 'lineage_names']);

# Query schema_map hash (taxon_id => schema dir name)
my $dir_name;
foreach my $lineage_id (reverse @{$res[0]->{lineage_ids} // [] }) {
    $lineage_id = int($lineage_id);
    if (exists $taxon_to_schema{$lineage_id}) {
        $dir_name = $taxon_to_schema{$lineage_id};
        # nb dev
        last;
    }
}

print STDERR "dir_name: $dir_name\n";

# based on the schema name, determine there is a cgmlst scheme to use
if ($dir_name) {
    my $tmp_dir = File::Temp->newdir(CLEANUP => 0);
    my $clean_fasta_dir = "$tmp_dir/clean_fastas";
    my $allele_call_out = "$tmp_dir/new_genomes_allele_call/";
    my $schema_path = application_backend_dir . "/CoreGenomeMLST/chewbbaca_schemas/" . $dir_name;
    
    mkdir $clean_fasta_dir or die "Cannot create clean fastas directory: $!\n";
    
    # rename to be certain it ends with 'fasta' for chewBBACA
    run_cmd(["cp", "$in_file", "$clean_fasta_dir/input.fasta"], "Add extension");

    # new genome allele call
    my @allele_call_cmd = (
        "chewBBACA.py", "AlleleCall",
        "--input-files", $clean_fasta_dir,
        "--schema-directory", $schema_path,
        "--output-directory", $allele_call_out,
        "--cpu", "4",
        "--output-unclassified",
        "--output-missing",
        "--output-novel",
        "--no-inferred"
    );
    run_cmd(\@allele_call_cmd, "Allele Call");
    
    # Check percent of loci have an exact allele match
    my $cluster_threshold  = 85;  # run clustering
    my $qc_threshold     = 70;  # god vs poor
    my $new_allele_call = "$tmp_dir/new_genomes_allele_call/results_alleles.tsv";
    my ($total, $exact, $pct) = (0, 0, 0);
    my $qc = "poor";              # good|poor
    my $do_cluster = 0;           # 1 => run clustering

    if ($opt->dry_run) {
        print STDERR "DRY RUN - would check allele call quality from $new_allele_call\n";
        $total = 100;
        $exact = 100;
        $pct = 100; 
    } else {
        open(my $check_fh, '<', $new_allele_call)
            or die "Cannot open $new_allele_call: $!";
        my $hdr = <$check_fh>;   # skip header
        my $data = <$check_fh>;  # data line
        close($check_fh);

        chomp $data;
        my @vals = split(/\t/, $data);
        shift @vals;  # remove first column (filename)

        $total = scalar @vals;
        $exact = grep { /^\d+$/ && $_ > 0 } @vals;
        $pct   = $total > 0 ? ($exact / $total) * 100 : 0;
        }

        $qc = ($pct >= $qc_threshold) ? "good" : "poor";
        $do_cluster = ($pct >= $cluster_threshold) ? 1 : 0;

        print STDERR sprintf(
        "Allele call quality: %d / %d exact matches (%.1f%%) => qc=%s, cluster=%s\n",
        $exact, $total, $pct, $qc, ($do_cluster ? "yes" : "no")
        );
    
        # Parse allele calls into a comma string ONCE (always saved)
        my $allele_call_string;
        if ($opt->dry_run) {
            print STDERR "DRY RUN - would parse allele calls from $new_allele_call\n";
            $allele_call_string = "DRY_RUN_ALLELE_CALLS";
        } else {
            open(my $allele_fh, '<', $new_allele_call) or die "Cannot open $new_allele_call: $!";
            my $header_line = <$allele_fh>;
            my $data_line   = <$allele_fh>;
            close($allele_fh);

            die "No allele-call data line found when parsing string from $new_allele_call\n"
                unless defined($data_line);

            chomp $data_line;
            $data_line =~ s/\t/,/g;        # tabs -> commas
            $data_line =~ s/^[^,]*,//;     # remove first column (filename)
            $allele_call_string = $data_line;
            print STDERR "Allele call string: $allele_call_string\n";
        }

         # Build analysis event once
        my $event = {
            tool_name      => "p3x-compute-cgmlst",
            parameters     => [map { "$_" } @allele_call_cmd],
            execution_time => scalar gettimeofday,
            hostname       => $hostname,
        };
        my $event_id = $gto_in->add_analysis_event($event);

        my $sequence_typing = {
            typing_method => "sequence_typing",
            allele_calls  => $allele_call_string,
            schema_name   => $dir_name,
            event_id      => $event_id,
            loci_total    => $total,
            loci_called   => $exact,
            loci_missing  => $total - $exact,
            pct_called    => $pct,
            qc            => $qc, 
        };

            if (!$do_cluster) {
        print STDERR "Skipping clustering — allele call quality below ${cluster_threshold}%.\n";
        push(@{$gto_in->{sequence_types}}, $sequence_typing);
        } else {

            # Join new allele call with master allele call using chewBBACA JoinProfiles
            my $master_table = application_backend_dir . "/CoreGenomeMLST/precomputed_clusters/refs/" . $dir_name . "_11_25_2025_joined.tsv";
            my $master_joined = "$tmp_dir/master_joined.tsv";
            my $copy_of_master = "$tmp_dir/master_copy.tsv";

            run_cmd(["cp",  $master_table, $copy_of_master], "Copy Master Table");

            run_cmd([
                "chewBBACA.py", "JoinProfiles",
                "--profiles", $copy_of_master, $new_allele_call,
                "--output-file", $master_joined
            ], "Join Profiles");

            run_cmd([
                "python3", "/home/nbowers/bvbrc-dev/dev_container/modules/bvbrc_CoreGenomeMLST/service-scripts/core-genome-mlst-utils.py",
                "clean-allelic-profile", $master_joined
            ], "Clean Allele Call");

            my $precomputed_clusters_dir  = application_backend_dir . "/CoreGenomeMLST/precomputed_clusters/";
            my $precomputed_clusters_path = $precomputed_clusters_dir . lc($dir_name) . ".cgMLSTv1.npz";
            my $local_clusters_path       = "$tmp_dir/precomputed_clusters.npz";
            my $heircc_out                = "$tmp_dir/cluster";
            run_cmd(["cp", $precomputed_clusters_path, $local_clusters_path], "Copy Master npz");

            run_cmd([
                "python3", "/home/nbowers/bvbrc-dev/dev_container/cgmlst_for_all/dev_clustering/testing_heirCC/testing_pHierCC/git_repo/pHierCC/pHierCC.py",
                "--profile", "$tmp_dir/master_joined_clean.tsv",
                "--output",  $heircc_out,
                "--append",  $precomputed_clusters_path
            ], "Cluster");

            my $heircc_tsv = "$tmp_dir/cluster.HierCC.gz";
            run_cmd(["gunzip", $heircc_tsv], "gunzip");

            # Parse HierCC output
            my $heircc_unzipped = "$tmp_dir/cluster.HierCC";
            my %cluster_kv;
                    if ($opt->dry_run) {
            %cluster_kv = (HC1 => "DRY_RUN", HC5 => "DRY_RUN");
        } else {
            open(my $heircc_fh, '<', $heircc_unzipped) or die "Cannot open $heircc_unzipped: $!";
            my $heircc_header = <$heircc_fh>;
            chomp $heircc_header;
            my @heircc_keys = split(/\t/, $heircc_header);

            my $heircc_data;
            while (my $line = <$heircc_fh>) {
                chomp $line;
                my @cols = split(/\t/, $line);
                if ($cols[0] eq 'input') {
                    $heircc_data = $line;
                    last;
                }
            }
            close($heircc_fh);

            die "Could not find input row in $heircc_unzipped\n" unless $heircc_data;
            my @heircc_values = split(/\t/, $heircc_data);

            # change all hc levels
            # shift @heircc_values;  # drop genome id col
            # shift @heircc_keys;    # drop genome id col
            # @cluster_kv{@heircc_keys} = @heircc_values;
            shift @heircc_values;
            shift @heircc_keys;

            my %cgmlst_hc;
            for my $i (0 .. $#heircc_keys) {
                my $k = $heircc_keys[$i];
                my $v = $heircc_values[$i];

                # keep only HC columns
                next unless $k =~ /^HC(\d+)$/i;

                my $level = $1;

                # optional: restrict to specific HC levels
                next unless $level =~ /^(0|2|5|10|20|50|100)$/;

                $cgmlst_hc{"cgmlst_hc$level"} = $v;
        }

        # Attach clustering results and push sequence_typing
        # if we want all HC levels
        #$sequence_typing->{cluster_key_value_pairs} = \%cluster_kv;
        $sequence_typing->{cgmlst_hc} = \%cgmlst_hc;
        push(@{$gto_in->{sequence_types}}, $sequence_typing);
        }
    }

} else {
    print STDERR "Schema does not exist for this species\n";
}

if ($opt->out)
{
    $gto_in->destroy_to_file($opt->out);
}
else
{
    $gto_in->destroy_to_file(\*STDOUT);
}