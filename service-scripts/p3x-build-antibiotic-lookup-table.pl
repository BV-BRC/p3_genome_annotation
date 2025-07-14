#
# Query antibiotics collection to build mapping from classifier name to
# BV-BRC antibiotic name
#

use strict;
use P3DataAPI;
use Data::Dumper;
use Text::CSV_XS qw(csv);

my $api = P3DataAPI->new;

my %antibiotic_names;
for my $model_dir (@ARGV)
{
    my $tab = "$model_dir/Antibiotic.abbr.sort2.tab";
    
    my $data = csv(file => $tab, sep_char => "\t", headers => "auto");
    $antibiotic_names{$_->{Antibiotic}} = 1 foreach @$data;
}

my @res = $api->query("antibiotics", ['eq', 'antibiotic_name', '*'], [select => "antibiotic_name"]);
for my $ent (@res)
{
    my $name = $ent->{antibiotic_name};
    my $mapped = $name;
    $mapped =~ s,[/ ],_,g;
    print "$name\t$mapped\n";
    if (!$antibiotic_names{$name})
    {
	print STDERR "No $name\n";
    }
    delete $antibiotic_names{$name};
}

#
# Remaining are the antibiotics from the models that are not in BV-BRC
#
for my $name (sort keys %antibiotic_names)
{
    my $revmap = $name;
    $revmap =~ s/_/ /g;
    print "$revmap\t$name\tmodel_only\n";
}
