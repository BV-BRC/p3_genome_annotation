#
# Query antibiotics collection to build mapping from classifier name to
# BV-BRC antibiotic name
#

use strict;
use P3DataAPI;

my $api = P3DataAPI->new;

print "Query\n";
my @res = $api->query("antibiotics", ['eq', 'antibiotic_name', '*'], [select => "antibiotic_name"]);
for my $ent (@res)
{
    my $name = $ent->{antibiotic_name};
    my $mapped = $name;
    $mapped =~ s,[/ ],_,g;
    print "$name\t$mapped\n";
}
