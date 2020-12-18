#!/usr/bin/env perl

#
# Build the VigorTaxonMap.pm module from a template based on the tabular text file.
#

use strict;
use Template;
use Data::Dumper;
@ARGV == 3 or die "Usage: $0 map-file template output\n";

my $map = shift;
my $template_file = shift;
my $out_file = shift;

-f $template_file or die "Template file $template_file is missing\n";
open(M, "<", $map) or die "Cannot open $map: $!";

my @map;
while (<M>)
{
    if (/^(\S+)_db\t([^\t]+)\t(\d+)$/)
    {
	push(@map, { db => $1, name => $2, taxon => $3 });
    }
}
close(M);

my %vars = ( map => \@map );

open(O, ">", $out_file) or die "Cannot write $out_file: $!";

my $templ = Template->new();
print Dumper(\%vars);
my $ok = $templ->process($template_file, \%vars, \*O);
$ok or die "Error processing $template_file: " . $templ->error();
