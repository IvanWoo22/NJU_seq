#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use AlignDB::IntSpan;

# Explicitly choose transcript-first/last/both; no inferred "patch" semantics.
use Getopt::Long;
my $which;
GetOptions('which=s' => \$which) or die "Invalid options\n";
defined($which) && $which =~ /^(?:first|last|both)$/
    or die "Specify --which first|last|both (in transcript orientation)\n";
while (<STDIN>) {
    chomp;
    my @f = split /\t/;
    @f >= 5 && $f[4] =~ /^[+-]$/ or die "Expected exon distance detail with strand in column 5\n";
    my $set = AlignDB::IntSpan->new($f[3]);
    my $first = $set->at_island($f[4] eq '+' ? 1 : -1);
    my $last = $set->at_island($f[4] eq '+' ? -1 : 1);
    my $island = $set->find_islands($f[0]);
    print "$_\n" if (($which eq 'first' || $which eq 'both') && $first eq $island)
        || (($which eq 'last' || $which eq 'both') && $last eq $island);
}
