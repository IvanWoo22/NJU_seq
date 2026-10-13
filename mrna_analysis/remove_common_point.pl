#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
use NJUSeq::SiteTable qw(read_sites);
my ($rows, $members, $by_file) = read_sites(@ARGV);
for my $i (0 .. $#ARGV) {
    my @points = grep { scalar(keys %{$members->{$_}}) < @ARGV } sort keys %{$by_file->[$i]};
    next unless @points;
    print "#$ARGV[$i]\n";
    print "$by_file->[$i]{$_}\n" for @points;
    print "\n";
}
