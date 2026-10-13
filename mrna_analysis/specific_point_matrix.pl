#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
use NJUSeq::SiteTable qw(read_sites);
my ($rows, $members) = read_sites(@ARGV);
my %groups;
for my $key (sort keys %$rows) {
    my $mask = join ",", sort { $a <=> $b } keys %{$members->{$key}};
    push @{$groups{$mask}}, $key;
}
for my $mask (sort keys %groups) {
    my @selected = split /,/, $mask;
    my %selected = map { $_ => 1 } @selected;
    my $matrix = join "\t", map { $selected{$_} // 0 } 0 .. $#ARGV;
    print "$rows->{$_}\t$matrix\n" for @{$groups{$mask}};
}
