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
    print "#", join(" ", @ARGV[@selected]), "\n";
    print "$rows->{$_}\n" for @{$groups{$mask}};
    print "\n";
}
