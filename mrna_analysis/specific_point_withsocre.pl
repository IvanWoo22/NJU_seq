#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
use NJUSeq::SiteTable qw(read_sites);
# Historical 26-column annotated score tables only; not current score.pl output.
my ($rows, $members, $by_file) = read_sites(@ARGV);
for my $table (@$by_file) {
    for my $row (values %$table) {
        scalar(split /\t/, $row, -1) == 26 or die "specific_point_withsocre.pl requires legacy 26-column input; use specific_point_matrix.pl for current score tables\n";
    }
}
my %groups;
for my $key (keys %$rows) {
    my @selected = sort { $a <=> $b } keys %{$members->{$key}};
    my @info = split /\t/, $rows->{$key}, -1;
    my (@combined, $score);
    $score = 0;
    for my $i (@selected) {
        my @f = split /\t/, $by_file->[$i]{$key}, -1;
        $f[19] =~ /^-?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$/ or die "Invalid legacy score\n";
        $score += $f[19];
        for my $j (0 .. 10) { push @{$combined[$j]}, $f[$j+9]; }
    }
    my $line = join "\t", @info[0..8,20..25], (map { join ";", @$_ } @combined), $score;
    push @{$groups{join ",", @selected}}, [$score, $line];
}
for my $mask (sort keys %groups) {
    print "#", join(";", @ARGV[split /,/, $mask]), "\n";
    print "$_->[1]\n" for sort { $b->[0] <=> $a->[0] || $a->[1] cmp $b->[1] } @{$groups{$mask}};
}
