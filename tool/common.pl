#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use FindBin;
use lib "$FindBin::Bin/../lib";
use NJUSeq::SiteTable qw(read_sites);
use Getopt::Long;
my $score_col;
GetOptions('score-col=i' => \$score_col) or die "Invalid options\n";
@ARGV == 3 or die "Usage: common.pl [--score-col N] A.tsv B.tsv C.tsv\n";
defined($score_col) && $score_col < 1 and die "Score column must be positive\n";
my ($rows, $members, $by_file) = read_sites(@ARGV);
for my $key (sort keys %$rows) {
    next unless keys(%{$members->{$key}}) == 3;
    my $sum = 0;
    for my $table (@$by_file) {
        my @f = split /\t/, $table->{$key};
        my $col = $score_col // (@f == 10 ? 10 : (@f >= 15 && (@f-12)%3 == 0 ? @f-1 : 0));
        $col && defined($f[$col-1]) && $f[$col-1] =~ /^-?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$/
            or die "Unknown score schema; provide --score-col N\n";
        $sum += $f[$col-1];
    }
    my @f = split /\t/, $rows->{$key};
    print join("\t", @f[0..6], $sum), "\n";
}
