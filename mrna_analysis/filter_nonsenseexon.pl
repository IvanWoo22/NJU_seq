#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long;

Getopt::Long::GetOptions(
    'help|h'         => sub { Getopt::Long::HelpMessage(0) },
    'transwording=s' => \my $mRNA,
    'transid=s'      => \my $trans_id,
    'exonid=s'       => \my $exon_id,
) or Getopt::Long::HelpMessage(1);

# Collect selected transcript IDs first; GFF child order need not be adjacent.
my @lines = <STDIN>;
my %transcripts;
for my $line (@lines) {
    next if $line =~ /^#/;
    my @f = split /\t/, $line;
    next unless @f == 9 && $f[2] eq $mRNA;
    $f[8] =~ /$trans_id([^;\s]+)/ or die "Missing transcript ID: $line";
    $transcripts{$1} = 1;
}
for my $line (@lines) {
    next if $line =~ /^#/;
    my @f = split /\t/, $line;
    next unless @f == 9;
    if ($f[2] eq "gene" || $f[2] eq $mRNA) { print $line; }
    elsif ($f[2] eq "exon") {
        $f[8] =~ /$exon_id([^;\s]+)/ or die "Missing exon parent: $line";
        print $line if grep { exists $transcripts{$_} } split /,/, $1;
    }
}
