#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use AlignDB::IntSpan;
use Getopt::Long;

# Six-column regions: chr,start,end,unused,gene,strand. Default is 1-based
# inclusive; --bed converts zero-based half-open coordinates explicitly.
my $bed = 0;
GetOptions('bed' => \$bed) or die "Invalid options\n";
@ARGV == 3 or die "Usage: judge_splice_exon.pl [--bed] ASE.tsv CSE.tsv sites.tsv\n";
my (%regions, %by_chr);
for my $kind (0, 1) {
    open my $fh, '<', $ARGV[$kind];
    while (<$fh>) {
        chomp;
        next if /^#/ || /^\s*$/;
        my ($chr,$a,$b,undef,$gene,$strand) = split /\t/;
        defined($strand) && $a =~ /^\d+$/ && $b =~ /^[1-9][0-9]*$/ && $strand =~ /^[+-]$/
            or die "Invalid six-column region in $ARGV[$kind]\n";
        $a++ if $bed;
        $a >= 1 && $b >= $a or die "Invalid region coordinates (use --bed for BED)\n";
        $chr =~ s/^chr//;
        my $group = "$chr\t$strand";
        my $key = "$group\t$gene";
        $regions{$key} //= [AlignDB::IntSpan->new, AlignDB::IntSpan->new];
        $regions{$key}[$kind]->add_range($a,$b);
        $by_chr{$group}{$key} = 1;
    }
    close $fh;
}
my (@count, @length);
@count = @length = (0,0);
my (%seen_site, %seen_gene);
open my $sites, '<', $ARGV[2];
while (<$sites>) {
    chomp;
    next if /^#/ || /^(?:Chromosome|Position)\t/ || /^\s*$/;
    my ($chr,$pos,$strand) = split /\t/;
    $pos =~ /^[1-9][0-9]*$/ && $strand =~ /^[+-]$/ or die "Invalid site row\n";
    $chr =~ s/^chr//;
    next if $seen_site{"$chr\t$pos\t$strand"}++;
    my @hits = (0,0);
    for my $key (keys %{$by_chr{"$chr\t$strand"} // {}}) {
        my @in = map { $regions{$key}[$_]->contains_all($pos) ? 1 : 0 } (0,1);
        next unless $in[0] || $in[1];
        $hits[$_] += $in[$_] for 0,1;
        unless ($seen_gene{$key}++) {
            $length[$_] += $regions{$key}[$_]->size for 0,1;
        }
    }
    if ($hits[0]+$hits[1]) { $count[$_] += $hits[$_]/($hits[0]+$hits[1]) for 0,1; }
}
print join("\n", @count, @length), "\n";
