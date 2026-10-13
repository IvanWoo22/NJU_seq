#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use IO::Zlib;

sub SEQ_REV_COMP {
    my $SEQ = reverse shift;
    $SEQ = uc $SEQ;
    $SEQ =~ tr/U/T/;
    $SEQ =~ tr/ACGTRYMKBDHVSWN/TGCAYRKMVHDBSWN/;
    $SEQ =~ tr/T/U/;
    return $SEQ;
}

sub SEQ_TR_TU {
    my $SEQ = shift;
    $SEQ = uc $SEQ;
    return ($SEQ =~ tr/T/U/r);
}

@ARGV == 4 && $ARGV[2] =~ /^\d+$/ && $ARGV[3] =~ /^\d+$/ or die "Usage: motif_nm.pl genome.fa sites.tsv left_nt right_nt\n";
open( my $FASTA, "<", $ARGV[0] );
my %fasta;
my $title_name;
while (<$FASTA>) {
    if (m/^>(\S+)/) {
        $title_name = $1;
        $title_name =~ s/^chr//i;
    }
    else {
        $_ =~ s/\r?\n//;
        $fasta{$title_name} .= $_;
    }
}
close($FASTA);

open( my $SEG, "<", $ARGV[1] );
my $l_len = $ARGV[2];
my $r_len = $ARGV[3];
while (<$SEG>) {
    s/\r?\n//;
    next if /^#/ || /^(?:Chromosome|Position)\t/ || /^\s*$/;
    my ( $chr, $position, $dir ) = split /\t/;
    $chr =~ s/^chr//i;
    $position =~ /^[1-9][0-9]*$/ && $dir =~ /^[+-]$/ or die "Invalid site row: $_\n";
    my $info = $_;
    my ( $start, $end );
    if ( $dir eq "+" ) {
        $start = $position - $l_len;
        $end   = $position + $r_len;
    }
    else {
        $start = $position - $r_len;
        $end   = $position + $l_len;
    }
    if ( exists( $fasta{$chr} ) ) {
        if ($start < 1 || $end > length($fasta{$chr})) {
            warn "Skipping out-of-bounds motif window: $chr:$start-$end\n";
            next;
        }
        my $length = $end - $start + 1;
        my $seq    = substr( $fasta{$chr}, $start - 1, $length );
        if ( $dir eq "-" ) {
            $seq = SEQ_REV_COMP($seq);
        }
        else {
            $seq = SEQ_TR_TU($seq);
        }
        print("$info\t$seq\n");
    }
    else {
        warn("Sorry, there is no such a segment: $_\n");
    }
}
close($SEG);

__END__
