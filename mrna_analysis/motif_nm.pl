#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use IO::Zlib;

sub SEQ_REV_COMP {
    my $SEQ = reverse shift;
    $SEQ =~ tr/Tt/Uu/;
    return ( $SEQ =~ tr/AGCUagcu/UCGAucga/r );
}

sub SEQ_TR_TU {
    my $SEQ = shift;
    return ( $SEQ =~ tr/Tt/Uu/r );
}

open( my $FASTA, "<", $ARGV[0] );
my %fasta;
my $title_name;
while (<$FASTA>) {
    if (m/^>(\S+)/) {
        $title_name = $1;
        $title_name =~ s/chr//;
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
    my ( $chr, $position, $dir ) = split /\t/;
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
        my $length = abs( $end - $start ) + 1;
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
