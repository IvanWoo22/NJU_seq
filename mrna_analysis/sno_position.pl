#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use PerlIO::gzip;

sub SEQ_REV_COMP {
    my $SEQ = reverse shift;
    $SEQ =~ tr/Uu/Tt/;
    return ( $SEQ =~ tr/AGCTagct/TCGAtcga/r );
}

sub SEQ_TR_TU {
    my $SEQ = shift;
    return ( $SEQ =~ tr/Uu/Tt/r );
}

open( my $FASTA, "<:gzip", $ARGV[0] );
open( my $SEG,   "<",      $ARGV[1] );

my %fasta;
my $title_name;
while (<$FASTA>) {
    if (/^>(\S+)/) {
        $title_name = $1;
    }
    else {
        $_ =~ s/\r?\n//;
        $fasta{$title_name} .= $_;
    }
}
close($FASTA);

while (<$SEG>) {
    s/\r?\n//;
    my ( $chr, $dir, $position, undef, $base, undef ) = split /\s+/;
    my $info = $_;
    my ( $start, $end );
    if ( $dir eq "+" ) {
        $start = $position - 2;
        $end   = $position + 5;
    }
    else {
        $start = $position - 5;
        $end   = $position + 2;
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
        $seq = SEQ_REV_COMP($seq);
        print("$info\t$seq\n");
    }
    else {
        warn("Sorry, there is no such a segment: $_\n");
    }
}
close($SEG);

__END__
