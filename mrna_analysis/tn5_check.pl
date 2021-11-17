#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open( my $FASTA, "<", $ARGV[0] );
open( my $SEG,   "<", $ARGV[1] );

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

sub SEQ_REV_COMP {
    my $SEQ = reverse shift;
    $SEQ =~ tr/Uu/Tt/;
    return ( $SEQ =~ tr/NAGCTagct/NTCGAtcga/r );
}

sub FIND {
    my $TITLE  = shift;
    my $START  = shift;
    my $LENGTH = shift;
    my $REFSEQ = shift;
    my $SEQ;
    if (    ( exists( $fasta{$TITLE} ) )
        and ( length( $fasta{$TITLE} ) >= $START + $LENGTH ) )
    {
        $SEQ = substr( $fasta{$TITLE}, $START - 1, $LENGTH );
    }
    else {
        $SEQ = "NNNNNN";
        #warn("Sorry, there is no such a segment: $TITLE\n");
    }
    my @SEQ       = split( //, $SEQ );
    my @REFSEQ    = split( //, $REFSEQ );
    my $MATCH_NUM = 0;
    foreach ( 0 .. 5 ) {
        $MATCH_NUM++ if $SEQ[$_] eq $REFSEQ[$_];
    }
    return ("$SEQ\t$MATCH_NUM");
}

while (<$SEG>) {
    s/\r?\n//;
    my (
        $qname, undef, $title, $start, undef,
        undef,  undef, undef,  undef,  $seq
    ) = split( /\t/, $_ );
    my ( $id, $umi ) = split( /_/, $qname );
    my ( $fasta_title, undef ) = split( /\|/, $title );
    my $end    = $start + length($seq);
    my $refseq = SEQ_REV_COMP($umi);
    my $output = FIND( $fasta_title, $end, 6, $refseq );
    print("$id\t$umi\t$output\n");
}
close($SEG);

__END__
