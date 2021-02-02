#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open( my $ANNO,  "<", $ARGV[0] );
open( my $SCORE, "<", $ARGV[1] );

my %geneanno;

readline($ANNO);
while (<$ANNO>) {
    chomp;
    my @tmp = split "\t";
    foreach my $i ( 0 .. 12 ) {
        if ( !exists( $tmp[$i] ) ) {
            $tmp[$i] = "";
        }
    }
    if ( !exists( $geneanno{ $tmp[1] } ) ) {
        $geneanno{ $tmp[1] } = join( "\t", @tmp[ 4 .. 12 ] );
    }
}
close($ANNO);

while (<$SCORE>) {
    chomp;
    my @tmp = split "\t";
    my ( $id, undef ) = split( "/", $tmp[6] );
    print("$_\t$geneanno{$id}\n");
}
close($SCORE);

__END__
