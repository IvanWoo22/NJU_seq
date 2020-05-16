#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Statistics::Basic qw(median);

open( my $IN1, "<", $ARGV[0] );
open( my $IN2, "<", $ARGV[1] );
open( my $IN3, "<", $ARGV[2] );

sub LIST_INPUT {
    my $IN = shift;
    my ( @ARRAY_SCORE, @ARRAY_LIST, %INFO, %SCORE );
    while (<$IN>) {
        chomp;
        my (
            $CHR,   $POS,     $DIR,  $BASE1, $BASE2,
            $BASE3, $GENE_ID, undef, undef,  $SCORE
        ) = split /\t/;
        push( @ARRAY_LIST,  $CHR . "\t" . $POS . "\t" . $DIR );
        push( @ARRAY_SCORE, $SCORE );
        $INFO{ $CHR . "\t" . $POS . "\t" . $DIR } =
            $CHR . "\t"
          . $POS . "\t"
          . $DIR . "\t"
          . $BASE1 . "\t"
          . $BASE2 . "\t"
          . $BASE3 . "\t"
          . $GENE_ID;
        $SCORE{ $CHR . "\t" . $POS . "\t" . $DIR } = $SCORE;
    }
    my $MEDIAN             = median(@ARRAY_SCORE);
    my @ARRAY_SCORE_MEDIAN = map( { abs($_ - $MEDIAN) } @ARRAY_SCORE );
    my $MAD                = median(@ARRAY_SCORE_MEDIAN);
    return ( $MAD );
}

my ($a_mad) = LIST_INPUT($IN1);
my ($b_mad) = LIST_INPUT($IN2);
my ($c_mad) = LIST_INPUT($IN3);

print("$a_mad\n$b_mad\n$c_mad\n");

