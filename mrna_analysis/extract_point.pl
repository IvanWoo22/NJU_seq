#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open( my $IN1, "<", $ARGV[0] );
open( my $IN2, "<", $ARGV[1] );
open( my $IN3, "<", $ARGV[2] );

sub LIST_INPUT {
    my $IN = shift;
    my ( @ARRAY_LIST, %INFO, %SCORE );
    while (<$IN>) {
        chomp;
        my (
            $CHR,   $POS,     $DIR,  $BASE1, $BASE2,
            $BASE3, $GENE_ID, undef, undef,  $SCORE
        ) = split /\t/;
        push( @ARRAY_LIST, $CHR . "\t" . $POS . "\t" . $DIR );
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
    return ( \@ARRAY_LIST, \%INFO, \%SCORE );
}

my ( $a_array, $a_info, $a_score ) = LIST_INPUT($IN1);
my ( $b_array, undef,   $b_score ) = LIST_INPUT($IN2);
my ( $c_array, undef,   $c_score ) = LIST_INPUT($IN3);

my (
    @a_only,     @b_only, @c_only,    @a_b_common, @a_c_common,
    @b_c_common, @common, @a_extract, @b_extract,  @c_extract,
    %a_list,     %b_list, %c_list
);

my $extract = $ARGV[3];
my $step    = $ARGV[4];

@a_extract =
  grep( { $a_score->{$_} >= $a_score->{ $a_array->[ $extract - 1 ] } }
    @{$a_array} );
@b_extract =
  grep( { $b_score->{$_} >= $b_score->{ $b_array->[ $extract - 1 ] } }
    @{$b_array} );
@c_extract =
  grep( { $c_score->{$_} >= $c_score->{ $c_array->[ $extract - 1 ] } }
    @{$c_array} );
%a_list = map( { $_ => 1 } @a_extract );
%b_list = map( { $_ => 1 } @b_extract );
%c_list = map( { $_ => 1 } @c_extract );

@common =
  grep ( { $a_list{$_} && $b_list{$_} } @c_extract );

while ( $#common + 1 >= 0.5 * $extract ) {
    $extract += $step;
    @a_extract =
      grep( { $a_score->{$_} >= $a_score->{ $a_array->[ $extract - 1 ] } }
        @{$a_array} );
    @b_extract =
      grep( { $b_score->{$_} >= $b_score->{ $b_array->[ $extract - 1 ] } }
        @{$b_array} );
    @c_extract =
      grep( { $c_score->{$_} >= $c_score->{ $c_array->[ $extract - 1 ] } }
        @{$c_array} );
    %a_list = map( { $_ => 1 } @a_extract );
    %b_list = map( { $_ => 1 } @b_extract );
    %c_list = map( { $_ => 1 } @c_extract );
    @common =
      grep ( { $a_list{$_} && $b_list{$_} } @c_extract );
}

$extract -= $step;

@a_extract =
  grep( { $a_score->{$_} >= $a_score->{ $a_array->[ $extract - 1 ] } }
    @{$a_array} );
@b_extract =
  grep( { $b_score->{$_} >= $b_score->{ $b_array->[ $extract - 1 ] } }
    @{$b_array} );
@c_extract =
  grep( { $c_score->{$_} >= $c_score->{ $c_array->[ $extract - 1 ] } }
    @{$c_array} );
%a_list = map( { $_ => 1 } @a_extract );
%b_list = map( { $_ => 1 } @b_extract );
%c_list = map( { $_ => 1 } @c_extract );

@common =
  grep ( { $a_list{$_} && $b_list{$_} } @c_extract );

foreach (@common) {
    my $score = $a_score->{$_} + $b_score->{$_} + $c_score->{$_};
    print( $a_info->{$_} . "\t" . $score . "\n" );
}

