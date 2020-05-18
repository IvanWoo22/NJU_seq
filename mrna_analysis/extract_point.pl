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
my ( $b_array, $b_info, $b_score ) = LIST_INPUT($IN2);
my ( $c_array, $c_info, $c_score ) = LIST_INPUT($IN3);

my ( @a_only, @b_only, @c_only, @a_b_common, @a_c_common, @b_c_common,
    @common );

my $extract = $ARGV[3];
my $step    = $ARGV[4];

while (( $#common >= $#a_only )
    || ( $#common >= $#b_only )
    || ( $#common >= $#b_only ) )
{
    $extract += $step;
    my @a_extract =
      grep( { $a_score->{$_} >= $a_score->{ $a_array->[ $extract - 1 ] } }
        @{$a_array} );
    my @b_extract =
      grep( { $b_score->{$_} >= $b_score->{ $b_array->[ $extract - 1 ] } }
        @{$b_array} );
    my @c_extract =
      grep( { $c_score->{$_} >= $c_score->{ $c_array->[ $extract - 1 ] } }
        @{$c_array} );
    my %a_list = map( { $_ => 1 } @a_extract );
    my %b_list = map( { $_ => 1 } @b_extract );
    my %c_list = map( { $_ => 1 } @c_extract );
    @a_only = grep ( { !$b_list{$_} && !$c_list{$_} } @a_extract );
    @b_only = grep ( { !$a_list{$_} && !$c_list{$_} } @b_extract );
    @c_only = grep ( { !$a_list{$_} && !$b_list{$_} } @c_extract );
    @common =
      grep ( { $a_list{$_} && $b_list{$_} } @c_extract );
}

$extract -= $step;

my @a_extract =
  grep( { $a_score->{$_} >= $a_score->{ $a_array->[ $extract - 1 ] } }
    @{$a_array} );
my @b_extract =
  grep( { $b_score->{$_} >= $b_score->{ $b_array->[ $extract - 1 ] } }
    @{$b_array} );
my @c_extract =
  grep( { $c_score->{$_} >= $c_score->{ $c_array->[ $extract - 1 ] } }
    @{$c_array} );
my %a_list = map( { $_ => 1 } @a_extract );
my %b_list = map( { $_ => 1 } @b_extract );
my %c_list = map( { $_ => 1 } @c_extract );

@a_b_common =
  grep ( { $b_list{$_} && !$c_list{$_} } @a_extract );
@a_c_common =
  grep ( { $a_list{$_} && !$b_list{$_} } @c_extract );
@b_c_common =
  grep ( { !$a_list{$_} && $c_list{$_} } @b_extract );
@common =
  grep ( { $a_list{$_} && $b_list{$_} } @c_extract );

foreach (@common) {
    my $score = $a_score->{$_} + $b_score->{$_} + $c_score->{$_};
    print( $a_info->{$_} . "\t" . $score . "\n" );
}
foreach (@a_b_common) {
    my $score = $a_score->{$_} + $b_score->{$_};
    print( $a_info->{$_} . "\t" . $score . "\n" );
}
foreach (@a_c_common) {
    my $score = $a_score->{$_} + $c_score->{$_};
    print( $c_info->{$_} . "\t" . $score . "\n" );
}
foreach (@b_c_common) {
    my $score = $b_score->{$_} + $c_score->{$_};
    print( $b_info->{$_} . "\t" . $score . "\n" );
}

