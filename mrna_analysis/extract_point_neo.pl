#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use List::Util qw/max min/;

open( my $IN1, "<", $ARGV[0] );
open( my $IN2, "<", $ARGV[1] );
open( my $IN3, "<", $ARGV[2] );

sub LIST_INPUT {
    my $IN = shift;
    my ( @ARRAY_LIST, %INFO, %SCORE, %COUNT );
    my $COUNT = 0;
    while (<$IN>) {
        $COUNT++;
        chomp;
        my (
            $CHR,     $POS,   $DIR, $BASE1, $BASE2, $BASE3,
            $GENE_ID, $TR_P1, $TR,  $NC_P1, $NC,    $SCORE
        ) = split /\t/;
        push( @ARRAY_LIST, $CHR . "\t" . $POS . "\t" . $DIR );
        $INFO{ $CHR . "\t" . $POS . "\t" . $DIR } =
            $CHR . "\t"
          . $POS . "\t"
          . $DIR . "\t"
          . $BASE1 . "\t"
          . $BASE2 . "\t"
          . $BASE3 . "\t"
          . $GENE_ID . "\t"
          . $NC_P1 . "\t"
          . $NC;
        $COUNT{ $CHR . "\t" . $POS . "\t" . $DIR } = $TR_P1 . "\t" . $TR;
        $SCORE{ $CHR . "\t" . $POS . "\t" . $DIR } = $SCORE;
    }
    return ( \@ARRAY_LIST, \%INFO, \%SCORE, \%COUNT, $COUNT );
}

my ( $a_array, $info, $a_score, $a_end, $a_count ) = LIST_INPUT($IN1);
my ( $b_array, undef, $b_score, $b_end, $b_count ) = LIST_INPUT($IN2);
my ( $c_array, undef, $c_score, $c_end, $c_count ) = LIST_INPUT($IN3);
my $count_max = List::Util::max( $a_count, $b_count, $c_count );

my ( @common, @a_extract, @b_extract, @c_extract, %a_list, %b_list, %c_list );

my $extract = $ARGV[3];
my $step    = $ARGV[4];

my ( $a_extract, $b_extract, $c_extract );
$a_extract = List::Util::min( $extract - 1, $a_count - 1 );
$b_extract = List::Util::min( $extract - 1, $b_count - 1 );
$c_extract = List::Util::min( $extract - 1, $c_count - 1 );
@a_extract =
  grep( { $a_score->{$_} >= $a_score->{ $a_array->[$a_extract] } }
    @{$a_array} );
@b_extract =
  grep( { $b_score->{$_} >= $b_score->{ $b_array->[$b_extract] } }
    @{$b_array} );
@c_extract =
  grep( { $c_score->{$_} >= $c_score->{ $c_array->[$c_extract] } }
    @{$c_array} );
%a_list = map( { $_ => 1 } @a_extract );
%b_list = map( { $_ => 1 } @b_extract );
%c_list = map( { $_ => 1 } @c_extract );

@common =
  grep ( { $a_list{$_} && $b_list{$_} } @c_extract );

while ( ( $#common <= $ARGV[3] ) and ( $extract < $count_max ) ) {
    $extract += $step;
    $a_extract = List::Util::min( $extract - 1, $a_count - 1 );
    $b_extract = List::Util::min( $extract - 1, $b_count - 1 );
    $c_extract = List::Util::min( $extract - 1, $c_count - 1 );
    @a_extract =
      grep( { $a_score->{$_} >= $a_score->{ $a_array->[$a_extract] } }
        @{$a_array} );
    @b_extract =
      grep( { $b_score->{$_} >= $b_score->{ $b_array->[$b_extract] } }
        @{$b_array} );
    @c_extract =
      grep( { $c_score->{$_} >= $c_score->{ $c_array->[$c_extract] } }
        @{$c_array} );
    %a_list = map( { $_ => 1 } @a_extract );
    %b_list = map( { $_ => 1 } @b_extract );
    %c_list = map( { $_ => 1 } @c_extract );
    @common =
      grep ( { $a_list{$_} && $b_list{$_} } @c_extract );
}

@a_extract =
  grep( { $a_score->{$_} >= $a_score->{ $a_array->[$a_extract] } }
    @{$a_array} );
@b_extract =
  grep( { $b_score->{$_} >= $b_score->{ $b_array->[$b_extract] } }
    @{$b_array} );
@c_extract =
  grep( { $c_score->{$_} >= $c_score->{ $c_array->[$c_extract] } }
    @{$c_array} );
%a_list = map( { $_ => 1 } @a_extract );
%b_list = map( { $_ => 1 } @b_extract );
%c_list = map( { $_ => 1 } @c_extract );

@common =
  grep ( { $a_list{$_} && $b_list{$_} } @c_extract );

foreach (@common) {
    my $score = $a_score->{$_} + $b_score->{$_} + $c_score->{$_};
    print(  $info->{$_} . "\t"
          . $a_end->{$_} . "\t"
          . $b_end->{$_} . "\t"
          . $c_end->{$_} . "\t"
          . $score
          . "\n" );
}

if ( $extract < $count_max ) {
    warn("$extract\n");
}
else {
    warn("$extract\texceeded!\n");
}

__END__
