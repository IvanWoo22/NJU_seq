#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

sub JUDGE_ZERO {
    my $IN = shift;
    if ( $IN == 0 ) {
        $IN = 1;
    }
    return ($IN);
}

sub SCORE {
    my $TR_START_COUNT = $_[0];
    my $TR_END_COUNT   = $_[1];
    my $NC_START_COUNT = $_[2];
    my $NC_END_COUNT   = $_[3];

    my @SCORE;
    $SCORE[0] = 0;
    for my $CURRENT ( 1 .. $#{$TR_END_COUNT} - 2 ) {
        my (
            $T_END,    $T_END_P1,   $T_END_P2, $N_END_M1,   $N_END,
            $N_END_P1, $N_END_P2,   $T_START,  $T_START_P1, $T_START_P2,
            $N_START,  $N_START_P1, $N_START_P2
          )
          = (
            JUDGE_ZERO( ${$TR_END_COUNT}[$CURRENT] ),
            JUDGE_ZERO( ${$TR_END_COUNT}[ $CURRENT + 1 ] ),
            JUDGE_ZERO( ${$TR_END_COUNT}[ $CURRENT + 2 ] ),
            JUDGE_ZERO( ${$NC_END_COUNT}[ $CURRENT - 1 ] ),
            JUDGE_ZERO( ${$NC_END_COUNT}[$CURRENT] ),
            JUDGE_ZERO( ${$NC_END_COUNT}[ $CURRENT + 1 ] ),
            JUDGE_ZERO( ${$NC_END_COUNT}[ $CURRENT + 2 ] ),
            JUDGE_ZERO( ${$TR_START_COUNT}[$CURRENT] ),
            JUDGE_ZERO( ${$TR_START_COUNT}[ $CURRENT + 1 ] ),
            JUDGE_ZERO( ${$TR_START_COUNT}[ $CURRENT + 2 ] ),
            JUDGE_ZERO( ${$NC_START_COUNT}[$CURRENT] ),
            JUDGE_ZERO( ${$NC_START_COUNT}[ $CURRENT + 1 ] ),
            JUDGE_ZERO( ${$NC_START_COUNT}[ $CURRENT + 2 ] )
          );

        if ( ( $T_END < 4 ) and ( $N_END < 4 ) ) {
            if ( $T_END > $N_END ) {
                $N_END = $T_END;
            }
            else {
                $T_END = $N_END;
            }
        }

        if ( ( $T_END_P2 < 4 ) and ( $N_END_P2 < 4 ) ) {
            if ( $T_END_P2 > $N_END_P2 ) {
                $N_END_P2 = $T_END_P2;
            }
            else {
                $T_END_P2 = $N_END_P2;
            }
        }

        if ( ( $T_START_P1 < 4 ) and ( $N_START_P1 < 4 ) ) {
            if ( $T_START_P1 > $N_START_P1 ) {
                $N_START_P1 = $T_START_P1;
            }
            else {
                $T_START_P1 = $N_START_P1;
            }
        }
        my ( $SCORE1, $SCORE2, $SCORE3, $SCORE4, $SCORE5, $SCORE );
        $SCORE1 = $T_END_P1 / $T_END;
        $SCORE2 =
          ( $T_START + $T_START_P1 + $T_START_P2 ) / ( $T_START_P1 * 3 );
        $SCORE3 = ( $N_END_M1 + $N_END + $N_END_P1 ) / ( $N_END * 3 );
        $SCORE4 =
          ( $N_START + $N_START_P1 + $N_START_P2 ) / ( $N_START_P1 * 3 );
        $SCORE5 = ( $T_END_P1 / $T_END_P2 ) / ( $N_END_P1 / $N_END_P2 );
        $SCORE  = $SCORE1 * $SCORE2 * $SCORE3 * $SCORE4 * $SCORE5;
        push( @SCORE, $SCORE );
    }
    push( @SCORE, 0 );
    push( @SCORE, 0 );
    return (@SCORE);
}

my @site;
my @base;
my @start_count;
my @end_count;
my @score;
open( my $IN_NC, "<", $ARGV[0] );
while (<$IN_NC>) {
    chomp;
    my ( $site, $base, $start_count, $end_count ) = split /\t/;
    push( @site,                $site );
    push( @base,                $base );
    push( @{ $start_count[0] }, $start_count );
    push( @{ $end_count[0] },   $end_count );
}
close($IN_NC);

foreach my $sample ( 1 .. $#ARGV ) {
    open( my $IN_TR, "<", $ARGV[$sample] );
    while (<$IN_TR>) {
        chomp;
        my ( undef, undef, $start_count, $end_count ) = split /\t/;
        push( @{ $start_count[$sample] }, $start_count );
        push( @{ $end_count[$sample] },   $end_count );
    }
    @{ $score[$sample] } = SCORE(
        \@{ $start_count[$sample] },
        \@{ $end_count[$sample] },
        \@{ $start_count[0] },
        \@{ $end_count[0] }
    );
    close($IN_TR);
}

foreach my $site ( 0 .. $#site ) {
    print(
"$site[$site]\t$base[$site]\t$start_count[0][$site]\t$end_count[0][$site]\t\t"
    );
    foreach my $sample ( 1 .. $#ARGV ) {
        print("$start_count[$sample][$site]\t$end_count[$sample][$site]\t");
    }
    print("\t");
    foreach my $sample ( 1 .. $#ARGV ) {
        print("$score[$sample][$site]\t");
    }
    print("\n");
}

__END__
