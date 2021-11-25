#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Statistics::Lite;

sub JUDGE_ZERO {
    my $IN = shift;
    if ( $IN == 0 ) {
        $IN = 1;
    }
    return ($IN);
}

sub RATIO {
    my $TR_END_COUNT = $_[0];
    my $NC_END_COUNT = $_[1];
    if ( $NC_END_COUNT == 0 ) {
        return ($TR_END_COUNT);
    }
    else {
        return ( $TR_END_COUNT / $NC_END_COUNT );
    }
}

sub SCORE {
    my $TR_START_COUNT = $_[0];
    my $TR_END_COUNT   = $_[1];
    my $NC_START_COUNT = $_[2];
    my $NC_END_COUNT   = $_[3];
    my $LEVEL          = $_[4];
    my @SCORE;
    $SCORE[0] = 0;
    $SCORE[1] = 0;

    for my $CURRENT ( 2 .. $#{$TR_END_COUNT} - 2 ) {
        my (
            $T_END,      $T_END_P1, $T_END_P2,   $N_END_M2, $N_END_M1,
            $N_END,      $N_END_P1, $N_END_P2,   $T_START,  $T_START_P1,
            $T_START_P2, $N_START,  $N_START_P1, $N_START_P2
          )
          = (
            JUDGE_ZERO( ${$TR_END_COUNT}[$CURRENT] ),
            JUDGE_ZERO( ${$TR_END_COUNT}[ $CURRENT + 1 ] ),
            JUDGE_ZERO( ${$TR_END_COUNT}[ $CURRENT + 2 ] ),
            JUDGE_ZERO( ${$NC_END_COUNT}[ $CURRENT - 2 ] ),
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
        my $SCORE;
        if ( ( $T_END_P1 < $N_END_P1 * $LEVEL ) or ( $T_END_P1 < 10 ) ) {
            $SCORE = 1;
        }
        else {
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
            my ( $SCORE1, $SCORE2, $SCORE3, $SCORE4, $SCORE5 );
            $SCORE1 = $T_END_P1 / $T_END - $N_END_P1 / $N_END;
            $SCORE1 = 1 if ( $SCORE1 < 1 );
            $SCORE2 = ( $T_START * $T_START_P2 )**0.5 / $T_START_P1;
            $SCORE2 = sqrt($SCORE1) if ( $SCORE2 > sqrt($SCORE1) );
            $SCORE2 = 5             if ( $SCORE2 > 5 );
            $SCORE2 = 0.2           if ( $SCORE2 < 0.2 );
            $SCORE3 =
              ( $N_END_M2 * $N_END_M1 * $N_END_P1 * $N_END_P2 )**0.25 / $N_END;
            $SCORE3 = sqrt($SCORE1) if ( $SCORE3 > sqrt($SCORE1) );
            $SCORE3 = 5             if ( $SCORE3 > 5 );
            $SCORE3 = 0.2           if ( $SCORE3 < 0.2 );
            $SCORE4 = ( $N_START * $N_START_P2 )**0.5 / $N_START_P1;
            $SCORE4 = sqrt($SCORE1) if ( $SCORE4 > sqrt($SCORE1) );
            $SCORE4 = 5             if ( $SCORE4 > 5 );
            $SCORE4 = 0.2           if ( $SCORE4 < 0.2 );
            $SCORE5 = ( $T_END_P1 / $T_END_P2 ) / ( $N_END_P1 / $N_END_P2 );
            $SCORE5 = sqrt($SCORE1) if ( $SCORE5 > sqrt($SCORE1) );
            $SCORE5 = 5             if ( $SCORE5 > 5 );
            $SCORE5 = 0.2           if ( $SCORE5 < 0.2 );
            $SCORE  = $SCORE1 * $SCORE2 * $SCORE3 * $SCORE4 * $SCORE5;
            $SCORE  = 1 if ( $SCORE < 1 );
        }
        push( @SCORE, $SCORE );
    }
    push( @SCORE, 0 );
    push( @SCORE, 0 );
    foreach my $HEAD ( 0 .. 19 ) {
        $SCORE[$HEAD] = 0;
    }
    foreach my $TAIL ( $#SCORE - 19 .. $#SCORE ) {
        $SCORE[$TAIL] = 0;
    }
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
    my @ratio;
    open( my $IN_TR, "<", $ARGV[$sample] );
    while (<$IN_TR>) {
        chomp;
        my ( undef, undef, $start_count, $end_count ) = split /\t/;
        push( @{ $start_count[$sample] }, $start_count );
        push( @{ $end_count[$sample] },   $end_count );
        if ( $end_count != 0 ) {
            my $ratio = RATIO( $end_count,
                ${ $end_count[0] }[ $#{ $end_count[$sample] } ] );
            push( @ratio, $ratio );
        }
    }
    my $mean  = Statistics::Lite::mean(@ratio);
    my $dev   = Statistics::Lite::stddev(@ratio);
    my $level = $mean + $dev / sqrt( $#ratio + 1 );
    @{ $score[$sample] } = SCORE(
        \@{ $start_count[$sample] },
        \@{ $end_count[$sample] },
        \@{ $start_count[0] },
        \@{ $end_count[0] }, $level
    );
    close($IN_TR);
    warn("MEAN=$mean\tSTDDEV=$level\n");
}

foreach my $site ( 0 .. $#site ) {
    print(
"$site[$site]\t$base[$site]\t\t$start_count[0][$site]\t$end_count[0][$site]\t\t"
    );
    foreach my $sample ( 1 .. $#ARGV ) {
        print("$start_count[$sample][$site]\t$end_count[$sample][$site]\t");
    }
    print("\t");
    my $score_sum = 0;
    foreach my $sample ( 1 .. $#ARGV ) {
        print("$score[$sample][$site]\t");
        $score_sum += $score[$sample][$site];
    }
    print("\t$score_sum\t");
    print("\n");
}

__END__
