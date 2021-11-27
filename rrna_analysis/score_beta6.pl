#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use POSIX;

sub JUDGE_ZERO {
    my $IN = shift;
    if ( $IN == 0 ) {
        $IN = 1;
    }
    return ($IN);
}

sub SCORE {
    my $TR_END_COUNT = $_[0];
    my $NC_END_COUNT = $_[1];
    my $TR_TOTAL     = $_[2];
    my $NC_TOTAL     = $_[3];
    my @SCORE;
    my @END_RATIO;
    $SCORE[0]     = 0;
    $SCORE[1]     = 0;
    $END_RATIO[0] = 0;
    $END_RATIO[1] = 0;

    for my $CURRENT ( 2 .. $#{$TR_END_COUNT} - 2 ) {
        my ( $T_END, $T_END_P1 ) = (
            JUDGE_ZERO( ${$TR_END_COUNT}[$CURRENT] ),
            JUDGE_ZERO( ${$TR_END_COUNT}[ $CURRENT + 1 ] )
        );
        my ( $N_END, $N_END_P1 );
        if ( ${$NC_END_COUNT}[$CURRENT] == 0 ) {
            $N_END = 1;
        }
        else {
            $N_END =
              POSIX::ceil( ${$NC_END_COUNT}[$CURRENT] * $TR_TOTAL / $NC_TOTAL );
        }
        if ( ${$NC_END_COUNT}[ $CURRENT + 1 ] == 0 ) {
            $N_END_P1 = 1;
        }
        else {
            $N_END_P1 =
              POSIX::ceil(
                ${$NC_END_COUNT}[ $CURRENT + 1 ] * $TR_TOTAL / $NC_TOTAL );
        }

        if ( ( $T_END < 4 ) and ( $N_END < 4 ) ) {
            if ( $T_END > $N_END ) {
                $N_END = $T_END;
            }
            else {
                $T_END = $N_END;
            }
        }

        my $SCORE;
        my $END_RATIO = 1;
        $SCORE     = $T_END_P1 / $T_END - $N_END_P1 / $N_END;
        $SCORE     = 0 if $SCORE < 0;
        $END_RATIO = $SCORE;
        $END_RATIO = 0 if $T_END_P1 < $N_END_P1 * 5;
        push( @SCORE,     $SCORE );
        push( @END_RATIO, $END_RATIO );
    }
    push( @SCORE,     0 );
    push( @SCORE,     0 );
    push( @END_RATIO, 0 );
    push( @END_RATIO, 0 );
    foreach my $HEAD ( 0 .. 19 ) {
        $SCORE[$HEAD]     = 0;
        $END_RATIO[$HEAD] = 0;
    }
    foreach my $TAIL ( $#SCORE - 19 .. $#SCORE ) {
        $SCORE[$TAIL]     = 0;
        $END_RATIO[$TAIL] = 0;
    }
    return ( \@SCORE, \@END_RATIO );
}

my @site;
my @base;
my @end_count;
my @score;
my @end_ratio;
my @total;
$total[0] = 0;
open( my $IN_NC, "<", $ARGV[0] );

while (<$IN_NC>) {
    chomp;
    my ( $site, $base, undef, $end_count ) = split /\t/;
    push( @site,              $site );
    push( @base,              $base );
    push( @{ $end_count[0] }, $end_count );
    $total[0] += $end_count;
}
close($IN_NC);

foreach my $sample ( 1 .. $#ARGV ) {
    open( my $IN_TR, "<", $ARGV[$sample] );
    $total[$sample] = 0;
    while (<$IN_TR>) {
        chomp;
        my ( undef, undef, undef, $end_count ) = split /\t/;
        push( @{ $end_count[$sample] }, $end_count );
        $total[$sample] += $end_count;
    }
    ( $score[$sample], $end_ratio[$sample] ) = SCORE(
        \@{ $end_count[$sample] },
        \@{ $end_count[0] },
        $total[$sample], $total[0]
    );
    close($IN_TR);
}

foreach my $site ( 0 .. $#site ) {
    print("$site[$site]\t$base[$site]\t\t$end_count[0][$site]\t\t");
    foreach my $sample ( 1 .. $#ARGV ) {
        print("$end_count[$sample][$site]\t");
    }
    my $score_sum     = 0;
    my $end_ratio_sum = 0;
    foreach my $sample ( 1 .. $#ARGV ) {
        print("${$score[$sample]}[$site]\t${$end_ratio[$sample]}[$site]\t");
        $end_ratio_sum += ${ $end_ratio[$sample] }[$site];
        $score_sum     += ${ $score[$sample] }[$site];
    }
    print("$score_sum\t$end_ratio_sum\n");
}

__END__
