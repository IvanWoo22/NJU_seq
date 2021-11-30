#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use POSIX;

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
        my ( $N_END, $N_END_P1, $T_END, $T_END_P1 );
        if ( ${$TR_END_COUNT}[$CURRENT] == 0 ) {
            $T_END = 1;
        }
        else {
            $T_END =
              POSIX::ceil( ${$TR_END_COUNT}[$CURRENT] * $NC_TOTAL / $TR_TOTAL );
        }
        $T_END_P1 = POSIX::ceil(
            ${$TR_END_COUNT}[ $CURRENT + 1 ] * $NC_TOTAL / $TR_TOTAL );
        if ( ${$NC_END_COUNT}[$CURRENT] == 0 ) {
            $N_END = 1;
        }
        else {
            $N_END = ${$NC_END_COUNT}[$CURRENT];
        }
        $N_END_P1 = ${$NC_END_COUNT}[ $CURRENT + 1 ];
        if ( ( $T_END < 4 ) and ( $N_END < 4 ) ) {
            if ( $T_END > $N_END ) {
                $N_END = $T_END;
            }
            else {
                $T_END = $N_END;
            }
        }
        my $SCORE;
        $SCORE = $T_END_P1 / $T_END - $N_END_P1 / $N_END;
        push( @SCORE,     $SCORE );
        push( @END_RATIO, $T_END_P1 );
    }
    push( @SCORE,     0 );
    push( @SCORE,     0 );
    push( @END_RATIO, 0 );
    push( @END_RATIO, 0 );
    foreach my $TERMIN ( 0 .. 19 ) {
        $SCORE[$TERMIN]                 = 0;
        $END_RATIO[$TERMIN]             = 0;
        $SCORE[ $#SCORE - $TERMIN ]     = 0;
        $END_RATIO[ $#SCORE - $TERMIN ] = 0;
    }
    return ( \@SCORE, \@END_RATIO );
}

my ( @site, @base, @end_count, @score, @end_count_cor, @total );
my ( @socre_cor1, @socre_cor2, @socre_cor3 );
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
    ( $score[$sample], $end_count_cor[$sample] ) = SCORE(
        \@{ $end_count[$sample] },
        \@{ $end_count[0] },
        $total[$sample], $total[0]
    );
    close($IN_TR);
}

foreach my $site ( 0 .. $#site - 1 ) {
    print("$site[$site]\t$base[$site]\t$end_count[0][$site]\t");
    my $fto = 0;    # fewer than one
    my $soa = 0;    # summary of all
    foreach my $sample ( 1 .. $#ARGV ) {
        $socre_cor1[$sample][$site] = $score[$sample][$site];
        $socre_cor2[$sample][$site] = $score[$sample][$site];
        $socre_cor3[$sample][$site] = $score[$sample][$site];
        $soa += $end_count_cor[$sample][$site];
        if ( $score[$sample][$site] < 20 ) {
            $fto++;
        }
    }
    if ( $fto > 0 ) {
        foreach my $sample ( 1 .. $#ARGV ) {
            $socre_cor1[$sample][$site] = 0;
            $socre_cor3[$sample][$site] = 0;
        }
    }
    if ( $soa < 15 * $end_count[0][ $site + 1 ] ) {
        foreach my $sample ( 1 .. $#ARGV ) {
            $socre_cor2[$sample][$site] = 0;
            $socre_cor3[$sample][$site] = 0;
        }
    }
    my $score_sum      = 0;
    my $score_cor1_sum = 0;
    my $score_cor2_sum = 0;
    my $score_cor3_sum = 0;
    foreach my $sample ( 1 .. $#ARGV ) {
        $score_cor1_sum += $socre_cor1[$sample][$site];
        $score_cor2_sum += $socre_cor2[$sample][$site];
        $score_cor3_sum += $socre_cor3[$sample][$site];
        $score_sum      += $score[$sample][$site];
        print(
"$end_count[$sample][$site]\t$end_count_cor[$sample][$site]\t$score[$sample][$site]\t$socre_cor1[$sample][$site]\t$socre_cor2[$sample][$site]\t$socre_cor3[$sample][$site]\t"
        );
    }
    print("$score_sum\t$score_cor1_sum\t$score_cor2_sum\t$score_cor3_sum\n");
}

__END__
