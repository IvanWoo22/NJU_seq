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
    my %SCORE;
    my %END_COR;
    for my $CURRENT ( keys( %{$TR_END_COUNT} ) ) {
        my ( $CHR, $DIR, $POS ) = split( /\t/, $CURRENT );
        my ($FORMAL_POS);
        if ( $DIR eq "+" ) {
            $FORMAL_POS = $POS - 1;
        }
        else {
            $FORMAL_POS = $POS + 1;
        }
        my ( $N_END, $N_END_P1, $T_END, $T_END_P1 );
        if (
            not exists(
                ${$TR_END_COUNT}{ $CHR . "\t" . $DIR . "\t" . $FORMAL_POS }
            )
          )
        {
            $T_END = 1;
        }
        else {
            $T_END = POSIX::ceil(
                ${$TR_END_COUNT}{ $CHR . "\t" . $DIR . "\t" . $FORMAL_POS } *
                  $NC_TOTAL /
                  $TR_TOTAL );
        }
        if (
            not exists(
                ${$NC_END_COUNT}{ $CHR . "\t" . $DIR . "\t" . $FORMAL_POS }
            )
          )
        {
            $N_END = 1;
        }
        else {
            $N_END =
              ${$NC_END_COUNT}{ $CHR . "\t" . $DIR . "\t" . $FORMAL_POS };
        }
        $T_END_P1 =
          POSIX::ceil( ${$TR_END_COUNT}{$CURRENT} * $NC_TOTAL / $TR_TOTAL );

        if ( not exists( ${$NC_END_COUNT}{$CURRENT} ) ) {
            $N_END_P1 = 1;
        }
        else {
            $N_END_P1 = ${$NC_END_COUNT}{$CURRENT};
        }

        my $SCORE = $T_END_P1 / $T_END - $N_END_P1 / $N_END;
        $SCORE{$CURRENT}   = $SCORE;
        $END_COR{$CURRENT} = $T_END_P1;
    }
    return ( \%SCORE, \%END_COR );
}

my ( @end_count, @score, @end_count_cor, @total, %info, %all_site_id );
open( my $IN_NC, "<", $ARGV[0] );
while (<$IN_NC>) {
    chomp;
    my @tmp = split /\t/;
    my $id  = $tmp[0] . "\t" . $tmp[2] . "\t" . $tmp[1];
    ${ $end_count[0] }{$id} = $tmp[8];
}
close($IN_NC);

foreach my $sample ( 1 .. $#ARGV ) {
    $total[0] = 0;
    open( my $IN_TR, "<", $ARGV[$sample] );
    while (<$IN_TR>) {
        chomp;
        my @tmp = split /\t/;
        my $id  = $tmp[0] . "\t" . $tmp[2] . "\t" . $tmp[1];
        if ( exists( $all_site_id{$id} ) ) {
            $all_site_id{$id}++;
        }
        else {
            $info{$id}        = join( "\t", @tmp[ 3 .. 6 ] );
            $all_site_id{$id} = 1;
        }
        ${ $end_count[$sample] }{$id} = $tmp[8];
        $total[$sample] += $tmp[8];
        if ( exists( ${ $end_count[0] }{$id} ) ) {
            $total[0] += ${ $end_count[0] }{$id};
        }
    }
    close($IN_TR);
    warn("$total[$sample], $total[0]\n");
    ( $score[$sample], $end_count_cor[$sample] ) = SCORE(
        \%{ $end_count[$sample] },
        \%{ $end_count[0] },
        $total[$sample], $total[0]
    );
}

foreach my $id ( keys(%all_site_id) ) {
    if ( $all_site_id{$id} == $#ARGV ) {
        my ( $chr, $dir, $pos ) = split( /\t/, $id );
        my $NC_END_COUNT;
        if ( exists( ${ $end_count[0] }{$id} ) ) {
            $NC_END_COUNT = ${ $end_count[0] }{$id};
        }
        else {
            $NC_END_COUNT = 0;
        }
        my $soas = 0;
        my $soac = 0;
        foreach my $sample ( 1 .. $#ARGV ) {
            $soac += ${ $end_count_cor[$sample] }{$id};
            $soas += ${ $score[$sample] }{$id};
        }
        if ( ( $soas >= 90 ) and ( $soac >= 9 * $NC_END_COUNT ) ) {
            print("$chr\t$dir\t$pos\t$info{$id}\t");
            print("$NC_END_COUNT\t");
            foreach my $sample ( 1 .. $#ARGV ) {
                print(
"${$end_count[$sample]}{$id}\t${$end_count_cor[$sample]}{$id}\t${$score[$sample]}{$id}\t"
                );
            }
            print("$soas\n");
        }
    }
}

__END__
