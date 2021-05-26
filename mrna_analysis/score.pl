#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use POSIX;

my %NC_count;
my $NC_total = 0;
open( my $IN_NC, "<", $ARGV[0] );
while (<$IN_NC>) {
    chomp;
    my @tmp = split /\t/;
    my $id  = $tmp[0] . "\t" . $tmp[2] . "\t" . $tmp[1];
    $NC_count{$id} = $tmp[8];
    $NC_total += $tmp[8];
}
close($IN_NC);

my ( %TR_count, %info );
my $TR_total = 0;
open( my $IN_TR, "<", $ARGV[1] );
while (<$IN_TR>) {
    chomp;
    my @tmp = split /\t/;
    my $id  = $tmp[0] . "\t" . $tmp[2] . "\t" . $tmp[1];
    $info{$id}     = join( "\t", @tmp[ 0 .. 6 ] );
    $TR_count{$id} = $tmp[8];
    $TR_total += $tmp[8];
}
close($IN_TR);

foreach my $id ( keys(%TR_count) ) {
    my ( $chr, $dir, $pos ) = split( /\t/, $id );
    my ( $formal_pos, $score );
    if ( $dir eq "+" ) {
        $formal_pos = $pos - 1;
    }
    else {
        $formal_pos = $pos + 1;
    }

    my ( $T_END, $T_END_P1, $N_END, $N_END_P1 );
    if ( not exists( $TR_count{ $chr . "\t" . $dir . "\t" . $formal_pos } ) ) {
        $T_END = 1;
    }
    else {
        $T_END = $TR_count{ $chr . "\t" . $dir . "\t" . $formal_pos };
    }
    if ( not exists( $NC_count{ $chr . "\t" . $dir . "\t" . $formal_pos } ) ) {
        $N_END = 1;
    }
    else {
        $N_END =
          POSIX::ceil( $NC_count{ $chr . "\t" . $dir . "\t" . $formal_pos } *
              $TR_total /
              $NC_total );
    }

    if ( ( $T_END < 4 ) and ( $N_END < 4 ) ) {
        if ( $T_END > $N_END ) {
            $N_END = $T_END;
        }
        else {
            $T_END = $N_END;
        }
    }
    if ( not exists( $NC_count{$id} ) ) {
        $N_END_P1 = 0;
    }
    else {
        $N_END_P1 = POSIX::ceil( $NC_count{$id} * $TR_total / $NC_total );
    }
    $T_END_P1 = $TR_count{$id};
    $score    = $T_END_P1 / $T_END - $N_END_P1 / $N_END;

    if ( ( $score >= 10 ) and ( $T_END_P1 > $N_END_P1 * 5 ) ) {
        print("$info{$id}\t$T_END_P1\t$T_END\t$N_END_P1\t$N_END\t$score\n");
    }
}

__END__
