#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use POSIX;
use Statistics::Lite;

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

my @start_count;
my @end_count;
my %info;
my @score;

open( my $IN_NC, "<", $ARGV[0] );
while (<$IN_NC>) {
    chomp;
    my @tmp = split /\t/;
    my $id  = $tmp[0] . "\t" . $tmp[2] . "\t" . $tmp[1];
    ${ $start_count[0] }{$id} = $tmp[7];
    ${ $end_count[0] }{$id}   = $tmp[8];
}
close($IN_NC);

foreach my $sample ( 1 .. $#ARGV ) {
    my @ratio;
    open( my $IN_TR, "<", $ARGV[$sample] );
    while (<$IN_TR>) {
        chomp;
        my @tmp = split /\t/;
        my $id  = $tmp[0] . "\t" . $tmp[2] . "\t" . $tmp[1];
        $info{$id} = join( "\t", @tmp[ 0 .. 6 ] );
        ${ $start_count[$sample] }{$id} = $tmp[7];
        ${ $end_count[$sample] }{$id}   = $tmp[8];
        if ( $tmp[8] > 0 ) {
            my $ratio;
            if ( exists( ${ $end_count[0] }{$id} ) ) {
                $ratio = RATIO( $tmp[8], ${ $end_count[0] }{$id} );
            }
            else {
                $ratio = RATIO( $tmp[8], 0 );
            }
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

# foreach my $id ( keys(%TR_count) ) {
#     my ( $chr, $dir, $pos ) = split( /\t/, $id );
#     my ( $formal_pos, $score );
#     if ( $dir eq "+" ) {
#         $formal_pos = $pos - 1;
#     }
#     else {
#         $formal_pos = $pos + 1;
#     }
#
#     my ( $T_END, $T_END_P1, $N_END, $N_END_P1 );
#     if ( not exists( $TR_count{ $chr . "\t" . $dir . "\t" . $formal_pos } ) ) {
#         $T_END = 1;
#     }
#     else {
#         $T_END = $TR_count{ $chr . "\t" . $dir . "\t" . $formal_pos };
#     }
#     if ( not exists( $NC_count{ $chr . "\t" . $dir . "\t" . $formal_pos } ) ) {
#         $N_END = 1;
#     }
#     else {
#         $N_END =
#           POSIX::ceil( $NC_count{ $chr . "\t" . $dir . "\t" . $formal_pos } *
#               $TR_total /
#               $NC_total );
#     }
#
#     if ( ( $T_END < 4 ) and ( $N_END < 4 ) ) {
#         if ( $T_END > $N_END ) {
#             $N_END = $T_END;
#         }
#         else {
#             $T_END = $N_END;
#         }
#     }
#     if ( not exists( $NC_count{$id} ) ) {
#         $N_END_P1 = 0;
#     }
#     else {
#         $N_END_P1 = POSIX::ceil( $NC_count{$id} * $TR_total / $NC_total );
#     }
#     $T_END_P1 = $TR_count{$id};
#     $score    = $T_END_P1 / $T_END - $N_END_P1 / $N_END;
#
#     if ( ( $score >= 10 ) and ( $T_END_P1 > $N_END_P1 * 5 ) ) {
#         print("$info{$id}\t$T_END_P1\t$T_END\t$N_END_P1\t$N_END\t$score\n");
#     }
# }

__END__
