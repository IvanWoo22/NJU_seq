#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my %NC_count;
my %NC_info;
open( my $IN_NC, "<", $ARGV[0] );
while (<$IN_NC>) {
    chomp;
    my ( $chr, $pos, $dir, undef, undef, undef, undef, undef, $end_count ) =
      split /\t/;
    my $id = $chr . "\t" . $dir . "\t" . $pos;
    $NC_info{$id}  = $_;
    $NC_count{$id} = $end_count;
}
close($IN_NC);

my %TR_count;
my %TR_info;
open( my $IN_TR, "<", $ARGV[1] );
while (<$IN_TR>) {
    chomp;
    my ( $chr, $pos, $dir, undef, undef, undef, undef, undef, $end_count ) =
      split /\t/;
    my $id = $chr . "\t" . $dir . "\t" . $pos;
    $TR_info{$id}  = $_;
    $TR_count{$id} = $end_count;
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
    unless ( exists( $TR_count{ $chr . "\t" . $dir . "\t" . $formal_pos } ) ) {
        $TR_count{ $chr . "\t" . $dir . "\t" . $formal_pos } = 1;
    }
    unless ( exists( $NC_count{ $chr . "\t" . $dir . "\t" . $formal_pos } ) ) {
        $NC_count{ $chr . "\t" . $dir . "\t" . $formal_pos } = 1;
    }
    $T_END = $TR_count{ $chr . "\t" . $dir . "\t" . $formal_pos };
    $N_END = $NC_count{ $chr . "\t" . $dir . "\t" . $formal_pos };
    if ( ( $T_END < 4 ) and ( $N_END < 4 ) ) {
        if ( $T_END > $N_END ) {
            $N_END = $T_END;
        }
        else {
            $T_END = $N_END;
        }
    }
    unless ( exists( $NC_count{$id} ) ) {
        $NC_count{$id} = 0;
    }
    $T_END_P1 = $TR_count{$id};
    $N_END_P1 = $NC_count{$id};
    $score    = $T_END_P1 / $T_END - $N_END_P1 / $N_END;

    if ( $score >= 10 ) {
        print("$TR_info{$id}\t$score\n");
    }
}

__END__
