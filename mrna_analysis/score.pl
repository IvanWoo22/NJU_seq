#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my %NC_count;
my %NC_info;
open( my $IN_NC, "<", $ARGV[0] );
while (<$IN_NC>) {
    chomp;
    my ( $chr, $pos, $dir, undef, undef, undef, undef, $end_count ) =
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
    my ( $chr, $pos, $dir, undef, undef, undef, undef, $end_count ) =
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
    if ( exists( $TR_count{ $chr . "\t" . $dir . "\t" . $formal_pos } ) ) {
        $score =
          $TR_count{$id} / $TR_count{ $chr . "\t" . $dir . "\t" . $formal_pos };
    }
    else {
        $score = $TR_count{$id} * 2;
    }
    if ( exists( $NC_count{$id} ) ) {
        if ( exists( $NC_count{ $chr . "\t" . $dir . "\t" . $formal_pos } ) ) {
            $score = $score - $NC_count{$id} /
              $NC_count{ $chr . "\t" . $dir . "\t" . $formal_pos };
        }
        else {
            $score = $score - $NC_count{$id} * 2;
        }
    }
    if ( $score > 20 ) {
        print("$TR_info{$id}\t$score\n");
    }
}

__END__
