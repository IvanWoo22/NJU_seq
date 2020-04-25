#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my $norm = $ARGV[0];

my %NC_end_count;
my %NC_info;
my $NC_sum   = 0;
my $NC_count = 0;
open( my $IN_NC, "<", $ARGV[1] );
while (<$IN_NC>) {
    $NC_count++;
    chomp;
    my ( $chr, $pos, $dir, undef, undef, undef, undef, undef, $end_count ) =
      split /\t/;
    my $id = $chr . "\t" . $dir . "\t" . $pos;
    $NC_info{$id}      = $_;
    $NC_end_count{$id} = $end_count;
    $NC_sum += $end_count;
}
close($IN_NC);

my %TR_end_count;
my %TR_info;
my $TR_sum   = 0;
my $TR_count = 0;
open( my $IN_TR, "<", $ARGV[2] );
while (<$IN_TR>) {
    $TR_count++;
    chomp;
    my ( $chr, $pos, $dir, undef, undef, undef, undef, undef, $end_count ) =
      split /\t/;
    my $id = $chr . "\t" . $dir . "\t" . $pos;
    $TR_info{$id}      = $_;
    $TR_end_count{$id} = $end_count;
    $TR_sum += $end_count;
}
close($IN_TR);

foreach my $id ( keys(%TR_end_count) ) {
    my ( $chr, $dir, $pos ) = split( /\t/, $id );
    my ( $formal_pos, $score );
    if ( $dir eq "+" ) {
        $formal_pos = $pos - 1;
    }
    else {
        $formal_pos = $pos + 1;
    }
    if ( exists( $TR_end_count{ $chr . "\t" . $dir . "\t" . $formal_pos } ) ) {
        $score =
          $TR_end_count{$id} /
          $TR_end_count{ $chr . "\t" . $dir . "\t" . $formal_pos };
    }
    else {
        $score = $TR_end_count{$id} * 2 * $norm * $TR_count / $TR_sum;
    }
    if ( exists( $NC_end_count{$id} ) ) {
        if (
            exists( $NC_end_count{ $chr . "\t" . $dir . "\t" . $formal_pos } ) )
        {
            $score = $score - $NC_end_count{$id} /
              $NC_end_count{ $chr . "\t" . $dir . "\t" . $formal_pos };
        }
        else {
            $score =
              $score - $NC_end_count{$id} * 2 * $norm * $NC_count / $NC_sum;
        }
    }
    if ( $score >= 10 ) {
        print("$TR_info{$id}\t$score\n");
    }
}

__END__
