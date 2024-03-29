#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open( my $IN1, "<", $ARGV[0] );
open( my $IN2, "<", $ARGV[1] );
open( my $IN3, "<", $ARGV[2] );

my ( %count, %info, %score );
while (<$IN1>) {
    chomp;
    my (
        $chr,   $pos,     $dir,  $base1, $base2,
        $base3, $gene_id, undef, undef,  $score
    ) = split /\t/;
    if ( exists( $count{ $chr . "\t" . $dir . "\t" . $pos } ) ) {
        $count{ $chr . "\t" . $dir . "\t" . $pos }++;
        $score{ $chr . "\t" . $dir . "\t" . $pos } += $score;
    }
    else {
        $count{ $chr . "\t" . $dir . "\t" . $pos } = 1;
        $info{ $chr . "\t" . $dir . "\t" . $pos } =
            $chr . "\t"
          . $dir . "\t"
          . $pos . "\t"
          . $base1 . "\t"
          . $base2 . "\t"
          . $base3 . "\t"
          . $gene_id;
        $score{ $chr . "\t" . $dir . "\t" . $pos } = $score;
    }
}
while (<$IN2>) {
    chomp;
    my (
        $chr,   $pos,     $dir,  $base1, $base2,
        $base3, $gene_id, undef, undef,  $score
    ) = split /\t/;
    if ( exists( $count{ $chr . "\t" . $dir . "\t" . $pos } ) ) {
        $count{ $chr . "\t" . $dir . "\t" . $pos }++;
        $score{ $chr . "\t" . $dir . "\t" . $pos } += $score;
    }
    else {
        $count{ $chr . "\t" . $dir . "\t" . $pos } = 1;
        $info{ $chr . "\t" . $dir . "\t" . $pos } =
            $chr . "\t"
          . $dir . "\t"
          . $pos . "\t"
          . $base1 . "\t"
          . $base2 . "\t"
          . $base3 . "\t"
          . $gene_id;
        $score{ $chr . "\t" . $dir . "\t" . $pos } = $score;
    }
}
while (<$IN3>) {
    chomp;
    my (
        $chr,   $pos,     $dir,  $base1, $base2,
        $base3, $gene_id, undef, undef,  $score
    ) = split /\t/;
    if ( exists( $count{ $chr . "\t" . $dir . "\t" . $pos } ) ) {
        $count{ $chr . "\t" . $dir . "\t" . $pos }++;
        $score{ $chr . "\t" . $dir . "\t" . $pos } += $score;
    }
    else {
        $count{ $chr . "\t" . $dir . "\t" . $pos } = 1;
        $info{ $chr . "\t" . $dir . "\t" . $pos } =
            $chr . "\t"
          . $dir . "\t"
          . $pos . "\t"
          . $base1 . "\t"
          . $base2 . "\t"
          . $base3 . "\t"
          . $gene_id;
        $score{ $chr . "\t" . $dir . "\t" . $pos } = $score;
    }
}

foreach ( keys(%count) ) {
    if ( $count{$_} == 3 ) {
        print("$info{$_}\t$score{$_}\n");
    }
}

__END__
