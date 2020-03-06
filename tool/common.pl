#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open( my $in1, "<", $ARGV[0] );
open( my $in2, "<", $ARGV[1] );
open( my $in3, "<", $ARGV[2] );

my %count;
my %info;
while (<$in1>) {
    chomp;
    my ( $chr, $pos, $dir, undef, $base2, undef, $gene ) = split /\t/;
    if ( exists( $count{ $chr . "\t" . $dir . "\t" . $pos } ) ) {
        $count{ $chr . "\t" . $dir . "\t" . $pos }++;
    }
    else {
        $count{ $chr . "\t" . $dir . "\t" . $pos } = 1;
        $info{ $chr . "\t" . $dir . "\t" . $pos } =
          $chr . "\t" . $dir . "\t" . $pos . "\t" . $base2 . "\t" . $gene;
    }
}
while (<$in2>) {
    chomp;
    my ( $chr, $pos, $dir, undef, $base2, undef, $gene ) = split /\t/;
    if ( exists( $count{ $chr . "\t" . $dir . "\t" . $pos } ) ) {
        $count{ $chr . "\t" . $dir . "\t" . $pos }++;
    }
    else {
        $count{ $chr . "\t" . $dir . "\t" . $pos } = 1;
        $info{ $chr . "\t" . $dir . "\t" . $pos } =
          $chr . "\t" . $dir . "\t" . $pos . "\t" . $base2 . "\t" . $gene;
    }
}
while (<$in3>) {
    chomp;
    my ( $chr, $pos, $dir, undef, $base2, undef, $gene ) = split /\t/;
    if ( exists( $count{ $chr . "\t" . $dir . "\t" . $pos } ) ) {
        $count{ $chr . "\t" . $dir . "\t" . $pos }++;
    }
    else {
        $count{ $chr . "\t" . $dir . "\t" . $pos } = 1;
        $info{ $chr . "\t" . $dir . "\t" . $pos } =
          $chr . "\t" . $dir . "\t" . $pos . "\t" . $base2 . "\t" . $gene;
    }
}

foreach ( keys(%count) ) {
    if ( $count{$_} == 3 ) {
        print("$info{$_}\n");
    }
}

__END__
