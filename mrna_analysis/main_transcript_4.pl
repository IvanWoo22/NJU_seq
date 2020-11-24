#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use POSIX;

my ( @line, @dis, @type, @weight );
my ( $fu, $cds, $tu ) = ( 0, 0, 0 );
while (<STDIN>) {
    chomp;
    push( @line, $_ );
    my @tmp = split "\t";
    $fu  += $tmp[5];
    $cds += $tmp[6];
    $tu  += $tmp[7];
    push( @dis,    $tmp[4] );
    push( @type,   $tmp[3] );
    push( @weight, $tmp[9] );
}

my ( @dis_fu, @dis_cds, @dis_tu );
@dis_fu  = (0) x 6;
@dis_cds = (0) x ( POSIX::ceil( 5 * $cds / $fu ) + 1 );
@dis_tu  = (0) x ( POSIX::ceil( 5 * $tu / $fu ) + 1 );
foreach ( 0 .. $#line ) {
    if ( $type[$_] eq "five_utr" ) {
        $dis_fu[ POSIX::ceil( $dis[$_] * 5 ) ] += $weight[$_];
    }
    elsif ( $type[$_] eq "cds" ) {
        $dis_cds[ POSIX::ceil( POSIX::ceil( 5 * $cds / $fu ) * $dis[$_] ) ] +=
          $weight[$_];
    }
    else {
        $dis_tu[ POSIX::ceil( POSIX::ceil( 5 * $tu / $fu ) * $dis[$_] ) ] +=
          $weight[$_];
    }
}
$dis_fu[1]  += $dis_fu[0];
$dis_cds[1] += $dis_cds[0];
$dis_tu[1]  += $dis_tu[0];

foreach my $bin ( 1 .. $#dis_fu ) {
    my $x = $bin - 0.5;
    my $y = $dis_fu[$bin];
    print "five_utr\t$x\t$y\n";
}

foreach my $bin ( 1 .. $#dis_cds ) {
    my $x = $bin - 0.5;
    my $y = $dis_cds[$bin];
    print "cds\t$x\t$y\n";
}

foreach my $bin ( 1 .. $#dis_tu ) {
    my $x = $bin - 0.5;
    my $y = $dis_tu[$bin];
    print "three_utr\t$x\t$y\n";
}
