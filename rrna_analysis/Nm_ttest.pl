#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Statistics::R;

sub TTEST {
    my ( $a1, $a2 ) = @_;
    my $R = Statistics::R->new();
    $R->set( 'a1', \@{$a1} );
    $R->set( 'a2', \@{$a2} );
    $R->run('c <- t.test(a1,a2)$p.value');
    if ( $R->get('c') < 0.05 ) {
        return (1);
    }
    else {
        return (0);
    }
}

open( my $IN1, "<", $ARGV[0] );
open( my $IN2, "<", $ARGV[1] );
while (<$IN1>) {
    chomp;
    my @tmp1 = split "\t";
    my @a    = @tmp1[ 13 .. 15 ];
    chomp( my $tmp = <$IN2> );
    my @tmp2 = split( "\t", $tmp );
    my @b    = @tmp2[ 13 .. 15 ];
    if (    ( ( $tmp1[17] > 1000 ) or ( $tmp2[17] > 1000 ) )
        and ( TTEST( \@a, \@b ) ) )
    {
        print "$tmp1[0]\t$tmp1[1]\n";
    }
}

__END__
