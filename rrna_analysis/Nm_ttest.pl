#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
# Historical unadjusted Welch test for three biological replicates per group.
# Not a candidate-calling stage of the current mRNA pipeline.
use Getopt::Long;
my $legacy;
GetOptions('legacy-three-replicates' => \$legacy) or die "Invalid options\n";
$legacy && @ARGV == 2 or die "Legacy 18+ column interface only: specify --legacy-three-replicates; columns 14-16 must be biological replicates, column 18 the original coverage gate\n";
require Statistics::R;

sub TTEST {
    my ( $a1, $a2 ) = @_;
    my $R = Statistics::R->new();
    $R->set( 'a1', \@{$a1} );
    $R->set( 'a2', \@{$a2} );
    $R->run('c <- tryCatch(t.test(a1,a2)$p.value, error=function(e) NA_real_)');
    my $p = $R->get("c");
    $R->stop();
    if ( defined($p) && $p =~ /^[0-9.eE+-]+$/ && $p < 0.05 ) {
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
    my $tmp = <$IN2>;
    defined($tmp) or die "Tables have different row counts\n";
    chomp $tmp;
    my @tmp2 = split( "\t", $tmp );
    @tmp1 >= 18 && @tmp2 >= 18 && $tmp1[0] eq $tmp2[0] && $tmp1[1] eq $tmp2[1]
        or die "Legacy table columns or coordinate order mismatch\n";
    for my $value (@tmp1[13..15,17], @tmp2[13..15,17]) {
        $value =~ /^-?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$/ or die "Non-numeric legacy measurement\n";
    }
    my @b = @tmp2[13..15];
    if (    ( ( $tmp1[17] > 1000 ) or ( $tmp2[17] > 1000 ) )
        and ( TTEST( \@a, \@b ) ) )
    {
        print "$tmp1[0]\t$tmp1[1]\n";
    }
}

eof($IN2) or die "Second table has extra rows\n";
__END__
