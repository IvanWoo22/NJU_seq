#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my $filedir = dirname(__FILE__);
my $envdir = $ENV{'PWD'};

open( my $IN, "<", $ARGV[0] );
open(my $OUT, ">", $envdir."/".$ARGV[0].".tmp");

my %sig;
while (<$IN>) {
    chomp;
    my ( undef, undef, undef, $lef, $nm, $rgt, undef ) = split;
    if ( exists( $sig{ $lef . $nm . $rgt } ) ) {
        $sig{ $lef . $nm . $rgt }++;
    }
    else {
        $sig{ $lef . $nm . $rgt } = 1;
    }
}

my $data_table;
foreach my $nm (qw(A G C T)) {
    foreach my $lef (qw(A G C T)) {
        foreach my $rgt (qw(A G C T)) {
            $data_table .= "$lef$nm$rgt\t$sig{$lef.$nm.$rgt}\n";
        }
    }
}

print $OUT ("$data_table");

system(
    "Rscript $filedir/picture_draw.R $envdir/$ARGV[0].tmp $envdir/$ARGV[1]"
);