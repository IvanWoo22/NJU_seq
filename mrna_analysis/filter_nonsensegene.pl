#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my ( @line, @gene );
my $i = 0;
while (<>) {
    chomp;
    push( @line, $_ );
    my ( undef, undef, $type, undef, undef, undef, undef, undef, undef ) =
      split /\t/;
    if ( $type eq "gene" ) {
        push( @gene, $i );
    }
    $i++;
}

foreach ( 0 .. $#gene - 1 ) {
    if ( $gene[ $_ + 1 ] - $gene[$_] > 1 ) {
        foreach my $t ( $gene[$_] .. $gene[ $_ + 1 ] - 1 ) {
            print("$line[$t]\n");
        }
    }
}
if ( $gene[-1] < $#line ) {
    foreach my $t ( $gene[-1] .. $#line ) {
        print("$line[$t]\n");
    }
}

__END__
