#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use AlignDB::IntSpan;

while (<STDIN>) {
    chomp;
    my @tmp = split "\t";
    my $set = AlignDB::IntSpan->new( $tmp[3] );
    if (   ( $set->at_island(1) eq $set->find_islands( $tmp[0] ) )
        or ( $set->at_island(1) eq $set->find_islands( $tmp[0] ) ) )
    {
        print("$_\n");
    }
}
