#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

my ( @line, @gene, @gene_range, @gene_position, @mark );
my $i = 0;
while (<>) {
    chomp;
    push( @line, $_ );
    my ( $chr, undef, $type, $start, $end, undef, $dir, undef, undef ) =
      split /\t/;
    if ( $type eq "gene" ) {
        push( @gene, $i );
        my $position = $chr . "\t" . $dir;
        push( @gene_position, $position );
        my $range = AlignDB::IntSpan->new;
        $range->add_range( $start, $end );
        push( @gene_range, $range );
    }
    $i++;
}

$mark[0] = 0;
foreach ( 1 .. $#gene ) {
    $mark[$_] = 0;
    if (
            ( $gene_position[ $_ - 1 ] eq $gene_position[$_] )
        and
        ( $gene_range[ $_ - 1 ]->AlignDB::IntSpan::overlap( $gene_range[$_] )
            > 0 )
      )
    {
        $mark[ $_ - 1 ] = 1;
        $mark[$_] = 1;
    }
}

foreach ( 0 .. $#gene ) {
    if ( ( $mark[$_] == 0 ) and ( $_ < $#gene ) ) {
        foreach my $t ( $gene[$_] .. $gene[ $_ + 1 ] - 1 ) {
            print("$line[$t]\n");
        }
    }
    elsif ( ( $mark[$_] == 0 ) and ( $_ == $#gene ) ) {
        foreach my $t ( $gene[$_] .. $#line ) {
            print("$line[$t]\n");
        }
    }
}

__END__
