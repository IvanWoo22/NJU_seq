#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

while (<>) {
    chomp;
    my (
        $id,            undef,        undef,           undef,
        $exon_length,   $exon_number, undef,           $intron_length,
        $intron_number, undef,        $variant_length, $variant_number
    ) = split /\t/;
    my ( $exon_cov, $intron_cov, $variant_cov );
    if ( $exon_length > 0 ) {
        $exon_cov = $exon_number / $exon_length;
    }
    else {
        $exon_cov = 0;
    }
    if ( $intron_length > 0 ) {
        $intron_cov = $intron_number / $intron_length;
    }
    else {
        $intron_cov = 0;
    }
    if ( $variant_length > 0 ) {
        $variant_cov = $variant_number / $variant_length;
    }
    else {
        $variant_cov = 0;
    }

    print "$id\t$exon_cov\t$intron_cov\t$variant_cov\n";
}

__END__
