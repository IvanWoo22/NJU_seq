#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

while (<STDIN>) {
    chomp;
    my (
        $id,   undef,             undef,
        undef, $five_utr_length,  $five_utr_number,
        undef, $cds_length,       $cds_number,
        undef, $three_utr_length, $three_utr_number,
        undef
    ) = split /\t/;
    my ( $five_utr_cov, $cds_cov, $three_utr_cov );
    if ( $five_utr_length > 0 ) {
        $five_utr_cov = $five_utr_number / $five_utr_length;
    }
    else {
        $five_utr_cov = 0;
    }
    if ( $cds_length > 0 ) {
        $cds_cov = $cds_number / $cds_length;
    }
    else {
        $cds_cov = 0;
    }
    if ( $three_utr_length > 0 ) {
        $three_utr_cov = $three_utr_number / $three_utr_length;
    }
    else {
        $three_utr_cov = 0;
    }
    print "$id\t$five_utr_cov\t$cds_cov\t$three_utr_cov\n";
}

__END__
