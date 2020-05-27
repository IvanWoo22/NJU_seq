#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open( my $in_sam1, "<", $ARGV[0] );
my %list;
while (<$in_sam1>) {
    chomp;
    my $id = split( /\t/, $_, 2 );
    if ( exists( $list{$id} ) ) {
        $list{$id}++;
    }
    else {
        $list{$id} = 1;
    }
}
close($in_sam1);

open( my $in_sam2,  "<", $ARGV[0] );
open( my $out_sam1, ">", $ARGV[1] );
open( my $out_sam2, ">", $ARGV[2] );
while (<$in_sam2>) {
    chomp;
    my $id = split( /\t/, $_, 2 );
    if ( $list{$id} == 1 ) {
        print $out_sam1 "$_\n";
    }
    else {
        print $out_sam2 "$_\n";
    }
}
close($in_sam2);
close($out_sam1);
close($out_sam2);

__END__
