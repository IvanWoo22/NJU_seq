#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

my %point_set;
while (<STDIN>) {
    chomp;
    my ( $chr, $point, $dir ) = split /\t/;
    my $chr_dir = $chr . "\t" . $dir;
    if ( $dir eq "+" ) {
        $point--;
    }
    else {
        $point++;
    }
    if ( exists( $point_set{$chr_dir} ) ) {
        $point_set{$chr_dir}->AlignDB::IntSpan::add($point);
    }
    else {
        $point_set{$chr_dir} = AlignDB::IntSpan->new;
        $point_set{$chr_dir}->AlignDB::IntSpan::add($point);
    }
}

open( my $GENE, "<", $ARGV[0] );
while (<$GENE>) {
    chomp;
    my (
        $chr,             undef,       $dir,
        undef,            undef,       undef,
        $five_utr_string, $cds_string, $three_utr_string
    ) = split /\t/;
    $chr =~ s/chr//;
    my $chr_dir = $chr . "\t" . $dir;

    my $five_utr_set     = AlignDB::IntSpan->new($five_utr_string);
    my $cds_set          = AlignDB::IntSpan->new($cds_string);
    my $three_utr_set    = AlignDB::IntSpan->new($three_utr_string);
    my $five_utr_length  = $five_utr_set->cardinality;
    my $cds_length       = $cds_set->cardinality;
    my $three_utr_length = $three_utr_set->cardinality;
    my ( @five_utr_point, @cds_point, @three_utr_point );
    if ( exists( $point_set{$chr_dir} ) ) {
        @five_utr_point =
          AlignDB::IntSpan::intersect( $point_set{$chr_dir}, $five_utr_set )
          ->as_array;
        @cds_point =
          AlignDB::IntSpan::intersect( $point_set{$chr_dir}, $cds_set )
          ->as_array;
        @three_utr_point =
          AlignDB::IntSpan::intersect( $point_set{$chr_dir}, $three_utr_set )
          ->as_array;
    }
    if ( $dir eq "+" ) {
        foreach (@five_utr_point) {
            my $index        = $five_utr_set->index($_);
            my $distribution = $index / $five_utr_length;
            print("$chr_dir\t$_\tfive_utr\t$distribution\n");
        }
        foreach (@cds_point) {
            my $index        = $cds_set->index($_);
            my $distribution = $index / $cds_length;
            print("$chr_dir\t$_\tcds\t$distribution\n");
        }
        foreach (@three_utr_point) {
            my $index        = $three_utr_set->index($_);
            my $distribution = $index / $three_utr_length;
            print("$chr_dir\t$_\tthree_utr\t$distribution\n");
        }
    }
    else {
        foreach (@five_utr_point) {
            my $index        = $five_utr_set->index($_);
            my $distribution = 1 - $index / $five_utr_length;
            print("$chr_dir\t$_\tfive_utr\t$distribution\n");
        }
        foreach (@cds_point) {
            my $index        = $cds_set->index($_);
            my $distribution = 1 - $index / $cds_length;
            print("$chr_dir\t$_\tcds\t$distribution\n");
        }
        foreach (@three_utr_point) {
            my $index        = $three_utr_set->index($_);
            my $distribution = 1 - $index / $three_utr_length;
            print("$chr_dir\t$_\tthree_utr\t$distribution\n");
        }
    }
}
close($GENE);

__END__
