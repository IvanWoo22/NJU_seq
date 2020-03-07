#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

my %point_set;
while (<STDIN>) {
    chomp;
    my ( $chr, $dir, $point ) = split /\t/;
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
        $chr,             $gene_id,    $dir,
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
    my ( $five_utr_number, $cds_number, $three_utr_number );
    if ( exists( $point_set{$chr_dir} ) ) {
        $five_utr_number =
          AlignDB::IntSpan::intersect( $point_set{$chr_dir}, $five_utr_set )
          ->cardinality;
        $cds_number =
          AlignDB::IntSpan::intersect( $point_set{$chr_dir}, $cds_set )
          ->cardinality;
        $three_utr_number =
          AlignDB::IntSpan::intersect( $point_set{$chr_dir}, $three_utr_set )
          ->cardinality;
    }
    else {
        $five_utr_number  = 0;
        $cds_number       = 0;
        $three_utr_number = 0;
    }
    my $stat_val;

    if ( $five_utr_number + $cds_number + $three_utr_number >= 5 ) {
        $stat_val = "T";
    }
    else {
        $stat_val = "F";
    }

    print(
"$gene_id\t$chr_dir\t$five_utr_string\t$five_utr_length\t$five_utr_number\t$cds_string\t$cds_length\t$cds_number\t$three_utr_string\t$three_utr_length\t$three_utr_number\t$stat_val\n"
    );
}
close($GENE);

__END__
