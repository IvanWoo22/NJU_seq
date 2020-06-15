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
    my ( $chr, $gene_id, $dir, $constant_exon, $constant_intron,
        $variable_area ) = split /\t/;
    $chr =~ s/chr//;
    my $exon_set       = AlignDB::IntSpan->new($constant_exon);
    my $intron_set     = AlignDB::IntSpan->new($constant_intron);
    my $variant_set    = AlignDB::IntSpan->new($variable_area);
    my $chr_dir        = $chr . "\t" . $dir;
    my $exon_length    = $exon_set->cardinality;
    my $intron_length  = $intron_set->cardinality;
    my $variant_length = $variant_set->cardinality;
    my ( $exon_number, $intron_number, $variant_number ) = ( 0, 0, 0 );

    if ( exists( $point_set{$chr_dir} ) ) {
        $exon_number =
          AlignDB::IntSpan::intersect( $point_set{$chr_dir}, $exon_set )
          ->cardinality;
        $intron_number =
          AlignDB::IntSpan::intersect( $point_set{$chr_dir}, $intron_set )
          ->cardinality;
        $variant_number =
          AlignDB::IntSpan::intersect( $point_set{$chr_dir}, $variant_set )
          ->cardinality;
    }
    my $stat_val;

    if ( $exon_number + $intron_number + $variant_number >= 5 ) {
        $stat_val = "T";
    }
    else {
        $stat_val = "F";
    }

# Gene_ID Chromosome Direction Exon_Range Exon_Length Exon_Point Intron_Range Intron_Length Intron_Point Variant_Range Variant_Length Variant_Point Stat_Value
    print(
"$gene_id\t$chr_dir\t$constant_exon\t$exon_length\t$exon_number\t$constant_intron\t$intron_length\t$intron_number\t$variable_area\t$variant_length\t$variant_number\t$stat_val\n"
    );

}

close($GENE);

__END__
