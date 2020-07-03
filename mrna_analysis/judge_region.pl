#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use AlignDB::IntSpan;

my ( %cse, %csi, %ase );
open( my $ALTERSPLICE_GENE, "<", $ARGV[0] );
while (<$ALTERSPLICE_GENE>) {
    chomp;
    my ( $chr, undef, $dir, $constant_exon, $constant_intron, $variable_area )
      = split /\t/;
    $chr =~ s/chr//;
    my $chr_dir = $chr . "\t" . $dir;
    if ( exists( $cse{$chr_dir} ) ) {
        $cse{$chr_dir}->AlignDB::IntSpan::add_runlist($constant_exon);
    }
    else {
        $cse{$chr_dir} = AlignDB::IntSpan->new($constant_exon);
    }
    if ( exists( $csi{$chr_dir} ) ) {
        $csi{$chr_dir}->AlignDB::IntSpan::add_runlist($constant_intron);
    }
    else {
        $csi{$chr_dir} = AlignDB::IntSpan->new($constant_intron);
    }
    if ( exists( $ase{$chr_dir} ) ) {
        $ase{$chr_dir}->AlignDB::IntSpan::add_runlist($variable_area);
    }
    else {
        $ase{$chr_dir} = AlignDB::IntSpan->new($variable_area);
    }
}
close($ALTERSPLICE_GENE);

my ( %utr5, %cds, %utr3 );
open( my $TRANSCRIPT_REGION_GENE, "<", $ARGV[1] );
while (<$TRANSCRIPT_REGION_GENE>) {
    chomp;
    my ( $chr, undef, $dir, undef, $utr5, $cds, $utr3 ) = split /\t/;
    $chr =~ s/chr//;
    my $chr_dir = $chr . "\t" . $dir;
    if ( exists( $utr5{$chr_dir} ) ) {
        $utr5{$chr_dir}->AlignDB::IntSpan::add_runlist($utr5);
    }
    else {
        $utr5{$chr_dir} = AlignDB::IntSpan->new($utr5);
    }
    if ( exists( $cds{$chr_dir} ) ) {
        $cds{$chr_dir}->AlignDB::IntSpan::add_runlist($cds);
    }
    else {
        $cds{$chr_dir} = AlignDB::IntSpan->new($cds);
    }
    if ( exists( $utr3{$chr_dir} ) ) {
        $utr3{$chr_dir}->AlignDB::IntSpan::add_runlist($utr3);
    }
    else {
        $utr3{$chr_dir} = AlignDB::IntSpan->new($utr3);
    }
}
close($TRANSCRIPT_REGION_GENE);

while (<STDIN>) {
    chomp;
    my ( $chr, $point, $dir ) = split /\t/;
    my $chr_dir = $chr . "\t" . $dir;
    my $set     = AlignDB::IntSpan->new($point);
    my $reg1    = 0;
    if ( ( exists( $cse{$chr_dir} ) ) and ( $set->subset( $cse{$chr_dir} ) ) ) {
        $reg1 = "CSE";
    }
    if ( ( exists( $csi{$chr_dir} ) ) and ( $set->subset( $csi{$chr_dir} ) ) ) {
        $reg1 = "CSI";
    }
    if ( ( exists( $ase{$chr_dir} ) ) and ( $set->subset( $ase{$chr_dir} ) ) ) {
        $reg1 = "ASE";
    }
    my $reg2 = 0;
    if ( ( exists( $utr5{$chr_dir} ) ) and ( $set->subset( $utr5{$chr_dir} ) ) )
    {
        $reg2 = "UTR5";
    }
    if ( ( exists( $utr3{$chr_dir} ) ) and ( $set->subset( $utr3{$chr_dir} ) ) )
    {
        $reg2 = "UTR3";
    }
    if ( ( exists( $cds{$chr_dir} ) ) and ( $set->subset( $cds{$chr_dir} ) ) ) {
        $reg2 = "CDS";
    }
    print("$_\t$reg1\t$reg2\n");
}

__END__
