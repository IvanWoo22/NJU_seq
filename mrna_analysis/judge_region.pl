#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

my ( %cse, %csi, %ase );
open( my $ALTERSPLICEGENE, "<", $ARGV[0] );
while (<$ALTERSPLICEGENE>) {
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
close($ALTERSPLICEGENE);

while (<STDIN>) {
    chomp;
    my ( $chr, $point, $dir ) = split /\t/;
    my $chr_dir = $chr . "\t" . $dir;
    my $set     = AlignDB::IntSpan->new($point);
    my $reg     = 0;
    if ( ( exists( $cse{$chr_dir} ) ) and ( $set->subset( $cse{$chr_dir} ) ) ) {
        $reg = "CSE";
    }
    if ( ( exists( $csi{$chr_dir} ) ) and ( $set->subset( $csi{$chr_dir} ) ) ) {
        $reg = "CSI";
    }
    if ( ( exists( $ase{$chr_dir} ) ) and ( $set->subset( $ase{$chr_dir} ) ) ) {
        $reg = "ASE";
    }
    print("$_\t$reg\n");
}

__END__
