#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

my ( %cer, %cir, %vr );
open( my $ALTERSPLICEGENE, "<", $ARGV[0] );
while (<$ALTERSPLICEGENE>) {
    chomp;
    my ( $chr, undef, $dir, $constant_exon, $constant_intron, $variable_area )
      = split /\t/;
    $chr =~ s/chr//;
    my $chr_dir = $chr . "\t" . $dir;
    if ( exists( $cer{$chr_dir} ) ) {
        $cer{$chr_dir}->AlignDB::IntSpan::add_runlist($constant_exon);
    }
    else {
        $cer{$chr_dir} = AlignDB::IntSpan->new($constant_exon);
    }
    if ( exists( $cir{$chr_dir} ) ) {
        $cir{$chr_dir}->AlignDB::IntSpan::add_runlist($constant_intron);
    }
    else {
        $cir{$chr_dir} = AlignDB::IntSpan->new($constant_intron);
    }
    if ( exists( $vr{$chr_dir} ) ) {
        $vr{$chr_dir}->AlignDB::IntSpan::add_runlist($variable_area);
    }
    else {
        $vr{$chr_dir} = AlignDB::IntSpan->new($variable_area);
    }
}
close($ALTERSPLICEGENE);

while (<STDIN>) {
    chomp;
    my ( $chr, $point, $dir ) = split /\t/;
    my $chr_dir = $chr . "\t" . $dir;
    my $set     = AlignDB::IntSpan->new($point);
    my $reg     = 0;
    if ( ( exists( $cer{$chr_dir} ) ) and ( $set->subset( $cer{$chr_dir} ) ) ) {
        $reg = "CER";
    }
    if ( ( exists( $cir{$chr_dir} ) ) and ( $set->subset( $cir{$chr_dir} ) ) ) {
        $reg = "CIR";
    }
    if ( ( exists( $vr{$chr_dir} ) ) and ( $set->subset( $vr{$chr_dir} ) ) ) {
        $reg = "VR";
    }
    print("$_\t$reg\n");
}

__END__
