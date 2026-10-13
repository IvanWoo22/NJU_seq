#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
my $motif_col;
GetOptions('motif-col=i' => \$motif_col) or die "Invalid options\n";
@ARGV == 2 && (!defined($motif_col) || $motif_col > 0) or die "Usage: motif_compare.pl [--motif-col N] mirna.tsv sites_with_motif.tsv (default: last column)\n";
sub NORMALIZE { my $s = uc shift; $s =~ tr/T/U/; return $s; }
open( my $MIRNA_FH, "<", $ARGV[0] );
open( my $NM_FH,    "<", $ARGV[1] );

my %mirna_motif;
while (<$MIRNA_FH>) {
    chomp;
    my ( $name, $motif ) = split /\t/;
    defined($motif) or die "Invalid miRNA motif row\n";
    $motif = NORMALIZE($motif);
    if ( exists( $mirna_motif{$motif} ) ) {
        $mirna_motif{$motif} .= ",$name";
    }
    else {
        $mirna_motif{$motif} = $name;
    }
}
close($MIRNA_FH);
while (<$NM_FH>) {
    chomp;
    next if /^#/ || /^(?:Chromosome|Position)\t/;
    my @f = split /\t/;
    my $motif = $f[defined($motif_col) ? $motif_col - 1 : -1];
    defined($motif) && $motif =~ /^[ACGTURYSWKMBDHVNacgturyswkmbdhvn]+$/ or die "Motif column does not contain a sequence; specify --motif-col\n";
    $motif = NORMALIZE($motif);
    if ( exists( $mirna_motif{$motif} ) ) {
        my $mirna_num = split( /,/, $mirna_motif{$motif} );
        print("$_\t$mirna_num\t$mirna_motif{$motif}\n");
    }
    else {
        print("$_\t0\tNULL\n");
    }
}
close($NM_FH);
__END__
