#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use AlignDB::IntSpan;

#---------------#
# GetOpt section
#---------------#

=head1 NAME
merge.pl -- Merge the different transcript sites which on the same position of the genome.
=head1 SYNOPSIS
    perl merge.pl --refstr "ID=exon:" --geneid "ENST" -i input.tmp -o output.tmp
        Options:
            --help\-h  Brief help message
            --refstr  The sign before the gene ID in the reference
            --geneid  Characteristics shared by genes
            --in\-i  Input file with path
            --out\-o  Output file with path
=cut

Getopt::Long::GetOptions(
    'help|h'   => sub { Getopt::Long::HelpMessage(0) },
    'refstr=s' => \my $refstr,
    'geneid=s' => \my $geneid,
    'in|i=s'   => \my $in_fq,
    'out|o=s'  => \my $out_fq,
) or Getopt::Long::HelpMessage(1);

my %trans_range;
my %trans_chr;
my %trans_dir;
while (<>) {
    chomp;
    my ( $chr, $start, $end, $dir, $info ) = split /\t/;
    $chr  =~ s/chr//;
    $info =~ /$refstr([A-Z,a-z,0-9]+\.[0-9]+)/;
    if ( exists( $trans_chr{$1} ) ) {
        $trans_range{$1}->add_range( $start, $end );
    }
    else {
        $trans_chr{$1}   = $chr;
        $trans_dir{$1}   = $dir;
        $trans_range{$1} = AlignDB::IntSpan->new();
        $trans_range{$1}->add_range( $start, $end );
    }
}

sub COORDINATE_POS {
    my $INDEX = $_[0];
    my $SITE  = $_[1];
    my $ISLAND;
    if ( $trans_dir{$INDEX} eq "+" ) {
        $ISLAND = $trans_range{$INDEX}->at($SITE);
    }
    else {
        $ISLAND = $trans_range{$INDEX}->at( -$SITE );
    }
    my $ABS_SITE =
      $trans_chr{$INDEX} . "\t" . $ISLAND . "\t" . $trans_dir{$INDEX};
    return ($ABS_SITE);
}

open( my $IN_FH, "<", $in_fq );
my ( %gene_id, %gene, %base, %t5, %t3 );
while (<$IN_FH>) {
    chomp;
    my ( $trans_id, $site, $base1, $base2, $base3, $t5, $t3 ) = split /\t/;
    $trans_id =~ /^($geneid[A-Z,a-z,0-9]+\.[0-9]+)/;
    $trans_id = $1;
    $trans_id =~ /^($geneid[A-Z,a-z,0-9]+)/;
    my $gene_id = $1;
    my $gene    = $gene_id;
    my $asite   = COORDINATE_POS( $trans_id, $site );

    if ( exists( $gene{$asite} ) ) {
        unless ( $gene_id{$asite} =~ /$gene_id/ ) {
            $gene_id{$asite} .= "/" . $gene_id;
            $gene{$asite}    .= "/" . $gene;
        }
        $t5{$asite} += $t5;
        $t3{$asite} += $t3;
    }
    else {
        $gene_id{$asite} = $gene_id;
        $gene{$asite}    = $gene;
        $t5{$asite}      = $t5;
        $t3{$asite}      = $t3;
        $base{$asite}    = $base1 . "\t" . $base2 . "\t" . $base3;
    }
}
close($IN_FH);

open( my $OUT_FH, ">", $out_fq );
foreach ( keys(%gene) ) {
    print $OUT_FH (
        "$_\t$base{$_}\t$gene_id{$_}\t$gene{$_}\t$t5{$_}\t$t3{$_}\n");
}
close($OUT_FH);

__END__
