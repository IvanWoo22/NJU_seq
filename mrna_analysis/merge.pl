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
            --transid Characteristics shared by transcripts
            --in\-i  Input file with path
            --out\-o  Output file with path
=cut

Getopt::Long::GetOptions(
    'help|h'    => sub { Getopt::Long::HelpMessage(0) },
    'refstr=s'  => \my $refstr,
    'geneid=s'  => \my $geneid,
    'transid=s' => \my $transid,
    'in|i=s'    => \my $in_fq,
    'out|o=s'   => \my $out_fq,
) or Getopt::Long::HelpMessage(1);

my %trans_range;
my %trans_chr;
my %trans_dir;
while (<>) {
    chomp;
    my ( $chr, $start, $end, $dir, $info ) = split /\t/;
    $chr  =~ s/chr//;
    $info =~ /$refstr(\w+\.[0-9]+)/;
    if ( exists( $trans_chr{$1} ) ) {
        $trans_range{$1}->AlignDB::IntSpan::add_range( $start, $end );
    }
    else {
        $trans_chr{$1}   = $chr;
        $trans_dir{$1}   = $dir;
        $trans_range{$1} = AlignDB::IntSpan->new();
        $trans_range{$1}->AlignDB::IntSpan::add_range( $start, $end );
    }
}

sub COORDINATE_POS {
    my $INDEX = $_[0];
    my $SITE  = $_[1];
    my $ISLAND;
    if ( $trans_dir{$INDEX} eq "+" ) {
        $ISLAND = $trans_range{$INDEX}->AlignDB::IntSpan::at($SITE);
    }
    else {
        $ISLAND = $trans_range{$INDEX}->AlignDB::IntSpan::at( -$SITE );
    }
    my $ABS_SITE =
      $trans_chr{$INDEX} . "\t" . $ISLAND . "\t" . $trans_dir{$INDEX};
    return ($ABS_SITE);
}

open( my $IN_FH, "<", $in_fq );
my ( %gene_id, %base, %t5, %t3 );
while (<$IN_FH>) {
    chomp;
    my ( $info, $site, $base1, $base2, $base3, $t5, $t3 ) = split /\t/;
    $info =~ /($transid\w+\.[0-9]+)/;
    my $trans_id = $1;
    $info =~ /($geneid\w+)/;
    my $gene_id  = $1;
    my $abs_site = COORDINATE_POS( $trans_id, $site );

    if ( exists( $gene_id{$abs_site} ) ) {
        unless ( $gene_id{$abs_site} =~ /$gene_id/ ) {
            $gene_id{$abs_site} .= "/" . $gene_id;
        }
        $t5{$abs_site} += $t5;
        $t3{$abs_site} += $t3;
    }
    else {
        $gene_id{$abs_site} = $gene_id;
        $t5{$abs_site}      = $t5;
        $t3{$abs_site}      = $t3;
        $base{$abs_site}    = $base1 . "\t" . $base2 . "\t" . $base3;
    }
}
close($IN_FH);

open( my $OUT_FH, ">", $out_fq );
foreach ( keys(%gene_id) ) {
    print $OUT_FH ("$_\t$base{$_}\t$gene_id{$_}\t$t5{$_}\t$t3{$_}\n");
}
close($OUT_FH);

__END__
