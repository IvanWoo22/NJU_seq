#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use AlignDB::IntSpan;

=head1 NAME
multimatch_handle.pl -- Handle the reads matched different sites.
=head1 SYNOPSIS
    perl multimatch_handle.pl --refstr "ID=exon:" --geneid "ENST" -i input.tmp -o output.tmp
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
    my ( $chr, $start, $end, $dir, $info ) = split( /\t/, $_ );
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
    my $ABS_SITE = $trans_chr{$INDEX} . "\t" . $ISLAND;
    return ($ABS_SITE);
}

open( my $IN_FH,  "<", $in_fq );
open( my $OUT_FH, ">", $out_fq );
my %exist;
while (<$IN_FH>) {
    chomp;
    my ( $read_name, $trans_info, $site, undef, $seq ) = split /\s+/;
    $trans_info =~ /^($geneid[A-Z,a-z,0-9]+\.[0-9]+)/;
    my $id            = $1;
    my $abs_site      = COORDINATE_POS( $id, $site );
    my $read_abs_site = $read_name . "\t" . $abs_site;
    unless ( exists( $exist{$read_abs_site} ) ) {
        $exist{$read_abs_site} = 1;
        print $OUT_FH ("$read_abs_site\t$seq\n");
    }
}
close($IN_FH);
close($OUT_FH);

__END__
