#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

# Diagnostic only: forward R1 flank versus the reverse-complement tag.
use Getopt::Long;
my $flank_length = 6;
GetOptions('length=i' => \$flank_length) or die "Invalid options\n";
@ARGV == 2 && $flank_length > 0 or die "Usage: tn5_check.pl [--length N] transcript.fa alignments.sam\n";
open( my $FASTA, "<", $ARGV[0] );
open( my $SEG,   "<", $ARGV[1] );

my %fasta;
my $title_name;
while (<$FASTA>) {
    if (/^>(\S+)/) {
        $title_name = $1;
    }
    else {
        $_ =~ s/\r?\n//;
        $fasta{$title_name} .= $_;
    }
}
close($FASTA);

sub SEQ_REV_COMP {
    my $SEQ = reverse shift;
    $SEQ =~ tr/Uu/Tt/;
    return ( $SEQ =~ tr/NAGCTagct/NTCGAtcga/r );
}

sub FIND {
    my $TITLE  = shift;
    my $START  = shift;
    my $LENGTH = shift;
    my $REFSEQ = shift;
    my $SEQ;
    if (    ( exists( $fasta{$TITLE} ) )
        and ( $START >= 1 && length( $fasta{$TITLE} ) >= $START + $LENGTH - 1 ) )
    {
        $SEQ = substr( $fasta{$TITLE}, $START - 1, $LENGTH );
    }
    else {
        $SEQ = "N" x $LENGTH;
        #warn("Sorry, there is no such a segment: $TITLE\n");
    }
    my @SEQ       = split( //, $SEQ );
    my @REFSEQ    = split( //, $REFSEQ );
    my $MATCH_NUM = 0;
    foreach ( 0 .. $LENGTH - 1 ) {
        $MATCH_NUM++ if $SEQ[$_] eq $REFSEQ[$_];
    }
    return ("$SEQ\t$MATCH_NUM");
}

while (<$SEG>) {
    next if /^@/;
    chomp;
    my @f = split /\t/;
    @f >= 11 or die "Expected SAM, not a five-column filtered table\n";
    next if $f[1] & (4 | 16 | 128); # unmapped, reverse or read 2 unsupported
    my ($qname, $title, $start, $cigar) = @f[0,2,3,5];
    $qname =~ /^(.*)_([ACGTUNacgtun]+)$/ or die "Missing terminal UMI: $qname\n";
    my ($id, $umi) = ($1, uc $2);
    length($umi) >= $flank_length or die "UMI shorter than --length\n";
    $cigar =~ /^(?:[1-9][0-9]*[=XMID])+$/ or die "Require end-to-end CIGAR\n";
    my $span = 0;
    $span += $1 while $cigar =~ /([0-9]+)[=XMD]/g;
    my ($fasta_title) = split /\|/, $title;
    my $output = FIND($fasta_title, $start + $span, $flank_length, SEQ_REV_COMP($umi));
    print "$id\t$umi\t$output\n";
}
close($SEG);

__END__
