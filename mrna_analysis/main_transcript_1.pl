#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Getopt::Long;

Getopt::Long::GetOptions(
    'help|h'    => sub { Getopt::Long::HelpMessage(0) },
    'geneid=s'  => \my $geneid,
    'transid=s' => \my $transid,
    'tag'       => \my $use_tag,
) or Getopt::Long::HelpMessage(1);

sub TAG_SCORE {
    my $INFO = shift;
    my $SCORE;
    my $TSL;
    if ( $INFO =~ m/transcript_support_level=([1-5]);/ ) {
        $TSL = 15 - 5 * $1;
    }
    else {
        $TSL = -10;
    }
    my @TMP  = ( $INFO =~ m/tag=([A-Z,a-z,0-9,_,\,]*)/ );
    my @TAG = ();
    if (defined $TMP[0]) {
        @TAG = split( /,/, $TMP[0] );
    }
    my $MANE = 0;
    $MANE = grep /^MANE_Select$/i, @TAG;
    $MANE *= 30;
    my $AP;
    if ( $INFO =~ m/appris_principal_([1-5])/ ) {
        $AP = 30 - 5 * $1;
    }
    else {
        $AP = 0;
    }
    my $ALT;
    if ( $INFO =~ m/appris_alternative_([1-2])/ ) {
        $ALT = 15 - 5 * $1;
    }
    else {
        $ALT = 0;
    }
    my ( $INCOMPLETE, $MRNA_S, $MRNA_E, $CDS_S, $CDS_E ) = ( 0, 0, 0, 0, 0 );
    $MRNA_S     = grep /^mRNA_start_NF$/i, @TAG;
    $MRNA_E     = grep /^mRNA_end_NF$/i,   @TAG;
    $CDS_S      = grep /^cds_start_NF$/i,  @TAG;
    $CDS_E      = grep /^cds_end_NF$/i,    @TAG;
    $INCOMPLETE = -( $CDS_S + $CDS_E ) * 100 - ( $MRNA_S + $MRNA_E ) * 10;
    my $BASIC = grep /^basic$/i, @TAG;
    $BASIC *= 60;
    $BASIC = -60 if $BASIC == 0;
    my $CCDS = grep /^basic$/i, @TAG;
    $CCDS *= 5;
    $SCORE = $TSL + $MANE + $BASIC + $INCOMPLETE + $AP + $ALT + $CCDS;
    return ($SCORE);
}

sub GET_ID {
    my ( $FEATURE, $INFO ) = @_;
    my $ID;
    if ( $INFO =~ m/$FEATURE(\w+\.?[0-9]*)/ ) {
        $ID = $1;
    }
    else {
        warn("There is a problem in $INFO;\n");
    }
    return $ID;
}

my %gene_data;
while (<STDIN>) {
    chomp;
    my ( undef, undef, undef, $start, $end, undef, undef, undef, $info ) =
      split /\t/;
    my $trans_name = GET_ID( $transid, $info );
    my $gene_name  = GET_ID( $geneid,  $info );
    
    if (defined $trans_name && defined $gene_name) {
        if ($use_tag) {
            my $score = TAG_SCORE($info);
            if ( $score > 20 ) {
                $gene_data{$gene_name}{$trans_name} = $score;
            }
        } else {
            my $length = abs($end - $start) + 1;
            $gene_data{$gene_name}{$trans_name} = $length;
        }
    }
}

foreach my $gene ( sort( keys(%gene_data) ) ) {
    my $number = 0;
    foreach my $trans (
        sort { ${ $gene_data{$gene} }{$b} <=> ${ $gene_data{$gene} }{$a} }
        keys %{ $gene_data{$gene} }
      )
    {
        $number++;
        print "$trans\t$gene\t$number\n";
    }
}

__END__
