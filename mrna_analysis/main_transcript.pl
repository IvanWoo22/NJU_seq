#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use Getopt::Long;

Getopt::Long::GetOptions(
    'help|h'    => sub { Getopt::Long::HelpMessage(0) },
    'geneid=s'  => \my $geneid,
    'transid=s' => \my $transid,
) or Getopt::Long::HelpMessage(1);

sub TAG_SCORE {
    my $INFO = shift;
    my $SCORE;
    my $TSL;
    if ( $INFO =~ m/transcript_support_level=([1-5]);/ ) {
        $TSL = 3 - $1;
    }
    else {
        $TSL = -3;
    }
    my @TMP  = ( $INFO =~ m/tag=([A-Z,a-z,0-9,_,\,]*)/ );
    my @TAG  = split( /,/, $TMP[0] );
    my $MANE = 0;
    $MANE = grep /^MANE_Select$/i, @TAG;
    $MANE *= 3;
    my $AP;
    if ( $INFO =~ m/appris_principal_([1-5])/ ) {
        $AP = 6 - $1;
    }
    else {
        $AP = 0;
    }
    my $ALT;
    if ( $INFO =~ m/appris_alternative_([1-2])/ ) {
        $ALT = 3 - $1;
    }
    else {
        $ALT = 0;
    }
    my ( $INCOMPLETE, $MRNA_S, $MRNA_E, $CDS_S, $CDS_E ) = ( 0, 0, 0, 0, 0 );
    $MRNA_S = grep /^mRNA_start_NF$/i, @TAG;
    $MRNA_E = grep /^mRNA_end_NF$/i,   @TAG;
    $CDS_S  = grep /^cds_start_NF$/i,  @TAG;
    $CDS_E  = grep /^cds_end_NF$/i,    @TAG;
    $INCOMPLETE = -( $MRNA_S + $MRNA_E + $CDS_S + $CDS_E ) * 2;
    my $BASIC = 0;
    $BASIC = grep /^basic$/i, @TAG;
    $BASIC *= 3;
    $SCORE = $TSL + $MANE + $BASIC + $INCOMPLETE + $AP + $ALT;
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

my %gene_info;
while (<STDIN>) {
    chomp;
    my ( undef, undef, undef, undef, undef, undef, undef, undef, $info ) =
      split /\t/;
    my $trans_name = GET_ID( $transid, $info );
    my $gene_name  = GET_ID( $geneid,  $info );
    $gene_info{$gene_name}{$trans_name} = 1 if TAG_SCORE($info) > 0;
}

foreach my $gene ( sort( keys(%gene_info) ) ) {
    foreach my $trans ( sort( keys( %{ $gene_info{$gene} } ) ) ) {
        print "$trans\t$gene\n";
    }
}

__END__
