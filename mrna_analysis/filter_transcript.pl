#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Getopt::Long;

Getopt::Long::GetOptions(
    'help|h'    => sub { Getopt::Long::HelpMessage(0) },
    'geneid=s'  => \my $geneid,
    'transid=s' => \my $transid,
) or Getopt::Long::HelpMessage(1);
defined($geneid) && length($geneid) && defined($transid) && length($transid)
    or die "Both --geneid and --transid attribute prefixes are required\n";

sub TAG_SCORE {
    my $INFO = shift;
    my $SCORE;
    my @TMP = ( $INFO =~ m/tag=([A-Z,a-z,0-9,_,\,]*)/ );
    my @TAG = split( /,/, $TMP[0] // "" );
    my ( $INCOMPLETE, $CDS_S, $CDS_E ) = ( 0, 0, 0, );
    $CDS_S      = grep /^cds_start_NF$/i, @TAG;
    $CDS_E      = grep /^cds_end_NF$/i,   @TAG;
    $INCOMPLETE = -( $CDS_S + $CDS_E ) * 10;
    my $BASIC = 0;
    $BASIC = grep /^basic$/i, @TAG;
    $BASIC *= 5;
    $SCORE = $BASIC + $INCOMPLETE;
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
    my ( undef, undef, $type, undef, undef, undef, undef, undef, $info ) =
      split /\t/;
    next unless defined($info) && $type =~ /^(?:mRNA|transcript)$/;
    my $trans_name = GET_ID( $transid, $info );
    my $gene_name  = GET_ID( $geneid,  $info );
    next unless defined($gene_name) && defined($trans_name);
    $gene_info{$gene_name}{$trans_name} = 1 if TAG_SCORE($info) >= 5;
}

foreach my $gene ( sort( keys(%gene_info) ) ) {
    foreach my $trans ( sort( keys( %{ $gene_info{$gene} } ) ) ) {
        print "$trans\n";
    }
}

__END__
