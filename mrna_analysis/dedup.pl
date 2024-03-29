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
dedup.pl -- Deduplication by finding out these two sites of transcripts located on the same position.
=head1 SYNOPSIS
    perl dedup.pl --refstr "Parent=" --transid "ENST" --info data/hsa_exon.info
=cut

Getopt::Long::GetOptions(
    'help|h'    => sub { Getopt::Long::HelpMessage(0) },
    'refstr=s'  => \my $refstr,
    'transid=s' => \my $transid,
    'info=s'    => \my $info_fh,
) or Getopt::Long::HelpMessage(1);

open( my $IN_FH, "<", $info_fh );
my %trans_range;
my %trans_chr;
my %trans_dir;
while (<$IN_FH>) {
    chomp;
    my ( $chr, $start, $end, $dir, $info ) = split( /\t/, $_ );
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
close($IN_FH);

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
    my $ABS_SITE = $trans_chr{$INDEX} . "\t" . $ISLAND;
    return ($ABS_SITE);
}

my %exist;
while (<>) {
    chomp;
    my ( $read_name, $trans_info, $site ) = split /\s+/;
    $trans_info =~ /^($transid\w+\.[0-9]+)/;
    my $id            = $1;
    my $abs_site      = COORDINATE_POS( $id, $site );
    my $read_abs_site = $read_name . "\t" . $abs_site;
    unless ( exists( $exist{$read_abs_site} ) ) {
        $exist{$read_abs_site} = 1;
        print("$_\n");
    }
}

__END__
