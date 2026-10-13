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

$transid //= "";
defined($refstr) && defined($info_fh) or die "Require --refstr and --info\n";
open( my $IN_FH, "<", $info_fh );
my %trans_range;
my %trans_chr;
my %trans_dir;
while (<$IN_FH>) {
    chomp;
    my ( $chr, $start, $end, $dir, $info ) = split( /\t/, $_ );
    $chr  =~ s/^chr//;
    # An unmatched annotation must not reuse the previous regex capture.
    next unless $info =~ /$refstr($transid\w+\.[0-9]+)/;
    my $transcript_id = $1;
    if ( exists( $trans_chr{$transcript_id} ) ) {
        $trans_range{$transcript_id}->AlignDB::IntSpan::add_range( $start, $end );
    }
    else {
        $trans_chr{$transcript_id}   = $chr;
        $trans_dir{$transcript_id}   = $dir;
        $trans_range{$transcript_id} = AlignDB::IntSpan->new();
        $trans_range{$transcript_id}->AlignDB::IntSpan::add_range( $start, $end );
    }
}
close($IN_FH);

sub COORDINATE_POS {
    my $INDEX = $_[0];
    my $SITE  = $_[1];
    exists $trans_range{$INDEX} or die "Transcript absent from annotation: $INDEX\n";
    $SITE =~ /^[1-9][0-9]*$/ && $SITE <= $trans_range{$INDEX}->size
        or die "Coordinate outside transcript $INDEX: $SITE\n";
    my $ISLAND;
    if ( $trans_dir{$INDEX} eq "+" ) {
        $ISLAND = $trans_range{$INDEX}->AlignDB::IntSpan::at($SITE);
    }
    else {
        $ISLAND = $trans_range{$INDEX}->AlignDB::IntSpan::at( -$SITE );
    }
    my $ABS_SITE = join "\t", $trans_chr{$INDEX}, $trans_dir{$INDEX}, $ISLAND;
    return ($ABS_SITE);
}

my %exist;
while (<>) {
    chomp;
    my ( $read_name, $trans_info, $site ) = split /\s+/;
    $trans_info =~ /^($transid\w+\.[0-9]+)/ or die "Invalid transcript: $trans_info\n";
    my $id            = $1;
    my $abs_site      = COORDINATE_POS( $id, $site );
    my $read_abs_site = $read_name . "\t" . $abs_site;
    unless ( exists( $exist{$read_abs_site} ) ) {
        $exist{$read_abs_site} = 1;
        print("$_\n");
    }
}

__END__
