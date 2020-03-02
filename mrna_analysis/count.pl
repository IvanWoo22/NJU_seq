#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

sub COUNT_ALIGNMENT_LENGTH {
    my $CIGAR  = shift;
    my $LENGTH = 0;
    while ( $CIGAR =~ /([0-9]+)=/g )
    {    # Alignment column containing two identical letters.
        $LENGTH += $1;
    }
    while ( $CIGAR =~ /([0-9]+)X/g )
    {    # Alignment column containing a mismatch, i.e. two different letters.
        $LENGTH += $1;
    }
    while ( $CIGAR =~ /([0-9]+)D/g ) {  # Deletion (gap in the target sequence).
        $LENGTH += $1;
    }
    return $LENGTH;
}

open( my $IN_SAM, "<", $ARGV[0] );

my %site_end;
my %site_start;
my %site_base;

while (<$IN_SAM>) {
    chomp;
    my ( undef, $rname, $site, $cigar, $seq ) = split /\t/;
    my $length  = COUNT_ALIGNMENT_LENGTH($cigar);
    my $end_pos = $site + $length - 1;
    my $end     = $rname . "\t" . $end_pos;
    my @string  = split( //, $seq );

    if ( exists $site_end{$end} ) {
        $site_end{$end}++;
    }
    else {
        $site_end{$end}  = 1;
        $site_base{$end} = $string[-2];
    }

    my $start = $rname . "\t" . $site;
    if ( exists $site_start{$start} ) {
        $site_start{$start}++;
    }
    else {
        $site_start{$start} = 1;
    }
}
close($IN_SAM);

foreach ( keys %site_end ) {
    if ( exists $site_start{$_} ) {
        print "$_\t$site_base{$_}\t$site_start{$_}\t$site_end{$_}\n";
        delete $site_start{$_};
    }
    else {
        print "$_\t$site_base{$_}\t0\t$site_end{$_}\n";
    }
}

__END__
