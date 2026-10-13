#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

# Return the observed bases at the final three reference positions, not
# the final three query positions (which differ when an indel is present).
sub ALIGNMENT_TAIL {
    my ($cigar, $seq) = @_;
    $cigar =~ /^(?:[1-9][0-9]*[=XID])+$/ or die "Invalid --xeq CIGAR: $cigar\n";
    my ($query, $length, $tail) = (0, 0, "");
    while ($cigar =~ /([0-9]+)([=XID])/g) {
        my ($n, $op) = ($1, $2);
        if ($op eq "I") { $query += $n; next; }
        $length += $n;
        $tail .= $op eq "D" ? "N" x $n : substr($seq, $query, $n);
        $query += $n unless $op eq "D";
        $tail = substr($tail, -3) if length($tail) > 3;
    }
    $query == length($seq) && $length >= 3 or die "CIGAR/sequence length mismatch or fragment too short\n";
    return ($length, join "\t", split //, $tail);
}

@ARGV == 1 or die "Usage: count.pl five_column_alignments.tsv\n";

open( my $IN_SAM, "<", $ARGV[0] );

my %site_end;
my %site_start;
my %site_base;

while (<$IN_SAM>) {
    chomp;
    my ( undef, $rname, $site, $cigar, $seq ) = split /\t/;
    defined($seq) && $site =~ /^[1-9][0-9]*$/ or die "Invalid five-column alignment at line $.\n";
    my ($length, $tail) = ALIGNMENT_TAIL($cigar, $seq);
    my $end_pos = $site + $length - 2;
    my $end     = $rname . "\t" . $end_pos;

    if ( exists $site_end{$end} ) {
        $site_end{$end}++;
    }
    else {
        $site_end{$end} = 1;
        $site_base{$end} =
          $tail;
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
