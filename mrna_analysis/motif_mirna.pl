#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my %fasta;
my $title_name;
while (<>) {
    if (m/^>(\S+)/) {
        $title_name = $1;
    }
    else {
        $_ =~ s/\r?\n//;
        $fasta{$title_name} .= $_;
    }
}

foreach ( keys(%fasta) ) {
    my @seq   = split( //, $fasta{$_} );
    my $motif = join( "", @seq[ -8 .. -1 ] );
    print "$_\t$motif\n";
}

__END__
