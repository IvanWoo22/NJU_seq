#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my $trans;
while (<>) {
    chomp;
    my $line = $_;
    my ( undef, undef, $type, undef, undef, undef, undef, undef, $info ) =
      split( /\t/, $line );
    if ( $type eq "gene" ) {
        print "$line\n";
    }
    elsif ( $type eq "mRNA" ) {
        print "$line\n";
        $info =~ /ID=transcript:([A-Z,0-9,a-z]+\.*[0-9]+);/;
        $trans = $1;
    }
    elsif ( $type eq "exon" ) {
        $info =~ /Parent=transcript:([A-Z,0-9,a-z]+\.*[0-9]+);/;
        if ( $trans eq $1 ) {
            print "$line\n";
        }
    }
}

__END__
