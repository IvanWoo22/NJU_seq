#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long;

Getopt::Long::GetOptions(
    'help|h'         => sub { Getopt::Long::HelpMessage(0) },
    'transwording=s' => \my $mRNA,
    'transid=s'      => \my $trans_id,
    'exonid=s'       => \my $exon_id,
) or Getopt::Long::HelpMessage(1);

my $trans;
while (<STDIN>) {
    chomp;
    my $line = $_;
    my ( undef, undef, $type, undef, undef, undef, undef, undef, $info ) =
      split( /\t/, $line );
    if ( $type eq "gene" ) {
        print "$line\n";
    }
    elsif ( $type eq $mRNA ) {
        print "$line\n";
        $info =~ /$trans_id(\w+\.*[0-9]+);/;
        $trans = $1;
    }
    elsif ( $type eq "exon" ) {
        $info =~ /$exon_id(\w+\.*[0-9]+);/;
        if ( ( defined($trans) ) and ( $trans eq $1 ) ) {
            print "$line\n";
        }
    }
}

__END__
