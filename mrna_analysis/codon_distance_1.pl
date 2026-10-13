#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

sub GET_ID {
    my ( $FEATURE, $INFO ) = @_;
    my $ID;
    if ( $INFO =~ m/(?:transcript_id=|(?:Parent|ID)=transcript:)([^;\s,]+)/ ) {
        $ID = $1;
    }
    else {
        warn("There is a problem in $INFO;\n");
    }
    return $ID;
}

my %trans;
open( my $IN_FH, "<", $ARGV[0] );
while (<$IN_FH>) {
    chomp;
    my ( $trans, $gene, $no ) = split "\t";
    $trans{$trans} = $gene . "\t" . $no;
}
close($IN_FH);

while (<STDIN>) {
    chomp;
    my ( $chr, $start, $end, $strand, $info ) = split /\t/;
    defined($info) && $start =~ /^[1-9][0-9]*$/ && $end >= $start && $strand =~ /^[+-]$/
        or die "Expected five-column codon regions: chr,start,end,strand,attributes\n";
    $chr =~ s/chr//;
    my $trans = GET_ID( "transcript_id=", $info );
    if ( defined($trans) && exists( $trans{$trans} ) ) {
        print "$trans\t$trans{$trans}\t$chr\t$strand\t$start\t$end\n";
    }
}

__END__
