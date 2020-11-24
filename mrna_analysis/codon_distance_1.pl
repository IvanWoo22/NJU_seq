#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

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
    $chr =~ s/chr//;
    my $trans = GET_ID( "transcript_id=", $info );
    if ( exists( $trans{$trans} ) ) {
        print "$trans\t$trans{$trans}\t$chr\t$strand\t$start\t$end\n";
    }
}

__END__
