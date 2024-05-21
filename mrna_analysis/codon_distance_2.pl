#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use AlignDB::IntSpan;

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

my ( %trans, %gene, %no, %codon, %chr_dir, %exon, %intron, %mrna );

open( my $IN_FH1, "<", $ARGV[0] );
while (<$IN_FH1>) {
    chomp;
    my ( $trans, $gene, $no, $chr, $dir, $start, $end ) = split "\t";
    $trans{$trans} = 1;
    $gene{$trans}  = $gene;
    unless ( exists( $codon{$trans} ) ) {
        $codon{$trans} = AlignDB::IntSpan->new;
    }
    $codon{$trans}->AlignDB::IntSpan::add_range( $start, $end );
    $chr_dir{$trans} = "$chr\t$dir";
    $no{$trans}      = $no;
}
close($IN_FH1);

while (<STDIN>) {
    chomp;
    my ( undef, undef, $type, $start, $end, undef, undef, undef, $info ) =
      split /\t/;
    my $trans_id = GET_ID( "transcript_id=", $info );
    if ( exists( $trans{$trans_id} ) ) {
        if ( $type eq "transcript" ) {
            $mrna{$trans_id} = AlignDB::IntSpan->new;
            $mrna{$trans_id}->AlignDB::IntSpan::add_range( $start, $end );
        }
        elsif ( $type eq "exon" ) {
            if ( exists( $exon{$trans_id} ) ) {
                $exon{$trans_id}->AlignDB::IntSpan::add_range( $start, $end );
            }
            else {
                $exon{$trans_id} = AlignDB::IntSpan->new;
                $exon{$trans_id}->AlignDB::IntSpan::add_range( $start, $end );
            }
        }
    }
}

foreach my $trans ( keys(%trans) ) {
    $intron{$trans} = $mrna{$trans}->AlignDB::IntSpan::diff( $exon{$trans} );
    print
"$trans\t$gene{$trans}\t$no{$trans}\t$mrna{$trans}\t$exon{$trans}\t$intron{$trans}\t$codon{$trans}\n";
}

__END__
