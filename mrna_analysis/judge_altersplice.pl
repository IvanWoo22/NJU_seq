#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

sub ALTER_SPLICE {
    my ( $GENE, $TRANSCRIPT, $EXON ) = @_;
    my @INTRON;
    foreach ( 0 .. $#{$TRANSCRIPT} ) {
        $INTRON[$_] =
          ${$TRANSCRIPT}[$_]->AlignDB::IntSpan::diff( ${$EXON}[$_] );
    }
    my $ALL_EXON   = AlignDB::IntSpan::intersect( @{$EXON} );
    my $ALL_INTRON = AlignDB::IntSpan::intersect(@INTRON);
    my $ALTER      = $GENE->AlignDB::IntSpan::diff($ALL_EXON)
      ->AlignDB::IntSpan::diff($ALL_INTRON);
    return ( $ALL_EXON->as_string, $ALL_INTRON->as_string, $ALTER->as_string );
}

open( my $ALTER,  ">", $ARGV[0] );
open( my $UNIQUE, ">", $ARGV[1] );

my ( @temp_exon, @temp_transcript, $temp_gene, $temp_gene_info );
while (<STDIN>) {
    chomp;
    my ( $chr, undef, $type, $start, $end, undef, $dir, undef, $info ) =
      split /\t/;
    if ( $type eq "gene" ) {
        if ( @temp_transcript > 1 ) {
            my ( $exon, $intron, $alter ) =
              ALTER_SPLICE( $temp_gene, \@temp_transcript, \@temp_exon );
            print $ALTER ("$temp_gene_info\t$exon\t$intron\t$alter\n");
        }
        elsif ( @temp_transcript == 1 ) {
            my $exon   = $temp_exon[0]->AlignDB::IntSpan::as_string;
            my $intron = $temp_gene->AlignDB::IntSpan::diff( $temp_exon[0] );
            print $UNIQUE ("$temp_gene_info\t$exon\t$intron\t-\n");
        }
        splice( @temp_transcript, 0 );
        splice( @temp_exon,       0 );
        $temp_gene = AlignDB::IntSpan->new;
        $temp_gene->add_range( $start, $end );
        $temp_gene_info = "$chr\t$info\t$dir";
    }
    elsif ( $type eq "mRNA" ) {
        my $set1 = AlignDB::IntSpan->new;
        $set1->add_range( $start, $end );
        push( @temp_transcript, $set1 );
        my $set2 = AlignDB::IntSpan->new;
        push( @temp_exon, $set2 );
    }
    elsif ( $type eq "exon" ) {
        $temp_exon[-1]->AlignDB::IntSpan::add_range( $start, $end );
    }
}
close($ALTER);
close($UNIQUE);

__END__
