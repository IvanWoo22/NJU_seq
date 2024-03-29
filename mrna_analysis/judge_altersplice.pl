#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use AlignDB::IntSpan;

Getopt::Long::GetOptions(
    'help|h'         => sub { Getopt::Long::HelpMessage(0) },
    'transwording=s' => \my $mRNA,
    'geneid=s'       => \my $gene_id,
    'alter=s'        => \my $ALTER_FH,
    'unique=s'       => \my $UNIQUE_FH,
) or Getopt::Long::HelpMessage(1);

sub GET_INFO {
    my $INFO = shift;
    my $GENE_ID;
    if ( $INFO =~ m/$gene_id(\w+)\.?[0-9]*;/ ) {
        $GENE_ID = $1;
    }
    else {
        warn("There is a problem in $INFO;\n");
    }
    return $GENE_ID;
}

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

open( my $ALTER,  ">", $ALTER_FH );
open( my $UNIQUE, ">", $UNIQUE_FH );

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
        my $gene_info = GET_INFO($info);
        $temp_gene_info = "$chr\t$gene_info\t$dir";
    }
    elsif ( $type eq $mRNA ) {
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
