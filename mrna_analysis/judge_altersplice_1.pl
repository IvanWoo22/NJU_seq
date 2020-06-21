#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use AlignDB::IntSpan;

Getopt::Long::GetOptions(
    'help|h'      => sub { Getopt::Long::HelpMessage(0) },
    'transid=s'   => \my $transid,
    'rep_trans=s' => \my $IN_TRANS,
) or Getopt::Long::HelpMessage(1);

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

my ( %trans, %gene_info, %chr, %dir, %exon, %intron, %mrna );

open( my $IN_FH, "<", $IN_TRANS );
while (<$IN_FH>) {
    chomp;
    my ( $trans, $gene ) = split "\t";
    $trans{$trans} = 1;
    $gene_info{$gene}{$trans} = 1;
}
close($IN_FH);

foreach my $gene ( sort( keys(%gene_info) ) ) {
    if ( keys( %{ $gene_info{$gene} } ) == 1 ) {
        foreach my $trans ( sort( keys( %{ $gene_info{$gene} } ) ) ) {
            delete( $trans{$trans} );
        }
        delete( $gene_info{$gene} );
    }
}

while (<STDIN>) {
    chomp;
    my ( $chr, undef, $type, $start, $end, undef, $dir, undef, $info ) =
      split /\t/;
    my $name = GET_ID( $transid, $info );
    if ( exists( $trans{$name} ) ) {
        $chr{$name} = $chr unless exists( $chr{$name} );
        $dir{$name} = $dir unless exists( $dir{$name} );
        if ( $type eq "transcript" ) {

            $mrna{$name} = AlignDB::IntSpan->new;
            $mrna{$name}->AlignDB::IntSpan::add_range( $start, $end );
        }
        elsif ( $type eq "exon" ) {
            if ( exists( $exon{$name} ) ) {
                $exon{$name}->AlignDB::IntSpan::add_range( $start, $end );
            }
            else {
                $exon{$name} = AlignDB::IntSpan->new;
                $exon{$name}->AlignDB::IntSpan::add_range( $start, $end );
            }
        }
    }
}
