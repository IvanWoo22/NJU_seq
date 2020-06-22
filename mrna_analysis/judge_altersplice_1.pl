#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use AlignDB::IntSpan;

Getopt::Long::GetOptions(
    'help|h'      => sub { Getopt::Long::HelpMessage(0) },
    'geneid=s'    => \my $gene_id,
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

my (
    %trans,    %gene_info,  %chr,  %dir,
    %exon,     %intron,     %mrna, %gene_range,
    %all_exon, %all_intron, %alter_splice
);

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
    my $trans_name = GET_ID( $transid, $info );
    my $gene_name  = GET_ID( $gene_id, $info );
    if ( exists( $trans{$trans_name} ) ) {
        $chr{$gene_name} = $chr unless exists( $chr{$gene_name} );
        $dir{$gene_name} = $dir unless exists( $dir{$gene_name} );
        if ( $type eq "transcript" ) {
            $mrna{$trans_name} = AlignDB::IntSpan->new;
            $mrna{$trans_name}->AlignDB::IntSpan::add_range( $start, $end );
        }
        elsif ( $type eq "exon" ) {
            if ( exists( $exon{$trans_name} ) ) {
                $exon{$trans_name}->AlignDB::IntSpan::add_range( $start, $end );
            }
            else {
                $exon{$trans_name} = AlignDB::IntSpan->new;
                $exon{$trans_name}->AlignDB::IntSpan::add_range( $start, $end );
            }
        }
    }
}

foreach my $gene ( sort( keys(%gene_info) ) ) {
    $gene_range{$gene}   = AlignDB::IntSpan->new;
    $all_exon{$gene}     = AlignDB::IntSpan->new;
    $all_intron{$gene}   = AlignDB::IntSpan->new;
    $alter_splice{$gene} = AlignDB::IntSpan->new;
    foreach my $trans ( sort( keys( %{ $gene_info{$gene} } ) ) ) {
        $gene_range{$gene}->AlignDB::IntSpan::add( $mrna{$trans} );
        $intron{$trans} =
          $mrna{$trans}->AlignDB::IntSpan::diff( $exon{$trans} );
    }
    foreach my $point ( $gene_range{$gene}->AlignDB::IntSpan::as_array ) {
        my $exon_count   = 0;
        my $intron_count = 0;
        foreach my $trans ( sort( keys( %{ $gene_info{$gene} } ) ) ) {
            $exon_count++
              if $exon{$trans}->AlignDB::IntSpan::contains_all($point);
            $intron_count++
              if $intron{$trans}->AlignDB::IntSpan::contains_all($point);
        }
        $all_exon{$gene}->AlignDB::IntSpan::add($point)
          if ( $exon_count > 0 and $intron_count == 0 );
        $all_intron{$gene}->AlignDB::IntSpan::add($point)
          if ( $exon_count == 0 and $intron_count > 0 );
        $alter_splice{$gene}->AlignDB::IntSpan::add($point)
          if ( $exon_count > 0 and $intron_count > 0 );
    }
    my $all_exon     = $all_exon{$gene}->AlignDB::IntSpan::as_string;
    my $all_intron   = $all_intron{$gene}->AlignDB::IntSpan::as_string;
    my $alter_splice = $alter_splice{$gene}->AlignDB::IntSpan::as_string;
    print(
"$chr{$gene}\t$gene\t$dir{$gene}\t$all_exon\t$all_intron\t$alter_splice\n"
    );
}

__END__
