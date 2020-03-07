#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

sub GET_INFO {
    my $INFO = shift;
    my $GENE_ID;
    if ( $INFO =~ m/ID=gene:([A-Z,a-z,0-9]+);/ ) {
        $GENE_ID = $1;
    }
    else {
        warn("There is a problem in $INFO;\n");
    }
    return $GENE_ID;
}

my ( @line, @gene, %gene_info );
my $i = 0;
while (<STDIN>) {
    chomp;
    push( @line, $_ );
    my ( undef, undef, $type, undef, undef, undef, undef, undef, $info ) =
      split /\t/;
    if ( $type eq "gene" ) {
        push( @gene, $i );
        $gene_info{$info} = $#gene;
    }
    $i++;
}

open( my $GENE, "<", $ARGV[0] );
while (<$GENE>) {
    chomp;
    my ( $chr, $info, $dir, $exon, $intron, undef ) =
      split /\t/;
    my $transcript;
    my $trans_cds       = AlignDB::IntSpan->new;
    my $trans_five_utr  = AlignDB::IntSpan->new;
    my $trans_three_utr = AlignDB::IntSpan->new;
    foreach
      my $t ( $gene[ $gene_info{$info} ] .. $gene[ $gene_info{$info} + 1 ] - 1 )
    {
        my ( undef, undef, $type, $start, $end, undef, undef, undef,
            $line_info ) = split( /\t/, $line[$t] );
        if ( $type eq "mRNA" ) {
            if ( $line_info =~ /^ID=transcript:([0-9,a-z,A-Z,\.,\_]+);/ ) {
                $transcript = $1;
            }
        }
        elsif ( $type eq "CDS" ) {
            $trans_cds->add_range( $start, $end );
        }
        elsif ( $type eq "five_prime_UTR" ) {
            $trans_five_utr->add_range( $start, $end );
        }
        elsif ( $type eq "three_prime_UTR" ) {
            $trans_three_utr->add_range( $start, $end );
        }
    }
    my $trans_five_utr_string  = $trans_five_utr->as_string;
    my $trans_cds_string       = $trans_cds->as_string;
    my $trans_three_utr_string = $trans_three_utr->as_string;
    my $name = GET_INFO($info);
    print(
"$chr\t$name\t$dir\t$exon\t$intron\t$transcript\t$trans_five_utr_string\t$trans_cds_string\t$trans_three_utr_string\n"
    );
}
close($GENE);

__END__
