#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use AlignDB::IntSpan;

Getopt::Long::GetOptions(
    'help|h'         => sub { Getopt::Long::HelpMessage(0) },
    'transwording=s' => \my $mRNA,
    'geneid=s'       => \my $geneid,
    'transid=s'      => \my $transid,
    'ymlinput=s'     => \my $IN_FH,
) or Getopt::Long::HelpMessage(1);

sub GET_INFO {
    my $INFO = shift;
    my $GENE_ID;
    if ( $INFO =~ m/$geneid(\w+)\.?[0-9]*;/ ) {
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
        my $name = GET_INFO($info);
        push( @gene, $i );
        $gene_info{$name} = $#gene;
    }
    $i++;
}

open( my $GENE, "<", $IN_FH );
while (<$GENE>) {
    chomp;
    my ( $chr, $name, $dir, $exon, $intron, undef ) =
      split /\t/;
    my $transcript;
    my $trans_cds       = AlignDB::IntSpan->new;
    my $trans_five_utr  = AlignDB::IntSpan->new;
    my $trans_three_utr = AlignDB::IntSpan->new;
    foreach
      my $t ( $gene[ $gene_info{$name} ] .. $gene[ $gene_info{$name} + 1 ] - 1 )
    {
        my ( undef, undef, $type, $start, $end, undef, undef, undef,
            $line_info ) = split( /\t/, $line[$t] );
        if ( $type eq $mRNA ) {
            if ( $line_info =~ /$transid([0-9,a-z,A-Z,\.,\_]+);/ ) {
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
    print(
"$chr\t$name\t$dir\t$exon\t$intron\t$transcript\t$trans_five_utr_string\t$trans_cds_string\t$trans_three_utr_string\n"
    );
}
close($GENE);

__END__
