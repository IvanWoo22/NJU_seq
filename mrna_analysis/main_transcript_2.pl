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

my ( %rep_trans, %rep_gene, %chr, %dir, %cds, %utr5, %utr3, %no );

open( my $IN_FH, "<", $IN_TRANS );
while (<$IN_FH>) {
    chomp;
    my ( $trans, $gene, $no ) = split "\t";
    $rep_trans{$trans} = 1;
    $rep_gene{$trans}  = $gene;
    $no{$trans}        = $no;
}
close($IN_FH);

while (<STDIN>) {
    chomp;
    my ( $chr, undef, $type, $start, $end, undef, $dir, undef, $info ) =
      split /\t/;
    my $name = GET_ID( $transid, $info );
    if ( exists( $rep_trans{$name} ) ) {
        $chr{$name} = $chr unless exists( $chr{$name} );
        $dir{$name} = $dir unless exists( $dir{$name} );
        if ( $type eq "CDS" ) {
            if ( exists( $cds{$name} ) ) {
                $cds{$name}->AlignDB::IntSpan::add_range( $start, $end );
            }
            else {
                $cds{$name} = AlignDB::IntSpan->new;
                $cds{$name}->AlignDB::IntSpan::add_range( $start, $end );
            }
        }
        elsif ( $type eq "five_prime_UTR" ) {
            if ( exists( $utr5{$name} ) ) {
                $utr5{$name}->AlignDB::IntSpan::add_range( $start, $end );
            }
            else {
                $utr5{$name} = AlignDB::IntSpan->new;
                $utr5{$name}->AlignDB::IntSpan::add_range( $start, $end );
            }
        }
        elsif ( $type eq "three_prime_UTR" ) {
            if ( exists( $utr3{$name} ) ) {
                $utr3{$name}->AlignDB::IntSpan::add_range( $start, $end );
            }
            else {
                $utr3{$name} = AlignDB::IntSpan->new;
                $utr3{$name}->AlignDB::IntSpan::add_range( $start, $end );
            }
        }
    }
}

foreach my $trans ( keys(%rep_trans) ) {
    my ( $trans_five_utr_string, $trans_cds_string, $trans_three_utr_string );
    if ( exists( $utr5{$trans} ) ) {
        $trans_five_utr_string = $utr5{$trans}->AlignDB::IntSpan::as_string;
    }
    else {
        $trans_five_utr_string = "-";
    }
    if ( exists( $cds{$trans} ) ) {
        $trans_cds_string = $cds{$trans}->AlignDB::IntSpan::as_string;
    }
    else {
        $trans_cds_string = "-";
    }
    if ( exists( $utr3{$trans} ) ) {
        $trans_three_utr_string = $utr3{$trans}->AlignDB::IntSpan::as_string;
    }
    else {
        $trans_three_utr_string = "-";
    }
    print(
"$chr{$trans}\t$trans\t$rep_gene{$trans}\t$no{$trans}\t$dir{$trans}\t$trans_five_utr_string\t$trans_cds_string\t$trans_three_utr_string\n"
    );
}

__END__
