#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use AlignDB::IntSpan;

my ( %five_utr_set, %cds_set, %three_utr_set, %five_utr_length, %cds_length,
    %three_utr_length, %trans, %chr_dir );
open( my $GENE, "<", $ARGV[0] );
while (<$GENE>) {
    chomp;
    my ( $chr, $trans, $gene_tmp, $no, $dir,
        $five_utr_string, $cds_string, $three_utr_string )
      = split /\t/;
    $chr =~ s/chr//;
    my ( $gene, undef ) = split( /\./, $gene_tmp );
    ${ $trans{$gene} }[ $no - 1 ] = $trans;
    $chr_dir{$trans}       = $chr . "\t" . $dir;
    $five_utr_set{$trans}  = AlignDB::IntSpan->new($five_utr_string);
    $cds_set{$trans}       = AlignDB::IntSpan->new($cds_string);
    $three_utr_set{$trans} = AlignDB::IntSpan->new($three_utr_string);
    $five_utr_length{$trans} =
      $five_utr_set{$trans}->AlignDB::IntSpan::cardinality;
    $cds_length{$trans} = $cds_set{$trans}->AlignDB::IntSpan::cardinality;
    $three_utr_length{$trans} =
      $three_utr_set{$trans}->AlignDB::IntSpan::cardinality;
}
close($GENE);

open( my $POINT, "<", $ARGV[1] );
open( my $OUT, ">", $ARGV[2]);
while (<$POINT>) {
    chomp;
    my $line = $_;
    my %type;
    my ( undef, $site, $dir, undef, undef, undef, $gene_list, undef ) =
      split "\t";
    my @gene   = split /\//, $gene_list;
    my $weight = 1 / @gene;
    if ( $dir eq "+" ) {
        foreach my $gene (@gene) {
            foreach my $trans ( @{ $trans{$gene} } ) {
                if (
                    $five_utr_set{$trans}->AlignDB::IntSpan::contains_all($site)
                  )
                {
                    my $index =
                      $five_utr_set{$trans}->AlignDB::IntSpan::index($site);
                    my $distribution = $index / $five_utr_length{$trans};
                    print(
"$chr_dir{$trans}\t$site\tfive_utr\t$distribution\t$five_utr_length{$trans}\t$cds_length{$trans}\t$three_utr_length{$trans}\t$index\t$weight\n"
                    );
                    $type{"five_utr"} = 1;
                    last;
                }
                elsif (
                    $cds_set{$trans}->AlignDB::IntSpan::contains_all($site) )
                {
                    my $index =
                      $cds_set{$trans}->AlignDB::IntSpan::index($site);
                    my $distribution = $index / $cds_length{$trans};
                    print(
"$chr_dir{$trans}\t$site\tcds\t$distribution\t$five_utr_length{$trans}\t$cds_length{$trans}\t$three_utr_length{$trans}\t$index\t$weight\n"
                    );
                    $type{"cds"} = 1;
                    last;
                }
                elsif ( $three_utr_set{$trans}
                    ->AlignDB::IntSpan::contains_all($site) )
                {
                    my $index =
                      $three_utr_set{$trans}->AlignDB::IntSpan::index($site);
                    my $distribution = $index / $three_utr_length{$trans};
                    print(
"$chr_dir{$trans}\t$site\tthree_utr\t$distribution\t$five_utr_length{$trans}\t$cds_length{$trans}\t$three_utr_length{$trans}\t$index\t$weight\n"
                    );
                    $type{"three_utr"} = 1;
                    last;
                }
            }
        }
    }
    else {
        foreach my $gene (@gene) {
            foreach my $trans ( @{ $trans{$gene} } ) {
                if (
                    $five_utr_set{$trans}->AlignDB::IntSpan::contains_all($site)
                  )
                {
                    my $index =
                      $five_utr_set{$trans}->AlignDB::IntSpan::index($site);
                    my $distribution = 1 - $index / $five_utr_length{$trans};
                    print(
"$chr_dir{$trans}\t$site\tfive_utr\t$distribution\t$five_utr_length{$trans}\t$cds_length{$trans}\t$three_utr_length{$trans}\t$index\t$weight\n"
                    );
                    $type{"five_utr"} = 1;
                    last;
                }
                elsif (
                    $cds_set{$trans}->AlignDB::IntSpan::contains_all($site) )
                {
                    my $index =
                      $cds_set{$trans}->AlignDB::IntSpan::index($site);
                    my $distribution = 1 - $index / $cds_length{$trans};
                    print(
"$chr_dir{$trans}\t$site\tcds\t$distribution\t$five_utr_length{$trans}\t$cds_length{$trans}\t$three_utr_length{$trans}\t$index\t$weight\n"
                    );
                    $type{"cds"} = 1;
                    last;
                }
                elsif ( $three_utr_set{$trans}
                    ->AlignDB::IntSpan::contains_all($site) )
                {
                    my $index =
                      $three_utr_set{$trans}->AlignDB::IntSpan::index($site);
                    my $distribution = 1 - $index / $three_utr_length{$trans};
                    print(
"$chr_dir{$trans}\t$site\tthree_utr\t$distribution\t$five_utr_length{$trans}\t$cds_length{$trans}\t$three_utr_length{$trans}\t$index\t$weight\n"
                    );
                    $type{"three_utr"} = 1;
                    last;
                }
            }
        }
    }
    my ($fu,$cds,$tu) = (0,0,0);

    $fu = "five_utr" if exists($type{"five_utr"});
    $cds = "cds" if exists($type{"cds"});
    $tu = "three_utr" if exists($type{"three_utr"});
    print $OUT "$line\t$fu\t$cds\t$tu\n";
}
close($POINT);
close($OUT);

__END__
