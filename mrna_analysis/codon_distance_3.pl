use strict;
use warnings;
use autodie;
use AlignDB::IntSpan;
use POSIX;

sub DISTANCE {
    my ( $SITE, $SET_STRING ) = @_;
    my $SET = AlignDB::IntSpan->new;
    $SET->add($SET_STRING);
    my ( $ISLAND, $DIS, $DIR );
    if ( $SET->find_islands($SITE) eq "-" ) {
        $ISLAND = $SET->nearest_island($SITE);
        if ( $ISLAND =~ /,/ ) {
            ( $ISLAND, undef ) = split( ",", $ISLAND );
        }
        if ( $ISLAND =~ /-/ ) {
            my ( $START, $END ) = split( "-", $ISLAND );
            if ( $START > $SITE ) {
                $DIS = "$SITE-$START";
                $DIR = -1;
            }
            else {
                $DIS = "$END-$SITE";
                $DIR = 1;
            }
        }
        else {
            if ( $SITE > $ISLAND->as_string ) {
                $DIS = "$ISLAND-$SITE";
                $DIR = 1;
            }
            else {
                $DIS = "$SITE-$ISLAND";
                $DIR = -1;
            }
        }
    }
    else {
        $DIS = 0;
        $DIR = 0;
    }
    return ( $DIS, $DIR );
}

my @window;
foreach ( 0 .. 20 ) {
    $window[$_] = 0;
}

my ( %gene, %exon, %intron, %codon );
open( my $IN1, "<", $ARGV[0] );
while (<$IN1>) {
    chomp;
    my ( $trans, $gene_tmp, $no, undef, $exon, $intron, $codon ) =
      split "\t";
    my ( $gene, undef ) = split( /\./, $gene_tmp );
    $gene{$gene}[ $no - 1 ] = $trans;
    $exon{$trans} = AlignDB::IntSpan->new;
    $exon{$trans}->AlignDB::IntSpan::add($exon);
    $intron{$trans} = AlignDB::IntSpan->new;
    $intron{$trans}->AlignDB::IntSpan::add($intron);
    $codon{$trans} = AlignDB::IntSpan->new;
    $codon{$trans}->AlignDB::IntSpan::add_runlist($codon);
}
close($IN1);

open( my $IN2, "<", $ARGV[1] );
while (<$IN2>) {
    chomp;
    my ( undef, $site, $strand, undef, undef, undef, $gene_list, undef ) =
      split "\t";
    my @gene   = split /\//, $gene_list;
    my $weight = 1 / @gene;
    foreach my $gene (@gene) {
        foreach my $no ( 0 .. $#{ $gene{$gene} } ) {
            if (
                ( defined( $gene{$gene}[$no] ) )
                and ( $exon{ $gene{$gene}[$no] }
                    ->AlignDB::IntSpan::contains_all($site) )
              )
            {
                my $distance;
                my ( $distance_neo, $direction ) =
                  DISTANCE( $site, $codon{ $gene{$gene}[$no] } );
                if ( $distance_neo ne "0" ) {
                    my $set1 = AlignDB::IntSpan->new;
                    $set1->add($distance_neo);
                    my $set2 = $set1->diff( $intron{ $gene{$gene}[$no] } );
                    if ( $strand eq "+" ) {
                        $distance = ( $set2->cardinality - 1 ) * $direction * 1;
                    }
                    else {
                        $distance =
                          ( $set2->cardinality - 1 ) * $direction * (-1);
                    }
                }
                else {
                    $distance = 0;
                }

                my ( $left_length, $site_1, $right_length, $site_2 );
                if ( $strand eq "+" ) {
                    $site_1 = $exon{ $gene{$gene}[$no] }->AlignDB::IntSpan::min;
                    $site_2 = $exon{ $gene{$gene}[$no] }->AlignDB::IntSpan::max;
                }
                else {
                    $site_1 = $exon{ $gene{$gene}[$no] }->AlignDB::IntSpan::max;
                    $site_2 = $exon{ $gene{$gene}[$no] }->AlignDB::IntSpan::min;
                }
                my ( $distance_neo_1, $direction_1 ) =
                  DISTANCE( $site_1, $codon{ $gene{$gene}[$no] } );
                my ( $distance_neo_2, $direction_2 ) =
                  DISTANCE( $site_2, $codon{ $gene{$gene}[$no] } );
                if ( $distance_neo_1 ne "0" ) {
                    my $set1 = AlignDB::IntSpan->new;
                    $set1->add($distance_neo_1);
                    my $set2 = $set1->diff( $intron{ $gene{$gene}[$no] } );
                    if ( $strand eq "+" ) {
                        $left_length =
                          ( $set2->cardinality - 1 ) * $direction_1 * 1;
                    }
                    else {
                        $left_length =
                          ( $set2->cardinality - 1 ) * $direction_1 * (-1);
                    }
                }
                else {
                    $left_length = 0;
                }
                if ( $distance_neo_2 ne "0" ) {
                    my $set1 = AlignDB::IntSpan->new;
                    $set1->add($distance_neo_2);
                    my $set2 = $set1->diff( $intron{ $gene{$gene}[$no] } );
                    if ( $strand eq "+" ) {
                        $right_length =
                          ( $set2->cardinality - 1 ) * $direction_2 * 1;
                    }
                    else {
                        $right_length =
                          ( $set2->cardinality - 1 ) * $direction_2 * (-1);
                    }
                }
                else {
                    $right_length = 0;
                }

                if ( $distance_neo ne "0" ) {
                    my $set1 = AlignDB::IntSpan->new;
                    $set1->add($distance_neo);
                    my $set2 = $set1->diff( $intron{ $gene{$gene}[$no] } );
                    if ( $strand eq "+" ) {
                        $distance = ( $set2->cardinality - 1 ) * $direction * 1;
                    }
                    else {
                        $distance =
                          ( $set2->cardinality - 1 ) * $direction * (-1);
                    }
                }
                else {
                    $distance = 0;
                }

                if (    ( $distance < 500 )
                    and ( $distance >= -500 )
                    and ( abs($left_length) > $ARGV[3] )
                    and ( abs($right_length) > $ARGV[4] ) )
                {
                    $window[ POSIX::floor( ( $distance + 500 ) / 50 ) ] +=
                      $weight;
                }
                print
"$site\t$gene\t$gene{$gene}[$no]\t$exon{$gene{$gene}[$no]}\t$codon{$gene{$gene}[$no]}\t$distance\t$weight\n";
                last;
            }
        }
    }
}
close($IN2);

open( my $BAR, ">", $ARGV[2] );
foreach ( 0 .. 19 ) {
    my $x = $_ * 50 - 475;
    print $BAR "$x\t$window[$_]\n";
}

__END__
