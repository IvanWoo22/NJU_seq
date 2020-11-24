#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use AlignDB::IntSpan;
use POSIX;

sub DISTANCE {
    my ( $SITE, $SET_STRING ) = @_;
    my $SET = AlignDB::IntSpan->new;
    $SET->add($SET_STRING);
    my ( $ISLAND, $DIS );
    if ( $SET->find_islands($SITE) eq "-" ) {
        $ISLAND = $SET->nearest_island($SITE)->as_string;
        if ( $ISLAND =~ /,/ ) {
            ( $ISLAND, undef ) = split( ",", $ISLAND );
        }
        if ( $ISLAND =~ /\-/ ) {
            my ( $START, $END ) = split( "-", $ISLAND );
            if ( $START > $SITE ) {
                $DIS = $SITE - $START;
            }
            else {
                $DIS = $SITE - $END;
            }
        }
        else {
            $DIS = $SITE - $ISLAND;
        }
    }
    else {
        $DIS = 0;
    }
    return $DIS;
}

sub PORTA {
    my $LENGTH = shift;
    $LENGTH = POSIX::ceil( $LENGTH / 2 );
    my $PORTA = shift;
    my $CELL  = POSIX::floor( $LENGTH / 10 );
    foreach my $K ( 0 .. $CELL - 1 ) {
        ${$PORTA}[$K] += 5 / $LENGTH;
    }
    ${$PORTA}[$CELL] += ( $LENGTH % 10 ) / ( 2 * $LENGTH );
    return @{$PORTA};
}

my ( @porta, @near5, @near3 );
foreach ( 0 .. 30 ) {
    $porta[$_] = 0;
    $near5[$_] = 0;
    $near3[$_] = 0;
}

my ( %gene, %exon, %exon_site );
open( my $IN1, "<", $ARGV[0] );
while (<$IN1>) {
    chomp;
    my ( $trans, $gene_tmp, $no, $exon, $exon_site ) = split "\t";
    my ( $gene, undef ) = split( /\./, $gene_tmp );
    $gene{$gene}[ $no - 1 ] = $trans;
    $exon{$trans} = AlignDB::IntSpan->new;
    $exon{$trans}->AlignDB::IntSpan::add($exon);
    $exon_site{$trans} = AlignDB::IntSpan->new;
    $exon_site{$trans}->AlignDB::IntSpan::add($exon_site);
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
            if ( $exon{ $gene{$gene}[$no] }
                ->AlignDB::IntSpan::contains_all($site) )
            {
                my $distance = DISTANCE( $site,
                    $exon_site{ $gene{$gene}[$no] }
                      ->AlignDB::IntSpan::as_string );
                my $size =
                  $exon{ $gene{$gene}[$no] }
                  ->AlignDB::IntSpan::find_islands($site)->size;
                if ( $strand eq "-" ) {
                    $distance = -$distance;
                }
                if ( $distance < 0 ) {
                    $near3[ POSIX::ceil( abs($distance) / 10 ) ] += $weight
                      if POSIX::ceil( abs($distance) / 10 ) < 31;
                }
                else {
                    $near5[ POSIX::ceil( abs($distance) / 10 ) ] += $weight
                      if POSIX::ceil( abs($distance) / 10 ) < 31;
                }
                print
"$site\t$gene\t$gene{$gene}[$no]\t$exon{$gene{$gene}[$no]}\t$strand\t$exon_site{$gene{$gene}[$no]}\t$distance\t$size\t$weight\n";
                @porta = PORTA( $size, \@porta );
                last;
            }
        }
    }
}
close($IN2);

open( my $BAR, ">", $ARGV[2] );
foreach ( 1 .. 30 ) {
    my $x = -$_ * 10 + 5;
    print $BAR "near3\t$x\t$near3[$_]\n";
}
foreach ( 1 .. 30 ) {
    my $x = $_ * 10 - 5;
    print $BAR "near5\t$x\t$near5[$_]\n";
}

open( my $POR, ">", $ARGV[3] );
foreach ( 0 .. 29 ) {
    my $x = $_ * 10 + 5;
    print $POR "$x\t$porta[$_]\n";
}

__END__
