#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Algorithm::Combinatorics qw(combinations);

my ( %point, %point_info1, %point_info2, @all_set, %combine, %all_point );
foreach my $sample (@ARGV) {
    open( my $FH, "<", $sample );
    while (<$FH>) {
        chomp;
        my @tmp = split( "\t", $_ );
        push( @{ $point{$sample} }, $tmp[0] . $tmp[2] . $tmp[1] );
        $point_info1{ $tmp[0] . $tmp[2] . $tmp[1] } =
          join( "\t", @tmp[ 0 .. 8 ] );
        $point_info2{ $tmp[0] . $tmp[2] . $tmp[1] } =
          join( "\t", @tmp[ 20 .. 25 ] );
        @{ ${ $combine{$sample} }{ $tmp[0] . $tmp[2] . $tmp[1] } } =
          @tmp[ 9 .. 19 ];
    }
    push( @all_set, @{ $point{$sample} } );
}

foreach my $ELEMENT (@all_set) {
    $all_point{$ELEMENT}++;
}
my @all_point = keys %all_point;

foreach my $number ( 1 .. $#ARGV + 1 ) {
    my $iter = combinations( \@ARGV, $number );
    while ( my $choose = $iter->next ) {
        my ( @in, @out, %choose, @in_point, @out_point, %in_point, %out_point );
        foreach my $element ( @{$choose} ) {
            $choose{$element} = 1;
        }
        foreach my $sample ( 0 .. $#ARGV ) {
            my $file_name = $ARGV[$sample];
            if ( exists( $choose{$file_name} ) ) {
                push( @in, @{ $point{$file_name} } );
            }
            else {
                push( @out, @{ $point{$file_name} } );
            }
        }
        foreach my $ELEMENT (@in) {
            $in_point{$ELEMENT}++;
        }
        @in_point = grep { $in_point{$_} == $number; } ( keys %in_point );
        foreach my $ELEMENT (@out) {
            $out_point{$ELEMENT}++;
        }
        @out_point = keys %out_point;
        my @point =
          grep {
                  ( exists( $in_point{$_} ) )
              and ( $in_point{$_} == $number )
              and ( not( exists( $out_point{$_} ) ) )
          } @all_point;
        if ( @point > 0 ) {
            print "#" . join( ";", @$choose ) . "\n";
            my %output;
            foreach my $point (@point) {
                my @out_col = @{ ${ $combine{ ${$choose}[0] } }{$point} };
                my $score   = ${ ${ $combine{ ${$choose}[0] } }{$point} }[-1];
                foreach my $sample_id ( 1 .. $#{$choose} ) {
                    foreach my $col ( 0 .. 10 ) {
                        $out_col[$col] .=
";${ ${ $combine{${$choose}[$sample_id]} }{ $point } }[$col]";
                    }
                    $score +=
                      ${ ${ $combine{ ${$choose}[$sample_id] } }{$point} }[-1];
                }
                $output{ "$point_info1{$point}\t$point_info2{$point}" . "\t"
                      . join( "\t", @out_col ) } = $score;
            }
            foreach my $point_info (
                sort { $output{$b} <=> $output{$a} }
                keys %output
              )
            {
                print("$point_info\t$output{$point_info}\n");
            }
        }
    }
}

__END__
