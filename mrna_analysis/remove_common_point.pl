#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my ( %point, %point_info, @all_set );
foreach my $file_name (@ARGV) {
    open( my $FH, "<", $file_name );
    while (<$FH>) {
        chomp;
        my ( $chr, $pos, $dir ) = split( "\t", $_ );
        push( @{ $point{$file_name} }, $chr . $dir . $pos );
        $point_info{ $chr . $dir . $pos } = $_;
    }
    push( @all_set, @{ $point{$file_name} } );
}
my $number = $#ARGV + 1;
my %all_point;
foreach my $element (@all_set) {
    $all_point{$element}++;
}

my @common_point;
@common_point = grep { $all_point{$_} == $number; } ( keys %all_point );

foreach my $element (@ARGV) {
    my @point =
      grep { ( $all_point{$_} < $number ) } @{ $point{$element} };
    if ( @point > 0 ) {
        print "$element\n";
        foreach my $point (@point) {
            print "$point_info{$point}\n";
        }
        print "\n";
    }
}
__END__
