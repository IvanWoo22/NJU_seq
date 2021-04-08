#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open my $REF,    "<", $ARGV[0];
open my $IN_SAM, "<", $ARGV[1];

sub CIGAR_CONVERT {
    my $CIGAR = shift;
    my $SEQ   = shift;
    my @SEQ   = split( "", $SEQ );
    my @ALIGNMENT;
    while ( $CIGAR =~ /^([0-9]+)/ ) {
        my $LENGTH = $1;
        $CIGAR =~ s/^([0-9]+)//;
        if ( $CIGAR =~ /^[=|X]/ ) {
            $CIGAR =~ s/^[=|X]//;
            foreach ( 1 .. $LENGTH ) {
                push @ALIGNMENT, shift(@SEQ);
            }
        }
        elsif ( $CIGAR =~ /^I/ ) {
            $CIGAR =~ s/^I//;
            foreach ( 1 .. $LENGTH ) {
                shift(@SEQ);
            }
        }
        elsif ( $CIGAR =~ /^D/ ) {
            $CIGAR =~ s/^D//;
            foreach ( 1 .. $LENGTH ) {
                push @ALIGNMENT, "";
            }
        }
    }
    return \@ALIGNMENT;
}

my @site;
my @site_base;

readline($REF);
my $ref;
while (<$REF>) {
    s/\r?\n//;
    $ref .= $_;
}
close($REF);

@site_base = split( "", $ref );

foreach ( 0 .. $#site_base ) {
    ${ $site[$_] }{"A"} = 0;
    ${ $site[$_] }{"G"} = 0;
    ${ $site[$_] }{"C"} = 0;
    ${ $site[$_] }{"T"} = 0;
    ${ $site[$_] }{"N"} = 0;
}

while (<$IN_SAM>) {
    chomp;
    if ( $_ =~ /\t$ARGV[2]\t/ ) {
        my ( undef, undef, $position, $cigar, $seq ) = split( "\t", $_ );
        my $alignment = CIGAR_CONVERT( $cigar, $seq );
        my @alignment = @{$alignment};
        foreach ( 0 .. $#alignment ) {
            ${ $site[ $position - 1 + $_ ] }{ $alignment[$_] } += 1;
        }
    }
}
close($IN_SAM);

foreach ( 0 .. $#site_base ) {
    my $site_number = $_ + 1;
    my @base_sort =
      sort { ${ $site[$_] }{$b} <=> ${ $site[$_] }{$a} } keys %{ $site[$_] };
    my $max = shift(@base_sort);
    print(
"$site_number\t$site_base[$_]\t${ $site[$_] }{\"A\"}\t${ $site[$_] }{\"G\"}\t${ $site[$_] }{\"C\"}\t${ $site[$_] }{\"T\"}\t${ $site[$_] }{\"N\"}\t$max"
    );
    if ( $site_base[$_] ne $max ) {
        print("\tmut\n");
    }
    else {
        print("\n");
    }
}

__END__
