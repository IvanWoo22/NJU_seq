#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my @info;
chomp( my $first_read = <> );
my ( $first_qname, $first_rname, $first_site ) = split( /\t/, $first_read );
$info[0] = $first_read;

while (<>) {
    chomp( my $read = $_ );
    if ( $read =~ /^$first_qname\t/ ) {
        push( @info, $read );
    }
    else {
        if ( $#info > 0 ) {
            my %name;
            $name{ $first_rname . $first_site } = $info[0];
            print "$info[0]\n";
            foreach ( 1 .. $#info ) {
                my ( undef, $rname, $site ) = split( /\t/, $info[$_] );
                unless ( exists $name{ $rname . $site } ) {
                    $name{ $rname . $site } = $info[$_];
                    print "$info[$_]\n";
                }
            }
        }
        elsif ( $#info == 0 ) {
            print "$info[0]\n";
        }
        splice(@info);
        $info[0] = $read;
        ( $first_qname, $first_rname, $first_site ) = split( /\t/, $read );
    }
}

if ( $#info > 0 ) {
    my %name;
    $name{ $first_rname . $first_site } = $info[0];
    print "$info[0]\n";
    foreach ( 1 .. $#info ) {
        my ( undef, $rname, $site ) = split( /\t/, $info[$_] );
        unless ( exists $name{ $rname . $site } ) {
            $name{ $rname . $site } = $info[$_];
            print "$info[$_]\n";
        }
    }
}
elsif ( $#info == 0 ) {
    print "$info[0]\n";
}

__END__
