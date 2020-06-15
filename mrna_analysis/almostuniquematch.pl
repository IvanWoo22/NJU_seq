#!/usr/bin/env perl
use strict;
use warnings;
use PerlIO::gzip;

open( my $in_fq, "<:gzip", $ARGV[0] );
my %fq;
while (<$in_fq>) {
    chomp;
    chomp( my $seq = <$in_fq> );
    my ( $id, undef ) = split /\s/;
    $id =~ s/^@//;
    $fq{$id} = $seq;
    readline($in_fq);
    readline($in_fq);
}
close($in_fq);

open( my $in_sam1, "<", $ARGV[1] );
my %list;
while (<$in_sam1>) {
    chomp;
    my ( $id, undef, undef, undef, $seq ) = split /\t/;
    if ( $fq{$id} eq $seq ) {
        if ( exists( $list{$id} ) ) {
            $list{$id}++;
        }
        else {
            $list{$id} = 1;
        }
    }
}
close($in_sam1);

my %uniq_id;
foreach my $k ( keys(%list) ) {
    if ( $list{$k} == 1 ) {
        $uniq_id{$k} = 1;
    }
}

open( my $in_sam2, "<", $ARGV[1] );
open( my $out_sam, ">", $ARGV[2] );
while (<$in_sam2>) {
    chomp;
    my ( $id, undef, undef, undef, $seq ) = split /\t/;
    if ( ( exists( $uniq_id{$id} ) ) and ( $fq{$id} eq $seq ) ) {
        print $out_sam ("$_\n");
    }
}
close($in_sam2);
close($out_sam);

__END__
