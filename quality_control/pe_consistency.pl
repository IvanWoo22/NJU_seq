use strict;
use warnings;
use autodie;
use PerlIO::gzip;

sub REV_COMP {
    my $SEQ     = reverse(shift);
    my $R_C_SEQ = $SEQ =~ tr/AGTCagtc/TCAGtcag/r;
    return $R_C_SEQ;
}

open( my $in_fh1, "<:gzip", $ARGV[0] );
open( my $in_fh2, "<:gzip", $ARGV[1] );
open( my $out_fh, ">",      $ARGV[2] );

my $all    = 0;
my $proper = 0;
while (<$in_fh1>) {
    $all++;
    chomp( my $qname   = $_ );
    chomp( my $seq1    = <$in_fh1> );
    chomp( my $info    = <$in_fh1> );
    chomp( my $quality = <$in_fh1> );
    readline($in_fh2);
    chomp( my $seq2 = <$in_fh2> );
    readline($in_fh2);
    readline($in_fh2);
    my $r_c_seq2 = REV_COMP($seq2);
    my @seq1     = split( //, $seq1 );
    my @seq2     = split( //, $r_c_seq2 );

    if ( $#seq1 == $#seq2 ) {
        foreach ( 0 .. $#seq1 ) {
            if ( $seq1[$_] eq "N" ) {
                $seq1[$_] = $seq2[$_];
            }
            if ( $seq2[$_] eq "N" ) {
                $seq2[$_] = $seq1[$_];
            }
        }
        my $out_seq1 = join( "", @seq1 );
        my $out_seq2 = join( "", @seq2 );
        if ( $out_seq1 eq $out_seq2 ) {
            $proper++;
            print $out_fh ("$qname\n$seq1\n$info\n$quality\n");
        }
    }
}
my $proportion = $proper / $all * 100;
$proportion = sprintf( "%.2f", $proportion );
print "Total:\t$all\nConsistency:\t$proper\nProportion:\t$proportion";

close($in_fh1);
close($in_fh2);

__END__
