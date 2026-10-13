use strict;
use warnings;
use autodie;
use PerlIO::gzip;

sub REV_COMP {
    my $SEQ     = reverse(shift);
    my $R_C_SEQ = $SEQ =~ tr/AGTCagtc/TCAGtcag/r;
    return $R_C_SEQ;
}

@ARGV == 3 or die "Usage: pe_consistency.pl R1.fq.gz R2.fq.gz matching_R1.fq.gz\n";
if (-e $ARGV[2]) {
    my @out = stat($ARGV[2]);
    for my $in (@ARGV[0,1]) {
        my @s = stat($in);
        die "Input and output must differ\n" if @s && $s[0] == $out[0] && $s[1] == $out[1];
    }
}
open( my $in_fh1, "<:gzip", $ARGV[0] );
open( my $in_fh2, "<:gzip", $ARGV[1] );
open( my $out_fh, ">:gzip", $ARGV[2] );

my $all    = 0;
my $proper = 0;
while (<$in_fh1>) {
    $all++;
    chomp( my $qname   = $_ );
    chomp( my $seq1    = <$in_fh1> );
    chomp( my $info    = <$in_fh1> );
    chomp( my $quality = <$in_fh1> );
    my $qname2 = <$in_fh2>;
    my $seq2 = <$in_fh2>;
    my $plus2 = <$in_fh2>;
    my $qual2 = <$in_fh2>;
    defined($quality) && defined($qual2) or die "Truncated or unequal FASTQ pair\n";
    chomp($seq2, $qname2, $plus2, $qual2);
    my ($id1) = $qname =~ /^@(\S+)/;
    my ($id2) = $qname2 =~ /^@(\S+)/;
    defined($id1) && defined($id2) && $info =~ /^\+/ && $plus2 =~ /^\+/ or die "Malformed FASTQ pair\n";
    $id1 =~ s{/1$}{}; $id2 =~ s{/2$}{};
    $id1 eq $id2 or die "FASTQ IDs are out of sync: $id1 / $id2\n";
    length($seq1) == length($quality) && length($seq2) == length($qual2) or die "FASTQ sequence/quality mismatch\n";
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
eof($in_fh2) or die "R2 has extra records\n";
my $proportion = $all ? $proper / $all * 100 : 0;
$proportion = sprintf( "%.2f", $proportion );
print("Total:\t$all\nConsistency:\t$proper\nProportion:\t$proportion%\n");

close($in_fh1);
close($in_fh2);
close($out_fh);

__END__
