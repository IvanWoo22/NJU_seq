#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

foreach my $filename (@ARGV) {
    open( my $IN, "<", $filename );
    my $name = $filename;
    $name =~ s/\.bowtie2\.log$//;
    my ( $total_reads, $aligned_reads ) = ( 0, 0 );
    my ( $time, $rate ) = ("NA", "NA");
    while (<$IN>) {
        chomp;
        my $in = $_;
        if ( $in =~ m/(\d+) reads; of these:/ ) {
            $total_reads += $1;
        }
        if (
            (
                $in =~ m/(\d+) \([.0-9]+%\) aligned concordantly exactly 1 time/
            )
            or ( $in =~ m/(\d+) \([.0-9]+%\) aligned concordantly >1 times/ )
            or ( $in =~ m/(\d+) \([.0-9]+%\) aligned discordantly 1 time/ )
          )
        {
            $aligned_reads += $1;
        }
        if (   ( $in =~ m/(\d+) \([.0-9]+%\) aligned exactly 1 time/ )
            or ( $in =~ m/(\d+) \([.0-9]+%\) aligned >1 times/ ) )
        {
            $aligned_reads += $1 / 2;
        }
        if ( $in =~ m/Overall time: (\S+)/ ) {
            $time = $1;
        }
    }
    if ($total_reads > 0) { $rate = sprintf("%.2f%%", $aligned_reads / $total_reads * 100); }
    else { warn "No reads reported in $filename\n"; }
    print "$name\t$total_reads\t$aligned_reads\t$rate\t$time\n";
    close($IN);
}

__END__
