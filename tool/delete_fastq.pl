#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use PerlIO::gzip;
use Getopt::Long;

#---------------#
# GetOpt section
#---------------#

=head1 NAME
delete_fastq.pl -- Delete reads with specified IDs in a FastQ file.
=head1 SYNOPSIS
    perl delete_fastq.pl -n name_list.txt -i input.fq -o output.fq
        Options:
            --help\-h  Brief help message
            --name\-n  The sequences we want to fetch
            --in\-i  The FastQ file with path
            --out\-o  The FastQ file with path
=cut

Getopt::Long::GetOptions(
    'help|h'   => sub { Getopt::Long::HelpMessage(0) },
    'name|n=s' => \my $name_list,
    'in|i=s'   => \my $in_fq,
    'out|o=s'  => \my $out_fq,
) or Getopt::Long::HelpMessage(1);

open( my $name, "<", $name_list );
my %target_name;
while (<$name>) {
    chomp;
    $target_name{$_} = 0;
}

defined($name_list) && defined($in_fq) && defined($out_fq) or die "Require -n, -i and -o\n";
if (-e $out_fq) { my @i = stat($in_fq); my @o = stat($out_fq);
    die "Input and output must differ\n" if @i && $i[0] == $o[0] && $i[1] == $o[1]; }
my $in_fh;
if ( $in_fq =~ /\.gz$/ ) {
    open $in_fh, "<:gzip", $in_fq;
}
else {
    open( $in_fh, "<", $in_fq );
}

my $out_fh;
if ( $out_fq =~ /\.gz$/ ) {
    open $out_fh, ">:gzip", $out_fq;
}
else {
    open( $out_fh, ">", $out_fq );
}

while (<$in_fh>) {
    my $qname    = $_;
    my $sequence = <$in_fh>;
    my $t        = <$in_fh>;
    my $quality  = <$in_fh>;
    my ($read_id) = $qname =~ /^@(\S+)/;
    defined($quality) && defined($read_id) && $t =~ /^\+/
        or die "Malformed or truncated FASTQ record at line $.\n";
    (my $seq_check = $sequence) =~ s/\r?\n$//;
    (my $qual_check = $quality) =~ s/\r?\n$//;
    length($seq_check) == length($qual_check) or die "FASTQ sequence/quality length mismatch for $read_id\n";
    unless ( exists $target_name{$read_id} ) {
        print $out_fh "$qname$sequence$t$quality";
    }
}

close($in_fh);
close($out_fh);

__END__
