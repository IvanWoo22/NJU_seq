package NJUSeq::SiteTable;
use strict;
use warnings;
use autodie;
use Exporter 'import';
our @EXPORT_OK = qw(read_sites);

# Shared only by the legacy site-set utilities. Presence is per file, not row.
sub read_sites {
    my @files = @_;
    @files or die "Provide at least one site table\n";
    my (%rows, %members, @by_file);
    for my $i (0 .. $#files) {
        open my $fh, '<', $files[$i];
        $by_file[$i] = {};
        while (<$fh>) {
            chomp;
            next if /^\s*$/ || /^#/ || /^(?:Chromosome|Position)\t/;
            my @f = split /\t/, $_, -1;
            @f >= 3 && $f[1] =~ /^[1-9][0-9]*$/ && $f[2] =~ /^[+-]$/
                or die "Invalid site row in $files[$i] line $.\n";
            my $key = join "\t", @f[0..2];
            if (exists $by_file[$i]{$key} && $by_file[$i]{$key} ne $_) {
                die "Conflicting duplicate site in $files[$i]: $key\n";
            }
            $by_file[$i]{$key} = $_;
            $rows{$key} = $_; # preserve the old last-table metadata convention
            $members{$key}{$i} = 1;
        }
        close $fh;
    }
    return (\%rows, \%members, \@by_file);
}
1;
