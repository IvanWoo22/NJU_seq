#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Getopt::Long;

Getopt::Long::GetOptions(
    'help|h' => sub { Getopt::Long::HelpMessage(0) },
    'id=s'   => \my $id_feature,
    'name=s' => \my $name_feature,
    'file=s' => \my $file,
    'col=s'  => \my $col,
) or Getopt::Long::HelpMessage(1);

defined($id_feature) && defined($name_feature) && defined($file) && defined($col) && $col =~ /^[1-9][0-9]*$/
    or die "Require --id --name --file and positive --col (insertion column)\n";
sub GET_INFO {
    my ( $ID_FEATURE, $NAME_FEATURE, $INFO ) = @_;
    my ( $ID, $NAME );
    if ( $INFO =~ m/$ID_FEATURE(\w+)\.?[0-9]*/ ) {
        $ID = $1;
    }
    else {
        warn("There is a problem in $INFO;\n");
    }
    if ( $INFO =~
        m/$NAME_FEATURE([A-Z,a-z,0-9,\-,\.,\(,\), ,\[,\],%,:,\/,\+,',`]+);/ )
    {
        $NAME = $1;
    }
    else {
        warn("There is a problem in $INFO;\n");
        $NAME = "";
    }
    return ( $ID, $NAME );
}

my %gene_info;

while (<STDIN>) {
    chomp;
    next if /^#/;
    my ( undef, undef, $type, undef, undef, undef, undef, undef, $info ) =
      split /\t/;
    next unless defined($type) && $type eq "gene";
    my ( $id, $name ) = GET_INFO( $id_feature, $name_feature, $info );
    $gene_info{$id} = $name if defined $id;
}

open( my $REF_FH, "<", $file );
while (<$REF_FH>) {
    chomp;
    my @line = split "\t";
    $col <= @line + 1 or die "Insertion column exceeds table width\n";
    if (/^(?:#?Chromosome|#?Position)\t/) { splice @line, $col-1, 0, "GeneName"; print join("\t", @line), "\n"; next; }
    @line >= 7 or die "Site table needs gene IDs in column 7\n";
    my @id = split "/", $line[6];
    my $name;
    if ( $#id > 0 ) {
        my @name = map { $gene_info{$_} // "NA" } @id;
        $name = join( "/", @name );
    }
    else {
        $name = $gene_info{ $id[0] } // "NA";
    }
    splice @line, $col-1, 0, $name;
    my $out = join "\t", @line;
    print "$out\n";
}
close($REF_FH);
__END__
