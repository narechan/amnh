#!/usr/bin/perl -w

=head1 NAME

rmChars.pl

=head1 SYNOPSIS

rmChars.pl

Options:

--matrix is your alignfile
--type is uninf or constant
--outdir is your outdir

The matrix must be in nexus format.

Requires the bioperl libs. 

=head1 DESCRIPTION

This program removes uninf or invariant chars from a matrix
and recomputes the char partitions in order to compress large
datasets for subsequent parsimony analysis.

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2008 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use TreeSupports;
use Bio::AlignIO;
use Bio::SeqIO;

my ($help, $matrixfile, $type, $outdir);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    't|type=s'        => \$type,
    'o|outdir=s'      => \$outdir,   
    ) or pod2usage;
pod2usage if $help;

for my $option ($matrixfile, $type){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir/cmds`;
`mkdir -p $outdir/rmChars`;
`mkdir -p $outdir/logs`;

#####MAIN#####

# instantiate the object and load stuff we need
my $supportobj = TreeSupports->new;
$supportobj->load_aln    ($matrixfile);

# get the charsets
my $charsets = $supportobj->get_charsets;

# store alignment information
my $partitions = {};
my $lengths    = {};
my $alnin = Bio::AlignIO->new(-file   => "$matrixfile",
			      -format => "nexus");

# only one aln there
my $aln = $alnin->next_aln();
foreach my $seq ($aln->each_seq){
    my $id        = $seq->display_id;
    
    foreach my $charset (sort keys %$charsets){
	print STDERR "Storing\t$id\t$charset\n";

	my $coords = $charsets->{$charset};
	my ($start, $end) = split (/\-/, $coords);
	
	my $partition = $seq->subseq($start, $end);
	my $partlen   = length ($partition);
	
	$partitions->{$charset}->{$id} = $partition;
	$lengths->{$charset} = $partlen;
    }
}

# build nexus files for every charset
# and compress for selected chars
foreach my $charset (sort keys %$charsets){
    print STDERR "Compress $charset\n";

    open (C, ">$outdir/cmds/$charset.nxs");
    print C "#NEXUS\n";
    print C "BEGIN DATA;\n";
    print C "DIMENSIONS NTAX=$supportobj->{'alignment'}->{'ntax'} NCHAR=$lengths->{$charset};\n";
    print C "FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" DATATYPE=PROTEIN MISSING=? GAP=-;\n";
    print C "MATRIX\n\n";
    
    foreach my $id (sort keys %{$partitions->{$charset}}){
	print C "$id\t$partitions->{$charset}->{$id}\n";
    }
    
    print C ";\n";
    print C "END;\n\n";
    
    print C "BEGIN PAUP;\n";
    print C "log start file=$outdir/logs/$charset.log replace=yes;\n";
    print C "exclude $type;\n";
    print C "export file=$outdir/rmChars/$charset.rmChars format=text;\n";
    print C "END;\n";
    close (C);
    
    `paup -n $outdir/cmds/$charset.nxs`;
}
    
# calculate the the char partitions based on compression
# and report those charpars that have no informative chars
opendir (R, "$outdir/rmChars");
my @cfiles = grep (/^.+\..+$/, readdir(R));
closedir (R);

my $newcharsets = {};
my $newseqs = {};
my $nchars = 0;
my $start = 1;

open (NS, ">$outdir/omitted.charpars");
foreach my $cfile (sort @cfiles){
    print STDERR "Calc $cfile\n";

    open (F, "$outdir/rmChars/$cfile");
    my ($charset, $file) = split (/\./, $cfile);

    my $seqlen;
    my $end;
    my @seqs;
    my $signal = 0;
    while (my $line = <F>){
        chomp $line;

        my ($acc, $seq) = split (/\s+/, $line);
	push (@seqs, $line);
	$seqlen = length ($seq);
	$end = $start + ($seqlen - 1);

	($signal = 1) if ($seqlen == 0);
	($newcharsets->{$charset} = "$start-$end") unless ($signal == 1);
    }
    close (F);
    
    print STDERR "$charset\t$seqlen\n";
    
    $start = $end + 1;
    $nchars += $seqlen;
    ($newseqs->{$charset} = [@seqs]) unless ($signal == 1);
    (print NS "$charset\n") if ($signal == 1);
}


# print the new matrix
open (C, ">$outdir/$matrixfile.cmprss");
print C "#NEXUS\n";
print C "BEGIN DATA;\n";
print C "DIMENSIONS NTAX=$supportobj->{'alignment'}->{'ntax'} NCHAR=$nchars;\n";
print C "FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" DATATYPE=PROTEIN MISSING=? GAP=-;\n";
print C "MATRIX\n";

foreach my $charset (sort keys %$newseqs){
    print C "[Partition $charset chars $newcharsets->{$charset}]\n";
    
    my $seqblock = join ("\n", @{$newseqs->{$charset}});
    print C "$seqblock\n\n";
}

print C ";\n";
print C "END;\n\n";

print C "BEGIN SETS;\n";
foreach my $charset (sort keys %$newcharsets){
    print C "CHARSET $charset = $newcharsets->{$charset};\n";
}
print C "END;\n";

close (C);
