#! /usr/bin/perl -w

=head1 NAME

get_genbank_data_wgs.pl

=head1 SYNOPSIS

  get_genbank_data.pl

Options:

 --help     Show brief help and exit
 --list     List of accessions to query
 --outdir   Is your outdir

=head1 DESCRIPTION

Given a list of genbank assembly records, finds and writes the sequence

=head1 SEE ALSO

perl.

=head1 AUTHOR

Apurva Narechania

=head1 COPYRIGHT

Copyright (c) 2007 Cold Spring Harbor Laboratory

This library is free software;  you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::GenBank;
use Bio::SeqIO;

my ($help, $list, $outdir);
GetOptions(
	   'h|help'       => \$help,
	   'l|list=s'     => \$list,
	   'o|outdir=s'   => \$outdir,
	   );

pod2usage(2) if $help;

`mkdir -p $outdir/report`;
`mkdir -p $outdir/contigs`;
`mkdir -p $outdir/prots`;

my $queries = {};
open (L, "$list");
while (my $line = <L>){
    chomp $line;
    print STDERR "Working on $line\n";
    
    # get the assembly report
    `curl ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/$line.assembly.txt > $outdir/report/$line.report`;
    
    # parse the assmbly report and get the sequence for each contig
    open (A, "$outdir/report/$line.report");
    open (FO, ">$outdir/contigs/$line.fa");
    open (PO, ">$outdir/prots/$line.fa");
    while (my $line = <A>){
	chomp $line;
	next if ($line =~m/^\#/);
	
	my @line = split (/\t/, $line);
	my $gb = new Bio::DB::GenBank;
	my $seq_obj = $gb->get_Seq_by_acc($line[4]);
	
	# get the contigs
	my $id = $seq_obj->id();
	my $seq = $seq_obj->seq();
	print FO ">$id\n$seq\n";

	# get the proteins
	my @feat_ary = $seq_obj->top_SeqFeatures();
	foreach my $feat (@feat_ary) {
	    if ($feat->primary_tag() eq "CDS") {
		if (($feat->has_tag('translation')) and $feat->has_tag('locus_tag')){
		    my @locus = $feat->each_tag_value('locus_tag');
		    my @trans = $feat->each_tag_value('translation');
		    for (my $i = 0; $i < @locus; $i++) {
			print PO ">$locus[$i]\n$trans[$i]\n";
		    }
		}
		else{
		    next;
		}
	    }
	    else{
		next;
	    }
	}
    }
    close (FO);
    close (PO);
    close (A);
}
=cut

=head
foreach my $query (keys %$queries){
    next if ($query eq "AACK01000001:AACK01000140[accn]");
    my $q = Bio::DB::Query::GenBank->new
	(-query   => $query,
	 -db      => 'nucleotide');
    
    my $seqio = $gb->get_Stream_by_query($q);
 
    open (FO, ">$outdir/$queries->{$query}.fa");
    while(my $seq = $seqio->next_seq){
	my $id = $seq->id();
	my $sequence = $seq->seq();
	
	print FO ">$id\n$sequence\n";
    }
    close (FO);
}
=cut
