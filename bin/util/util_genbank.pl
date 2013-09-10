#! /usr/bin/perl -w

=head1 NAME

util_genbank.pl

=head1 SYNOPSIS

  util_genbank.pl

Options:

 --help     Show brief help and exit
 --list     List of accessions to query
 --tag      The tag that you want to process
 --outdir   Is your output dir

=head1 DESCRIPTION

Given a list of genbank accessions, finds and writes the sequence

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

####SETUP####

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::GenBank;
use Bio::SeqIO;

my ($help, $list, $primtag, $outdir);
GetOptions(
	   'h|help'       => \$help,
	   'l|list=s'     => \$list,
	   'p|primtag=s'  => \$primtag,
	   'o|outdir=s'   => \$outdir,
	   );

pod2usage(2) if $help;

`mkdir -p $outdir`;


####MAIN####

#my $gb = Bio::DB::GenBank->new(-format => 'Fasta');
my $gb = Bio::DB::GenBank->new();

open (LIST, "$list");
my $sequence = {};
while (my $entry = <LIST>){
    chomp $entry;
    my ($short, $line) = split (/\t/, $entry);
    
    print STDERR "Working on $entry\n";

    # get the record from remote genbank
    my $seq_obj = $gb->get_Seq_by_acc($line);

    # cycle through the features
    foreach my $feat ($seq_obj->get_SeqFeatures){
	next unless ($feat->primary_tag eq $primtag);

	# assign a name
	my @names;
	my $name;
#	if ($feat->has_tag("gene")){
#	    @names = $feat->get_tag_values("gene");
#	    $name = join "_", @names;
#	}
#	elsif ($feat->has_tag("product")){
#	    @names = $feat->get_tag_values("product");
#	    $name = join "_", @names;
#	}
#	else {
#	    $name = "no name found";
#	}
#	$name =~s/\s/\_/g;

	# get the locus
	if ($feat->has_tag("locus_tag")){
	    @names = $feat->get_tag_values("locus_tag");
	    $name = join "_", @names;
	}
	else {
	    $name = "NONAME";
	}
	
	# get the pt translation, first translation only
	my @pt;
	my $pt;
	if ($feat->has_tag("translation")){
	    @pt = $feat->get_tag_values("translation");
	    $pt = $pt[0];
	    
	    # get the genomic sequence                                                                 
	    my $seq = $feat->seq->seq;
	    
	    # store                                                                                         
	    $sequence->{$short}->{$name} = [$pt, $seq];
	}
	else {
	    print STDERR "$name\tno translation\n";
	}
    }
}
close (LIST);

# print out the material
foreach my $line (sort keys %$sequence){
    open (PT, ">$outdir/$line.pt.fasta");
    open (NT, ">$outdir/$line.nt.fasta");

    foreach my $acc (sort keys %{$sequence->{$line}}){
#	print PT ">$name-$acc\n$sequence->{$name}->{$acc}[0]\n";
#	print NT ">$name-$acc\n$sequence->{$name}->{$acc}[1]\n";
	print PT ">$acc\n$sequence->{$line}->{$acc}[0]\n";
        print NT ">$acc\n$sequence->{$line}->{$acc}[1]\n";
    }
    close (PT);
    close (NT);
}
