#!/usr/bin/perl -w

=head1 NAME

gutenberg_harvest.pl

=head1 SYNOPSIS

  gutenberg_harvest.pl -- this program mines gutenberg for desired texts

Options:

 --help        Show brief help and exit
 --infile      Author\tTitle infile
 --catalog     Is the RDF XML file of the entire gutenberg catalog
 --outdir      Is your output dir

=head1 DESCRIPTION

Mine gutenberg respectfully

Usage examp:

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

#####
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use XML::LibXML;
use LWP::Simple;

my ($help, $infile, $catalog, $outdir);
GetOptions(
    'h|help'          => \$help,
    'i|infile=s'      => \$infile,
    'c|catalog=s'     => \$catalog,
    'o|outdir=s'      => \$outdir,
    ) or pod2usage;

pod2usage if $help;

for my $option ($infile, $catalog, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# parse the catalog and store authors and titles
#binmode (STDOUT, ':utf8');
my $parser=XML::LibXML->new();
$parser->keep_blanks(0);
 
my $doc = $parser->parse_file('catalog.rdf'); #using DOM

my $books = {};
my @booknodes = $doc->findnodes('/rdf:RDF/pgterms:etext');

foreach my $booknode (@booknodes){
    my $etext_no = $booknode->getAttribute('rdf:ID');
    $etext_no =~s/^etext//;

    my $author;
    my $book;
    foreach my $creator ($booknode->findnodes('dc:creator//text()')) {
	$author = $creator->textContent;
    }
    foreach my $title ($booknode->findnodes('dc:title//text()')){
	$book = $title->textContent;
	$book =~s/\n/ /g;
    }

#    print "$etext_no\t$author\t$book\n";
    ($books->{$author}->{$book} = $etext_no) unless (exists ($books->{$author}->{$book}));
}

# cycle through authors and works to find matches in the catalog
my $matches = {};
open (MK, ">$outdir/matchkey");
open (IN, "$infile");
while (my $line = <IN>){
    chomp $line;
    
    my ($author, $work) = split (/\t/, $line);

    my $string = $author . "_" . $work;
    $string =~s/\s/_/g;
    $string =~s/[^\w\d]//g;

    my $match = 0;
    my $matchcode;
    foreach my $creator (keys %$books){
	if ($creator =~m/$author/i){
	    foreach my $book (keys %{$books->{$creator}}){
		if ($book =~m/$work/i){
		    $matchcode = $books->{$creator}->{$book};
		    $match = 1;
		}
		else {
		    next;
		}
	    }
	}
	else {
	    next;
	}
    }

    $matches->{$string} = $matchcode;
    print MK "$line\t$match\t$matchcode\n";
}
close (IN);
close (MK);	  

foreach my $work (keys %$matches){
    print STDERR "Getting $work\n";
    my $content = get("http://www.gutenberg.org/cache/epub/$matches->{$work}/pg$matches->{$work}.txt");
    next unless ($content);
    open (FH, ">$outdir/$work");
    print FH "$content\n";
    close (FH);
    sleep 2;
}
