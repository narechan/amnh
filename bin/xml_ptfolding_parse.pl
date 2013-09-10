#!/usr/bin/perl -w

=head1 NAME

xml_ptfolding_parse.pl

=head1 SYNOPSIS

  xml_ptfolding_parse.pl -- mines pt folding database for desired data

Options:

 --help        Show brief help and exit
 --infile      xml data file

=head1 DESCRIPTION

Mine protein folding xml

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
use Data::Dumper;

my ($help, $infile);
GetOptions(
    'h|help'          => \$help,
    'i|infile=s'      => \$infile,
    ) or pod2usage;

pod2usage if $help;

for my $option ($infile){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# parse the infile and store
my $parser=XML::LibXML->new();
$parser->keep_blanks(0);
 
# parse the file
my $doc = $parser->parse_file($infile); #using DOM

# get the family info
my $family = get_family($doc);
print Dumper ($family);

# get the protein info
my $proteins = get_proteins($doc);
print Dumper ($proteins);

# get the selection info
my $selection = get_selection($doc);
print Dumper ($selection);

print yellow;

#####SUBS#####

sub get_family{
    my $doc = shift;

    my $family = {};
    my @families = $doc->getElementsByTagName('family');
    my $id = $families[0]->getAttribute('id');
    my @name = $families[0]->getElementsByTagName('name');
    my $name = $name[0]->getFirstChild->getData;
    
    $family->{'id'} = $id;
    $family->{'name'} = $name;

    return ($family);
}

sub get_proteins{
    my $doc = shift;
    
    my $proteins = {};
    my @families = $doc->getElementsByTagName('family');
    foreach my $family (@families){
	my $id = $family->getAttribute('id');
	my @name = $family->getElementsByTagName('name');
	my $name = $name[0]->getFirstChild->getData;
	
	# get protein information
	my @proteins = $family->getElementsByTagName('proteins');
	my @pts = $proteins[0]->getElementsByTagName('protein');
	foreach my $protein (@pts){
	    my $id = $protein->getAttribute('id');
	    
	    my @domains = $protein->getElementsByTagName('domain');
	    my $dcounter = -1;
	    foreach my $domain (@domains){
		$dcounter++;
		my $dstart = $domain->getAttribute('start');
		my $dend   = $domain->getAttribute('stop');
		
		my @dtype = $domain->getElementsByTagName('domain_type');
		my $dtype = $dtype[0]->getFirstChild->getData;
		
		$proteins->{$id}->{$dcounter}->{'type'} = $dtype;
		$proteins->{$id}->{$dcounter}->{'start'} = $dstart;
		$proteins->{$id}->{$dcounter}->{'end'} = $dend;
		
		my @structs = $domain->getElementsByTagName('structures');
		my @structures = $structs[0]->getElementsByTagName('structure');
		my @strucdata;
		my $scounter = -1;
		foreach my $structure (@structures){
		    $scounter++;
		    my $stype = $structure->getAttribute('type');
		    
		    my @struc_id = $structure->getElementsByTagName('id');
		    my $struc_id = $struc_id[0]->getFirstChild->getData;
		    
		    my @chain = $structure->getElementsByTagName('chain');
		    (my $chain = $chain[0]->getFirstChild->getData) if (@chain);
		    
		    $proteins->{$id}->{$dcounter}->{'structure'}->{$scounter}->{'type'} = $stype;
		    $proteins->{$id}->{$dcounter}->{'structure'}->{$scounter}->{'id'} = $struc_id;
		    $proteins->{$id}->{$dcounter}->{'structure'}->{$scounter}->{'chain'} = $chain;
		    
		    my @scops = $structure->getElementsByTagName('scop');
		    my @sccs = $scops[0]->getElementsByTagName('sccs');
		    foreach my $sccs (@sccs){
			my $sccs_id = $sccs->getAttribute('id');
			
			my @prob = $sccs->getElementsByTagName('probability');
			my $prob = $prob[0]->getFirstChild->getData;
			
			$proteins->{$id}->{$dcounter}->{'structure'}->{$scounter}->{'scop'}->{$sccs_id} = $prob;
		    }
		}
	    }
	}
    }
    return ($proteins);
}

sub get_selection{
    my $doc = shift;

    $selection = {};
    my @models = $doc->getElementsByTagName('model');
    foreach my $model (@models){
	my $id = $model->getAttribute('id');
	
	next unless ($id == 8);
	
	my @sites = $model->getElementsByTagName('site');
	foreach my $site (@sites){
	    my @column = $site->getElementsByTagName('column');
	    my $column = $column[0]->getFirstChild->getData;
	    
	    my @prob = $site->getElementsByTagName('probability');
	    my $prob = $prob[0]->getFirstChild->getData;
	
	    $selection->{$column} = $prob;
	}
    }
    return ($selection);
}
    
