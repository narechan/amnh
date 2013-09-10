#!/usr/bin/perl -w

=head1 NAME

blast_parse.pl

=head1 SYNOPSIS

  blast_parse.pl -- 
              

Options:

 --help        Show brief help and exit
 --result      Is your blast result report

=head1 DESCRIPTION

Given a dir of blast reports
    parse the blast.

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
use Parallel::ForkManager;
use Bio::SearchIO;

my ($help, $result);
GetOptions(
    'h|help'          => \$help,
    'r|result=s'     => \$result,
    ) or pod2usage;

pod2usage if $help;

for my $option ($result){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

print STDERR "Parsing $result\n";
parse ($result);
    
#####SUBS#####

sub parse {
    my $infile = shift;
    my $in = Bio::SearchIO->new(-file   => "$infile",
				-format => 'blast');
    
    while( my $result = $in->next_result ) {
        my $hit_count = 0;  
        unless ($result->hits) {
            print join(
			   "\t",
			   $result->query_name,        #1
			   $result->query_length,      #2
			   'No hits found',            #3
			   '0',                        #4
			   '1',                        #5
			   '1',                        #6
			   '1000',                     #7
			   '0',                        #8
			   '0',                        #9
			   '0',                        #10
			   '0',                        #11
			   '0',                        #12
			   '0',                        #13
			   '0',                        #14
			   '0',                        #15
			   '0',                        #16
			   '0',                        #17
			   '0',                        #18
			   '0',                        #19
			   $result->query_description, #20
			   'NULL',                     #21
			   ), "\n";
        }
	
	while( my $hit = $result->next_hit ) {
	    $hit_count++;
#	    last if ($hit_count > 1); # get one hit
	    
	    my $hsp_count = 0;
	    while (my $hsp = $hit->next_hsp ) {
		$hsp_count++;
#		last if ($hsp_count > 1); # get one HSP 

		# length requirement if desired
#		my $hsplen = $hsp->hsp_length;
#		my $qrylen = $result->query_length;
#		my $lenreq = $hsplen / $qrylen;
#		if (($lenreq < 0.95) or ($lenreq > 1.05)){
#		    print OUT join(
#				   "\t",
#				   $result->query_name,        #1                                        
#				   $result->query_length,      #2                                            
#				   'Length insufficient',            #3                             
#				   '0',                        #4                                            
#				   '1',                        #5                                              
#				   '1',                        #6                                             
#				   '1000',                     #7                                              
#				   '0',                        #8                                           
#				   '0',                        #9                                            
#				   '0',                        #10                                           
#				   '0',                        #11                                           
#				   '0',                        #12                                        
#				   '0',                        #13                                           
#				   '0',                        #14                                           
#				   '0',                        #15                                           
#				   '0',                        #16                                             
#				   '0',                        #17                                         
#				   '0',                        #18                                         
#				   '0',                        #19                                          
#				   $result->query_description, #20                                           
#				   'NULL',                     #21                                            
#				   ), "\n";
#		    next;
#		}
		    
		print join(
			       "\t",
			       $result->query_name,           #1
			       $result->query_length,         #2
			       $hit->name,                    #3
			       $hit->length(),                #4
			       $hit_count,                    #5
			       $hsp->rank,                    #6
			       $hsp->evalue(),                #7
			       $hsp->score,                   #8
			       $hsp->frac_identical('total'), #9
			       $hsp->start('query'),          #10
			       $hsp->end('query'),            #11
			       $hsp->gaps('query'),           #12
			       $hsp->frac_identical('query'), #13 won't be accurate for blastx
			       $hsp->strand('query'),         #14
			       $hsp->start('hit'),            #15
			       $hsp->end('hit'),              #16
			       $hsp->gaps('hit'),             #17
			       $hsp->frac_identical('hit'),   #18
			       $hsp->strand('hit'),           #19
			       $result->query_description,    #20
			       $hit->description,             #21
			       ), "\n";
		
	    }
	    
	}
	
    }
}
