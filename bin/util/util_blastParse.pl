#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SearchIO;
use Bio::Search::Tiling::MapTiling;

my ($help, $result);
GetOptions(
    'h|help'          => \$help,
    'r|result=s'     => \$result,
    ) or pod2usage;

pod2usage if $help;

#####MAIN#####

print STDERR "Parsing $result\n";
parse ($result);
    
#####SUBS#####

sub parse {
    my $infile = shift;
    my $in = Bio::SearchIO->new(-file   => "$infile",
				-format => 'blast');
    
    # foreach query
    while (my $result = $in->next_result){
        my $hit_count = 0;  
	my $qname = $result->query_name;
	
	# foreach hit
	my @hits;
	my $tophit;
	my $topdesc;
	my $bachits = 0;
	while (my $hit = $result->next_hit){
	    $hit_count++;
#	    last if ($hit_count > 10); # get top ten hits
#	    last if ($hit_count > 30);

	    # get the hit data                                                                   
	    my $hname = $hit->name;
            my $hdesc = $hit->description;
            my $evalue = $hit->significance;

	    # see of the hit has a match to our genus
#	    if ($hdesc =~m/Rickettsia/){
#	    if ($hdesc =~m/Rickettsia|Anaplasma|Wolbachia|Orientia|Ehrlichia|Neorickettsia|Odyssella|Midichloria/){
#		$bachits++;
#	    }
	    
	    # get hit stats for the result
	    if ($hit_count == 1){
		$tophit = $hname;
		$topdesc = $hdesc;

#		if ($topdesc =~m/Rickettsia/){
		if ($hdesc =~m/Rickettsia|Anaplasma|Wolbachia|Orientia|Ehrlichia|Neorickettsia|Odyssella|Midichloria/){
		    $bachits++;
		}
		
	    }
	    
#	    unless ($hname =~m/plant\d+/){
#		$bachits++;
#	    }
	    
	    # tile the hsps
#	    my $tiling = Bio::Search::Tiling::MapTiling->new($hit);
#	    my $qiden  = $tiling->frac_identical(-type=>'query', -action=>'exact');
	    
	    my $joined = $qname . "\t" . $hname . "\t" . $hdesc . "\t" . $evalue;
	    push (@hits, $joined);
#	    print "$qname\t$hit_count\t$hname\t$evalue\t$qiden\n";
#	    print "$qname\t$hit_count\t$hname\t$hdesc\t$evalue\n";
	    
	}
	if ($bachits >= 1){
	    if ($tophit){
#	        print "==>$qname\t$tophit\t$topdesc\t$bachits\n";
		print "==>$qname\t$tophit\t$topdesc\n";
	    }
	    else {
#	        print "==>$qname\tNONE\tNONE\t$bachits\n";
		print "==>$qname\tNONE\tNONE\n";
	    }
	    
	    foreach my $hit (@hits){
		print "$hit\n";
	    }
	}
    }
}
