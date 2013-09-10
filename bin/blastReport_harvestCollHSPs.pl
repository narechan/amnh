#!/usr/bin/perl

# Purpose: This program uses standard blast output to harvest 
#          hit sequences.

#          Options:
#          -b is the blast report
#          -o is the outdir

# use the following modules                                                 
use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Location::Simple;
use Bio::DB::GenBank;
use strict;

# get inputs                                                                    
my %opts = ();
getopts ('b:o:', \%opts);
my $blastreport    = $opts{'b'};
my $outdir         = $opts{'o'};

`mkdir $outdir`;

# parse the blast report
print STDERR "Parsing blast report\n";
my $blastparse = parse ($blastreport, $outdir);

# mine data from the parse
foreach my $hit (sort keys %$blastparse){
    
    print STDERR "Getting $hit\n";
    my $start = $blastparse->{$hit}[0];
    my $end = $blastparse->{$hit}[1];
    my $ori = $blastparse->{$hit}[2];
    
    unless ($ori){
	print STDERR "No data for $hit\n";
	next;
    }

    # get the hit sequence
    my $gb = new Bio::DB::GenBank(-format => 'Fasta');
    my $seq_obj = $gb->get_Seq_by_acc($hit);
    my $id = $seq_obj->id();
    my $desc = $seq_obj->desc();
    my $seq = $seq_obj->seq();
    
    # create a location object and extract the subseq
    print STDERR "Extracting hit subseq\n";
    my $location = Bio::Location::Simple->new(-start  => $start,
					      -end    => $end,
					      -strand => $ori);
    
    my $subseq   = $seq_obj->subseq($location);

    print STDERR "$hit\t$start\t$end\t$ori\n";
    print ">$id\n$subseq\n";
}
close (F);



#####SUBS#####

sub parse {
    my $infile = shift;
    my $outdir = shift;
    my $in = Bio::SearchIO->new(-file   => "$infile",
				-format => 'blast');
    
    my $parse = {};
    open (OUT, ">$outdir/$infile.parse");
    while( my $result = $in->next_result ) {
        my $hit_count = 0;  
        unless ($result->hits) {
            print OUT join(
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
	    
	    my $minstart = 1000000000000000000000000000000;
	    my $maxend   = 0;
	    my $ori;
	    my $hitname = $hit->name;

	    my $hsp_count = 0;
	    my @pair;
	    while (my $hsp = $hit->next_hsp ) {
		$hsp_count++;
#		last if ($hsp_count > 1); # get one HSP 

		# parse out min start and max end among multiple hsps
		my $hspstart = $hsp->start('hit');
		my $hspend   = $hsp->end('hit');
		my $hspori   = $hsp->strand('hit');
		($minstart = $hspstart) if ($hspstart < $minstart);
		($maxend   = $hspend) if ($hspend > $maxend);
		$ori = $hspori;

		print OUT join(
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
	    push (@pair, $minstart, $maxend, $ori);
	    $parse->{$hitname} = \@pair;
	    
	}
	
    }
    return ($parse);
}
