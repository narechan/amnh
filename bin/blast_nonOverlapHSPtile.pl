#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SearchIO;
use Bio::AlignIO;

my ($help, $result);
GetOptions(
    'h|help'         => \$help,
    'r|result=s'     => \$result,
    ) or pod2usage;

pod2usage if $help;

#####MAIN#####

print STDERR "Parsing $result\n";
my $infile = shift;
my $in = Bio::SearchIO->new(-file   => "$result",
			    -format => "blast");

# foreach query                                                                                         
while (my $result = $in->next_result){
    my $hit_count = 0;
    my $qname = $result->query_name;
    
    # foreach hit                                                                                    
    while (my $hit = $result->next_hit){
	$hit_count++;

	# get the hit data                                                                                 
	my $hname = $hit->name;
	my $hdesc = $hit->description;
	my $evalue = $hit->significance;
	my ($alns, $keeps, $qstarts) = non_overlapping_hsps ($hit);

	my $concatq;
	my $concath;
	my $deflineq = $qname;
	my $deflineh = $hname;
#	foreach my $a (sort {$a <=> $b} keys %$alns){
	foreach my $a (sort {$qstarts->{$a} <=> $qstarts->{$b}} keys %$qstarts){
	    my $counter = 0;
	    my $alignment = $alns->{$a};
	    foreach my $seq ($alignment->each_seq()){
		$counter++;
		
		if ($counter == 1){
		    $concatq .= $seq->seq;
		    my $ints = $keeps->{$a}->{'query'}[0];
		    my $inte = $keeps->{$a}->{'query'}[1];
		    $deflineq .= " " . $ints . "-" . $inte;
		}
		else {
                    $concath .= $seq->seq;
		    my $ints = $keeps->{$a}->{'hit'}[0];
                    my $inte = $keeps->{$a}->{'hit'}[1];
                    $deflineh .= " " . $ints . "-" . $inte;
		}
	    }
	}
	print ">$deflineq\n$concatq\n>$deflineh\n$concath\n";
    }
}
    
#####SUBS#####

sub non_overlapping_hsps {
    my $hit = shift;

    my $hsp_count = 0;
    
    # get the hit data                                                                   
    my $hname = $hit->name;
    my $hdesc = $hit->description;
    my $evalue = $hit->significance;
    
    # get the first/best hsp
    my $hspalnkeep = {};
    my $hspkeep = {};
    my $hspqstart = {};
    my $firsthsp = $hit->next_hsp;
    my $fqstart   = $firsthsp->start('query');
    my $fqend     = $firsthsp->end('query');
    my $fhstart   = $firsthsp->start('hit');
    my $fhend     = $firsthsp->end('hit');
    my $faln      = $firsthsp->get_aln;
    
    $hspkeep->{$hsp_count}->{'query'} = [$fqstart, $fqend];
    $hspkeep->{$hsp_count}->{'hit'}   = [$fhstart, $fhend];
    $hspqstart->{$hsp_count} = $fqstart;
    $hspalnkeep->{$hsp_count} = $faln;

    # cycle through the hsps
    while (my $hsp = $hit->next_hsp){
	$hsp_count++;
	my $qstart   = $hsp->start('query');
	my $qend     = $hsp->end('query');
	my $hstart   = $hsp->start('hit');
	my $hend     = $hsp->end('hit');
	my $hspaln = $hsp->get_aln;

	# cycle through the kept hsps and see if 
	# new one overlaps with any kept ones
	# if so, skip
	my $signal = 0;
	foreach my $khsp (keys %$hspkeep){
	    
	    my $qkhspstart = $hspkeep->{$khsp}->{'query'}[0];
	    my $qkhspend   = $hspkeep->{$khsp}->{'query'}[1];
	    my $hkhspstart = $hspkeep->{$khsp}->{'hit'}[0];
	    my $hkhspend   = $hspkeep->{$khsp}->{'hit'}[1];
	    
	    if ( ( (($qstart < $qkhspstart) and ($qend < $qkhspstart)) or
		   (($qstart > $qkhspend) and ($qend > $qkhspend)) ) and
		 ( (($hstart < $hkhspstart) and ($hend < $hkhspstart)) or 
		   (($hstart > $hkhspend) and ($hend > $hkhspend)) ) ){
		next;
	    }
	    else {
		$signal++;
	    }	}
	
	# if non-overlapping, include in $hspkeep
	if ($signal == 0){
	    $hspalnkeep->{$hsp_count} = $hspaln;
	    $hspkeep->{$hsp_count}->{'query'} = [$qstart, $qend];
	    $hspkeep->{$hsp_count}->{'hit'}   = [$hstart, $hend];
	    $hspqstart->{$hsp_count} = $qstart;
	}
    }
    return ($hspalnkeep, $hspkeep, $hspqstart);
}

