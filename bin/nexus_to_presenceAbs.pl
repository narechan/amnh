#!/usr/bin/perl -w

## Builds a presence absence matrix

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use TreeSupports;
use Bio::AlignIO;
use Bio::SeqIO;

my ($help, $matrixfile, $listfile1, $listfile2);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    'a|list1=s'       => \$listfile1,
    'b|list2=s'       => \$listfile2,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($matrixfile){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# instantiate the object and load stuff we need
my $supportobj = TreeSupports->new;
$supportobj->load_aln ($matrixfile);

# get the charsets
my $charsets = $supportobj->get_charsets;
my $chars    = $supportobj->get_nchar;

# store alignment information
my $pa = {};
my $partdist = {};
my $lengths    = {};
my $alnin = Bio::AlignIO->new(-file   => "$matrixfile",
			      -format => "nexus");

# store the lists
my $list1 = {};
open (L1, "$listfile1");
while (my $line = <L1>){
    chomp $line;
    my ($acc, $cat) = split (/\t/, $line);
    $list1->{$acc} = 1;
}
close (L1);

my $list2 = {};
open (L2, "$listfile2");
while (my $line = <L2>){
    chomp $line;
    my ($acc, $cat) = split (/\t/, $line);
    $list2->{$acc} = 1;
}
close (L2);


# set up taxa counter for the header partition array
my @taxa;
my $tcounter = 0;

# only one aln there
my $aln = $alnin->next_aln();
foreach my $seq ($aln->each_seq){
    my $id = $seq->display_id;
    my $tcounter++;

    foreach my $charset (sort keys %$charsets){
#	print STDERR "Processsing\t$id\t$charset\n";
	(push (@taxa, $charset)) if ($tcounter == 1);

	# get the partition coordinates
	my $coords = $charsets->{$charset};
	my ($start, $end) = split (/\-/, $coords);
	
	# harvest the sequence and get its length
	my $partition = $seq->subseq($start, $end);
	my $partlen   = length ($partition);

	# check of the sequence is all missing data
	# and assign as present or absent
	if ($partition =~m/\?{$partlen}/){
#	    push (@{$pa->{$id}}, 0);
	    $pa->{$id}->{$charset} = 0;
	    $partdist->{$charset}->{'missing'}++;
	}
	else{
#	    push (@{$pa->{$id}}, 1);
	    $pa->{$id}->{$charset} = 1;

	    # check if it's in the either list
	    if (exists($list1->{$id})){
		$partdist->{$charset}->{$listfile1}++;
	    }
	    elsif (exists($list2->{$id})){
                $partdist->{$charset}->{$listfile2}++;
            }
	    else{
		$partdist->{$charset}->{'other'}++;
	    }
	}
    }
}


#	    if (($id eq "PE001") or ($id eq "PE002") or ($id eq "PE003") or ($id eq "PE004") or ($id eq "PE005")){
#		$partdist->{$charset}->{'low'}++;
#	    }
#	    elsif (($id eq "PE006") or ($id eq "PE007") or ($id eq "PE008") or ($id eq "PE009") or ($id eq "PE010")){
#		$partdist->{$charset}->{'high'}++;
#	    }
#	    else {
#		$partdist->{$charset}->{'other'}++;
#	    }	
#	}
#    }
#}

# header
#my $taxastring = join "\t", @taxa;
#print "\t$taxastring\n";

my $charsetstring = join "\t", sort keys %$charsets;
print "\t$charsetstring\n";

# spit out pa matrix
#foreach my $tax (sort keys %$pa){
#    my @paarray = @{$pa->{$tax}};
#    my $pastring = join "\t", @paarray;
#    print "$tax\t$pastring\n";
#}
foreach my $tax (sort keys %$pa){
    my @paa;
    foreach my $char (sort keys %{$pa->{$tax}}){
	push (@paa, $pa->{$tax}->{$char})
    }
    my $paast = join "\t", @paa;
    print "$tax\t$paast\n";
}


# print the stderr all the list bins
foreach my $gene (sort keys %$partdist){                                                                      
    my @cats;                                                                                                 
    foreach my $cat ($listfile1, $listfile2, "other", "missing"){                                   
        if ($partdist->{$gene}->{$cat}){                                                                       
            push (@cats, $partdist->{$gene}->{$cat});                                                          
        }                                                                                                    
        else {                                                                                                
            push (@cats, 0);                                                                                
        }                                                                                                    
    }                                                                                                        
    my $cats = join "\t", @cats;                                                                            
    print STDERR "$gene\t$cats\n";                                                                        
}                                                                                                            
=cut
