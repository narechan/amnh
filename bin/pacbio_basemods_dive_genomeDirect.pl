#!/usr/bin/perl -w

#####SETUP#####
use strict;
use warnings;
use Chart::Gnuplot;
use Getopt::Long;
use Pod::Usage;

my ($help, $motif, $outdir, $slope, $covcut, $csvfile, $methpos, $cutoff, $fastafile);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fastafile, #your assembly as a fasta
    'm|motif=s'       => \$motif, #only forward string required
    'o|outdir=s'      => \$outdir, #output lives here
    's|slope=s'       => \$slope, #indexing the cutoff to converage
    'c|coverage=s'    => \$covcut, #minimum coverage required
    'x|csv=s'         => \$csvfile, #pacbio ipdr raw file
    'y|methpos=s'     => \$methpos, #comma delimited string of meth positions: (f,r); one sided: (f).
    'z|cutoff=s'      => \$cutoff, #the ipdr rev/fw cutoff at which a position is designated unmethylated
	   ) or pod2usage;
pod2usage if $help;

`mkdir -p $outdir`;

#####MAIN#####

### parse the csv file
print STDERR "Parsing CSV file\n";
my $csv = {};
open (C, "$csvfile");
while (my $line = <C>){
    chomp $line;
        my ($chrom,
            $start,
            $strand,
            $base,
            $score,
            $tmean,
            $terr,
            $modpred,
            $ipdr,
            $coverage) = split (/\,/, $line);

    $csv->{$start}->{$strand}->{'score'} = $score;
    $csv->{$start}->{$strand}->{'ipdr'} = $ipdr;
    $csv->{$start}->{$strand}->{'coverage'} = $coverage;
}

### mine the genome for the motif of interest (only forward string required)
### and populate the ipdr datastruc for f/r ipdr plots in double motifs
### (if passes coverage screen)
print STDERR "Mining the genome for $motif\n";
my $methyls = {};
my $ipdrplot = {};
my $fpos;
my $rpos;
my $onesided = 0;
if ($methpos =~m/\,/){
    ($fpos, $rpos) = split (/\,/, $methpos);
}
else{
    $fpos = $methpos;
    $onesided = 1;
}

# get the genome into a string
my $genomestring;
open (F, "$fastafile");
while (my $line = <F>){
    chomp $line;
    next if ($line =~m/^>/);
    $genomestring .= $line;
}

# get it's reverse complement into a string
my $rcgenomestring = revdnacomp ($genomestring);

# index the RC position with respect to the forward position
my $rcposindex = {};
my $lenrcgenomestring = length ($rcgenomestring);
my $leniterator = $lenrcgenomestring;
for (my $i = 1; $i <= $lenrcgenomestring; $i++){
    $rcposindex->{$i} = $leniterator;
    $leniterator--;
}

# match to the forward strand
while ($genomestring =~m/($motif)/g) {
    my $pattern  = $1;
    my $prematch = $`;
    my $pos      = length($prematch);
    $methyls->{$pos+$fpos}->{1} = $pattern;
    ($methyls->{$pos+$rpos}->{0} = revdnacomp ($pattern)) unless ($onesided == 1);

    # consolidate the ipdrs for a motif and describe as f or r
    # if forward and reverse pass the coverage filter
    unless ($onesided == 1){
	if (($csv->{$pos+$fpos}->{1}->{'coverage'}>=$covcut) and ($csv->{$pos+$rpos}->{0}->{'coverage'}>=$covcut)){
	    $ipdrplot->{$pos+$fpos}->{1} = $csv->{$pos+$fpos}->{0}->{'ipdr'}; #rev bc fw std is 0 in csv!!!
	    $ipdrplot->{$pos+$fpos}->{0} = $csv->{$pos+$rpos}->{1}->{'ipdr'};
	}
	else{
	    next;
	}
    }
}

# match to the RC strand
while ($rcgenomestring =~m/($motif)/g) {
    my $pattern  = $1;
    my $prematch = $`;
    my $pos      = $rcposindex->{length($prematch)};
    $methyls->{$pos-$fpos}->{0} = $pattern;
    ($methyls->{$pos-$rpos}->{1} = revdnacomp ($pattern)) unless ($onesided == 1);
    
    # consolidate and filter
    unless ($onesided == 1){
	if (($csv->{$pos-$rpos}->{1}->{'coverage'}>=$covcut) and ($csv->{$pos-$fpos}->{0}->{'coverage'}>=$covcut)){
	    $ipdrplot->{$pos-$rpos}->{1} = $csv->{$pos-$rpos}->{0}->{'ipdr'};
	    $ipdrplot->{$pos-$rpos}->{0} = $csv->{$pos-$fpos}->{1}->{'ipdr'};
	}
	else {
	    next;
	}
    }
}

### Analyze basemods
print STDERR "Analyzing $motif basemods\n";

# globals for basemods and plots
my $table = {};
my @covspass;
my @covsfail;
my $covmax = 0;
my @ipdrspass;
my @ipdrsfail;
my $ipdrmax = 0;
my @qualspass;
my @qualsfail;
my $qualmax = 0;
foreach my $pos (sort {$a <=> $b} keys %$methyls){
    foreach my $trustr (sort keys %{$methyls->{$pos}}){
#	print STDERR "$pos\t$trustr\n";
	
	# to deal with the fact that 0/1 is flipped in the csv file
	my $str;
	if ($trustr == 1){
	    $str = 0;
	}
	else{
	    $str = 1;
	}
	
	# check to see if there's any coverage
	unless ($csv->{$pos}->{$str}->{'coverage'}){
	    $csv->{$pos}->{$str}->{'coverage'} = 0;
	    $csv->{$pos}->{$str}->{'score'} = 0;
	    $csv->{$pos}->{$str}->{'ipdr'} = 0;
	}
	
	# get all the data to output in one large table
	$table->{$pos}->{$trustr} = [$csv->{$pos}->{$str}->{'coverage'},
				     $csv->{$pos}->{$str}->{'score'},
				     $csv->{$pos}->{$str}->{'ipdr'}];
	
	
	# update shared maxima for pass and fail
	if ($covmax < $csv->{$pos}->{$str}->{'coverage'}){
	    $covmax = $csv->{$pos}->{$str}->{'coverage'};
	}
	if ($qualmax < $csv->{$pos}->{$str}->{'score'}){
	    $qualmax = $csv->{$pos}->{$str}->{'score'};
	}
	if ($ipdrmax < $csv->{$pos}->{$str}->{'ipdr'}){
	    $ipdrmax = $csv->{$pos}->{$str}->{'ipdr'};
	}

	# coverage screen
	unless ($csv->{$pos}->{$str}->{'coverage'} >= $covcut){
	    push (@covsfail, $csv->{$pos}->{$str}->{'coverage'});
	    push (@qualsfail, $csv->{$pos}->{$str}->{'score'});
	    push (@ipdrsfail, $csv->{$pos}->{$str}->{'ipdr'});
	    next;
	}
	    
	# slope based quality screen
	my $expectedScore = $slope * $csv->{$pos}->{$str}->{'coverage'};
	if ($csv->{$pos}->{$str}->{'score'} >= $expectedScore){
	    push (@covspass, $csv->{$pos}->{$str}->{'coverage'});
	    push (@qualspass, $csv->{$pos}->{$str}->{'score'});
	    push (@ipdrspass, $csv->{$pos}->{$str}->{'ipdr'});
	}
	else {
	    push (@covsfail, $csv->{$pos}->{$str}->{'coverage'});
	    push (@qualsfail, $csv->{$pos}->{$str}->{'score'});
	    push (@ipdrsfail, $csv->{$pos}->{$str}->{'ipdr'});
	}
    }
}

# print bulk output table
open (T, ">$outdir/motifs.data.all");
foreach my $p (sort {$a <=> $b} keys %$table){
    foreach my $s (sort {$a <=> $b} keys %{$table->{$p}}){
	my @data = @{$table->{$p}->{$s}};
	my $data = join "\t", @data;
	print T "$p\t$s\t$data\n";
    }
}
close (T);

# print out the coverage vs qv cloud for this motif
my $chart = Chart::Gnuplot->new(output => "$outdir/$motif.covVSqual.gif",
				xrange=>[0, $covmax + 1],
				yrange=>[0, $qualmax + 1],
				xlabel=>"Coverage",
				ylabel=>"QV",
				title=>"$csvfile: $motif",
				bg=>"white");
my $datapass = Chart::Gnuplot::DataSet->new(xdata=>\@covspass,
					    ydata=>\@qualspass,
					    style=>"points",
					    pointtype=>"fill-circle",
					    pointsize=>0.5);
my $datafail = Chart::Gnuplot::DataSet->new(xdata=>\@covsfail,
					    ydata=>\@qualsfail,
					    style=>"points",
					    pointtype=>"fill-circle",
					    pointsize=>0.5);

# print out the coverage vs ipdr cloud for this motif                                                     
my $chartipdr = Chart::Gnuplot->new(output => "$outdir/$motif.covVSipdr.gif",
				    xrange=>[0, $covmax + 1],
				    yrange=>[0, $ipdrmax + 1],
				    xlabel=>"Coverage",
				    ylabel=>"IPDRatio",
				    title=>"$csvfile: $motif",
				    bg=>"white");
my $datapassipdr = Chart::Gnuplot::DataSet->new(xdata=>\@covspass,
						ydata=>\@ipdrspass,
						style=>"points",
						pointtype=>"fill-circle",
						pointsize=>0.5);
my $datafailipdr = Chart::Gnuplot::DataSet->new(xdata=>\@covsfail,
						ydata=>\@ipdrsfail,
						style=>"points",
						pointtype=>"fill-circle",
						pointsize=>0.5);

# quals charts
if ((@qualspass) and (@qualsfail)){
    $chart->plot2d($datapass, $datafail);
}
elsif (@qualspass){
    $chart->plot2d($datapass);
}
elsif (@qualsfail){
    $chart->plot2d($datafail);
}
else {
    print STDERR "No data for QV charts\n";
}

# ipdrs charts
if ((@ipdrspass) and (@ipdrsfail)){
    $chartipdr->plot2d($datapassipdr, $datafailipdr);
}
elsif (@ipdrspass){
    $chartipdr->plot2d($datapassipdr);
}
elsif (@ipdrsfail){
    $chartipdr->plot2d($datafailipdr);
}
else {
    print STDERR "No data for IPDR charts\n";
}

# print stats on pass and fail observations to stderr
my $covspass = @covspass;
my $covsfail = @covsfail;

print STDERR "$motif pass: $covspass\n";
print STDERR "$motif fail: $covsfail\n";

# print out the coverage f/r ipdr cloud for this motif pair
# if required
if ($rpos){
    open (I, ">$outdir/motifs.data.low");
    my @forward;
    my $forwardmax = 0;
    my @reverse;
    my $reversemax = 0;
    foreach my $ipdr (keys %$ipdrplot){
	my $signal = 0;
	if ($ipdrplot->{$ipdr}->{1}){
	    push (@forward, $ipdrplot->{$ipdr}->{1});
	    if ($ipdrplot->{$ipdr}->{1} > $forwardmax){
		$forwardmax = $ipdrplot->{$ipdr}->{1};
	    }
	    if ($ipdrplot->{$ipdr}->{1} <= $cutoff){
		$signal++;
	    }
	}
	else {
	    print STDERR "IPDR 1 $ipdr error\n";
	}
	if ($ipdrplot->{$ipdr}->{0}){
	    push (@reverse, $ipdrplot->{$ipdr}->{0});
	    if ($ipdrplot->{$ipdr}->{0} > $reversemax){
		$reversemax = $ipdrplot->{$ipdr}->{0};
	    }
	    if ($ipdrplot->{$ipdr}->{0} <= $cutoff){
		$signal++;
	    }
	}
	else {
	    print STDERR "IPDR 0 $ipdr error\n";
	}
	(print I "$ipdr\n") if ($signal == 2);
    }
    
    my $ipdrp = Chart::Gnuplot->new(output => "$outdir/ipdr_fr_plot.gif",
				    xrange=>[0, $forwardmax + 1],
				    yrange=>[0, $reversemax + 1],
				    xlabel=>"IPDratio Forward Strand",
				    ylabel=>"IPDratio Reverse Strand",
				    title=>"$csvfile: $motif",
				    bg=>"white");
    my $dataipdrp = Chart::Gnuplot::DataSet->new(xdata=>\@forward,
						 ydata=>\@reverse,
						 style=>"points",
						 pointtype=>"fill-circle",
						 pointsize=>0.5);
    $ipdrp->plot2d($dataipdrp);
    close (I);
}


sub revdnacomp {
    my $dna = shift;
    my $revcomp = reverse($dna);

    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}
