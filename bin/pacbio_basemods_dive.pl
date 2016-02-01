#!/usr/bin/perl -w

#####SETUP#####
use strict;
use warnings;
use Chart::Gnuplot;
use Getopt::Long;
use Pod::Usage;

my ($help, $gff3, $motifst, $outdir, $slope, $covcut, $csvfile, $separation, $cutoff);
GetOptions(
    'h|help'          => \$help,
    'g|gff3=s'        => \$gff3,
    'm|motif=s'       => \$motifst, #you can enter the motif as a comma delimited pair if you want f/r IPD graph
    'o|outdir=s'      => \$outdir,
    's|slope=s'       => \$slope,
    'c|coverage=s'    => \$covcut,
    'x|csv=s'         => \$csvfile,
    'y|separation=s'  => \$separation, #the distance between modified bases of a pair (if you have 2 motifs)
    'z|cutoff=s'      => \$cutoff, #the ipdr rev/fw cutoff at which a position is designated unmethylated
	   ) or pod2usage;
pod2usage if $help;

`mkdir -p $outdir`;

#####MAIN#####

# var for IPD ratio f/r plots
my $ipdrplot = {};

my @motif = split (/\,/, $motifst);
foreach my $motif (@motif){
    print STDERR "Working on $motif\n";
    
    # parse the csv file                                                                                        
    print STDERR "Parsing $motif CSV file\n";
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
    
    # parse the gff file
    print STDERR "Parsing $motif GFF file\n";
    my $lines   = {};
    my $gff3str = {};
    open (G, "$gff3");
    while (my $line = <G>){
	next if ($line =~m/^\#/);
	chomp $line;
	
	my ($chrom,
	    $source,
	    $type,
	    $start,
	    $end,
	    $score,
	    $strand,
	    $phase,
	    $attr) = split (/\t/, $line);
	
	my $str;
	if ($strand eq "+"){
	    $str = 0;
	}
	else {
	    $str = 1;
	}
	
	$lines->{$start}->{$str} = $line;
	$gff3str->{$start}->{$str}->{'type'}   = $type;
	$gff3str->{$start}->{$str}->{'score'}  = $score;
	
	my @attrs = split (/\;/, $attr);
	foreach my $att (@attrs){
	    my ($key, $val) = split (/\=/, $att);
	    $gff3str->{$start}->{$str}->{$key} = $val;
	}
	unless (exists ($gff3str->{$start}->{$str}->{'identificationQv'})){
	    $gff3str->{$start}->{$str}->{'identificationQv'} = 0;
	}
    }
    close (G);
    
    # analyze the basemods
    print STDERR "Analyzing $motif basemods\n";
    my @covspass;
    my @covsfail;
    my $covmax = 0;
    my @ipdrspass;
    my @ipdrsfail;
    my $ipdrmax = 0;
    my @qualspass;
    my @qualsfail;
    my $qualmax = 0;
    foreach my $pos (sort {$a <=> $b} keys %$gff3str){
	foreach my $str (sort keys %{$gff3str->{$pos}}){
	    
	    # screen everything but our motif of interest
	    next unless (exists ($gff3str->{$pos}->{$str}->{'motif'}));
	    next unless ($gff3str->{$pos}->{$str}->{'motif'} eq $motif);
	    
#	print "$lines->{$pos}->{$str}\t$csv->{$pos}->{$str}->{'score'}\t$csv->{$pos}->{$str}->{'coverage'}\n";
	    
	    # coverage screen
	    unless ($gff3str->{$pos}->{$str}->{'coverage'} >= $covcut){
		push (@covsfail, $gff3str->{$pos}->{$str}->{'coverage'});
		push (@qualsfail, $gff3str->{$pos}->{$str}->{'score'});
		push (@ipdrsfail, $csv->{$pos}->{$str}->{'ipdr'});
		next;
	    }
	    
	    if ($motifst =~m/\,/){
		# after coverage screen and before quality screen,
		# compile ipdrs for forward and reverse motifs
		# given the separation and the motif identities
		my $addpos = $pos + $separation;
		my $subpos = $pos - $separation;
		my $oppstr;
		if ($str == 0){
		    $oppstr = 1;
		}
		else{
		    $oppstr = 0;
		}
		
		my $addposmotif;
		my $subposmotif;
		my $posmotif;
		my $signal = 0;
		if (exists ($ipdrplot->{$addpos})){
		    if (exists ($gff3str->{$addpos}->{$oppstr}->{'motif'})){
			$addposmotif = $gff3str->{$addpos}->{$oppstr}->{'motif'};
			$posmotif    = $gff3str->{$pos}->{$str}->{'motif'};
			if (($addposmotif ne $posmotif) and ($motifst =~m/$addposmotif/) and ($motifst =~m/$posmotif/)){
			    $ipdrplot->{$addpos}->{$str} = $csv->{$pos}->{$str}->{'ipdr'};
			    $signal++;
			}
		    }
		}
		if (exists ($ipdrplot->{$subpos})){
		    if (exists ($gff3str->{$subpos}->{$oppstr}->{'motif'})){
			$subposmotif = $gff3str->{$subpos}->{$oppstr}->{'motif'};
			$posmotif    = $gff3str->{$pos}->{$str}->{'motif'};
			if (($subposmotif ne $posmotif) and ($motifst =~m/$subposmotif/) and ($motifst =~m/$posmotif/)){
			    $ipdrplot->{$subpos}->{$str} = $csv->{$pos}->{$str}->{'ipdr'};
			    $signal++;
			}
		    }
		}
		if ($signal == 0){
		    $ipdrplot->{$pos}->{$str} = $csv->{$pos}->{$str}->{'ipdr'};
		}
	    }

	    # slope based quality screen
	    my $expectedScore = $slope * $gff3str->{$pos}->{$str}->{'coverage'};
	    if ($gff3str->{$pos}->{$str}->{'score'} >= $expectedScore){
#           if ($gff3str->{$pos}->{'score'} >= 30){
		push (@covspass, $gff3str->{$pos}->{$str}->{'coverage'});
		if ($covmax < $gff3str->{$pos}->{$str}->{'coverage'}){
		    $covmax = $gff3str->{$pos}->{$str}->{'coverage'};
		}
		push (@qualspass, $gff3str->{$pos}->{$str}->{'score'});
		if ($qualmax < $gff3str->{$pos}->{$str}->{'score'}){
		    $qualmax = $gff3str->{$pos}->{$str}->{'score'};
		}
		push (@ipdrspass, $csv->{$pos}->{$str}->{'ipdr'});
		if ($ipdrmax < $csv->{$pos}->{$str}->{'ipdr'}){
		    $ipdrmax = $csv->{$pos}->{$str}->{'ipdr'};
		}
	    }
	    else {
		push (@covsfail, $gff3str->{$pos}->{$str}->{'coverage'});
		push (@qualsfail, $gff3str->{$pos}->{$str}->{'score'});
		push (@ipdrsfail, $csv->{$pos}->{$str}->{'ipdr'});
	    }
	}
    }
    
    # print out the coverage vs qv cloud for this motif
    my $chart = Chart::Gnuplot->new(output => "$outdir/$motif.covVSqual.gif",
				    xrange=>[0, $covmax + 1],
				    yrange=>[0, $qualmax + 1],
				    xlabel=>"Coverage",
				    ylabel=>"QV",
				    title=>"$gff3: $motif",
				    bg=>"white");
#				    size=>'square 2');
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
					title=>"$gff3: $motif",
					bg=>"white");
#					size=>'square 2');
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

}

# print out the coverage f/r ipdr cloud for this motif pair
# if required
if ($motifst =~m/\,/){
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
	(print "$ipdr\n") if ($signal == 2);
    }
    
    my $ipdrp = Chart::Gnuplot->new(output => "$outdir/ipdr_fr_plot.gif",
				    xrange=>[0, $forwardmax + 1],
				    yrange=>[0, $reversemax + 1],
				    xlabel=>"IPDratio Forward Strand",
				    ylabel=>"IPDratio Reverse Strand",
				    title=>"$gff3: $motifst",
				    bg=>"white");
#				size=>'square 2');
    my $dataipdrp = Chart::Gnuplot::DataSet->new(xdata=>\@forward,
						 ydata=>\@reverse,
						 style=>"points",
						 pointtype=>"fill-circle",
						 pointsize=>0.5);
    $ipdrp->plot2d($dataipdrp);
}
