#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('f:', \%opts);
my $file  = $opts{'f'};

my $counter = 0;
my $variants = {};
my $paralogs = {};
open (F, "$file");
while (my $line = <F>){
    $counter++;
    chomp $line;

    # parse the line                                                                           
    my @line = split (/\t/, $line);
    my $comp  = $line[0];
    my $query = $line[2];
    my $hit   = $line[5];
    my $desig = $line[9];
    my $string = $query . ";" . $hit;
    
    # sort the variants and paralogs
    if ($desig eq "N"){
	$variants->{$counter} = $string;
    }
    else {
	$paralogs->{$counter} = $string;
    }
}
close (F);

# grow the variant groups; note that these must be sorted
# in the input file to allow for a single pass through
my $groups = {};
my $gcount = 0;
foreach my $vpair (sort {$a <=> $b} keys %$variants){
#foreach my $vpair (sort keys %$variants){
    my ($query, $hit) = split (/\;/, $variants->{$vpair});
    
    # cycle through groups and add to a group if one of a pair
    # exists in the group
    my $reporter = 0;
    foreach my $group (sort {$a <=> $b} keys %$groups){
	if ( (exists ($groups->{$group}->{$query})) or (exists ($groups->{$group}->{$hit})) ){
	    $groups->{$group}->{$query} = 1;
	    $groups->{$group}->{$hit} = 1;
	    $reporter++;

	    print STDERR "V: $vpair\t$query\t$hit\t$group existing\n";
	    last;
	}
	else {
	    next;
	}
    }

    # if the pair does not exist in the groups so far,
    # create a new one
    if ($reporter > 0){
	next;
    }
    else {
	$gcount++;
	$groups->{$gcount}->{$query} = 1;
	$groups->{$gcount}->{$hit} = 1;
	print STDERR "V: $vpair\t$query\t$hit\t$gcount created\n";
	
    }
}

# check the paralogs for consistency with variant groups
# and to add singletons
foreach my $para (sort {$a <=> $b} keys %$paralogs){
    my ($query, $hit) = split (/\;/, $paralogs->{$para});

    # find the group numbers of the paralogs, if any
    my $grpcacheqry;
    my $grpcachehit;
    foreach my $group (sort {$a <=> $b} keys %$groups){
	foreach my $member (keys %{$groups->{$group}}){
	    if ($member eq $query){
		$grpcacheqry = $group;
	    }
	    elsif ($member eq $hit){
		$grpcachehit = $group;
	    }
	    else {
		next;
	    }
	}
    }

    # if neither of the paralogs are defined in any group,
    # create a group for both
    if ( (!$grpcacheqry) and (!$grpcachehit) ){
	$gcount++;
        $groups->{$gcount}->{$query} = 1;
	$gcount++;
        $groups->{$gcount}->{$hit} = 1;
	print STDERR "P: $para\t$query\t$hit\tAdd group $query Add group $hit\n";
    }
    
    # if the paralogs are in the same group --> violation
    elsif ($grpcacheqry == $grpcachehit){
	print STDERR "P: $para\t$query\t$hit\tViolation\n";
    }
    
    # if the paralogs are in different groups --> pass
    elsif ( ($grpcacheqry) and ($grpcachehit) and ($grpcacheqry != $grpcachehit) ){
	print STDERR "P: $para\t$query\t$hit\tPass\n";
    }
    
    # if one of the paralogs doesn't exist in a group,
    # create a group for it
    elsif (!$grpcacheqry){
	$gcount++;
	$groups->{$gcount}->{$query} = 1;
	print STDERR "P: $para\t$query\t$hit\tAdd group $query\n";
    }    
    elsif (!$grpcachehit){
        $gcount++;
        $groups->{$gcount}->{$hit} = 1;
	print STDERR "P: $para\t$query\t$hit\tAdd group $hit\n";
    }
    
    # die if there is an unexpected scenario
    else {
	print STDERR "P: What the hell just happened?\n";
	die;
    }
}

# print out the groups
foreach my $group (sort {$a <=> $b} keys %$groups){
    print "$group\t";
    
    my @members;
    foreach my $member (keys %{$groups->{$group}}){
	push (@members, $member);
    }
    my $members_string = join ";", @members;
    
    print "$members_string\n";
}
















    
=head
    # first check to see whether paralogs are in the same                                                        
    # group ==> violation                                                                                        
    foreach my $group (sort {$a <=> $b} keys %$groups){
        if ( (exists ($groups->{$group}->{$query})) and (exists ($groups->{$group}->{$hit})) ){
            print "Violation\n";
        }
    }




    push (@lines, $line);
}
close (F);

# cluster vars
my $clusters = {};
my $clusterct = 1;

# prime the clustering from the first line
my $firstline = unshift @lines;

# parse the line
my @line = split (/\t/, $line);
my $comp  = $line[0];
my $query = $line[2];
my $hit   = $line[5];
my $desig = $line[9];

# if they are paralogs "Y", split out into sep clusters                                                     
# or put into same cluster if they are splice variants "N"
if ($desig eq "Y"){
    $clusters->{$clusterct} = $query;
    $clusters->{$clusterct++} = $hit;
}
else {
    $clusters->{$clusterct} = [$query, $hit];
}

# now cycle through the rest of the comp group
foreach my $line (@lines){
    my @line = split (/\t/, $line);
    my $comp  = $line[0];
    my $query = $line[2];
    my $hit   = $line[5];
    my $desig = $line[9];

    if ($desig eq "Y"){
	$clusters->{$clusterct++} = $query;
	$clusters->{$clusterct++} = $hit;
    }
    else {
	$clusters->{$clusterct} = [$query, $hit];
    }

my $clusters = {};
my $clusterct = 1;
open (F, "$file");
while (my $line = <F>){
    chomp $line;

    my @line = split (/\t/, $line);
    my $comp  = $line[0];
    my $query = $line[2];
    my $hit   = $line[5];
    my $desig = $line[9];

    # if they are paralogs, split out into sep clusters
    if ($desig eq "Y"){
	$clusters->{$clusterct} = $query;
	$clusters->{$clusterct++} = $hit;
	
    
}
close (F);
=cut
