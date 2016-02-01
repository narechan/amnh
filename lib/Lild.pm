package Lild;

=head1 NAME

Lild.pm - contains methods for the LILD pipeline

=head1 SYNOPSIS



=head1 DESCRIPTION



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

#####TODO: Look into constraint that contains the root: disallow?
#####      Look into multiple outgroups
#####      Document

use strict;
use warnings;
use Bio::TreeIO;
use Bio::AlignIO;
use Bio::SeqIO;

our $VERSION = "0.01";

sub new {
    my ($class, %args) = @_;
    
    # create the object with room for
    # a tree, and aln, and the config
    my $self = bless {
	'tree'          => undef,
	'configuration' => undef,
	'alignment'     => undef,
	
    }, $class;
 
    # if any args, store those
    foreach my $arg (keys %args){
	$self->{$arg} = $args{$arg};
    }
    
    return $self;
}
 
sub load_config{
    my $self = shift;
    my $file = shift;

    # get and store key/val pairs in config
    my $config = {};
    open (F, "$file") or
	die "Unable to open config file\n";
    while (my $line = <F>){
        chomp $line;
	my ($key, $value) = 
	    ($line =~ m/^\s*(\w+)\s*=\s*(.*)$/);
        $config->{$key} = $value;
    }
    close (F);

    # check for PAUP in user's path if not in config                                                 
    if (! exists($config->{'PAUP'})){
	my $paup = `which paup`;
	chomp $paup;

        if (defined ($paup)){
            $config->{'PAUP'} = $paup;
        }
        else {
            die "No access to PAUP!\n";
	}
    }

    $self->{'configuration'} = $config;
}

sub load_aln{
    my $self      = shift;
    my $alignfile = shift;

    open (NEX, "$alignfile");
    my $charset = {};
    my $nexus   = {};
    while (my $line = <NEX>){
	chomp $line;

        # take only first instances of all of these things                                     
        # header only                                                                 
        ($nexus->{'nchar'}    = $1) if (($line =~m/nchar\s*=\s*(\d+)/i) and (!$nexus->{'nchar'}));
        ($nexus->{'ntax'}     = $1) if (($line =~m/ntax\s*=\s*(\d+)/i) and (!$nexus->{'ntax'}));
        ($nexus->{'datatype'} = $1) if (($line =~m/datatype\s*=\s*(\w+)/i) and (!$nexus->{'datatype'}));
        ($nexus->{'missing'}  = $1) if (($line =~m/missing\s*=\s*(.{1})/i) and (!$nexus->{'missing'}));
        ($nexus->{'gap'}      = $1) if (($line =~m/gap\s*=\s*(.{1})/i) and (!$nexus->{'gap'}));

	if ($line =~m/outgroup/i){
	    $line =~s/outgroup//ig;
	    $line =~s/\s+//g;
	    $line =~s/\;//g;
	    
	    # any instances of more than one outgroup???? <====FIX
	    $nexus->{'outgroup'} = $line; 
	}
	
	if ($line =~m/charset/i){
	    my ($partition, $coords) =
                split (/\=/, $line);

#	    print STDERR "$line\n";
	    $partition =~s/charset\s*//ig;
	    $partition =~s/\s+//g;

	    $coords =~s/\;//g;
	    $coords =~s/^\s//;

	    $charset->{$partition} = $coords;
	}
	$nexus->{'charset'} = $charset;
    }
    close (NEX);
    
    $self->{'alignment'} = $nexus;
}

sub load_tree{
    my $self     = shift;
    my $treefile = shift;

    # generate the tree object                                                            
    my $in_tree  = Bio::TreeIO->new(-file   => "$treefile",
                                    -format => "nexus");
    my $tree     = $in_tree->next_tree;

    # get all nodes and internal nodes                                                 
    my @allnodes      = $tree->get_nodes;
    my @internalnodes = grep {!$_->is_Leaf} @allnodes;

    my $lineages = {};
    my $internal_counter = 0;
    foreach my $node (@internalnodes){
        $internal_counter++;
        my @descendents = $node->get_all_Descendents;
        my @children    = grep{$_->is_Leaf} @descendents;

	foreach my $leaf (@children){
            my $id = $leaf->id;
            push @{$lineages->{$internal_counter}}, $id;
        }
    }

    $self->{'tree'} = $lineages;
}

sub get_charsets{
    my $self = shift;
    return $self->{'alignment'}->{'charset'};
}

sub get_leaves{
    my $self = shift;
    return $self->{'tree'};
}

sub get_nchar{
    my $self = shift;
    return $self->{'alignment'}->{'nchar'};
}

sub build_lild{
    my $self   = shift;
    my $expt   = shift;
    my $matrix = shift;
    my $outdir_cmds = shift;
    my $outdir_logs = shift;
    
    open (E, "$expt");
    my $partition;
    my $coords;
    my $node;
    my $constraintstring;
    while (my $line = <E>){
	chomp $line;
	
	($partition, $coords, $node, $constraintstring) = 
	    split (/\t/, $line);
    }
    close (E);

    my @parts;
    open (LILD, ">$outdir_cmds/$partition.$node.nxs");
	    
    # header
    print LILD "#NEXUS\n";
    print LILD "BEGIN PAUP;\n";
    print LILD "SET maxtrees=100;\n";
    print LILD "SET autoinc=100;\n";
    print LILD "SET increase=auto;\n";
	    
    # execute the matrix
    print LILD "EXECUTE $matrix;\n";

    # print the constraints
    print LILD "CONSTRAINTS Node$node=(($constraintstring));\n";
	    
    # print the test charset
    print LILD "CHARSET TEST=$coords;\n";
    push (@parts, "TEST");
	    
    # generate and print the random partitions
    # right now LILD does not support complex charsets
    # with multiple ranges or combinations of ranges and
    # individual characers. Must be a single range, or a
    # space delimited list of characters

    # if the coords are a range, do this:
    my $partLen;
    if ($coords =~m/\-/){
	my ($start, $end) = split (/\-/, $coords);
	$partLen       = $end - $start + 1;
    }
    
    # otherwise they are a list of positions 
    else {
	my @coords = split (/\s/, $coords);
	$partLen = @coords;
    }

    for (my $sample = 1; $sample <= $self->{'configuration'}->{SAMPLES}; $sample++) {
		
	# inititalize for this partition
	my %population;
	my @sample;
	my $sample_idx = 0;

	# Get each element of this sample.
	for (my $sample_el = 1; $sample_el <= $partLen; $sample_el++) {

	    # Get random integer among all chars in the matrix
	    my $random_idx = int(rand($self->{'alignment'}->{'nchar'})) + 1;

	    #see if the random number chosen overlaps with an existing member. if so,
	    #retry until we get a good random number we haven't seen before
	    if (exists ($population{$random_idx})){
		($random_idx = int(rand($self->{'alignment'}->{'nchar'})) + 1) until
		    (! exists($population{$random_idx}));
	    }
	    $population{$random_idx} = 1;
	    push (@sample, $random_idx);
	}

	@sample = sort {$a <=> $b} @sample;
	my $randstring = join (" ", @sample);
	print LILD "CHARSET RANDPART$sample=$randstring;\n";
	push (@parts, "RANDPART$sample");

    }
	    
    # start the log
    print LILD "LOG start file=$outdir_logs/$partition.$node.log;\n";

    # searches foreach partition
    foreach my $part (@parts){
	print LILD "[!==>partition: $part]\n";
	print LILD "EXCLUDE all;\n";
	print LILD "INCLUDE $part;\n";
#	print LILD "EXCLUDE uninf;\n";
	
	if ($part =~m/TEST/){
	    print LILD "$self->{'configuration'}->{TESTSEARCH} enforce=no;\n";
	    print LILD "[!==>end_record]\n";
	    print LILD "$self->{'configuration'}->{TESTSEARCH} enforce=yes constraints=Node$node;\n";
	    print LILD "[!==>end_record]\n";
	}
	else {
	    print LILD "$self->{'configuration'}->{RANDSEARCH} enforce=no;\n";
	    print LILD "[!==>end_record]\n";
	    print LILD "$self->{'configuration'}->{RANDSEARCH} enforce=yes constraints=Node$node;\n";
	    print LILD "[!==>end_record]\n";
	}
    }
    
    # finish up
    print LILD "LOG stop;\n";
    print LILD "END;\n";
    close (LILD);
    
    return ($partition, $node);
}

sub parse_lild{
    my $self   = shift;
    my $file   = shift;
    
    my $constraint;
    my $length;
    my $partition;
    my $data = {};
    open (F, "$file");
    while (my $line = <F>){
	chomp $line;

        # mine the data we need                                                         
        ($constraint = 0) if ($line =~m/not\s*enforced/i);
	($constraint = 1) if ($line =~m/constraint-tree/i);
	($length     = $1) if ($line =~m/found\s*=\s*(\d+)/i);
	($partition  = $1) if ($line =~m/==>partition:\s*(\w+)/i);

	# look for end of record marker
	if ($line =~m/==>end_record/){
	    $data->{$partition}->{$constraint} = $length;
	}
	
	# check for PAUP errors
	if ($line =~m/Execution terminated due to errors/){
	    $data = {};
	    last;
	}
    }

    close (F);
    return ($data);
}

sub run_paup{
    my $self   = shift;
    my $outdir_cmds = shift;
    my $part   = shift;
    my $node   = shift;
    
    `$self->{'configuration'}->{PAUP} -n $outdir_cmds/$part.$node.nxs`;
    
}

sub generate_partitions{
    my $self   = shift;
    my $window = shift;
    my $motion = shift;
    my $st     = shift;
    my $nchar  = shift;

    my %parts;
    my $counter = 0;
    my $i = $st;
    while (1){
        $counter++;
        my $start = $i;
        my $end   = $i + ($window - 1);

        my $coords = $start . "-" . $end;

        if ($end >= $nchar){
            $end = $nchar;
            $coords = $start . "-" . $end;
            $parts{$counter} = $coords;
            last;
        }

        $parts{$counter} = $coords;
        $i += $motion;
    }

    return (\%parts);
}


1;
