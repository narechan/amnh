package Ild;

=head1 NAME

Ild.pm - contains methods for the ILD pipeline (pairwise and slideRule)

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

#####TODO:
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
	
	if ($line =~m/^charset/i){
	    $line =~s/charset//ig;
	    $line =~s/\s+//g;
	    $line =~s/\;//g;
	    
            my ($partition, $coords) =
                split (/\=/, $line);
	    $charset->{$partition} = $coords;
	}
	$nexus->{'charset'} = $charset;
    }
    close (NEX);
    
    $self->{'alignment'} = $nexus;
}

sub store_alignment_seq {
    my $self = shift;
    my $afile    = shift;
    my $charsets = shift;

    # store alignment information                                                                    
    my $partitions = {};
    my $lengths    = {};
    my $alnin = Bio::AlignIO->new(-file   => "$afile",
                                  -format => "nexus");

   # only one aln there                                                                               
    my $aln = $alnin->next_aln();
    foreach my $seq ($aln->each_seq){
        my $id        = $seq->display_id;

        foreach my $charset (sort keys %$charsets){
            my $coords = $charsets->{$charset};
            my ($start, $end) = split (/\-/, $coords);

            my $partition = $seq->subseq($start, $end);
            my $partlen   = length ($partition);

            $partitions->{$charset}->{$id} = $partition;
            $lengths->{$charset} = $partlen;
        }
    }
    
    return ($partitions, $lengths);
}

sub generate_nxs_pairwise{
    my $self = shift;
    my $seq1 = shift;
    my $seq2 = shift;
    my $len1 = shift;
    my $len2 = shift;
    my $part1 = shift;
    my $part2 = shift;
    my $outnxs = shift;

    my $start = 1;
    my $end   = 0;
    my $sets  = {};
    my $lcounter = 0;
    for my $len ($len1, $len2){
	$lcounter++;
        $end = $start + ($len - 1);
	$sets->{$lcounter} = "$start-$end";
        $start = $end + 1;
    }

    my $party = $part1 . "v" . $part2;
    open (C, ">$outnxs/$party.nxs");
    print C "#NEXUS\n";
    print C "BEGIN DATA;\n";
    print C "DIMENSIONS NTAX=$self->{'alignment'}->{'ntax'} NCHAR=$end;\n";
    print C "FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" DATATYPE=PROTEIN MISSING=? GAP=-;\n";
    print C "MATRIX\n\n";

    for my $seq ($seq1, $seq2){
	foreach my $id (sort keys %$seq){
	    print C "$id\t$seq->{$id}\n";
	}
	print C "\n";
    }
    
    print C ";\n";
    print C "END;\n\n";

    my $pcounter = 0;
    print C "BEGIN SETS;\n";
    for my $part ($part1, $part2){
	$pcounter++;
	print C "CHARSET $part = $sets->{$pcounter};\n";
    }
    print C "END;\n";

    close (C);
}

sub get_charsets{
    my $self = shift;
    return $self->{'alignment'}->{'charset'};
}

sub get_nchar{
    my $self = shift;
    return $self->{'alignment'}->{'nchar'};
}

sub run_paup{
    my $self   = shift;
    my $exptfile   = shift;
    
    `paup -n $exptfile`;
    
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

sub parse_treelens{
    my $self = shift;
    my $file = shift;
    
    my $lens = {};
    my $numtrees = {};

    my $part;
    my $ntrees;
    my $length;

    open (F, "$file");
    while (my $line = <F>){
	chomp $line;
	($part       = $1) if ($line =~m/==>ILD:\s*(\w+)/i);
	($ntrees     = $1) if ($line =~m/trees retained\s*=\s*(\d+)/i);
        ($length     = $1) if ($line =~m/found\s*=\s*(\d+)/i);

	# look for record ender
	if ($line =~m/==>end_record/){
            $lens->{$part}     = $length;
	    $numtrees->{$part} = $ntrees;
        }
	else {
	    next;
	}
    }
    close (F);
    return ($lens, $numtrees);
}

	
sub generate_ildstat_singlegene{
    my $self       = shift;
    my $counter    = shift;
    my $charset    = shift;
    my $outdir_cmd = shift;
    my $outdir_log = shift;
    my $alignfile  = shift;

    open (ILD, ">$outdir_cmd/$counter.nex");

    # print out the header                                                    
    print ILD "#NEXUS\n";
    print ILD "SET increase = auto;\n";
    print ILD "EXECUTE $alignfile;\n";
    print ILD "LOG START file = $outdir_log/$counter.log replace = yes;\n";
    
    # create the expt
    print ILD "EXCLUDE all;\n";
    print ILD "INCLUDE $charset;\n";

    print ILD "[!==>ILD: $charset]\n";
    print ILD "$self->{'configuration'}->{'TREECOMMAND'};\n";
    print ILD "[!==>end_record]\n";
    
    print ILD "LOG STOP;\n";
    print ILD "QUIT /warnTsave=no;\n";
}
    
sub generate_ildstat_batch_pw{
    my $self       = shift;
    my $counter    = shift;
    my $char1      = shift;
    my $charfile      = shift;
    my $outdir_cmd = shift;
    my $outdir_log = shift;
    my $alignfile  = shift;

    open (ILD, ">$outdir_cmd/$counter.nex");

    # print out the header                                            
    print ILD "#NEXUS\n";
    print ILD "SET increase = auto;\n";
    print ILD "EXECUTE $alignfile;\n";
    print ILD "LOG START file = $outdir_log/$counter.log replace = yes;\n";

    # cycle through the partitions to batch
    foreach my $char2 (sort keys %$charfile){
	next if ($char2 eq $char1);
	
	my $charpart = $char1 . "v" . $char2;
	
	# create the experiment                                            
	print ILD "EXCLUDE all;\n";
	print ILD "INCLUDE $char1 $char2;\n";
	
	print ILD "[!==>ILD: $charpart]\n";
	print ILD "$self->{'configuration'}->{'TREECOMMAND'};\n";
	print ILD "[!==>end_record]\n";
    }
    print ILD "LOG STOP;\n";
    print ILD "QUIT /warnTsave=no;\n";
}


sub generate_ild_part{
    my $self       = shift;
    my $counter    = shift;
    my $char1      = shift;
    my $char2      = shift;
    my $outdir_cmd = shift;
    my $outdir_log = shift;
    my $alignfile  = shift;

#    open (ILD, ">$outdir_cmd/$char1-$char2.nex");
    open (ILD, ">$outdir_cmd/$counter.nex");
    my $charpart = $char1 . "v" . $char2;

    # print out the header                                                                           
    print ILD "#NEXUS\n";
    print ILD "SET increase = auto;\n";
    print ILD "EXECUTE $alignfile;\n";
#    print ILD "LOG START file = $outdir_log/$char1-$char2.log replace = yes;\n";
    print ILD "LOG START file = $outdir_log/$counter.log replace = yes;\n";

    # redefining in case there are partitions with mult charsets
    my $ch1 = $char1;
    my $ch2 = $char2;
    $ch1 =~s/\_/ /g;
    $ch2 =~s/\_/ /g;

    # create the experiment              
    print ILD "EXCLUDE all;\n";
    print ILD "INCLUDE $ch1 $ch2;\n";
    print ILD "BEGIN PAUP;\n";
    print ILD "[!==>ILD: $charpart]\n";
    print ILD "CHARPARTITION $charpart = $char1:$ch1, $char2:$ch2;\n";
    print ILD "$self->{'configuration'}->{'HOMPART'} ";
    print ILD "partition = $charpart";

    if ($self->{'configuration'}->{'SEARCH'}){
        print ILD " / $self->{'configuration'}->{'SEARCH'};\n";
    }
    else{
        print ILD ";\n";
    }

    print ILD "[!==>end_record]\n";

    print ILD "END;\n";
    print ILD "LOG STOP;\n";
    print ILD "QUIT /warnTsave=no;\n";
}

sub generate_ild_sr{
    my $self       = shift;
    my $start      = shift;
    my $end        = shift;
    my $outdir_cmd = shift;
    my $outdir_log = shift;
    my $nchar      = shift;
    my $partition  = shift;
    my $alignfile  = shift;

    open (ILD, ">$outdir_cmd/$partition.nex");

    # print out the header                                       
    print ILD "#NEXUS\n";
    print ILD "SET increase = auto;\n";
    print ILD "EXECUTE $alignfile;\n";
    print ILD "LOG START file = $outdir_log/$partition.log replace = yes;\n";

    # create the partition set                                                    
    print ILD "BEGIN SETS;\n";
    print ILD "CHARSET C$partition = ";
    print ILD "$start-$end;\n";
    print ILD "END;\n";

    # create the experiment                                                       
    print ILD "BEGIN PAUP;\n";

    my $fivestart  = 1;
    my $fiveend    = $start - 1;

    my $threestart = $end + 1;
    my $threeend   = $nchar;

    print ILD "[!==>ILD: C$partition]\n";
    if ($fiveend == 0){
        print ILD "CHARPARTITION CharPart$partition = ";
        print ILD "1:C$partition, 2:$threestart-$threeend;\n";
    }
    elsif ($threestart > $nchar){
	print ILD "CHARPARTITION CharPart$partition = ";
        print ILD "1:C$partition, 2:$fivestart-$fiveend;\n";
    }
    else{
        print ILD "CHARPARTITION CharPart$partition = ";
        print ILD "1:C$partition, 2:$fivestart-$fiveend $threestart-$threeend;\n";
    }

    print ILD "$self->{'configuration'}->{HOMPART} ";
    print ILD "partition = CharPart$partition";

    if ($self->{'configuration'}->{SEARCH}){
        print ILD " / $self->{'configuration'}->{SEARCH};\n";
    }
    else{
        print ILD ";\n";
    }

    print ILD "[!==>end_record]\n";

    print ILD "END;\n";
    print ILD "LOG STOP;\n";
    print ILD "QUIT /warnTsave=no;\n";
}

sub parse_ild{
    my $self    = shift;
    my $file    = shift;

    my $pval;
    my $len;
    my $ildchars;
    my $stump;
    my $signal = 0;

#    print STDERR "$file\n";
    
    open (ILD, "$file");
    while (my $line = <ILD>){
        chomp $line;
	
	if (($signal == 1) and ($line)){ #at end of $ildchars have blank line
	    $ildchars .= $line;
	}
	elsif (($signal == 1) and !($line)){ #shut off signal when blank
	    $signal = 0;
	}
	else {
	    if ($line =~m/\s*P\s{1}value\s*\=.*\=\s*(.+)/){
		$pval = $1;
	    }
	    if ($line =~m/\s*(\d+)\*\s*\d+/){
		$len  = $1;
	    }
	    if ($line =~m/\=\=\>ILD\:/){
		$line =~s/\s//g;
		($stump, $ildchars) = split (/\:/, $line);
		unless ($ildchars){ #hack for PAUP putting shit on the next line!
		    $signal = 1;
		}
	    }
	}
    }
    close (ILD);
    
    $len =~s/\*//;
    return ($pval, $len, $ildchars);
}


1;
