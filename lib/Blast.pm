package Blast;

=head1 NAME                                                                       
                                                                                  
Blast.pm - contains methods for running blast
                                                                                  
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

use strict;
use warnings;
use Bio::SearchIO;

sub new {
    my ($class, %args) = @_;

    # create the object with room for                                             
    # a tree, and aln, and the config                                             
    my $self = bless {
        'configuration' => undef,

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
            ($line =~ m/^\s*(.*)\s*=\s*(.*)$/);
        $config->{$key} = $value;
    }
    close (F);

    $self->{'configuration'} = $config;
}

sub generate_blastdb{
    my $self = shift;
    my $file = shift;
    my $logdir = shift;
    
    my $formatdb ="formatdb";
    ($formatdb .= " $self->{'configuration'}->{'FORMATDB'}")
	if ($self->{'configuration'}->{'FORMATDB'});

    $formatdb .= " -i $file";
    $formatdb .= " -l $logdir/formatdb.log";
    `$formatdb`;
}

sub run_blast{
    my $self = shift;
    my $fasta  = shift;
    my $db     = shift;
    my $out    = shift;

    my $fasta_name;
    if ($fasta =~/\//g){
        $fasta =~m/.*\/(.*)$/;
        $fasta_name = $1;
    }

    else {
        $fasta_name = $fasta;
    }

    my $cmd = "blastall";
    ($cmd .= " $self->{'configuration'}->{'BLAST'}")
	if ($self->{'configuration'}->{'BLAST'});
    $cmd .= " -i $fasta";
    $cmd .= " -d $db";
    $cmd .= " -o $out/$fasta_name.out";
    `$cmd`;

    return ("$fasta_name.out");
}

sub parse_blast {
    my $self = shift;
    my $indir = shift;
    my $infile = shift;
    my $outdir = shift;
    my $identity = shift;

    my $in = Bio::SearchIO->new(-file   => "$indir/$infile",
                                -format => 'blast');

    open (OUT, ">$outdir/$infile.parse");
    my $qrydata = {};
    my @refdata;
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
#           last if ($hit_count > 1); # get one hit                               

            my $hsp_count = 0;
            while (my $hsp = $hit->next_hsp ) {
                $hsp_count++;
#               last if ($hsp_count > 1); # get one HSP                           
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
                if ($hsp->frac_identical('total') >= $identity){
                    push @refdata,
                    {
                        'chrom' => $hit->name,
                        'start' => $hsp->start('hit'),
                        'end'   => $hsp->end('hit'),
                    };

                    $qrydata->{$result->query_name}->{$hit_count}++
			
		    }

            }

        }

    }

    return (\@refdata, $qrydata);
}


sub genome_coverage {
    my $self    = shift;
    my $refdata = shift;

    # sort ref data by chrom and by start                                        
    my @refdata = sort {$a->{chrom} cmp $b->{chrom} ||
                            $a->{start} <=> $b->{start}
                    } @$refdata;

    # get the first element                                                      
    my $first  = shift @refdata;
    my $chrom1 = $first->{chrom};
    my $start1 = $first->{start};
    my $end1   = $first->{end};

    # do joins on overlapping hits to get the                                    
    # non redundant coverage                                                     
    my $coverage = {};
    my $counter  = 0;
    foreach my $point (@refdata){

        # get data                                                              
	my $chrom2 = $point->{chrom};
        my $start2 = $point->{start};
        my $end2   = $point->{end};

        if ($chrom1 eq $chrom2){
            if (($start1 > $end2) or ($end1 < $start2)) {
		$counter++;

                # they don't overlap at all                                     
                $coverage->{$chrom1}->{$counter} = [$start1, $end1];
                $chrom1 = $chrom2;
                $start1 = $start2;
		$end1   = $end2;
            }
            else {
		my @coor = sort {$a <=> $b} ($start1, $start2, $end1, $end2);

#		print STDERR "Overlap\t$start1\t$end1\t$start2\t$end2\n";
		
                # get the extremes as the new range                         
                $start1 = shift @coor;
                $end1   = pop @coor;
            }
	}
        else {
	    # put in the last values                                
	    $coverage->{$chrom1}->{$counter + 1} = [$start1, $end1];
            
	    # redefine for new chrom
	    $counter = 0;
	    $chrom1 = $point->{chrom};
	    $start1 = $point->{start};
	    $end1   = $point->{end};
        }

    }
    # put in the last values                                                   
    $coverage->{$chrom1}->{$counter + 1} = [$start1, $end1];
    
    return ($coverage);

}

1;
