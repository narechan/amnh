package Bio::Phylo::Util::IDPool;
use strict;
{
    my @reclaim;
    my $obj_counter = 1;

    sub _initialize {
        my $obj_ID = 0;
        if (@reclaim) {
            $obj_ID = shift(@reclaim);
        }
        else {
            $obj_ID = $obj_counter;
            $obj_counter++;
        }
        return \$obj_ID;
    }

    sub _reclaim {
        my ( $class, $obj ) = @_;

        #        push @reclaim, $obj->get_id;
    }
}
1;
__END__

=head1 NAME

Bio::Phylo::Util::IDPool - Utility class for generating object IDs. No serviceable parts inside.

=head1 DESCRIPTION

This package defines utility functions for generating and reclaiming object
IDs. These functions are called by object constructors and destructors,
respectively. There is no direct usage.

=head1 SEE ALSO

=over

=item L<Bio::Phylo::Manual>

Also see the manual: L<Bio::Phylo::Manual> and L<http://rutgervos.blogspot.com>.

=back

=head1 CITATION

If you use Bio::Phylo in published research, please cite it:

B<Rutger A Vos>, B<Jason Caravas>, B<Klaas Hartmann>, B<Mark A Jensen>
and B<Chase Miller>, 2011. Bio::Phylo - phyloinformatic analysis using Perl.
I<BMC Bioinformatics> B<12>:63.
L<http://dx.doi.org/10.1186/1471-2105-12-63>

=head1 REVISION

 $Id: IDPool.pm 1593 2011-02-27 15:26:04Z rvos $

=cut
