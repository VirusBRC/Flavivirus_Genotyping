package Genotype_Draw;
# This file draws a bar graph according to the recombination analysis

use strict;
use warnings;
use POSIX qw(strftime ceil floor);
use Scalar::Util 'refaddr';

use English;
use Carp;
use Data::Dumper;
require GD::Graph::bars;
use version; our $VERSION = qv('1.1');
use Getopt::Long;

my $debug_all = 0;
my $debug = 0;

############## Subroutines ################


sub setDebugAll {
    my ($debugAll) = @_;
    $debug_all = $debugAll;
} # sub getShortSpecies


=head2 draw_bar
Draws a bar graph
=cut

sub draw_bar {
    my ( $graph_datas, $dir_path, $ALGO) = @_;

    my $debug = 0 || $debug_all;
    my $subname = 'Genotype_Draw::draw_bar';

#    $debug && print STDERR "$subname: \$graph_datas=\n".Dumper($graph_datas)."End of \$graph_datas\n";
    my $BIfileNames = [];
  for my $graph_data (@$graph_datas) {
    my $data = $graph_data->{data};
    my $legends = $graph_data->{legends};
    my $labels = $graph_data->{labels};

    my $BIfileName = '';
    my @data = ();
#    $debug && print STDERR "$subname: \$data=\n".Dumper($data)."End of \$data\n";
    if ($debug && $legends->[$#{$legends}] eq 'NA') {
        $legends = [ (@$legends)[0 .. $#{$legends}-1] ];
        @data = @data[0 .. $#{$legends}-1];
    }
    if ( $debug ) {
        print STDERR "\n$subname:\n";
        print STDERR "position\t";
        for my $legend (@$legends) {
            print STDERR "$legend\t";
        }
        print STDERR "\n";
        for my $j (0 .. $#{$data->[0]}) {
          for my $i (-1 .. $#{$legends}) {
            if (defined($data->[$i+1]->[$j])) {
                printf STDERR ($i==-1) ? "%8d\t" : "%.3f\t", "$data->[$i+1]->[$j]";
            } else {
                print STDERR "undef\t";
            }
          }
          print STDERR "\n";
        }
    }

    Genotype_Draw::trim_data( $data, $legends);

    for my $i (0 .. $#{$data}) {
        push @data, [ @{$data->[$i]} ];
    }
    if (!$legends || $#{$legends}<0) {
        for my $i (1 .. $#{$data}) {
            push @$legends, "Data$i";
        }
    }
    return if ($#data<0);

    my $my_graph = GD::Graph::bars->new(800,400);
    $my_graph->set(
        x_label => (exists($labels->{x_label})) ? $labels->{x_label} : 'Alignment position (nt)',
        y_label => (exists($labels->{y_label})) ? $labels->{y_label} : 'Branching Index score',
        title => (exists($labels->{title})) ? $labels->{title} : 'Genotype profile',
        y_max_value => 1.05,
        y_min_value => 0,
        y_tick_number => 21,
        y_label_skip => 2,
        'x_labels_vertical'=> 1,
        cumulate => 2,
        borderclrs => $my_graph->{dclrs},
#        cycle_clrs => 2,
        bar_spacing => 2,
#        shadow_depth => 4,

        transparent => 0,
    );

    $my_graph->set_legend( @$legends );
    #$my_graph->plot(\@data);

    #my $myimage = $my_graph->plot(\@data) or die $my_graph->error;
    my $myimage = $my_graph->plot($data) or die $my_graph->error;

    my $format = $my_graph->export_format;
    my $image_fn = (exists($labels->{job_name}) && $labels->{job_name}) ? "$labels->{job_name}_" : '';
    $image_fn .= "$graph_data->{region}_" if ($graph_data->{region} !~ /^\s*NA\s*$/i); # Only include region is necessary
    $image_fn .= 'BI_profile';
    $debug && print STDERR "$subname: BI profile exported to file '$image_fn'\n";
    $image_fn = "$image_fn.$format";
    $image_fn = ($dir_path) ? "$dir_path/$image_fn" : "$image_fn";
    my $pwd = `pwd`;
    chomp $pwd;
    $debug && print STDERR "$subname: \$pwd=$pwd\n";
    chomp($pwd);
    $debug && print STDERR "$subname: BI profile exported to file '$pwd/$image_fn'\n";
    $BIfileName = $image_fn;
    push @$BIfileNames, $BIfileName;
    open(IMG, ">$image_fn") or die $!;
    binmode IMG;
    print IMG $myimage->$format();
    close IMG;
    $debug && print STDERR "$subname: BI profile exported to file='$image_fn'\n";
  }

    $debug && print STDERR "$subname: leaving subroutine\n";
    return $BIfileNames;
} # sub draw_bar


=head2
  sub trim_data calculates the sum for each column in $data, then removes those beyond the first 7 columns

=cut

sub trim_data {
    my ($data, $legends) = @_;

    my $debug = 0 || $debug_all;
    my $subname = 'Genotype_Draw::trim_data';

#    $debug && print STDERR "$subname: \$data=\n".Dumper($data)."End of \$data\n";

#    return if ($#{$data} <= 6);
if ($#{$data} > 6) {

    my $sum;
    for my $i (0 .. $#{$data}) {
        my $t = 0;
        for my $j (0 .. $#{$data->[$i]}) {$t += $data->[$i]->[$j] if ($data->[$i]->[$j]); }
        $sum->[$i] = $t;
    }
    $debug && print STDERR "$subname: \$sum='@$sum'\n";
    my $threshold = [ sort {$b <=> $a} @$sum ];
    $debug && print STDERR "$subname: \$threshold='@$threshold'\n";
    $threshold = $threshold->[6];

    my $d;
    my $l;
    for my $i (0 ..  $#{$data}) {
        next if ($sum->[$i] < $threshold);
        my $col;
        for my $j (0 .. $#{$data->[$i]}) {
            $col->[$j] = $data->[$i]->[$j];
        }
        push @$d, $col;
        push @$l, $legends->[$i-1] if ($i>=1);
    }

    $_[0] = $d;
    $_[1] = $l;
    $data = $d;
    $legends = $l;
}

    # print the data in nice format
    if ($debug) {
        print STDERR "\n$subname:\n";
        print STDERR "position\t";
        for my $legend (@$legends) {
            print STDERR "$legend\t";
        }
        print STDERR "\n";
        for my $j (0 .. $#{$data->[0]}) {
          for my $i (-1 .. $#{$legends}) {
            if (defined($data->[$i+1]->[$j])) {
                printf STDERR ($i==-1)?"%8d\t":"%.3f\t", "$data->[$i+1]->[$j]";
            } else {
                print STDERR "undef\t";
            }
          }
          print STDERR "\n";
        }
    }

    $debug && print STDERR "$subname: leaving subroutine\n";
    return $debug;
} # sub trim_data

1;
