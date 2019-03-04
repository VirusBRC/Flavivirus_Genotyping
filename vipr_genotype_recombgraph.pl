#!/usr/bin/perl
# This script takes the output file from a recombination analysis in tab-delimited text file,
# and draws a graph according to the BI values
# Requirements:
# 1. The GD library
# 2. Perl module GD::Graph;
use strict;
use warnings;
use File::Spec::Functions;
#use Bio::SeqIO;
use POSIX qw(strftime);

use English;
use Carp;
use Data::Dumper;
#use File::Temp qw/ tempfile tempdir /;

use version; our $VERSION = qv('1.2.2'); # Feb. 08, 2013
use Cwd 'abs_path'; # needed to find absolute path for input file
use File::Basename;
my $libPath = './';
BEGIN { # The modules need to exist in same dir as the script
    $libPath = (-l __FILE__) ? dirname(readlink(__FILE__)) : dirname(__FILE__);
    $libPath = abs_path($libPath);
}
use lib ("$libPath");

use Getopt::Long;
use Genotype_Def;
use Genotype_recomb;
use Genotype_Draw;

my $debug   = 0;

# Get user-defined options
my $infile        = '';
my $dir_path      = '';

my $exe_dir  = './';
my $exe_name = $0;
if ($exe_name =~ /^(.*[\/])([^\/]+[.]pl)$/i) {
    ($exe_dir, $exe_name)  = ($1, $2);
} else {
    $exe_dir  = `pwd`;
}
print STDERR "\n";
print STDERR "$exe_name: $0 $VERSION executing...\n";
print STDERR "$exe_name: command='$0 @ARGV'\n";
my $useropts = GetOptions(
                 "i=s"    => \$infile,        # [inputFile.gbk] in full path
                 "d=s"    => \$dir_path,      # directory of the input file
                 "debug"  => \$debug,         # Turn on debug in all subroutines
#                 "o=s"    => \$outfile,       # outfile in full path
                 );
print STDERR "$exe_name: \$infile='$dir_path/$infile' \$debug='$debug'\n";
$dir_path = './' if (!$dir_path && $infile);
$dir_path = abs_path($dir_path) . '/';
#$outfile = $dir_path.'/'.$outfile if ($dir_path && $outfile && $outfile !~ /[\/]/);
#print STDERR "$exe_name: \$outfile='$outfile'\n";

if ($debug) {
    Genotype_Draw::setDebugAll( $debug);
}

############# Subroutines #############

sub Usage{
    my $msg = "Usage: ./vipr_genotype_recombgraph.pl -d ./ -i AB119282_seq001_recombination_summary.tsv\n";
    return $msg;
}


############# Main Program #############

    # First, look for the .bgk or fasta files to be processed
    my $accs = []; # Stores the list of files that need processing, not actually processing them
    if ($infile) {
        # If there is a simple requested input file
        if (-f "$dir_path/$infile") {
            push @$accs, [$#{$accs}+1, "$infile"];
#            $debug && print STDERR "$exe_name: \$accs = \n".Dumper($accs)."End of \$accs\n\n";
        } else {
            print STDERR "$exe_name: \$dir_path=$dir_path, \$infile=$infile, can't find such file, abort.\n";
        }
    }
    $debug && print STDERR "$exe_name: \$accs = \n".Dumper($accs)."End of \$accs\n\n";
    for (@$accs) { print STDERR "$exe_name: \$accs='$_->[0]'\t'$_->[1]'\n"; }

    # Input genbank file is required; Quit if absent
    if ($#{$accs}<0) {
        print &Usage();
        exit(0);
    }

    # Process each gbk or fasta file
    my $ALGO = '';
    for my $i (0 .. $#{$accs}) {
        my $msgs_all = [];
        my $infile = $accs->[$i]->[1];
        my $job_name = $accs->[$i]->[1];
        $job_name =~ s/[.]\w{1,7}$//; # remove the file name extension
        my $pos = index($job_name, '_');
        $pos = index($job_name, '_seq'); # underscore alone '_' causes problem in NC_004102
        if ($pos>=0) {
            $debug && print STDERR "$exe_name: \$job_name=$job_name \$pos=$pos\n";
            $pos = index($job_name, '_', $pos+1);
            $debug && print STDERR "$exe_name: \$job_name=$job_name \$pos=$pos\n";
            $job_name = substr($job_name, 0, $pos) if ($pos>0);
        }
        my $str = ($infile) ? $infile : 'undef';
        print STDERR "===================================================\n";
        print STDERR "$exe_name: Processing \$i=$i genome \$infile='$str'\n";

       # Process each sequence
#        my $genotypes = Genotype::genotype( $taxon, $faa, $tempdir, $dir_path, $progs, $no_rerun, $recomb, $debug);
        my $genotypes = [];
        open my $INF, '<', "$dir_path/$infile" or croak "Can't open infile '$dir_path/$infile': $OS_ERROR";
        while (my $line = <$INF>) {
            if ($line =~ m/^# Algorithm=(\w+)/) {
                $ALGO = $1;
                $debug && print STDERR "$exe_name: ALGO=$ALGO\n";
                next;
            } elsif ($line =~ m/^#/) {
                next;
            }
            $debug && print STDERR "$exe_name: \$line=$line\n";
            my $genotype = Genotype_Def::parseGenotypeString( $line, $ALGO);
            $debug && print STDERR "$exe_name: \$genotype=\n".Dumper($genotype)."End of \$genotype\n";
            push(@$genotypes, $genotype) if ($genotype);
        }
        close $INF or croak "Can't close infile '$dir_path/$infile': $OS_ERROR";
        print STDERR "$exe_name: loaded $#{$genotypes} genotypes from '$dir_path/$infile'\n";
        $debug && print STDERR "$exe_name: \$genotypes=\n".Dumper([@$genotypes[0..5]])."End of \$genotypes\n";

        # Draw a bar graph
        my $graph_data = Genotype_recomb::get_graph_data( $genotypes, $job_name, $dir_path, $ALGO);
        $debug && print STDERR "$exe_name: \$job_name='$job_name'\n";
        my $BIfileNames = Genotype_Draw::draw_bar( $graph_data, $dir_path);

        for my $BIfileName (@$BIfileNames) {
            print STDERR "$exe_name: BI profile exported to file='$BIfileName'\n";
        }

    } # for my $i (0 .. $#{$accs})

    print STDERR "\n$exe_name: Finished.\n";

exit;

1;
