#!/usr/bin/perl
# This file analyzes a single HCV sequence for genotype and recombination
# Requires access to:
# 1. muscle, to generate multiple sequence alignment with refseqs of all genotypes,
# 2. dnsdist, to generate distance matrix from phylip-formatted alignment (DNA only), and
# 3. fastme, to generate phylogeny tree from distance matrix
#use strict;
use warnings;
use File::Spec::Functions;
use Bio::SeqIO;
use POSIX qw(strftime);

use English;
use Carp;
use Data::Dumper;
use File::Temp qw/ tempfile tempdir /;
use Parallel::ForkManager; # To install this module: 1)perl -MCPAN -e shell; 2)force install Parallel::ForkManager

use version; our $VERSION = qv('2.0.0'); # Sep. 09, 2016
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
use Genotype;
use Genotype_util;

my $debug   = 0;

# Program locations, may leave blank if the programs are accessible from the prompt
# Default for fastme is fastme_linux32
my $progs = Genotype_Def::getProgs();

# Get user-defined options
my $infile        = '';
#my $outfile       = '';
my $dir_path      = '';
my $recomb        = '';
my $taxon         = '';
my $no_rerun      = '';
my $nproc         = 1;

my $exe_dir  = './';
my $exe_name = $0;
if ($exe_name =~ /^(.*[\/])([^\/]+[.]pl)$/i) {
    $exe_dir  = $1;
    $exe_name = $2;
} else {
    $exe_dir  = `pwd`;
}
if (1) {
print STDERR "$exe_name: $0 $VERSION executing... ".POSIX::strftime("%m/%d/%Y %H:%M:%S", localtime)."\n";
} else {
print STDERR "$exe_name: $0 $VERSION executing...\n"; # Need to keep the time off for comparison with earlier result
}
my $argv = [ @ARGV ];
print STDERR "$exe_name: command='$0 @$argv'\n";
my $useropts = GetOptions(
                 "recomb" => \$recomb,        # Recombination job
                 "no_rerun"  => \$no_rerun,   # use existing alignment, and generate newick file
                 "i=s"    => \$infile,        # [inputFile.gbk] in full path
                 "d=s"    => \$dir_path,      # directory of the input file
                 "t=s"    => \$taxon,         # taxon id of the species, HCV, DENG, or other
                 "debug"  => \$debug,         # Turn on debug in all subroutines
                 "nproc=s"  => \$nproc,         # number of processes, for parallel processing
#                 "o=s"    => \$outfile,       # outfile in full path
                 );
print STDERR "$exe_name: \$recomb='$recomb'  \$taxon='$taxon' \$infile='$infile' \$debug='$debug'\n";
$nproc = 0 if ($nproc<0);
$nproc = 20 if ($nproc>20);
print STDERR "$exe_name: \$nproc=$nproc\n";
$dir_path = './' if (!$dir_path && $infile);
$dir_path = abs_path($dir_path) . '/';
#$outfile = $dir_path.'/'.$outfile if ($dir_path && $outfile && $outfile !~ /[\/]/);
#print STDERR "$exe_name: \$outfile='$outfile'\n";

if ($debug) {
    Genotype_util::setDebugAll( $debug);
    my $msgs = Genotype_util::getPerlVersions;
    for (@$msgs) {
        print STDERR "$exe_name: $_\n";
    }
    print STDERR "\n";
}
$debug   = 0;

############# Subroutines #############

sub Usage{
    print STDERR "Usage: ./vipr_genotype.pl -t HCV -d ./ -i AB119282.faa\n";
    print STDERR "Usage: ./vipr_genotype.pl -t HCV -d ./\n";
    print STDERR "Usage: ./vipr_genotype.pl -t HCV -recomb -d ./\n";
    print STDERR "Usage: ./vipr_genotype.pl -recomb -d test -i EU155347.gb\n";
    my $taxons = Genotype_Def::get_taxon_allowed();
    print STDERR "Usage: Allowed taxon are:\n @$taxons\n";
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
            &Usage() && exit(0);
        }

    } elsif ($dir_path) {
        # If there is no specific input file, but has a directory
        $dir_path = "$dir_path";
        if (!-d $dir_path) {
            print STDERR "$exe_name: \$dir_path=$dir_path, can't find the directory, abort.\n";
            &Usage() && exit(0);
        }
        # First, look for any genbank file
        my $ptn = '^\s*([^\s]+)(\.(gb|gbk|genbank))\s*$';
        $accs = Genotype_util::list_dir_files( "$dir_path", $ptn);
        my $num = $#{$accs}+1;
        print STDERR "$exe_name: from directory: $dir_path, found $num gbk files.\n";
        if ($#{$accs}<0) {
            # If there is no genbank file, look for any fasta file
            $ptn = '^\s*([^\s]+)(\.(fa|faa|fasta))\s*$';
            $accs = Genotype_util::list_dir_files( "$dir_path", $ptn);
            $num = $#{$accs}+1;
            print STDERR "$exe_name: from directory: $dir_path, found $num fasta files.\n";
        }
        for my $j (0 .. $#{$accs}) {
            print STDERR "$exe_name: \$acc[$j]='$accs->[$j]->[0]' '$accs->[$j]->[1]'\n";
        }
        print STDERR "\n";
#        $debug && print STDERR "$exe_name: \$accs = \n".Dumper($accs)."End of \$accs\n\n";
    }
    $debug && print STDERR "$exe_name: \$accs = \n".Dumper($accs)."End of \$accs\n\n";
    for (@$accs) { print STDERR "$exe_name: \$accs='$_->[0]'\t'$_->[1]'\n"; }

    # Input genbank/fasta file is required; Quit if absent
    ($#{$accs}<0) && &Usage() && exit(0);

    # Set up parallel process, each input gbk/fasta file is one fork
    # Set $max_procs = 0; for debugging, ie. not forking
    # Set $max_procs to (intended number of process)-1 for processing
    my $max_procs = ($nproc>0) ? $nproc : 0;
    my $pm =  Parallel::ForkManager->new( $max_procs);
    $debug && print STDERR "$exe_name: ForkManager started with \$max_procs=$max_procs\n";
    $pm->run_on_start(
        sub { my ($pid, $ident) = @_;
            print STDERR "** $ident started, pid: $pid\n";
        }
    );
    $pm->run_on_finish(
        sub { my ($pid, $exit_code, $ident) = @_;
          if (!$exit_code) {
            print STDERR "** $ident finished PID $pid and exit_code: $exit_code\n";
          } else {
            print STDERR "** $ident ERROR: finished PID $pid and exit_code: $exit_code\n";
            print STDERR "** $ident ERROR: needs checking\n";
          }
        }
    );
    my $tally = { 0 => 0  }; # Key is the number of the file, value is the number of seqs in that file

    # Process each input gbk/fasta file
    for my $i (0 .. $#{$accs}) {
        # starts the fork in child, and returns to for loop in main
        my $pid = $pm->start("$i:$accs->[$i]->[1]") and next;

        my $infile = $accs->[$i]->[1];
        my $job_name = $accs->[$i]->[1];
        $job_name =~ s/[.]\w{1,7}$//; # remove the file name extension

        my $msgs_all = [];
        my $str = ($infile) ? $infile : 'undef';
        print STDERR "\n===================================================\n";
        print STDERR "$exe_name: Processing \$i=$i genome \$infile='$str'\n";

        my $finish_fn = "${job_name}_finished.txt";
        if (-f "$dir_path/$finish_fn") {
            print STDERR "$exe_name: Found existing file $dir_path/$finish_fn before starting, deleting it\n";
            `rm $dir_path/$finish_fn`;
        }

        # Separate multiple fasta records in 1 file into separate files
        my $faas = Genotype_util::get_faas( $infile, $job_name, $dir_path, $taxon);
        for (@$faas) {
            my $s1 = ($_->{accession}) ? $_->{accession} : 'undef';
            my $s2 = ($_->{file}) ? $_->{file} : 'undef';
            print STDERR "$exe_name: accession=$_->{long_name}\tname=$_->{job_name} file=$s1 long_name=$s2\n";
        }

#        my $tempdir = File::Temp->newdir();
        my $tempdir = tempdir( CLEANUP => 1);
        $debug && print STDERR "$exe_name: \$tempdir=$tempdir\n";
        # save result to $outfile
        my $OUTF;

        # Process each sequence within each file
        my $index_fn = '';
        my $INDEXF;
        my $iresult = 0;
        for my $j (0 .. $#{$faas}) {
            my $msgs = [];
            my $faa = $faas->[$j];
            $debug && print STDERR "$exe_name: \$j=$j \$faa=\n".Dumper($faa)."End of \$faa\n";
            $job_name = $faa->{job_name};
            $job_name = 'test' if (!$job_name);
            my $acc = $faa->{accession};
            print STDERR "\n$exe_name: Processing #$j accession=$acc long_name=$faa->{long_name} infile=$faa->{file}\n";
            if ($j==0) {
                # Open a file for summary of results, called <>_index.tsv, requested by NG web team
                $index_fn = "${job_name}_index.tsv";
                print STDERR "$exe_name: Processing \$i=$i \$index_fn=$index_fn\n";
                open $INDEXF, '>', "$dir_path/$index_fn" or croak "Can't open outfile '$index_fn': $OS_ERROR";
                print $INDEXF "#\tJob_name\tshort_name\tdefline\tspecies\tregion";
                print $INDEXF "\tstatus(geno)\tBI(geno)\tgenotype\tsubtype";
                print $INDEXF "\tstatus(recomb)\tBI(recomb)\trecombinant\tsubtype\tcomment(geno)\n";
            }

            my $genotypes = Genotype::genotype( $faa, $tempdir, $dir_path, $progs, $no_rerun, $recomb, 0);

            # Convert resuling array of hashes to text
            for my $g (@$genotypes) { my $msg=Genotype_Def::to_String( $g); print STDERR "'$msg'\n"; push @$msgs, $msg; }

            $debug && print STDERR "$exe_name: \$genotypes\n".Dumper($genotypes)."End of \$genotypes\n";
            if ($INDEXF) {
              for my $ii (0 .. $#{$genotypes}) {
                next if ($genotypes->[$ii]->{type} !~ m/^\s*genotype\s*$/i );
                $iresult++;
                printf $INDEXF ("%d\t$job_name\t$faa->{accession}\t", $iresult);
                print $INDEXF "$genotypes->[$ii]->{long_name}\t";
                print $INDEXF $genotypes->[$ii]->{taxon} ? "$genotypes->[$ii]->{taxon}\t" : "\t";
                print $INDEXF "$genotypes->[$ii]->{region}\t";
                print $INDEXF "$genotypes->[$ii]->{status}\t";
                print $INDEXF "$genotypes->[$ii]->{BI}\t";
                print $INDEXF ($genotypes->[$ii+1] && $genotypes->[$ii+1]->{type} eq 'recomb_summ') ? "$genotypes->[$ii+1]->{status}\t" : "Skip\t";
                print $INDEXF ($genotypes->[$ii+1] && $genotypes->[$ii+1]->{type} eq 'recomb_summ') ? "$genotypes->[$ii+1]->{BI}" : "BI=\tgenotype=NA\tsub=NA";
                print $INDEXF "\t$genotypes->[$ii]->{comment}";
                print $INDEXF "\n";
              }
            } else {
                print STDERR "$exe_name: \$i=$i \$j=$j index file not open\n";
            }
            # save result to outfile
            my $outfilename;
            # Save genotype result to <job_name>_genotype_summary.tsv
            if ($genotypes && $#{$genotypes}>=0) {
                $outfilename = "$dir_path/${job_name}${acc}_genotype_summary.tsv";
                print STDERR "$exe_name: genotype result: \$outfilename=$outfilename\n";
                $OUTF = undef;
                open $OUTF, '>', $outfilename or croak "Can't open outfile '$outfilename': $OS_ERROR";
#                print $OUTF "$0 @$argv\n";
#                print $OUTF "Summary:\n";
                  print $OUTF "Job_name\tJob_type\tregion\tstart..end\tBranch_Index\tgenotype\tsubtype\ttaxon";
                  #print $OUTF "\tDate(ymdhms)";
                  print $OUTF "\tStatus\tComment\n";
                for my $ii (0 .. $#{$msgs}) {
                    next if ($msgs->[$ii] !~ /\tgenotype\t/i);
                    print $OUTF "$msgs->[$ii]\n";
                }
                print $OUTF "\n";
                close $OUTF or croak "Can't close outfile '$outfilename': $OS_ERROR";
            }

            if ($j==$#{$faas}) {
                # close summary of results, called <>_index.tsv, requested by NG web team
                print STDERR "$exe_name: Closing \$index_fn=$index_fn\n";
                close $INDEXF or croak "Can't close outfile '$index_fn': $OS_ERROR";
                $debug && Genotype_util::check_file("$dir_path/$index_fn");
            }
            $tally->{$i}++;
        } # for my $j (0 .. $#{$faas})

        $tally->{files}++;
        # Create a file <>_finished.txt to indicated the finish of the script
        open my $FINISHF, '>', "$dir_path/$finish_fn" or croak "Can't open outfile '$finish_fn': $OS_ERROR";
        printf $FINISHF ("$exe_name: Processed $tally->{$i} seq(s) in File #%d: $infile.\n", $i+1);
        print $FINISHF "$exe_name: Finished for file $infile.\n";
        close $FINISHF or croak "Can't close outfile '$finish_fn': $OS_ERROR";
        print STDERR "$exe_name: Created \$finish_fn=$finish_fn\n";
        $debug && Genotype_util::check_file("$dir_path/$finish_fn");

        $pm->finish($job_name); # do the exit in the child process, pass an exit code to finish
    } # for my $i (0 .. $#{$accs})
    $pm->wait_all_children; # Wait for all child processes to finish

    $debug && print STDERR "$exe_name: \$tally\n".Dumper($tally)."End of \$tally\n";
    print STDERR "\n$exe_name: Finished ".POSIX::strftime("%m/%d/%Y %H:%M:%S", localtime).".\n";

exit;

1;
