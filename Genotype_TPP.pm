package Genotype_TPP;

# This file analyzes a single HCV sequence for genotype and recombination
# Requires access to:
# 1. MUSCLE, to generate multiple sequence alignment with refseqs of all genotypes,
# 2. dnsdist, to generate distance matrix from phylip-formatted alignment (DNA only), and
# 3. FastME to generate phylogeny tree from distance matrix
#use strict;
use warnings;
use File::Spec::Functions;
use IO::String;
use Bio::SeqIO;
use POSIX qw(strftime ceil floor);
use Bio::TreeIO;
use Bio::AlignIO;
#use Bio::SimpleAlign;
use Scalar::Util 'refaddr';
use Bio::Tools::Run::StandAloneBlast;

use English;
use Carp;
use Data::Dumper;
use version; our $VERSION = qv('2.0.0'); # Jun 01, 2017
#use Getopt::Long;

use Annotate_Def;
use Genotype_Def;
use Genotype_newickBI;
use Genotype_util;
use Genotype_calctree;
use Genotype_recomb;
use config::Config_Reader;
use util::Utils;

my $debug_all = Genotype_util::getDebugAll();

############## Constants ################
# To add a new species:
# 0. get all refseqs for the genotype in fasta file, make sure defline is >1a.AB123456
# 1. create refseq MSA: muscle -in refseq_stlouisenvephalitis_33.fasta -out refseq_stlouisenvephalitis_33.aln -stable
# 2. add the $RMSA_species to following 2 lists: $RMSA_<species> and $RMSA hash
# 3. add any refseq to Genotype_refseq_def.txt, and define the shorthand
# 4. add any long refseq to refseq_Flaviviridae<nn>.faa, and change the file name.
#    To create database needed by blast, run: formatdb -i refseq_Flaviviridae84.faa -p F
# Note: The automatic detection of species is for convenience, but not a fool proof method, user input is preferred.
# The St Louis and Bovine diarrhea are especially difficult.

# MUSCLE/ClustalW/MAFFT
    my $ALN_PROG = 
#                   "muscle";
                   "clustalw";
#                   "mafft"; # Actual command: mafft-profile
    my $RMSA = Genotype_Def::getRMSA();

#    my $newHCV_ICTV = 1; # Sep 30, 2013
#    my $fmtGenotype = '\d[-._\dA-Za-z]*';
#    $fmtGenotype = '([0-9]\w+|\d-[A-Za-z]+|NA|.+)' if (!$newHCV_ICTV);
    my $fmtGenotype = Genotype_Def::getfmtGenotype;

############## Subroutines ################

sub getProgs {
    my ($progs) = @_;
    my $debug = 0 || $debug_all || Genotype_util::getDebugAll();
    my $subn = "Genotype_TPP::getProgs";

#    $progs = {};
    $progs->{mafft}    = "/usr/bin/mafft-profile" if (!exists($progs->{mafft}) || !$progs->{mafft});
    $progs->{taxit}    = "/usr/bin/taxit" if (!exists($progs->{taxit}) || !$progs->{taxit});
    $progs->{pplacer}  = "/home/gsun/prog/pplacer-Linux-v1.1.alpha19/pplacer" if (!exists($progs->{pplacer}) || !$progs->{pplacer});
    $progs->{guppy}    = "/home/gsun/prog/pplacer-Linux-v1.1.alpha19/guppy" if (!exists($progs->{guppy})|| !$progs->{guppy});
    $progs->{clustalw} = '/usr/bin/clustalw' if (!exists($progs->{clustalw})|| !$progs->{clustalw});
    $progs->{cladinator}= "";
    my $GENOTYPE_HOME = Genotype_Def::getGenotypeHome;
    $debug && print STDERR "$subn: cladinator='".$GENOTYPE_HOME.'/bin/forester.jar'."'\n";
    if (-e $GENOTYPE_HOME.'/bin/forester.jar') {
      $progs->{cladinator} = "java -Xmx1024m -cp $GENOTYPE_HOME/bin/forester.jar org.forester.application.cladinator";
      $progs->{cladinator} .= " -rs -c=0.7 -m=$GENOTYPE_HOME/rmsa/HCV_ICTV_mapping.tsv -S='(\\d+)([a-z?]*)_.+'";
    } else {
      print STDERR "$subn: pwd='".`pwd`."'\n";
      print STDERR "$subn: Genotype_TPP needs cladinator program, but not found\n";
      exit(1);
    }
    $debug && print STDERR "$subn: \$progs=\n".Dumper($progs)."End of \$progs\n";

    return $progs;
} # sub getProgs

##################################################################
## algo w/ Taxit/Pplacer/guppy/cladinator
sub TPP {
  my ($faa, $genotype, $classifier, $postMapping, $reportUnifiedTree, $skipTaxit) = @_;
  my $debug = 0 || $debug_all || Genotype_util::getDebugAll();
  my $subn = "Genotype_TPP::TPP";

  $debug && print STDERR "$subn: \$classifier=".ref($classifier)."\n";
#  $debug && print STDERR "$subn: \$classifier=\n".Dumper($classifier)."End of \$classifier\n";
  my $config = $classifier->{config};
  my $accession = $faa->{accession};
  my $seq = $genotype->{seqAligned};
  my $class = 'ND';
  my $probability = 0;

  # 1. taxit
  $skipTaxit = 0 if (!defined($skipTaxit));
  $debug && print STDERR "$subn: \$skipTaxit='$skipTaxit'\n";
  $debug && print STDERR "$subn: taxit w/ ref MSA='output.profile.afa', seq file='seqs.afa'\n";
  if (!$skipTaxit || -e 'HCV.refpkg') {
    $classifier->taxit($accession);
    $debug && print STDERR "$subn: After classifier->taxit 'ls -l '=". `ls -l ` ."\n";
  } else {
    print STDERR "$subn: Skipped classifier->taxit 'ls -l '=". `ls -l ` ."\n";
  }
  
  # 2 pplacer/guppy
  $debug && print STDERR "$subn: tree: accession=$accession seq=$seq\n";
#  $debug && print STDERR "$subn: \$classifier=\n".Dumper($classifier)."End of \$classifier\n";
  my ($treeFile, $treeData) = $classifier->tree( $accession, $seq );
  $debug && print STDERR "$subn: After classifier->tree 'ls -l '=". `ls -l ` ."\n";

  # 2a Change short ID to original ID
  $debug && print STDERR "$subn: before restoreSeqID tree: accession=$accession tree=$treeFile\n";
  $treeData = &restoreSeqID( $faa, $treeFile );
  

  # 3 Get classification according to lookup file
  my $fileOutput = '';
  if ( !defined($treeFile) ) {
    $debug && print STDERR "$subn: After classifier->tree for $accession, can't get tree=$treeFile\n";
    $genotype->{comment} = "Empty result from classifier->tree";
  } else {
    $debug && print STDERR "$subn: before getClassification tree: accession=$accession tree=$treeFile\n";
    $debug && print STDERR "$subn: lookupfile: ".$config->getValue("TempDir")
         . "/" . $config->getValue("ClassifierDir") . $config->getValue("lookupfile")."\n";
    $debug && print STDERR "$subn: treeFile=$treeFile\n";
if (0) {
    ($class, $probability) = $classifier->getClassification( $treeFile, $accession );
} elsif (0) {
    ($class, $probability) = $classifier->findGenotype( $treeFile, $accession );
} elsif (1) {
    ($class, $probability, $fileOutput) = $classifier->callCladinator( $treeFile, $faa->{job_name} . $faa->{accession} );
}
    $genotype->{cladinatorOutput} = $fileOutput;
    $thirdColumn = $treeData;
    $debug && print STDERR "$subn: After classifier->getClassification\n";
    $debug && print STDERR "$subn: accession:class:thirdColumn= $accession\t$class\t$thirdColumn\n";
#    print LOG "$accession\t$class\n";
#    print $OUTFILE "$accession\t$class\t$thirdColumn\n" if ( defined($OUTFILE) );

  # 4 Post mapping
    if ( scalar keys %{$postMapping} > 0 && defined( $postMapping->{$class} ) ) {
        $debug && print STDERR "Post mapping class = $class -> $postMapping->{$class}";
        $class = $postMapping->{$class};
    }

  # 5 Output classification
    my @data = ( $faa->{jobName0} );
    push (@data, ( $class, $faa->{accession} ));
    if ($reportUnifiedTree) {
      push( @data, $tree );
    }
    print STDOUT ( join( "\t", ("$subn:", @data) ) . "\n" );

    # Save TPP result to $genotype
    if ($class eq 'ND' || $class eq '?') {
#      $genotype->{status} = 'Fail';
      $genotype->{comment} = "$class returned by pplacer";
#      $genotype->{BI} = "BI=1.000\tgenotype=$class";
#      $genotype->{BI} .= "\tsub=NA";
    } else {
      $genotype->{status} = 'Success';
      $genotype->{comment} = "";
      $genotype->{BI} = sprintf("BI=%5.3f\tgenotype=$class", $probability);
      $genotype->{BI} .= "\tsub=NA";
    }
  }

  # Replace ref accession with gt|accession in final tree
  if (1 && !$skipTaxit) {
      my $temp_tre = "$genotype->{accession}.tog.tre";
      if (!-e "$temp_tre" && -e "$genotype->{accession}.sing.tre") {
        $temp_tre = "$genotype->{accession}.sing.tre";
      }
      my $lookup = $classifier->{lookup};
      $debug && print STDERR "$subn: faa=\n".Dumper($faa)."End of faa\n";
      $debug && print STDERR "$subn: lookup=\n".Dumper($lookup)."End of lookup\n";
      if (-e "$temp_tre" && scalar(keys %{$lookup})>0) {
        `mv ${temp_tre} ${temp_tre}0`; # Backup original tree
        open my $fileIn, '<', "${temp_tre}0" or croak "Can'f open file '${temp_tre}0'";
        open my $fileOut, '>', "${temp_tre}";
        while (my $inline=<$fileIn>) {
          foreach my $k (sort {$lookup->{$a} cmp $lookup->{$b}} keys %{$lookup}) {
            $inline =~ s/$k/$lookup->{$k}|$k/;
#            $debug && print STDERR "$subn: $fileOut='$inline'\n";
          }
          $inline =~ s/$faa->{accession}/$faa->{long_name}/g; # Replace internal accession w/ original defline
          print $fileOut "$inline\n";
        }
        $debug && print STDERR "$subn: ${temp_tre}=".`cat ${temp_tre}`;
        close $fileIn;
        close $fileOut;
      }
  }
#exit;

  $debug && print STDERR "$subn: leaving subroutine\n";
  return $genotype;
} # sub TPP

##################################################################
## restoreSeqID, replaces the short ID of a sequence to original ID
## Returns the new data
sub restoreSeqID {
  my ($faa, $treeFile, $lookup) = @_;
  my $debug = 0 || $debug_all || Genotype_util::getDebugAll();
  my $subn = "Genotype_TPP::restoreSeqID";

  $debug && print STDERR "$subn: \$treeFile=$treeFile\n";
  $debug && print STDERR "$subn: \$faa=\n".Dumper($faa)."End of \$faa\n";

  my $data = '';
  $lookup = {} if (!$lookup);
  if ( -e "$treeFile" ) {
        `mv ${treeFile} ${treeFile}0`; # Backup original file
        open my $fileIn, '<', "${treeFile}0" or croak "Can'f open file '${treeFile}0'";
        open my $fileOut, '>', "${treeFile}";
        while (my $inline=<$fileIn>) {
          foreach my $k (sort {$lookup->{$a} cmp $lookup->{$b}} keys %{$lookup}) {
            $inline =~ s/$k/$lookup->{$k}|$k/;
#            $debug && print STDERR "$subn: $fileOut='$inline'\n";
          }
          $inline =~ s/$faa->{accession}/$faa->{long_name}/g; # Replace internal accession w/ original defline
          $data .= $inline;
          print $fileOut "$inline\n";
        }
        $debug && print STDERR "$subn: ${treeFile}=".`cat ${treeFile}`;
        close $fileIn;
        close $fileOut;
  }

#exit;
  return $data;
} # sub restoreSeqID


1;
