package Genotype_recomb;

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
use Genotype_util;
use Genotype_newickBI;

my $debug_all = 0;

# Hardcoded constants
my $GAP_SIZE = 400;
my $STEP_SIZE = 100;

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

#    my $newHCV_ICTV = 1; # Sep 30, 2013
#    my $fmtGenotype = '\d[-._\dA-Za-z]*';
#    $fmtGenotype = '([0-9]\w+|\d-[A-Za-z]+|NA|.+)' if (!$newHCV_ICTV);
    my $fmtGenotype = Genotype_Def::getfmtGenotype;

############## Subroutines ################

=head1
sub recombination, takes a input alignment file in phylip format, determines the recombination profile,
 and returns a reference of an array of hashes.
=cut

sub recombination {
#    my ($job_name, $taxon, $seq_name, $long_name, $temp_plp, $tempdir, $dir_path, $progs, $no_rerun) = @_;
    my ($faa, $seq_name, $taxon, $region, $temp_plp, $deflines, $tempdir, $dir_path, $progs, $no_rerun, $ALGO) = @_;
    $debug_all = Genotype_util::getDebugAll();
    my $debug = 0 || $debug_all;
    my $subn = "Genotype_recomb::recombination";
#my $GAP_SIZE = 400;
#my $STEP_SIZE = 100;

    $debug && print STDERR "$subn: Start\n";
    $debug && print STDERR "$subn: faa=\n".Dumper($faa)."End of faa\n";
#   push @$faas, {job_name=>$job_name, accession=>$acc, file=>$temp_faa, long_name=>$long_id};
    my $job_name = $faa->{job_name};
    $job_name = 'test' if (!$job_name);
    my $acc      = $faa->{accession};
    my $infile   = $faa->{file};
    my $job_name1 = $infile;
    $job_name1 =~ s/[.]\w{1,7}//i;
    my $long_name = $faa->{long_name};

    if (!$temp_plp) {
        print STDERR "$subn: \$temp_plp=null, abort.\n";
        exit(1);
    } elsif (!-e "$tempdir/$temp_plp") {
        print STDERR "$subn: file=$tempdir/$temp_plp" ." can't be found, abort.\n";
        exit(1);
    }

    $debug && print STDERR "$subn: Creating a new alignment object\n";
    my $recombs = [];
    my $alnio = Bio::AlignIO->new('-file' => "$tempdir/$temp_plp");
#    $debug && print STDERR "$subn: \$alnio = \n".Dumper($alnio)."End of \$alnio\n\n";
    my $aln = $alnio->next_aln();
    $debug && print STDERR "$subn: Creating a new alignment object, finished\n";
#    $debug && print STDERR "$subn: \$aln=\n".Dumper($aln)."End of \$aln\n\n";

    my $iter_min = 0;
    my $iter_max = floor(($aln->length() - 50) /$STEP_SIZE - $GAP_SIZE/$STEP_SIZE +1);
    if ($ALGO eq 'TPP' && defined($faa->{seqAligned})) {
        $iter_max = floor((length($faa->{seqAligned}) - 50) /$STEP_SIZE - $GAP_SIZE/$STEP_SIZE +1);
    }
    if ($debug) {
#        $iter_min = $iter_max -20 if ($debug && $iter_max>20);
#        $iter_min = 35;
#        $iter_max = 35 if ($debug && $iter_max>50);
    }
    $debug && print STDERR "$subn: \$iter_min=$iter_min \$iter_max=$iter_max length=".$aln->length()."\n";

    my $oldDebug = $debug;
    my $oldDebugAll = $debug_all; # Needed to turn off debug in other subroutines
    for my $iseg ($iter_min .. $iter_max) {
      # removes most files from recombination even with $debug=1, except 1st and last
      $debug = ($iseg==$iter_min || $iseg==$iter_max) ? $oldDebug : 0;
      $debug_all = ($iseg==$iter_min || $iseg==$iter_max) ? $oldDebugAll : 0;
      my $errcode = '';
      my $start = $iseg*$STEP_SIZE;
      my $end = $iseg*$STEP_SIZE + $GAP_SIZE -1;
      $end = $aln->length()-1 if ($end>$aln->length());
      ($debug && $iseg<5) && print STDERR "$subn: \$iseg=$iseg \$start=$start \$end=$end\n";

      # create a hash to hold the result
      my $recomb;
      $recomb = Genotype_Def::newGenotype( $faa->{job_name}, $acc, $long_name);
      $recomb->{type}   = "recomb_window";
      $recomb->{range1} = sprintf("%4s", $start+1);
      $recomb->{range2} = sprintf("%s", $end+1);
      $recomb->{taxon}  = $taxon;
      $recomb->{region} = $region;
      $recomb->{status} = 'Failed';
      $recomb->{comment}= "Starting sub recombination";
      $recomb->{BI} = "BI=\tgenotype=NA\tsub=NA";
      $recomb->{algo} = $ALGO;

      # Take a slice of the MSA and write to a new fasta file
      my $seg_afa = sprintf("%s_%d_%d.afa", "${job_name}_$acc", $start+1, $end+1);
      my $seg_plp = sprintf("%s_%d_%d.phylip", "${job_name}_$acc", $start+1, $end+1);

      my $msgs = '';
      my $seqExcluded;
#################### Genotype Analysis by MSA ####################
      if ($ALGO eq 'MSA') {
        $debug && print STDERR "$subn: getting MSA file for: \$iseg=$iseg $start..$end length=".$aln->length()."\n";
        ($seg_afa, $seqExcluded, $errcode) = Genotype_util::writeSliceFasta( $aln, $job_name, $acc, $seq_name, $start, $end, $tempdir);
        $debug && print STDERR "$subn: head $tempdir/$seg_afa=\n";
        $debug && print STDERR `head $tempdir/$seg_afa`;
        if ($errcode) {
            $debug && print STDERR "$subn: \$start=$start \$errcode=$errcode\n";
            $recomb->{comment} = $errcode;
            push @$recombs, $recomb;
            next;
        }
        if (0) {
          # Convert fasta MSA to phylip format
          my $in  = Bio::AlignIO->new('-file' => "$tempdir/$seg_afa", '-format' => 'fasta',);
          my $in1 = $in->next_aln;
#          $debug && print STDERR "$subn: \$in1=\n".Dumper($in1)."End of \$in1\n\n";
          my $out = Bio::AlignIO->new('-file' => ">$tempdir/$seg_plp", '-format' => 'phylip',);
          $out->write_aln($in1);
        } else {
          Genotype_util::fasta2phylip("$tempdir/$seg_afa", "$tempdir/$seg_plp");
        }
        $debug && print STDERR "$subn: head $tempdir/$seg_afa=\n";
        $debug && print STDERR `head $tempdir/$seg_afa`;
        $debug && print STDERR "$subn: head $tempdir/$seg_plp=\n";
        $debug && print STDERR `head $tempdir/$seg_plp`;

        # Generate phylogeny tree for each MSA slice
        ($temp_new, $errcode) = Genotype_calctree::calc_tree( $faa, $tempdir, $seg_plp, $progs, $no_rerun, $dir_path, $debug);
        $debug && print STDERR "$subn: \$temp_new=$temp_new \$errcode='$errcode'\n";

        if (scalar(keys %$deflines)>0) {
            my $defline_fail = Genotype_calctree::restoreDefline($deflines, $seqExcluded, "$tempdir/$temp_new");
            (scalar(keys %$defline_fail)>0) && print STDERR "$subn: \$defline_fail=\n".Dumper($defline_fail)."End of \$defline_fail\n";
        }

        # check the phylogeny tree result
        if ($errcode) {
            $debug && print STDERR "$subn: \$start=$start \$errcode=$errcode\n";
            $recomb->{comment} = $errcode;
            push @$recombs, $recomb;
            next;
        } else {
            $debug && system("cp $tempdir/$temp_new $dir_path");
            $debug && print STDERR `ls -l $tempdir/$temp_new`;
        }

# Run newickBI.pl to get branching index (BI) value and genotype
#        $cmd = "./phyloplace/newickBI.pl $temp_new \"AF011753\" ";
#        $cmd = "./phyloplace/newickBI.pl $tempdir/$temp_new \"$acc\" ";
#        $debug && print STDERR "Command to run newickBI.pl is:\n$cmd\n";
#        $msgs = `$cmd`;
        $msgs = Genotype_newickBI::newickBI("$tempdir/$temp_new", $seq_name, $taxon, 0);
      }

#################### Genotype Analysis by TPP ####################
      if ($ALGO eq 'TPP') { # 'TPP'
        my $seqAligned = $faa->{seqAligned};
        $debug && print STDERR "$subn: seqAligned=$seqAligned\n";
        # Contruct sequence by leaving only the window and replacing all others by -
        my $str = '';
        for my $i (0 .. $start-1) { $str .= '-'; }
        $debug && print STDERR "$subn: seqAligned=$str\n";
        for my $i ($start .. $end) { $str .= substr($seqAligned, $i, 1); }
        $debug && print STDERR "$subn: seqAligned=$str\n";
        for my $i ($end+1 .. length($seqAligned)-1) { $str .= '-'; }
        $debug && print STDERR "$subn: seqAligned=$str\n";

        my $sequence = $str;
        $sequence = substr($sequence, $start, $end-$start+1);
        my $GAP_PCT = 80;
        $debug && printf STDERR "$subn: sequence=$sequence\n";
        my $pct = ($sequence =~ tr/-Nn//) / length($sequence) *100;
        $debug && printf STDERR "$subn:      pct=$pct%%\n";
        $debug && printf STDERR "$subn: sequence=$sequence\n";
        if ( $pct > $GAP_PCT ) {
            $errcode = ($pct<=99.9) ? sprintf("%2d%% of $seq_name in $job_name is gap \@$start..$end", $pct)
                                    : "$seq_name in $job_name is empty \@$start..$end";
            $debug && printf STDERR "$subn: skip $errcode\n";
            $recomb->{comment} = $errcode;
#            next;

        } else {
            my $skipTaxit = 1;
            $recomb->{seqAligned} = $str;
            $recomb = Genotype_TPP::TPP($faa, $recomb, $faa->{classifier}, $faa->{postMapping}, $reportUnifiedTree, $skipTaxit);
        }
        push @$recombs, $recomb;
      }
      $debug && print STDERR "$subn: from newickBI: \$iseg=$iseg \$start=$start \$msgs='$msgs'\n";

      for my $msg (split /\n/, $msgs) {
          chomp($msg);
          $debug && print STDERR "$subn: \$msg='$msg'\n";
          ($recomb->{status}, $msg, $recomb->{subtype}, $recomb->{comment}) = Genotype_util::parseBI( $msg, $taxon);
          $recomb->{Date} = Genotype_util::getTime();
          $msg = Genotype_Def::convertGT( $msg, $recomb->{taxon});
          $debug && print STDERR "$subn: \$msg='$msg'\n";
          $recomb->{BI} = $msg;
#          if ($msg =~ /^(BI=([.0-9]+)\tgenotype=(([0-9]\w+|\d-[A-Za-z]+|NA))\s*(sub=[^\t]*))/) {
          if ($msg =~ /^(BI=([.0-9]+)\tgenotype=($fmtGenotype)\s*(sub=[^\t]*))/) {
              $recomb->{BI} = $1;
              $debug && print STDERR "$subn: \$1='$1' \$2='$2' \$3='$3'\n";
          }
          push @$recombs, $recomb;
      }
      $debug && print STDERR "$subn: \$recomb=\n".Dumper($recomb)."\n";

    } # for my $iseg (0 .. 10)

    # add subroutine to parse recombs results to determine the recombination status
    # output: accession, recomb, 1..n, BI=0.0, genotype=1a,1b, date, ...
    my $s = [ $aln->each_seq_with_id($seq_name) ];
    my $seq_length = 0;
    $seq_length = length($s->[0]->seq) - ($s->[0]->seq =~ tr/-Nn//);
    my $msa_length = []; # Can't recall this part's logic
    foreach my $seq ($aln->each_seq) {
       my $sequence = $seq->seq; 
       push @$msa_length, length($s->[0]->seq) - ($s->[0]->seq =~ tr/-Nn//);
#       $debug && print STDERR "$subn: \$msa_length=@$msa_length\n";
    }
    $msa_length = [ reverse sort @$msa_length ];
    $debug && print STDERR "$subn: \$msa_length=@$msa_length\n";
    my $eff_length = ($seq_length<$msa_length->[2]) ? $seq_length : $msa_length->[2];
    $recombs = Genotype_recomb::summ_recomb( $recombs, $seq_length);

    # Delete all files if no debugging
    if (1 && !$debug) {
        #print STDERR "$subn: Removing temporary files\n";
        `rm $tempdir/${job_name}_*_*.afa` if (-e "$tempdir/${job_name}_*_*.afa"); # Delete the MSA in fasta format
        `rm $tempdir/${job_name}_*_*.phylip` if (-e "$tempdir/${job_name}_*_*.phylip"); # Delete the MSA in phylip format
        `rm $tempdir/${job_name}_*_*.phylip` if (-e "$tempdir/${job_name}_*_*.phylip"); # Delete the MSA in phylip format
        `rm $tempdir/${job_name}_*_*.phylip_fastme_tree.new` if (-e "$tempdir/${job_name}_*_*.phylip_fastme_tree.new"); # Delete the phylogeny tree
        `rm $tempdir/${job_name}_*_*.phylip_fastme_tree.new` if (-e "$tempdir/${job_name}_*_*.phylip_fastme_tree.new"); # Delete the phylogeny tree
        `rm $tempdir/${job_name}_*_*.phylip_outfile` if (-e "$tempdir/${job_name}_*_*.phylip_outfile"); # Delete dnadist output
        `rm $tempdir/${job_name}_*_*.phylip_outfile` if (-e "$tempdir/${job_name}_*_*.phylip_outfile"); # Delete dnadist output
        `rm $dir_path/${job_name}_*_*.phylip_outfile` if (-e "$dir_path/${job_name}_*_*.phylip_outfile"); # Delete dnadist output
    }
    $debug = $oldDebug;
    $debug_all = $oldDebugAll; # Needed to turn off debug in other subroutines
    $debug && print STDERR "$subn: \$recomb[1]=\n".Dumper($recombs->[1])."\n";
    $debug && print STDERR "$subn: \$recomb[last]=".scalar(@$recombs)."\n";
    $debug && print STDERR "$subn: \$recomb[last]=\n".Dumper($recombs->[scalar(@$recombs)-1])."\n";
#    $debug && print STDERR "$subn: \n";
#    for (@$recombs) { print STDERR Genotype_util::to_String( $_)."\n"; }

    $debug && print STDERR "$subn: leaving subroutine\n";
    return $recombs;
} # sub recombination


=head1
sub summ_recomb looks through all recombs, and determines the dominant genotypes among the windows
  Prepends result to the input array
=cut

sub summ_recomb {
    my ($recombs, $seq_length) = @_;
    my $debug = 0 || $debug_all;
    my $subn = "Genotype_recomb::summ_recomb";

    my $TH_BI = 0.711;
    if ($debug) {
        $TH_BI = 0.4;
        print STDERR "$subn: threshold reduced to 0.4 (from 0.711) due to \$debug=$debug\n";
    }
    # take the required min. # of windows for a genotype to be considered as present as one less than square root of
    # total number of windows from the input sequence (not MSA), and the min. value is 1.
    #my $TH_min_num_windows = sqrt($seq_length/100 -2) -1 ;
    my $TH_min_num_windows = $seq_length/100 ;
    $TH_min_num_windows -= 2 if ($TH_min_num_windows>3);
    $TH_min_num_windows = sqrt($TH_min_num_windows) -1;
    $TH_min_num_windows = ($TH_min_num_windows<1) ? 1 : $TH_min_num_windows;
    $TH_min_num_windows = sprintf("%4.2f", $TH_min_num_windows);
    print STDERR "$subn: MSA_length=$seq_length \$TH_min_num_windows=$TH_min_num_windows\n";

    my $recomb_summary = [];
    my $sub_summary = {};
    my $recomb0 = Genotype_Def::newGenotype();
    $recomb0->{job_name}  = $recombs->[0]->{job_name};
    $recomb0->{accession} = $recombs->[0]->{accession};
    $recomb0->{long_name} = $recombs->[0]->{long_name};
    $recomb0->{algo} = $recombs->[0]->{algo};
    $recomb0->{type} = 'recomb_summ';
    $recomb0->{range1} = $recombs->[0]->{range1};
    $recomb0->{range2} = $recombs->[0]->{range2};
    $recomb0->{taxon}  = $recombs->[0]->{taxon};
    $recomb0->{region} = $recombs->[0]->{region};
    $recomb0->{Date} = Genotype_util::getTime();
    $recomb0->{status} = 'Failed';
    $recomb0->{comment} = '';
    unshift @$recombs, $recomb0;

    my $geno_old = '';
    my $BI_max = 0;
if (0) {
    for my $i (1 .. $#{$recombs}) {
        my $recomb = $recombs->[$i];
        $recomb0->{range1} = $recomb->{range1} if ($recomb->{range1} < $recomb0->{range1});
        $recomb0->{range2} = $recomb->{range2} if ($recomb0->{range2} < $recomb->{range2});

        my $BI = $recomb->{BI};
        $debug && print STDERR "$subn: range:$recomb->{range1}..$recomb->{range2} \$BI='$BI'\t";
        if (!($BI =~ /^BI=([.0-9]+)\tgenotype=($fmtGenotype)\tsub=(\w*)\t(.*)/
                   || $BI =~ /^BI=([.0-9]+)\tgenotype=($fmtGenotype)\tsub=(\w*)/)) {
            $debug && print STDERR ": Problem!\n";
            next;
        }

        $debug && print STDERR "\$1='$1' \$2='$2' \$3='$3'\n";
#        $BI_max = $1 if ($1>$BI_max);
      if (0) {
        $recomb_summary->{$2}++ if ($1>=$TH_BI);
      } else {
        my $seen = 0;
        for my $j (0 .. $#{$recomb_summary}) {
            next if ($recomb_summary->[$j]->[0] ne $2);
            $debug && print STDERR "$subn: \$geno_old='$geno_old' \$2=$2\n";
            $recomb_summary->[$j]->[1]++;
            $recomb_summary->[$j]->[2] = $1 if ($recomb_summary->[$j]->[2] < $1); # BI
            $geno_old = $2;
            $debug && print STDERR "$subn: \$recomb_summary->[$j]='@{$recomb_summary->[$j]}' \$2=$2\n";
            if ($recomb_summary->[$j]->[1]>=$TH_min_num_windows && $recomb_summary->[$j]->[1]<1+$TH_min_num_windows) {
                # list the first genotype w/ enough windows first, and so on, e.g. JX826592: 3a/2a
                for my $k (0 .. $j) {
                  $debug && print STDERR "$subn: j=$j:'@{$recomb_summary->[$j]}' k=$k:'@{$recomb_summary->[$k]}'\n";
                  next if ($k==$j);
                  next if ($recomb_summary->[$k]->[1] >= $recomb_summary->[$j]->[1]);
#                  for my $m (0 .. $#{$recomb_summary}) {
#                      $debug && print STDERR "$subn: m=$m:'@{$recomb_summary->[$m]}'\n";
#                  }
                  # insert the new genotype to where it belongs
                  my $new_summary = [ @{$recomb_summary}[0.. $k-1],
                                      $recomb_summary->[$j],
                                      @{$recomb_summary}[$k .. $j-1],
                                      @{$recomb_summary}[$j+1 .. $#{$recomb_summary}],
                                    ];
                  $recomb_summary = $new_summary;
                  for my $m (0 .. $#{$recomb_summary}) { $debug && print STDERR "$subn: m=$m:'@{$recomb_summary->[$m]}'\n"; }
                  last;
                }
            }
            $seen = 1;
            last;
        }
        if (!$seen) {
            push @$recomb_summary, [$2, 1, $1];
        }
        $geno_old = $2;
        #$debug && print STDERR "$subn: \$recomb_summary=\n".Dumper($recomb_summary)."\n";
      }
        $sub_summary->{$3}++ if ($1>=$TH_BI);
    }
} else {
    for my $i (1 .. $#{$recombs}) {
        my $recomb = $recombs->[$i];
        $recomb0->{range1} = $recomb->{range1} if ($recomb->{range1} < $recomb0->{range1});
        $recomb0->{range2} = $recomb->{range2} if ($recomb0->{range2} < $recomb->{range2});

        my $BI = $recomb->{BI};
        $debug && print STDERR "$subn: range:$recomb->{range1}..$recomb->{range2} \$BI='$BI'\t";
        if (!($BI =~ /^BI=([.0-9]+)\tgenotype=($fmtGenotype)\tsub=(\w*)\t(.*)/
                   || $BI =~ /^BI=([.0-9]+)\tgenotype=($fmtGenotype)\tsub=(\w*)/)) {
            $debug && print STDERR ": Problem!\n";
            next;
        }

        $debug && print STDERR "\$1='$1' \$2='$2' \$3='$3'\n";
        my $seen = 0;
        if ($#{$recomb_summary}<0 || ($recomb_summary->[$#{$recomb_summary}]->[0] ne $2 && $2!~m/^\sna*\s*$/i)) {
            push @$recomb_summary, [$2, 1, $1];
        } elsif ($2!~m/^\sna*\s*$/i) {
            $recomb_summary->[$#{$recomb_summary}]->[1]++;
            $recomb_summary->[$#{$recomb_summary}]->[2] = $1 if ($recomb_summary->[$#{$recomb_summary}]->[2] < $1); # BI
        }
        for my $m (0 .. $#{$recomb_summary}) { $debug && print STDERR "$subn: m0=$m:'@{$recomb_summary->[$m]}'\n"; }
        $sub_summary->{$3}++ if ($1>=$TH_BI);
    }

    # Exclude those segments not long enough
    my $new_summ = [];
    for my $j (0 .. $#{$recomb_summary}) {
        push @$new_summ, $recomb_summary->[$j] if ($recomb_summary->[$j]->[1] >= $TH_min_num_windows);
    }
    $recomb_summary = $new_summ;
    for my $m (0 .. $#{$recomb_summary}) { $debug && print STDERR "$subn: m1=$m:'@{$recomb_summary->[$m]}'\n"; }

    # get the longest segment for the same genotype
    $new_summ = [ $recomb_summary->[0] ];
    for my $j (1 .. $#{$recomb_summary}) {
        my $seen = 0;
        for my $k (0 .. $j-1) {
            next if ($recomb_summary->[$j]->[0] ne $recomb_summary->[$k]->[0]);
            if ($recomb_summary->[$k]->[1] < $recomb_summary->[$j]->[1]) {
                $recomb_summary->[$k] = $recomb_summary->[$j];
#                $recomb_summary->[$k]->[1] = $recomb_summary->[$j]->[1];
#                $recomb_summary->[$k]->[2] = $recomb_summary->[$j]->[2];
            }
            $seen = 1;
        }
        (!$seen) && push @$new_summ, $recomb_summary->[$j];
    }
    $recomb_summary = $new_summ;
    for my $m (0 .. $#{$recomb_summary}) { $debug && print STDERR "$subn: m2=$m:'@{$recomb_summary->[$m]}'\n"; }
}

    $debug && print STDERR "$subn: \$BI_max=$BI_max\n";
    $debug && print STDERR "$subn: \$recomb0=\n".Dumper($recomb0)."\n";
    $debug && print STDERR "$subn: \$recomb_summary=\n".Dumper($recomb_summary)."\n";
    $debug && print STDERR "$subn: \$recomb_summary=\n".Dumper($sub_summary)."\n";
    my $BI = sprintf("BI=%5.3f\tgenotype=\tsub=NA", $BI_max);
    my $gt_str = '';
    # find BI_max among the max of each genotype. This is used to exclude any genotype w/ too low BI, eg EU155250.
    for my $j (0 .. $#{$recomb_summary}) {
        $k = $recomb_summary->[$j]->[0];
        $debug && print STDERR "$subn: #$j:'@{$recomb_summary->[$j]}'\n";
        next if (!$k || $k eq 'NA');
        next if ($recomb_summary->[$j]->[1] < $TH_min_num_windows);
        $BI_max = $recomb_summary->[$j]->[2] if ($recomb_summary->[$j]->[2] >$BI_max); # take max BI from valid genotype
        $debug && print STDERR "$subn: \$BI_max=$BI_max \$recomb_summary='@{$recomb_summary->[$j]}'\n";
    }
    for my $j (0 .. $#{$recomb_summary}) {
        $k = $recomb_summary->[$j]->[0];
        $debug && print STDERR "$subn: #$j:'@{$recomb_summary->[$j]}'\n";
        if (!$k || $k eq 'NA' || $k !~ /$fmtGenotype/) {
            print STDERR "$subn: BI='@{$recomb_summary->[$j]}', does't match the required format '$fmtGenotype'. Skip\n";
            next;
        }
        if ($recomb_summary->[$j]->[1] < $TH_min_num_windows) {
            print STDERR "$subn: BI='@{$recomb_summary->[$j]}', less than TH_min_num_windows=$TH_min_num_windows. Skip\n";
            next;
        }
        if ($recomb_summary->[$j]->[2] < $BI_max*0.5) {
            print STDERR "$subn: BI='@{$recomb_summary->[$j]}', <50% of BI_max=$BI_max, though met min occurrence. Skip\n";
            next;
        }
      if ( 0 ) {
        $BI .= '/' if ($BI =~ m/genotype=$fmtGenotype/);
        $BI .= $k;
        $recomb0->{status} = 'Success';
        $debug && print STDERR "$subn: \$BI='$BI'\n";
      } else {
#        $BI_max = $recomb_summary->[$j]->[2] if ($recomb_summary->[$j]->[2] >$BI_max); # take max BI from valid genotype
#        $debug && print STDERR "$subn: \$BI_max=$BI_max\n";
        $gt_str .= '/' if ($gt_str);
        $gt_str .= $k;
        $recomb0->{status} = 'Success';
        $debug && print STDERR "$subn: \$gt_str='$gt_str'\n";
      }
    }
    $BI_max = 0.000 if (!$gt_str);
    $gt_str = 'NA' if (!$gt_str);
    $BI = sprintf("BI=%5.3f\tgenotype=$gt_str", $BI_max);
    $debug && print STDERR "$subn: \$BI_max=$BI_max\n";
    $debug && print STDERR "$subn: \$BI='$BI'\n";
    $BI .= "\tsub=";
    for $k (sort keys %$sub_summary) {
        $debug && print STDERR "$subn: \$k=$k:@{$sub_summary->{$k}} \$TH_min_num_windows=$TH_min_num_windows\n";
        next if ($k eq 'NA');
        next if ($sub_summary->{$k} < $TH_min_num_windows);
#        $BI .= ',' if ($BI =~ m/genotype=([0-9][A-Za-z]*|\d-[A-Za-z]+|.+)\s+sub=\w+$/);
        $BI .= '/' if ($BI =~ m/genotype=$fmtGenotype\s+sub=\w+$/);
        $BI .= $k;
#        $recomb0->{status} = 'Success';
        $debug && print STDERR "$subn: \$BI='$BI'\n";
    }
    $BI =~ s/sub=$/sub=NA/i;
    $recomb0->{BI} = $BI;
    $debug && print STDERR "$subn: \$recomb0=\n".Dumper($recomb0)."\n";

    $debug && print STDERR "$subn: leaving subroutine\n";
    return $recombs;
} # sub summ_recomb


=head1
sub get_graph_data, takes an array of hashed of the re-combination analysis,
  draw a bar graph based on the BI value vs location of the window.
=cut

sub get_graph_data {
    my ($genotypes, $job_name, $dir_path, $ALGO) = @_;
    my $subn = "Genotype_recomb::get_graph_data";

    my $debug = 0 || $debug_all;
    my $errcode = '';
    my $result = '';
    $job_name = '' if (!defined($job_name));

    $debug && print STDERR "$subn: \$genotypes=\n".Dumper($genotypes)."End of \$genotypes\n";
    if ( 0 ) {
        $debug && print STDERR "$subn: \n";
        for (@$genotypes) { print STDERR Genotype_Def::to_String( $_)."\n"; }
    }
    my $data;
#    my $data1;
    my $legends = {};
    # Gather the legends only, not the data
    for my $genotype (@$genotypes) {
        next if ($genotype->{type} ne 'recomb_window');
#        next if ($genotype->{BI} !~ /^BI=([0-9.\-+]*)\s*genotype=(\w*)/);
        next if ($genotype->{BI} !~ /^BI=([0-9.\-+]*)\s*genotype=(.*)\tsub=/);
        my $region = $genotype->{region};
        my $BI = $1;
        $BI = 0.000 if ($BI eq '');
        my $gt = $2;
        $gt = 'NA' if ($gt eq '');
        my $legend_seen = 0;
        for my $legend (@{$legends->{$region}}) {
#            $debug && print STDERR "$subn: \$legend=$legend \$gt=$gt \$legends='@$legends'\n";
            if ($legend eq $gt) {
                $legend_seen = 1;
                last;
            }
        }
        if (!$legend_seen) {
            $debug && print STDERR "$subn: \$gt=$gt region=$genotype->{region} \$legends=\n".Dumper($legends)."End of \$legends\n";
            if (exists($legends->{$region})
                && $#{$legends->{$region}}>=0
                && $legends->{$region}->[$#{$legends->{$region}}] eq 'NA') {
                my @a1 = (@{$legends->{$region}});
                $legends->{$region} = [ @a1[0 .. $#a1-1], $gt, $a1[$#a1] ];
            } else {
                push @{$legends->{$region}}, $gt;
            }
        }
#        $debug && print STDERR "$subn: \$BI=$BI \$gt=$gt \$legends='@$legends'\n";
    }
    $debug && print STDERR "$subn: \$legends=\n".Dumper($legends)."End of \$legends\n";

    # Gather the data, according to genotypes
    my $graph_data = [];
    for my $genotype (@$genotypes) {
#        $debug && print STDERR "$subn: \$genotype=\n".Dumper($genotype)."End of \$genotype\n";
        my $region = $genotype->{region};
        if ($genotype->{type} ne 'recomb_window') {
            next if ($genotype->{type} ne 'recomb_summ');
            my $data1 = [];
            $data = $data1;
            my $title = $genotypes->[0]->{long_name};
            $title = $1 if ($title =~ m/^\s*([^:]+):/i);
            $title = "Genotype profile for $title";
            $title .= ":$genotype->{region}" if ($region !~ m/^\s*na\s*$/i);
            $title .= " from the recombination analysis";
            my $labels = {
                   x_label => 'Alignment position (nt)',
                   y_label => ($ALGO eq 'TPP') ? 'Probability'
                            : 'Branching Index',
#                   title   => "Genotype profile for '$job_name'",
#                  title   => "Genotype profile for $genotypes->[0]->{long_name}:$genotype->{region} from the recombination analysis",
                   title   => $title,
                   job_name => $job_name,
                      };
#            $debug && print STDERR "$subn: \$labels=\n".Dumper($labels)."End of \$labels\n";
            push @$graph_data, {
                       region => $region,
                       taxon => $genotype->{taxon},
                       data => $data,
                       legends => $legends->{$region},
                       labels => $labels,
                     };
            next;
        }
#        next if ($genotype->{BI} !~ /^BI=([0-9.\-+]*)\s*genotype=(\w*)/);
        next if ($genotype->{BI} !~ /^BI=([0-9.\-+]*)\s*genotype=(.*)\tsub=/);
        push @{$data->[0]}, floor( ($genotype->{range1}+$genotype->{range2}-1)/2 );
        my $BI = ($1 eq '') ? 0.000 : $1;
        my $gt = ($2 eq '') ? 'NA'  : $2;

        for my $i (0 .. $#{$legends->{$region}}) {
            push @{$data->[$i+1]}, ($legends->{$region}->[$i] eq $gt) ? $BI : undef;
        }
#        $debug && print STDERR "$subn: \$graph_data=\n".Dumper($graph_data)."End of \$graph_data\n";
    }
    $debug && print STDERR "$subn: \$graph_data=\n".Dumper($graph_data)."End of \$graph_data\n";

    if (!$data || (scalar keys %$legends)<0) {
        return $errcode;
    }

# Comment out Draw_graph if no need for graph
#    Draw_graph::draw_bar( $data, $legends, $labels, $dir_path);

    $debug && print STDERR "$subn: leaving subroutine\n";

    return $graph_data;
} # sub get_graph_data

=head1
Save recombination result to <job_name>_recombination_summary.tsv
=cut
sub saveGraphData {
    my ($genotypes, $job_name, $acc, $dir_path, $argv) = @_;
    my $subn = "Genotype_recomb::saveGraphData";
    my $debug = 0 || $debug_all;

    my $fileRecomb = '';
    my $OUTF;
    for my $ii (1 .. $#{$genotypes}) {
        if ($ii==1) { # open the file if first record
            $fileRecomb = "${job_name}${acc}_recombination_summary.tsv";
            print STDERR "$subn: genotype result: \$fileRecomb=$fileRecomb\n";
            $OUTF = undef;
            open $OUTF, '>', "$dir_path/$fileRecomb" or croak "Can't open outfile '$dir_path/$fileRecomb': $OS_ERROR";
#           print $OUTF "$0 @$argv\n";
            print $OUTF "# Algorithm=$genotypes->[$ii]->{algo}\n";
            print $OUTF "# Summary:\n";
            print $OUTF "# Job_name\tJob_type\tregion\tstart..end\tBranch_Index\tgenotype\tsubtype\ttaxon";
            #print $OUTF "\tDate(ymdhms)";
            print $OUTF "\tStatus\tComment\n";
        }
        my $msg = Genotype_Def::to_String( $genotypes->[$ii]);
        print $OUTF "$msg\n";
        $debug && print STDERR "$subn: \$msg='$msg'\n";
#        print $OUTF "$msgs->[$ii]\n";
        if ($ii==$#{$genotypes}) { # close the file if last record
            print $OUTF "\n";
            close $OUTF or croak "Can't close outfile '$dir_path/$fileRecomb': $OS_ERROR";
        }
    }
    return $fileRecomb;
} # sub saveGraphData


=head1
=cut
sub drawRecombGraph {
    my ($genotypes, $job_name, $acc, $dir_path) = @_;
    my $subn = "Genotype_recomb::drawRecombGraph";
    my $debug = 0 || $debug_all;

    print STDERR "$subn: Drawing profile for the recombination result:\n";
    my $fileRecomb = '';
    $fileRecomb =Genotype_recomb::saveGraphData($genotypes, $job_name, $acc, $dir_path);
    # Draw a bar graph, using the recombination summary as input
#    Genotype_util::draw_graph_bar( $genotypes, "${job_name}${acc}", $dir_path);
    eval("require GD::Graph::bars;");
    if ($@) {
      print STDERR "$subn: The required Perl module GD::Graph::bars is missing, abort\n";
    } else {
      my $cmd = Genotype_Def::getGenotypeHome . "vipr_genotype_recombgraph.pl ";
      $cmd .= "-debug " if ($debug);
      $cmd .= "-d $dir_path -i $fileRecomb ";
      $cmd .= "2>&1";
      $debug && print STDERR "$subn: \$cmd='$cmd'\n";
      my $result = ` $cmd `;
      print STDERR "$result";
    }
}

1;
