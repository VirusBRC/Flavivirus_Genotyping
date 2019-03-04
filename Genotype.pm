package Genotype;

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
use Genotype_TPP;

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
                   "muscle";
#                   "clustalw";
#                   "mafft"; # Actual command: mafft-profile
    my $RMSA = Genotype_Def::getRMSA();

#    my $newHCV_ICTV = 1; # Sep 30, 2013
#    my $fmtGenotype = '\d[-._\dA-Za-z]*';
#    $fmtGenotype = '([0-9]\w+|\d-[A-Za-z]+|NA|.+)' if (!$newHCV_ICTV);
    my $fmtGenotype = Genotype_Def::getfmtGenotype;

############## Subroutines ################


=head1
sub genotype, takes a input genome object, runs the genotype process, and returns a reference of
 an array of hashes.
Input: Genotype::genotype: $faa=
$VAR1 = {
          'long_name' => 'EU155347',
          'file' => 'EU155347_seq001.faa',
          'job_name' => 'EU155347',
          'taxon' => 'HCV',
          'accession' => '_seq001'
        };
Output: Genotype::genotype: $genotypes=
$VAR1 = [
          {
            'status' => 'Success',
            'long_name' => 'EU155347',
            'Date' => '20130606113717',
            'taxon' => 'HCV',
            'accession' => '_seq001',
            'comment' => '',
            'range1' => '1',
            'job_name' => 'EU155347',
            'type' => 'genotype',
            'range2' => 9289,
            'BI' => 'BI=1.000   genotype=1a'
          }
        ];
=cut

##################################################################
sub genotype {
  my ($faa, $tempdir, $dir_path, $progs, $no_rerun, $run_recomb) = @_;
  my $debug = 0 || $debug_all || Genotype_util::getDebugAll();
  my $subn = "Genotype::genotype";

  my $taxon = $faa->{taxon};
  my $ALGO = 'MSA';
  my $USE_TPP = 1;
  if ($USE_TPP &&  $taxon =~ m/^HCV$/i) {
    $ALGO = 'TPP'; # taxit/pplacer/guppy
    $ALN_PROG = ( 1 ) ? "mafft"
              : "clustalw";
  }

  if ($ALGO eq 'MSA') {
    $progs->{muscle}  = "muscle"  if (!exists($progs->{muscle}) || !$progs->{muscle});
    $progs->{dnadist} = "dnadist" if (!exists($progs->{dnadist})|| !$progs->{dnadist});
    $progs->{fastme}  = "fastme"  if (!exists($progs->{fastme}) || !$progs->{fastme});
  } elsif ($ALGO eq 'TPP') {
    $progs = Genotype_TPP::getProgs($progs);
  }
  $debug && print STDERR "$subn: \$progs=\n".Dumper($progs)."End of \$progs\n";

  Genotype_util::check_programAll( $progs);

#   push @$faas, {job_name=>$jobName, accession=>$acc, file=>$temp_faa, long_name=>$long_id};
  $debug && print STDERR "$subn: \$faa=\n".Dumper($faa)."End of \$faa\n";
  my $jobName = $faa->{job_name};
  $jobName = 'test' if (!$jobName);
  my $acc      = $faa->{accession};
  my $jobName0 = $faa->{file};
  $jobName0 =~ s/[.]\w{1,7}//i;
  $faa->{jobName0} = $jobName0;
  $faa->{tempdir}  = $tempdir;
  my $long_name = $faa->{long_name};
  $debug && print STDERR "$subn: NAMES: jobName=$jobName acc=$acc jobName0=$jobName0\n";

  $debug && print STDERR "$subn: \$faa=\n".Dumper($faa);
  my $genotypes = [];
  my $errcode = '';
  # a hash to hold the result
  my $genotype = Genotype_Def::newGenotype( $faa->{job_name}, $acc, $long_name);
  if ($USE_TPP &&  $taxon =~ m/^HCV$/i) {
    $genotype->{BI} = "BI=1.000\tgenotype=NA\tsub=NA";
  }
  $genotype->{algo} = $ALGO;
#  push @$genotypes, $genotype;

  if (!$faa->{file}) {
      $genotype->{comment} = "Input genome file is null";
      return [ $genotype ];
  }

  my $READ_GBK_GENOTYPE = 0;
  if ( $READ_GBK_GENOTYPE ) {
      return [ $genotype ]; # testing the detection of genotype from genbank
  }

  # setup the sequence object
  my $seqName;
  my $in  = Bio::SeqIO->new( -file   => "<$dir_path/$faa->{file}",
                             -format => 'fasta',
                           );
  my $inseq = $in->next_seq();
  $debug && print STDERR "$subn: \$inseq=\n".Dumper($inseq)."\nEnd of \$inseq\n";
  if (!$inseq) {
      $genotype->{comment} = "No sequence found in input genome";
      return [ $genotype ];
  } elsif ($in->next_seq()) {
      $genotype->{comment} = "fasta $faa->{file} has >1 sequences, please break into separate files";
      return [ $genotype ];
  } else {
      $seqName = $inseq->primary_id;
      $seqName = substr($seqName, 0, 10);
      $faa->{seqName} = $seqName;
      $faa->{seq} = $inseq->seq;
      $debug && print STDERR "$subn: \$inseq=\n".Dumper($inseq)."\nEnd of \$inseq\n";
  }
  $genotype->{range1} = sprintf("%4s", 1);
  $genotype->{range2} = length($faa->{seq});

  # determine the taxon based on blast against all refseqs
#  if (!$taxon || $taxon =~ m/^\s*NOV[^I]*\s*$/i) {
  if (!$taxon) {
      $taxon = uc Genotype_util::determine_taxon( $faa->{file}, $dir_path, $tempdir);
      print STDERR "$subn: sub determine_taxon returned: infile=$faa->{file} taxon=$taxon\n";
  }

  $genotype->{taxon} = $taxon;
  $genotype->{accession} = $seqName;

  my $QUERY_SEQ_LENGTH_MIN = 200; # Per request by Yun/Richard, July 14, 2017
  $QUERY_SEQ_LENGTH_MIN = 400; # Per request by Yun/Richard, Nov 14, 2017, BRC-11704
  if (!Genotype_Def::check_taxon_allowed( $taxon)) {
      $genotype->{comment} = "$faa->{file} has invalid taxon: '$taxon'";
      return [ $genotype ];
#  } elsif ($genotype->{range2}<100) {
  } elsif ($genotype->{range2} < $QUERY_SEQ_LENGTH_MIN) {
      $genotype->{comment} = "Input genome is too short: $genotype->{range2}";
      return [ $genotype ];
  }

  my $temp_plp = '';
  my $region = 'NA';
  my $deflines;
  print STDERR "$subn: taxon=$taxon \tlong_name=$long_name\tacc=$acc\tALGO=$ALGO\n";


    my $timeExe = time;
#################### Genotype Analysis by TPP ####################
    if ($ALGO eq 'TPP') {
      $genotypes = Genotype::genotypeTPP($faa, $genotype, $tempdir, $dir_path, $progs, $no_rerun, $run_recomb, $region);
      $temp_plp = $genotype->{MSA_phylip};

    } # if ($ALGO eq 'TPP')

#################### Genotype Analysis by MSA ####################
    if ($ALGO eq 'MSA') {
      for my $rn (0 .. $#{$RMSA}) {
        # The current approach uses regions of MSA defined in genotype1, not different ref MSA for different regions
        my $r = $RMSA->[$rn];
        $debug && print STDERR "$subn: \$rn=$rn\n";
        next if (!exists($r->{uc($taxon)}));
        $debug && print STDERR "$subn: \$rn=$rn \$r=".$r->{uc($taxon)}."\n";
        $region = 'NA';
        my $ref_fasta = $r->{uc($taxon)};
#        if ($r->{uc($taxon)} ne $RMSA->[0]->{uc($taxon)}) {
        if ($rn > 0) {
            # create a new record of genotype to store new result
            $genotype = Genotype_Def::newGenotype( $faa->{job_name}, $acc, $long_name);
            $genotype->{taxon} = $taxon;
            $genotype->{accession} = $seqName;
            $genotype->{range1} = sprintf("%4s", 1);
            $genotype->{range2} = length($inseq->seq);
#            push @$genotypes, $genotype;
        }
        if ($taxon =~ /^NOV/i) { # See if the seq covers ORF1/ORF2 if taxon=NOV
          if ( 1 ) {
            $region = Genotype_util::determine_region( $faa->{file}, $taxon, $ref_fasta, $dir_path, $tempdir);
            $debug && print STDERR "$subn: sub determine_region returned \$region='$region' for \$taxon='$taxon'\n";
            if ($region !~ /^(ORF[12])$/i) {
                print STDERR "$subn: \$taxon='$taxon' \$region='$region', ref='$r' region not recognized, skip.\n";
#                $genotype->{comment} = 'Region for NOV sequence can\'t be determined.';
#                next;
            }
          } else {
            $region = ($ref_fasta eq $RMSA->[0]->{NOVI}) ? 'ORF1'
                    : ($ref_fasta eq $RMSA->[1]->{NOVI}) ? 'ORF2'
                    :                                      'NA' ;
          }
        }
        $genotype->{region} = $region;
        print STDERR "$subn: reference: \$taxon=$taxon \$region='$region' \$ref_fasta[$rn]='$ref_fasta'\n";

        # replace defline in phylip file with only accession, and only first 10 letters
        # save the mapping in order to restore the deflines in tree file -Apr, 2017
        # This change causes result change, since some refseqs were obscured in earlier version
        # 0: (V1.3.2) no cleanup for seq ids; 1: cleanup
        my $doIdCleanup = 1;
        if ($doIdCleanup) {
            $deflines = Genotype_calctree::mapDefline($ref_fasta);
        }


        # Perform 1 genotype analysis based on given rmsa
        $temp_plp = '';
#if (!($rn>1 && $taxon =~ /^NOV/i && $genotypes->[$rn-2]->{BI} !~ m/II[.]P*4/)) {
if (!($taxon =~ /^NOV/i && $genotype->{region} =~ m/NA/)) {
# skip the analysis for NOV cases other than II.3
        ($genotype, $temp_plp) = Genotype::genotype1( $faa, $genotype, $ref_fasta, $deflines, $jobName0, $tempdir, $dir_path, $progs, $no_rerun);
}
        $genotype->{rmsa} = $ref_fasta;
        push @$genotypes, $genotype;
        $debug && print STDERR "$subn: \$genotypes=\n".Dumper($genotypes)."\nEnd of \$genotypes\n";

        if ($taxon =~ /^NOV/i) {
            # for NOV, set the subtype for genotype II.4
            $debug && print STDERR "$subn: \$rn=$rn\n";
#            exit if ($rn==0);
            if (exists($RMSA->[3]->{NOV}) && $ref_fasta ne $RMSA->[1]->{NOV}) {
                next;
            } elsif ($#{$genotypes}>2) {
                # get subtype of II.4 for NOV
#                $genotypes = Genotype::summ_genotype_NOV($genotypes);
            }
        }

      } # for my $rn (0 .. $#{$RMSA}) {

    } # if ($ALGO eq 'MSA')
    print STDERR "$subn: time for genotype analysis: ".(time-$timeExe)." sec\n";


#################### Recombination Analysis ####################
#    $run_recomb = 1;
    if ($run_recomb) {
      if ( $taxon =~ m/^NOV/i && $region =~ m/^(ORF[12]|NA)$/i ) {
        print STDERR "$subn: taxon='$taxon', \$region='$region' recombination analysis is skipped b/c of taxon\n";
      } else {
        my $timeExe = time;
#        if ($run_recomb) {
        print STDERR "$subn: Entering Genotype_recomb::recombination \$taxon=$taxon\n";
#        my $recombs = &recombination($jobName, $taxon, $seqName, $long_name, $temp_plp, $tempdir, $dir_path, $progs, $no_rerun);
        my $recombs = Genotype_recomb::recombination( $faa, $seqName, $taxon, $region, $temp_plp, $deflines, $tempdir, $dir_path, $progs, $no_rerun, $ALGO);
        print STDERR "$subn: time for recombination analysis: ".(time-$timeExe)." sec\n";
        push @$genotypes, @$recombs;

        # draw graph for recombination analysis
        Genotype_recomb::drawRecombGraph($genotypes, $jobName, $acc, $dir_path);
      }
    }
    $debug && print STDERR "$subn: \n";
    for (@$genotypes) { $debug && print STDERR Genotype_Def::to_String( $_)."\n"; }

    # Clean up
    (!$debug && -f "$tempdir/$temp_plp") && `rm $tempdir/$temp_plp`;

    # Replace internal ID with user input defline, for both MSA and taxit/pplacer/guppy results
    $debug && print STDERR "$subn: faa=\n".Dumper($faa)."End of faa\n";
    if (-e "$dir_path/$faa->{jobName0}.fasta") {
        `mv $dir_path/$faa->{jobName0}.fasta $dir_path/$faa->{jobName0}.fasta0`; # Backup original tree
        open my $fileIn, '<', "$dir_path/$faa->{jobName0}.fasta0" or croak "Can'f open file '$dir_path/$faa->{jobName0}.fasta0'";
        open my $fileOut, '>', "$dir_path/$faa->{jobName0}.fasta";
        while (my $inline=<$fileIn>) {
          if ($inline =~ m/^>/) {
            $inline =~ s/$faa->{accession}/$faa->{long_name}/g; # Replace internal accession w/ original defline
#            $debug && print STDERR "$subn: $fileOut='$inline'\n";
          }
          print $fileOut "$inline";
        }
#        $debug && print STDERR "$subn: $dir_path/${temp_afa}=".`cat $dir_path/${temp_afa}`;
        close $fileIn;
        close $fileOut;
        if (!$debug) {
          `rm $dir_path/$faa->{jobName0}.fasta0`;
        } else {
          $debug && print STDERR "$subn: Kept $dir_path/$faa->{jobName0}.fasta0`\n";
        }

        my $deflineFail = Genotype_calctree::restoreDefline($deflines, {}, "$dir_path/$faa->{jobName0}.fasta");
        $debug && print STDERR "$subn: deflineFail=\n".Dumper($deflineFail)."End of deflineFail\n";

    }

    $debug && print STDERR "$subn: leaving subroutine\n";
    return $genotypes;
} # sub genotype

################################################################################
sub genotypeTPP {
  my ($faa, $genotype, $tempdir, $dir_path, $progs, $no_rerun, $run_recomb, $region) = @_;
  my $debug = 0 || $debug_all || Genotype_util::getDebugAll();
  my $subn = "Genotype::genotypeTPP";
  my $ALGO = 'TPP'; # taxit/pplacer/guppy

  my $taxon = $faa->{taxon};

  $debug && print STDERR "$subn: \$faa=\n".Dumper($faa);
  my $genotypes = [];
  my $errcode = '';

  # setup the sequence object
  my $seqName = $faa->{seqName};

  my $temp_plp = '';
  my $deflines;
  print STDERR "$subn: taxon=$taxon \tlong_name=$faa->{long_name}\tacc=$faa->{accession}\tALGO=$ALGO\n";


#################### Genotype Analysis by TPP ####################
  if ($ALGO eq 'TPP') {

      $debug && print STDERR "$subn: entering 'TPP'\n";
      my $genotypeHome = Genotype_Def::getGenotypeHome();

      # Create symbolic links to reference files
      my $ref_fasta = '';
      my $refFiles = [ 'CladeClassn.fasta', 'CladeClassn.node.lookup', 'CladeClassn.stats',
                       'CladeClassn.tree',  'CladeClassn.tree.lookup',
#                       '.refpkg', # Needed when taxit is not run
                     ];
      for my $refFile ( @$refFiles ) {
        $refFile = "${taxon}${refFile}";
        if ($refFile eq $taxon.'CladeClassn.fasta') {
          $ref_fasta = "$tempdir/${refFile}";
          $genotype->{rmsa} = $ref_fasta;
        }
#        $refFile = $genotypeHome . "referenceData/" . lc($taxon) . "/${taxon}${refFile}";
        if (!-e "$genotypeHome/referenceData/".lc($taxon)."/$refFile") {
          print STDERR "$subn: ERROR: Ref File not found: $refFile\n";
          next;
        } elsif (-e "$tempdir/${refFile}") {
          $debug && print STDERR "$subn: INFO: Ref File exits: $tempdir/$refFile\n";
          next;
        }
        my $cmd = " ln -s '$genotypeHome/referenceData/".lc($taxon)."/$refFile' '$tempdir' ";
        my $res = ` $cmd `;
        $debug && ($res) && print STDERR "$subn: \$cmd=$cmd\n";
        $debug && ($res) && print STDERR "$subn: \$res=$res\n";
      }

      my $config = config::Config_Reader->new( fileName => "$genotypeHome/referenceData/HCV-Advanced-IN-text-out-text.xml");
      $config->setValue( 'TempDir', $tempdir );
      $config->setValue( 'ClassifierDir', '' );
      $config->setValue('input', $faa->{file});
      $config->setValue('output', $faa->{file} . '.output');
      $config->setValue('cladinator', $progs->{cladinator});
#      $config->setValue('lookupfile', $config->getValue('TempDir') .'/'. $config->getValue('lookupfile'));
      $debug && print STDERR "$subn: \$config=\n".Dumper($config)."End of \$config\n";
      $debug && print STDERR "$subn: configMap=\n".Dumper($config->getConfigMap())."End of configMap\n";

      # Run ClustalW/MAFFT in profile mode
      my $jobName2 = $faa->{jobName0};
      $jobName2 .= "_$region" if ($region && $region !~ m/^\s*NA\s*$/i);
      my $temp_faa = "$jobName2.faa";
      my $temp_afa = "$jobName2.fasta";
      my $temp_plp = "$jobName2.phylip";
      my $temp_out = "$jobName2.phylip_outfile";
      my $temp_new = "$jobName2.phylip_phyml_tree.txt";
      $debug && print STDERR "$subn: \$temp_faa=$temp_faa \$temp_afa=$temp_afa \$temp_plp=$temp_plp \$temp_new=$temp_new\n";
      $debug && print STDERR "$subn: \$ref_fasta=$ref_fasta\n";
      $debug && print STDERR "$subn: \$faa=\n".Dumper($faa)."End of \$faa\n";

      my $cmd = "cp $dir_path/$faa->{jobName0}.faa $tempdir/$temp_faa";
      system($cmd); # Copy input FASTA file to tempdir
      print STDERR "$subn: cmd='$cmd'\n";
      $errcode = Genotype_util::check_file("$tempdir/$temp_faa");
      if ($errcode) {
        $genotype->{comment} = "fasta input $errcode";
        print STDERR "$subn: \$errcode=$genotype->{comment}\n";
#        return ($genotype, $temp_plp);
      }
      print STDERR "$subn: input sequence: file='$dir_path/$temp_faa'\n";

      $debug && print STDERR "$subn: files in $tempdir: \n".` ls -l $tempdir `;
      my $timeExe = time;
      my $pwd1 = ` pwd `;
      chomp($pwd1);
      chdir($tempdir);
      ($temp_plp, $errcode)
          = Genotype_calctree::run_MSA( $tempdir, $dir_path, $progs, $ref_fasta, $temp_faa, $temp_afa, $temp_plp, $deflines, $ALN_PROG);
      $debug && print STDERR "$subn: files in $tempdir: \n".` ls -l $tempdir `;
      chdir($pwd1);
      $genotype->{MSA_phylip} = $temp_plp;
      $debug && print STDERR "$subn: execution time of run_MSA is: ".(time-$timeExe)." sec\n";

      # Separate alignment to reference MSA and query sequence
      my $i = 0;
      my $printFreq = $debug ? 50 : 300;
      my $in  = Bio::SeqIO->new( -file   => "<$tempdir/$temp_afa", -format => 'fasta',);
      my $outRef  = Bio::SeqIO->new( -file   => ">$tempdir/output.profile.afa", -format => 'fasta',);
      my $outSeq  = Bio::SeqIO->new( -file   => ">$tempdir/seqs.afa", -format => 'fasta',);
      my $seqAligned = '';
      my $nextSeq = $in->next_seq();
      while (my $seq = $nextSeq) {
        $i++;
        $nextSeq = $in->next_seq();
        my $seqAcc = $seq->display_id;
        ($i==1 || $i % $printFreq ==0 || !defined($nextSeq) ) &&
            print STDERR "$subn: i=$i \tquery=$faa->{accession} \tseqAcc=$seqAcc\tlength=" . length($seq->seq) . "\n";
#        if ( $seqAcc eq $faa->{accession} ) {
        if ( !defined($nextSeq) ) { # Query sequence is the last one in MSA
          $outSeq->write_seq($seq);
          $seqAligned = $seq->seq;
          $genotype->{seqAligned} = $seqAligned;
          $faa->{seqAligned} = $seqAligned;
          $debug && print STDERR "$subn: i=$i \tFound query sequence: $faa->{accession}\n";
        } else {
          $outRef->write_seq($seq);
        }
      }
      $debug && print STDERR "$subn: After alignment by '$ALN_PROG', created ref MSA='output.profile.afa', seq file='seqs.afa'\n";

      # Setup classifier to run taxit/pplacer/guppy
      my $pwd = `pwd`;
      chomp($pwd);
      chdir($tempdir);
      my $type     = $config->getValue("Type");
      my $typeTree = $config->getValue("Type_tree");
      my $REFERENCE_POST_MAPPING_FILE = "$genotypeHome/referenceData/". lc($type). "/${type}CladeClassn.post.lookup";
      my $postMapping = _getPostMappingFile($REFERENCE_POST_MAPPING_FILE);
      $faa->{postMapping} = $postMapping;
      $debug && print STDERR "$subn: \$postMapping=\n".Dumper($postMapping)."End of \$postMapping\n";

      my $class    = "algo::$type";
      $class .= $typeTree if ( defined($typeTree) );
      $debug && print STDERR "$subn: class=$class\n";
      eval "use $class;"; # Load required module
      $debug && ($@) && print STDERR "$subn: \$\@=$@\n";
      my $classifier = $class->new(
        config   => $config,
        db_conn  => undef,
        blastout => undef,
        count    => 3
        );
      $faa->{classifier} = $classifier;
      $debug && print STDERR "$subn: pwd=".`pwd`."\n";

##############################################
    $debug && print STDERR "$subn: \$genotype->{taxon}=$taxon \$ref_fasta='$ref_fasta'\n";
    my $refAcc = 'DQ278892'; # The sequence used to locate the region, could be anyone in rmsa
    if (scalar(keys %$deflines)==0) {
        $refAcc = ($taxon eq 'HCV') ? '78' #'6w_DQ278892' # 'DQ278892'
#                : ($taxon eq 'ZIKA') ? '1b|AY632535'
                : '';
    } else {
        $refAcc = ($taxon eq 'HCV') ? '78' #'6w_DQ278892' # 'DQ278892'
                : '';
    }
#    $refAcc = '6w_DQ278892'; # Just for files from HCVI
    my $ref_region = {
           '78.disabled' => { # HCV:DQ278892: CE1=1..1155; NS5b=7744..9516, genotype-2.0.1
                     'CE1'  =>   {start =>  1, end   => 1155, },
                     'NS5b' =>   {start => 7744, end   => 9516, },
           },
           # The sampling at C/E1 and NS5B regions have been turned off per ticket BRC-11607, Nov 17, 2017
           '6w_DQ278892.disabled' => { # HCV:DQ278892: CE1=340..1488; NS5b=7621..9393, genotype-1.3.2
                     'CE1'  =>   {start =>  1, end   => 1155, }, #'CE1'  =>   {start =>  340, end   => 1488, },
                     'NS5b' =>   {start => 7744, end   => 9516, }, #'NS5b' =>   {start => 7621, end   => 9393, },
           },
           '6w|DQ278892' => { # HCV:DQ278892: CE1=340..1488; NS5b=7621..9393, genotype-1.3.2
                     'CE1'  =>   {start =>  340, end   => 1488, },
                     'NS5b' =>   {start => 7621, end   => 9393, },
           },
           'DQ278892' => { # HCV:DQ278892: CE1=340..1488; NS5b=7621..9393, genotype-1.3.3
                     'CE1'  =>   {start =>    1, end   => 1149 , }, # These values are from JCVI
                     'NS5b' =>   {start => 7282, end   => 9054 , },
           },
    };
    $debug && print STDERR "$subn: \$refAcc=$refAcc \$ref_region=\n".Dumper($ref_region)."End of \$ref_region\n\n";

    # process defined regions, and get a genotype for each region, eg. CE1 and NS5b for HCV
    if ( ($taxon eq 'HCV' || $taxon eq 'ZIKA')
       ) {
        my $alnio = Bio::AlignIO->new('-file' => "$tempdir/$temp_plp");
        #$debug && print STDERR "$subn: \$alnio = \n".Dumper($alnio)."End of \$alnio\n\n";
        my $aln = $alnio->next_aln();
        $debug && print STDERR "$subn: Creating a new alignment object to read full-length MSA, finished\n";
        #$debug && print STDERR "$subn: \$aln = \n".Dumper($aln)."End of \$aln\n\n";

        my $ref_loc = $ref_region->{$refAcc};
        foreach my $region (sort keys %$ref_loc) {
          my $acc      = $faa->{accession};
          my $long_name = $faa->{long_name};
          my $genotype1 = Genotype_Def::newGenotype( $faa->{job_name}, $acc, $long_name);
          $genotype1->{taxon} = $taxon;
          $genotype1->{accession} = $genotype->{accession};
          $genotype1->{rmsa} = $genotype->{rmsa};

          $genotype1->{region} = $region;
          $debug && print STDERR "$subn: \$refAcc=$refAcc region=$region \$loc=\n".Dumper($ref_loc->{$region})."End of \$loc\n\n";
#          my ($start, $end) = adjustMsaLocation( $aln, $jobName2, substr($refAcc, 0, 10), $ref_loc->{$region});
          my ($start, $end) = adjustMsaLocation( $aln, $faa->{accession}, substr($refAcc, 0, 10), $ref_loc->{$region});
          if ($region eq 'Trimmed') { # No need to adjust for gaps for Trimmed
              ($start, $end) = ($ref_loc->{$region}->{start}, $ref_loc->{$region}->{end});
          }
          $genotype1->{range1} = $start;
          $genotype1->{range2} = $end;
          $end += 1;
          $debug && print STDERR "$subn: region=$region \$start=$start \$end=$end\n";

          my $seqAligned = $faa->{seqAligned};
          $debug && print STDERR "$subn: seqAligned=$seqAligned\n";
          # Contruct sequence by leaving only the window and replacing all others by -
          my $str = '';
          for my $i (0 .. $start-1) { $str .= '-'; }
          $debug && print STDERR "$subn: seqAligned=$str\n";
          $debug && print STDERR "$subn: seqAligned=".length($seqAligned)." start=$start end=$end\n";
          for my $i ($start .. $end) {
            $debug && print STDERR "$subn: seqAligned=".length($seqAligned)." start=$start end=$end i=$i\n";
            $str .= substr($seqAligned, $i, 1);
          }
          $debug && print STDERR "$subn: seqAligned=$str\n";
          for my $i ($end+1 .. length($seqAligned)-1) { $str .= '-'; }
          $debug && print STDERR "$subn: seqAligned=$str\n";

          my $sequence = $str;
          $sequence = substr($sequence, $start, $end-$start+1);
          my $GAP_PCT = 80;
          $debug && printf STDERR "$subn: sequence=$sequence\n";
          my $pct = $sequence;
          $pct = ($pct =~ tr/-Nn//);
          $debug && printf STDERR "$subn:      pct=$pct\n";
          $pct = $pct / length($sequence) *100;
          $debug && printf STDERR "$subn:      pct=$pct%%\n";
          $debug && printf STDERR "$subn: sequence=$sequence\n";
          if ( $pct > $GAP_PCT ) {
            $errcode = ($pct<=99.9) ? sprintf("%2d", $pct+0.5)."%% of $seqName in $faa->{jobName0} is gap \@$start..$end"
                                    : "$seqName in $faa->{jobName0} is empty \@$start..$end";
            $debug && printf STDERR "$subn: skip $errcode\n";
            $genotype1->{comment} = $errcode;
#            next;

          } else {
            my $skipTaxit = 0;
            $genotype1->{seqAligned} = $str;
            $genotype1->{region} = $region;
            $genotype1 = Genotype_TPP::TPP($faa, $genotype1, $faa->{classifier}, $faa->{postMapping}, $reportUnifiedTree, $skipTaxit);
          }
          push @$genotypes, $genotype1;

        }
    }
    $debug && print STDERR "$subn: \$genotypes=\n".Dumper($genotypes)."End of \$genotypes\n";
##############################################

      # get genotype from the entire MSA
      $genotype = Genotype_TPP::TPP($faa, $genotype, $classifier, $postMapping, $reportUnifiedTree);
      push @$genotypes, $genotype;

      # Summarize the genotype from all regions and All
      $debug && print STDERR "$subn: Ready for summ_genotype_HCV taxon=$taxon\n";
      if ($taxon eq 'HCV' #&& $ref_fasta eq $RMSA->[0]->{HCV}
        #|| $taxon eq 'ZIKA' && $ref_fasta eq $RMSA->[0]->{ZIKA}
       ) {
        $genotype->{region} = 'All';
#        push @$genotypes, $genotype;
        $debug && print STDERR "$subn: \$genotypes=\n".Dumper($genotypes)."End of \$genotypes\n";

        $genotypes = summ_genotype_HCV($genotypes);
        $genotype = $genotypes->[0];
        $debug && print STDERR "$subn: \$genotypes=\n".Dumper($genotypes)."End of \$genotypes\n";
      }

      # Save tree from genotype analysis
      my $treeFile = "$jobName2.pplacer.tre";
      if (-e "$tempdir/$genotype->{accession}.tog.tre") {
        `cp $tempdir/$genotype->{accession}.tog.tre $dir_path/$treeFile`;
      } else {
        if (-e "$tempdir/$genotype->{accession}.sing.tre") {
          `cp $tempdir/$genotype->{accession}.sing.tre $dir_path/$treeFile`;
          # Copy pplacer tree to phylip_fastme_tree.new, so that it can be accessed from web download
          `cp $tempdir/$genotype->{accession}.sing.tre $dir_path/$jobName2.phylip_fastme_tree.new`;
        }
      }
      # Save Cladinator output
      if ($genotype->{cladinatorOutput} && -e "$tempdir/$genotype->{cladinatorOutput}") {
        `cp $tempdir/$genotype->{cladinatorOutput} $dir_path/$genotype->{cladinatorOutput}`;
        print STDERR "$subn: INFO: Cladinator output saved to $dir_path/$genotype->{cladinatorOutput}\n";
      } else {
        print STDERR "$subn: ERROR: no Cladinator output was found while in TPP\n";
      }

  } # if ($ALGO eq 'TPP')

  $debug && print STDERR "$subn: \n";
  for (@$genotypes) { $debug && print STDERR Genotype_Def::to_String( $_)."\n"; }

  # Clean up
  # No need to keep phylip file for TPP algorithm
  if (-e "$dir_path/$faa->{jobName0}.phylip") {
    `rm $dir_path/$faa->{jobName0}.phylip `;
    $debug && print STDERR "$subn: Removed $dir_path/$faa->{jobName0}.phylip \n";
  } else {
    $debug && print STDERR "$subn: Kept $dir_path/$faa->{jobName0}.phylip \n";
  }

  $debug && print STDERR "$subn: leaving subroutine\n";
  return $genotypes;
} # sub genotypeTPP


=head1
sub genotype1, takes an input genome object, runs the genotype process, and returns a reference of
 an array of hashes.
Input: Genotype::genotype: $faa=
$VAR1 = {
          'long_name' => 'EU155347',
          'file' => 'EU155347_seq001.faa',
          'job_name' => 'EU155347',
          'taxon' => 'HCV',
          'accession' => '_seq001'
        };
Output: Genotype::genotype: $genotypes=
$VAR1 = [
          {
            'status' => 'Success',
            'long_name' => 'EU155347',
            'Date' => '20130606113717',
            'taxon' => 'HCV',
            'accession' => '_seq001',
            'comment' => '',
            'range1' => '   1',
            'job_name' => 'EU155347',
            'type' => 'genotype',
            'range2' => 9289,
            'BI' => 'BI=1.000   genotype=1a'
            'MSA_phylip' => '',
            'rmsa' => '',
          }
        ];
=cut

sub genotype1 {
    my ($faa, $genotype, $ref_fasta, $deflines, $jobName0, $tempdir, $dir_path, $progs, $no_rerun) = @_;
    my $debug = 0 || Genotype_util::getDebugAll();
    my $subn = "Genotype::genotype1";

    $progs->{clustalw} = "clustalw" if (!exists($progs->{clustalw})|| !$progs->{clustalw});
    $progs->{muscle}   = "muscle"   if (!exists($progs->{muscle})  || !$progs->{muscle});
    $progs->{dnadist}  = "dnadist"  if (!exists($progs->{dnadist}) || !$progs->{dnadist});
    $progs->{fastme}   = "fastme"   if (!exists($progs->{fastme})  || !$progs->{fastme});
    $debug && print STDERR "$subn: \$progs=\n".Dumper($progs)."End of \$progs\n";
    foreach my $k (keys %$progs) {
        my $err = Genotype_util::check_program( $progs->{$k});
        if ($err) {
            print STDERR "$subn: \$progs->{$k}='$progs->{$k}': \$err='$err'\n";
            exit(1);
        }
    }

#   push @$faas, {job_name=>$jobName, accession=>$acc, file=>$temp_faa, long_name=>$long_id};
#    my $taxon = $faa->{taxon};
#    my $jobName = $faa->{job_name};
#    $jobName = 'test' if (!$jobName);
#    my $acc      = $faa->{accession};
#    my $jobName0 = $faa->{file};
#    $jobName0 =~ s/[.]\w{1,7}//i;
#    my $long_name = $faa->{long_name};

#    $debug && print STDERR "$subn: \$faa=\n".Dumper($faa);
#    my $genotypes = [];
    my $errcode = '';
#    # a hash to hold the result
#    my $genotype = Genotype_util::newGenotype( $faa->{job_name}, $acc, $long_name);

    my $seqName = $genotype->{accession};

#    $ref_fasta = '7-15_HCV_3ref_MAFFTaln.fas';
#    $ref_fasta = 'DENV1234-refseqs34_proper_type.afa';
    if (!$ref_fasta ) {
        $genotype->{comment} = "Refseq MSA '$ref_fasta' is empty for $faa->{file} of taxon=$genotype->{taxon}";
        return ($genotype, $temp_plp);
    } elsif (!-e $ref_fasta) {
        $genotype->{comment} = "Refseq MSA '$ref_fasta' not found for $faa->{file} of taxon=$genotype->{taxon}";
        return ($genotype, $temp_plp);
    }

    my $jobName2 = $jobName0;
    $jobName2 .= "_$region" if ($region && $region !~ m/^\s*NA\s*$/i);
    my $temp_faa = "$jobName2.faa";
    my $temp_afa = "$jobName2.fasta";
    my $temp_plp = "$jobName2.phylip";
    my $temp_out = "$jobName2.phylip_outfile";
    my $temp_new = "$jobName2.phylip_phyml_tree.txt";
    $debug && print STDERR "$subn: \$temp_faa=$temp_faa \$temp_afa=$temp_afa \$temp_plp=$temp_plp \$temp_new=$temp_new\n";

    system("cp $dir_path/$jobName0.faa $tempdir/$temp_faa"); # Copy input FASTA file to tempdir
    $errcode = Genotype_util::check_file("$tempdir/$temp_faa");
    if ($errcode) {
        $genotype->{comment} = "fasta input $errcode";
        print STDERR "$subn: \$errcode=$genotype->{comment}\n";
        return ($genotype, $temp_plp);
    }

    print STDERR "$subn: input sequence: file='$dir_path/$temp_faa'\n";

    my $timeExe = time;
    ($temp_plp, $errcode)
        = Genotype_calctree::run_MSA( $tempdir, $dir_path, $progs, $ref_fasta, $temp_faa, $temp_afa, $temp_plp, $deflines, $ALN_PROG);
    $genotype->{MSA_phylip} = $temp_plp;
    $debug && print STDERR "$subn: execution time of run_MSA is: ".(time-$timeExe)." sec\n";

    if ($errcode) {
        $errcode = "phylip alignment $errcode";
        print STDERR "$subn: \$errcode=$errcode\n";
        $genotype->{comment} = $errcode;
        return ($genotype, $temp_plp);
    }

=head1
    # If HCV, run 2 additional analyses using only the CE1 and NS5b regions, based on the MSA of full-length
    # Use DQ278892 and KC248197 as the ref for numbersing for CE1
> src=MSA|ACC=DQ278892|Ver=DQ278892.1|CDS=GI:89113921|ref=NC_004102.1|RM=342..914|Loc=340..912||AA=1..191|symbol=C|Partial=N|product=core protein|*new*
MSTLPKPQRITKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKASERSQPRG
...
LLSCLTVPTSA
> src=MSA|ACC=DQ278892|Ver=DQ278892.1|CDS=GI:89113921|ref=NC_004102.1|RM=915..1490|Loc=913..1488||AA=192..383|symbol=E1|Partial=N|product=E1 protein|*new*
LNYANKSGIYHLTNDCPNSSIVYEAETAILHLPGCVPCVKVANRSKCWVPATPTLAVQDE
...
LLVLLLFAGVDA
> src=MSA|ACC=DQ278892|Ver=DQ278892.1|CDS=GI:89113921|ref=NC_004102.1|RM=7602..9374|Loc=7621..9393||AA=2428..3018|symbol=NS5b|Partial=N|product=NS5B RNA-dependent RNA polymerase|*new*
SMSYSWTGAPITPCAAEEEKLPISPLSNGLLRHHNLVYSTTSRSAPLRQKKVTFDRLQVL
...
GAATLDLSGWFTSGYSGGDIFHSVSYARPRVLLLCLLLLTVGVGIFLLPAR
> src=MSA|ACC=KC248197|Ver=KC248197.1|CDS=GI:525546092|ref=NC_004102.1|RM=342..914|Loc=343..915||AA=1..191|symbol=C|Partial=N|product=core protein|*new*
MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRG
...
LLSCLTVPTSA
> src=MSA|ACC=KC248197|Ver=KC248197.1|CDS=GI:525546092|ref=NC_004102.1|RM=915..1490|Loc=916..1491||AA=192..383|symbol=E1|Partial=N|product=E1 protein|*new*
LEVRNVSGIYHVTNDCSNSSIVYEMDGMILHLPGCVPCVRENNSSRCWVSLTPTLAARDS
...
AVVLMLFAGVDA
> src=MSA|ACC=KC248197|Ver=KC248197.1|CDS=GI:525546092|ref=NC_004102.1|RM=7602..9374|Loc=7624..9396||AA=2428..3018|symbol=NS5b|Partial=N|product=NS5B RNA-dependent RNA polymerase|*new*
SMSYTWTGALITPCAAEETSLPINALSNSLLRHRNLVYSTTSRSAGQRQKKVTFDRLQVL
...
GAANLDLSNWFTAGYSGGDIYHSVSHARPRWFLVCLLLLSVGVGIYLLPNR
=cut

    my $taxon = $genotype->{taxon};
    my $genotypes = [];
    $debug && print STDERR "$subn: \$genotype->{taxon}=$taxon \$ref_fasta='$ref_fasta'\n";
    my $refAcc = 'DQ278892'; # The sequence used to locate the region, could be anyone in rmsa
    if (scalar(keys %$deflines)==0) {
        $refAcc = ($taxon eq 'HCV') ? '6w|DQ278892'
#                : ($taxon eq 'ZIKA') ? '1b|AY632535'
                : ($taxon eq 'ZIKA') ? '1a|AY535||'
                : '';
    } else {
        $refAcc = ($taxon eq 'HCV') ? '_DQ278892'
                : ($taxon eq 'ZIKA') ? '_AY535'
                : '';
    }
    my $ref_region = {
#                   'AY632535||Uganda' => {
                   '1a|AY535||' => { # ZIKA, genotype-1.3.2
                              'CME' =>   {start =>  107, end   => 2476, },
#                              '3UTR' =>   {start =>  10367, end   => 10794, },
                              'NS5' =>   {start => 7655, end   => 10363, },
                                    },
                   '_AY535' => { # ZIKA, genotype-1.3.3
                              'CME' =>   {start =>  107, end   => 2476, },
                              'NS5' =>   {start => 7655, end   => 10363, },
                                    },
                   '6w|DQ278892' => { # HCV:DQ278892: CE1=340..1488; NS5b=7621..9393, genotype-1.3.2
                              'CE1'  =>   {start =>  340, end   => 1488, },
                              'NS5b' =>   {start => 7621, end   => 9393, },
                                    },
                   '_DQ278892' => { # HCV:DQ278892: CE1=340..1488; NS5b=7621..9393, genotype-1.3.3
                              'CE1'  =>   {start =>  340, end   => 1488, },
                              'NS5b' =>   {start => 7621, end   => 9393, },
                                    },
                   '1l|KC248197' => { # HCV:KC248197: CE1=342..1491; NS5b=7602..9374, genotype-1.3.2
                              'CE1'  =>   {start =>  342, end   => 1491, },
                              'NS5b' =>   {start => 7624, end   => 9396, },
                                    },
                   '_KC248197' => { # HCV:KC248197: CE1=342..1491; NS5b=7602..9374, genotype-1.3.3
                              'CE1'  =>   {start =>  342, end   => 1491, },
                              'NS5b' =>   {start => 7624, end   => 9396, },
                                    },
                     };

    # Add a region for MSA after trimming leading/trailing gaps
    my $addTrimMSA = 0;
    if ($addTrimMSA) {
        getTrimRegion( $ref_region, $refAcc, $seqName, $tempdir, $temp_afa);
    }
    $debug && print STDERR "$subn: \$refAcc=$refAcc \$ref_region=\n".Dumper($ref_region)."End of \$ref_region\n\n";

    # process defined regions, and get a genotype for each region, eg. CE1 and NS5b for HCV
    if ( (#$taxon eq 'HCV' || 
          $taxon eq 'ZIKA')
         && $ref_fasta eq $RMSA->[0]->{$taxon}
       ) {
        my $alnio = Bio::AlignIO->new('-file' => "$tempdir/$temp_plp");
        #$debug && print STDERR "$subn: \$alnio = \n".Dumper($alnio)."End of \$alnio\n\n";
        my $aln = $alnio->next_aln();
        $debug && print STDERR "$subn: Creating a new alignment object to read full-length MSA, finished\n";
        #$debug && print STDERR "$subn: \$aln = \n".Dumper($aln)."End of \$aln\n\n";

        my $ref_loc = $ref_region->{$refAcc};
#        my $region = ($taxon eq 'HCV')  ? 'CE1'
#                   : ($taxon eq 'ZIKA') ? 'CME'
#                   : '';
        #foreach my $region (sort {$ref_loc->{$a}->{start}<=>$ref_loc->{$b}->{start}} keys %$ref_loc) {
        foreach my $region (sort keys %$ref_loc) {
          my $acc      = $faa->{accession};
          my $long_name = $faa->{long_name};
          my $genotype1 = Genotype_Def::newGenotype( $faa->{job_name}, $acc, $long_name);
          $genotype1->{taxon} = $taxon;
          $genotype1->{accession} = $genotype->{accession};
          $genotype1->{rmsa} = $genotype->{rmsa};

          $genotype1->{region} = $region;
          $debug && print STDERR "$subn: \$refAcc=$refAcc region=$region \$loc=\n".Dumper($ref_loc->{$region})."End of \$loc\n\n";
          my ($start, $end) = adjustMsaLocation( $aln, $jobName0, substr($refAcc, 0, 10), $ref_loc->{$region});
          if ($region eq 'Trimmed') { # No need to adjust for gaps for Trimmed
              ($start, $end) = ($ref_loc->{$region}->{start}, $ref_loc->{$region}->{end});
          }
          $genotype1->{range1} = $start;
          $genotype1->{range2} = $end;
          $end += 1;
          $debug && print STDERR "$subn: region=$region \$start=$start \$end=$end\n";

          my $jobName = $faa->{job_name};
          $jobName = 'test' if (!$jobName);
          my $errcode = '';
          my $seg_afa = '';
          my $seg_plp = '';
          my $seqExcluded = {};
          $debug && print STDERR "$subn: getting MSA for: \$jobName=$jobName region=$region $start..$end length=".$aln->length()."\n";
          ($seg_afa, $seqExcluded, $errcode) = Genotype_util::writeSliceFasta( $aln, $jobName, $acc, $seqName, $start, $end, $tempdir);
          $debug && print STDERR "$subn: head $tempdir/$seg_afa=\n" . `head -n100 $tempdir/$seg_afa`;
          if ($errcode) {
            $debug && print STDERR "$subn: \$start=$start \$errcode=$errcode\n";
            $genotype1->{BI} = "BI=0.000\tgenotype=NA\tsub=NA";
            $genotype1->{comment} = $errcode;
          } else {
            $seg_plp = sprintf("%s_%d_%d.phylip", "${jobName}_$acc", $start+1, $end+1);
            $seg_plp = sprintf("%s_%s.phylip", "${jobName}_$acc", $region); # saves the phylip file and tree file in more recognizable name
            Genotype_util::fasta2phylip("$tempdir/$seg_afa", "$tempdir/$seg_plp");

            my ($temp_plp1, $temp_new1);
            ($genotype1, $temp_plp1, $temp_new1) =
                Genotype::tree2genotype( $faa, $genotype1, $jobName0, $tempdir, $seg_plp,
                    $deflines, $seqExcluded, $progs, $no_rerun, $dir_path);
            $debug && print STDERR "cat $tempdir/$temp_new1=" . `cat $tempdir/$temp_new1`;
            (!$debug && -f "$dir_path/$temp_new1") && `rm $dir_path/$temp_new1`; # No need to keep the tree file
          }
          push @$genotypes, $genotype1;

        }
    }
    $debug && print STDERR "$subn: \$genotypes=\n".Dumper($genotypes)."End of \$genotypes\n";

    # get genotype from the entire MSA
if ( 1 ) {
    ($genotype, $temp_plp, $temp_new) = Genotype::tree2genotype( $faa, $genotype, $jobName0, $tempdir, $temp_plp, $deflines, $seqExcluded, $progs, $no_rerun, $dir_path);
        $debug && print STDERR "$subn: head $tempdir/$temp_afa=\n" . `head -n100 $tempdir/$temp_afa`;
} else {
=head1
=cut
}

    # summarize the genotype from all regions and All
    if ($taxon eq 'HCV' && $ref_fasta eq $RMSA->[0]->{HCV}
        || $taxon eq 'ZIKA' && $ref_fasta eq $RMSA->[0]->{ZIKA}
       ) {
        $genotype->{region} = 'All';
        push @$genotypes, $genotype;
        $debug && print STDERR "$subn: \$genotypes=\n".Dumper($genotypes)."End of \$genotypes\n";

        $genotypes = summ_genotype_HCV($genotypes);
        $genotype = $genotypes->[0];
        $debug && print STDERR "$subn: \$genotypes=\n".Dumper($genotypes)."End of \$genotypes\n";
    }


#   push @$faas, {job_name=>$jobName, accession=>$acc, file=>$temp_faa, long_name=>$long_id};
    # in the newick tree file, substitute the short seqName with long_name
    $debug && print STDERR "$subn: \$faa=\n".Dumper($faa)."End of \$faa\n";
    $debug && print STDERR "$subn: \$faa->{accession}=$faa->{accession} long_name=$faa->{long_name}\n";
    if (1 && $faa->{accession} ne $faa->{long_name}) {
        # replace the short temporary defline with real defline from input
        open my $file1, '<', "$dir_path/$temp_new" or croak "Can'f open file '$dir_path/$temp_new'";
        open my $file2, '>', "$dir_path/${temp_new}2";
        while (my $inline=<$file1>) {
            chomp($inline);
            next if (!$inline);
            $inline =~ s/$faa->{accession}/$faa->{long_name}/;
            print $file2 "$inline\n";
            $debug && print STDERR "$subn: 2 \$inline=$inline\n";
        }
        $debug && print STDERR "$subn: ${temp_new}2=".`cat $dir_path/${temp_new}2`;
        close $file1;
        close $file2;
        `mv $dir_path/${temp_new}2 $dir_path/${temp_new}`;
#        $debug && print STDERR "$subn: ${temp_new}2=".`cat $dir_path/${temp_new}2`;
#        $debug && print STDERR "$subn: $temp_new=".`cat $dir_path/$temp_new`;
    }
    $debug && print STDERR "$subn: \$genotype=\n".Dumper($genotype);

    $debug && print STDERR "$subn: \n";
    for ($genotype) { $debug && print STDERR Genotype_Def::to_String( $_)."\n"; }

    # For viprbrc.org backend processing only, remove the intermediate files
    if ( 0 ) {
        (-f "$dir_path/$temp_afa") && `rm $dir_path/$temp_afa`;
        (-f "$dir_path/$temp_new") && `rm $dir_path/$temp_new`;
    }
    if ( 1 ) { # no need to keep phylip file or CE1, NS5B files, or new2 file
#        (!$debug && -f "$tempdir/$temp_plp") && `rm $tempdir/$temp_plp`;
        (!$debug && -f "$tempdir/$temp_afa") && `rm $tempdir/$temp_afa`;
        (!$debug && -f "$tempdir/${temp_afa}0") && `rm $tempdir/${temp_afa}0`;
        (!$debug && -f "$tempdir/$temp_faa") && `rm $tempdir/$temp_faa`;
        (!$debug && -f "$tempdir/$temp_out") && `rm $tempdir/$temp_out`;
        (!$debug && -f "$tempdir/$temp_new") && `rm $tempdir/$temp_new`;
        (!$debug && -f "$dir_path/$temp_plp") && `rm $dir_path/$temp_plp`;
        (!$debug && -f "$dir_path/${temp_new}2") && `rm $dir_path/${temp_new}2`;
    }

    $debug && print STDERR "$subn: leaving subroutine\n";
    return ($genotype, $temp_plp);
} # sub genotype1


=head1
sub adjustMsaLocation, takes a start/end of the CE1 or NS5b region, account for the gaps in the MSA,
  returns the location of the start/end of a region in MSA.
=cut

sub adjustMsaLocation {
    my ($aln, $jobName0, $ref_accession, $loc) = @_;
    my $debug = 0 || $debug_all;
    my $subn = "Genotype::adjustMsaLocation";

    $debug && print STDERR "$subn: \$aln=\n".Dumper($aln)."End of \$aln\n";
#    $debug && print STDERR "$subn: \$loc=\n".Dumper($loc)."End of \$loc\n";
    $debug && print STDERR "$subn: \$ref_accession=$ref_accession jobName0=$jobName0\n";
    die "$subn: \$loc seems empty" if (!$loc);
    die "$subn: \$aln seems empty" if (!$aln);

    my $errcode = '';
    my $aln_q; 
    if (0) {
      $aln_q = msa_get_aln( $aln, $ref_accession);
      $debug && print STDERR "$subn: \$aln_q=\n".Dumper($aln_q)."End of \$aln_q\n";
      die "$subn: [ERROR] problem with \$aln_q as \$ref_accession='$ref_accession'" if (!$aln_q);
    } else {
      $aln_q = msa_get_aln( $aln, $jobName0);
      $debug && print STDERR "$subn: \$aln_q=\n".Dumper($aln_q)."End of \$aln_q\n";
      $aln->remove_seq($aln_q) if ($aln_q);
      $debug && print STDERR "$subn: after removing '$jobName0' \$aln=\n".Dumper($aln)."End of \$aln\n";
    }

    $debug && print STDERR "$subn: \$consensus=\n".Dumper($aln->consensus_string(0))."End of \$consensus\n";

    my $qseq;
#    $qseq = $aln_q->seq;
    $qseq = $aln->consensus_string(0);
    $debug && print STDERR "$subn: \$qseq='$qseq'\n";

    my ($istart, $iend) = ($loc->{start}, $loc->{end});
    $debug && print STDERR "$subn: initial values: \$istart=$istart \$iend=$iend\n";
    my $gap = 0;
    my @cs = split(//, $qseq);
    my $c;
    for (my $i=0; $i<$#cs; $i++ ) {
        $debug && print STDERR "$subn: \$i=$i \$gap=$gap \$istart=$istart \$cs[$i]=$cs[$i]\n";
        $gap ++ if ($cs[$i] eq '-');
        last if ($i == $gap + $istart-1);
    }
    $istart += $gap -1;

    $gap = 0;
    for (my $i=0; $i<$#cs; $i++ ) {
        $debug && print STDERR "$subn: \$i=$i \$gap=$gap \$iend=$iend \$cs[$i]=$cs[$i]\n";
        $gap ++ if ($cs[$i] eq '-');
        last if ($i == $gap + $iend-1);
    }
    $iend += $gap -1;

    my ($start, $end) = ($istart, $iend);
    return ($start, $end);
} # sub adjustMsaLocation


=head1
sub getTrimRegion, takes a MSA, look for the gaps of the query sequence in the MSA,
  returns the location of the start/end of query in MSA, excluding leading/trailing gaps
=cut

sub getTrimRegion {
    my ($ref_region, $refAcc, $seqName, $tempdir, $temp_afa) = @_;
    my $debug = 0 || Genotype_util::getDebugAll();
    my $subn = "Genotype::getTrimRegion";

    my $alnio = Bio::AlignIO->new('-file' => "$tempdir/$temp_afa");
    my $aln = $alnio->next_aln();
    $debug && print STDERR "$subn: afa='$tempdir/$temp_afa'\n";
    $debug && print STDERR "$subn: refAcc='$refAcc' seqName='$seqName'\n";
    $debug && print STDERR "$subn: \$aln=".ref($aln) ." has ".$aln->num_sequences." seqs, length=".$aln->length."\n";
    die "$subn: \$aln seems empty" if (!$aln);

    my $errcode = '';
    my $seq_q = msa_get_aln( $aln, $seqName);
    $debug && print STDERR "$subn: \$seq_q=\n".Dumper($seq_q)."End of \$seq_q\n";
    die "$subn: [ERROR] problem with \$seq_q as \$seqName='$seqName'" if (!$seq_q);

    my $qseq = $seq_q->seq;
    $debug && print STDERR "$subn: id=".$seq_q->id()."\tlength=".$seq_q->length()."\tseq='$qseq'\n";

    my ($istart, $iend) = (0, $aln->length-1);
    $debug && print STDERR "$subn: initial values: \$istart=$istart \$iend=$iend\n";
    my @cs = split(//, $qseq);
    my $c;

    # Check leading gaps, also any long gap within 100 bp of start
    my $THRESHOLD = 100; # Look within 100 bp wrt start/end; if there is long gap, move start/end
    my $MINIMUM_LENGTH = 20; # Ignore short (<=20) seqs
    my ($i, $e) = ($istart, $iend);
    my ($lenSeq, $lenGap, $lastGap) = (0, 0, $istart);
    while ($lenSeq<$THRESHOLD) {
        my $ls = 0;
        my $lg = 0;
        while (substr($qseq, $i, 1) ne '-' && $ls<$THRESHOLD) {
            $debug && print STDERR "$subn: lastGap=$lastGap i=$i seq=$ls gap=$lg \$cs[$i]=$cs[$i]\n";
            ++$ls;
            ++$i;
            last if ($i>=$e);
        }
        last if ($lenSeq>=$THRESHOLD);
        while (substr($qseq, $i, 1) eq '-') {
            $debug && print STDERR "$subn: lastGap=$lastGap i=$i seq=$ls gap=$lg \$cs[$i]=$cs[$i]\n";
            ++$lg;
            ++$i;
            last if ($i>=$e);
        }
        if ($lenSeq==0 && $ls==0) { # Leading gap
            $lastGap = $i;
        } elsif (($ls>0 && $lg>$ls*5)) { # only consider long gaps
            $lastGap = $i;
        }
        $lenSeq += $ls;
        $lenGap += $lg;
    }
    $istart = $lastGap;
    $debug && print STDERR "$subn: i=$i lastGap=$lastGap seq=$lenSeq gap=$lenGap \$e=$e \$cs[$i]=$cs[$i]\n";

    # Check trailing gap, and any long gap within 100 bp of end
    $lenSeq = 0;
    $lenGap = 0;
    $lastGap = $iend;
    while ($lenSeq<$THRESHOLD) {
        my $ls = 0;
        my $lg = 0;
        while (substr($qseq, $e, 1) ne '-' && $ls<$THRESHOLD) {
            $debug && print STDERR "$subn: lastGap=$lastGap e=$e seq=$ls gap=$lg \$cs[$e]=$cs[$e]\n";
            ++$ls;
            --$e;
            last if ($i>=$e);
        }
        last if ($ls>=$THRESHOLD);
        while (substr($qseq, $e, 1) eq '-') {
            $debug && print STDERR "$subn: lastGap=$lastGap e=$e seq=$ls gap=$lg \$cs[$e]=$cs[$e]\n";
            ++$lg;
            --$e;
            last if ($i>=$e);
        }
        if ($lenSeq==0 && $ls==0) { # Trailing gap
            $lastGap = $e;
        } elsif (($ls>0 && $lg>$ls*5)) { # only consider long gaps
            $lastGap = $e;
        }
        $lenSeq += $ls;
        $lenGap += $lg;
    }
    $iend = $lastGap;
    $debug && print STDERR "$subn: i=$i lastGap=$lastGap seq=$lenSeq gap=$lenGap \$e=$e \$cs[$e]=$cs[$e]\n";

    my $ok = 0;
    ++$istart;
    ++$iend;
    if ($istart>=1 && $iend<=$aln->length && $istart<=$iend) {
        $ref_region->{$refAcc}->{'Trimmed'}->{start} = $istart;
        $ref_region->{$refAcc}->{'Trimmed'}->{end}   = $iend;
        $ok = 1;
    }
    $debug && print STDERR "$subn: ok=$ok refAcc='$refAcc' \$ref_region=\n".Dumper($ref_region)."End of \$ref_region\n";

    return $ok;
} # sub getTrimRegion


=head1
sub tree2genotype, takes an input MSA in phylip format, builds a phylogeny tree and finds the BI value for one sequence,
 and returns a reference of an array of hashes.
Output: Genotype::genotype: $genotypes=
$VAR1 = [
          {
            'status' => 'Success',
            'long_name' => 'EU155347',
            'Date' => '20130606113717',
            'taxon' => 'HCV',
            'accession' => '_seq001',
            'comment' => '',
            'range1' => '   1',
            'job_name' => 'EU155347',
            'type' => 'genotype',
            'range2' => 9289,
            'BI' => 'BI=1.000   genotype=1a'
            'MSA_phylip' => '',
            'rmsa' => '',
          }
        ];
=cut

sub tree2genotype {
#($faa, $genotype, $ref_fasta, $jobName0, $tempdir, $dir_path, $progs, $no_rerun) = @_;
    my ($faa, $genotype, $jobName0, $tempdir, $temp_plp, $deflines, $seqExcluded, $progs, $no_rerun, $dir_path) = @_;
    my $debug = 0 || $debug_all || Genotype_util::getDebugAll();
    my $subn = "Genotype::tree2genotype";

    my $errcode = '';
    my $taxon = $genotype->{taxon};
    my $seqName = $genotype->{accession};

    my $jobName2 = $jobName0;
    $jobName2 .= "_$region" if ($region && $region !~ m/^\s*NA\s*$/i);
    my $temp_new = "$jobName2.phylip_phyml_tree.txt";
    $debug && print STDERR "$subn: \$temp_plp=$temp_plp \$temp_new=$temp_new\n";
    #$debug && print STDERR "$subn: \$deflines=\n".Dumper($deflines)."End of \$deflines\n";

    $debug && print STDERR "$subn: input sequence: file='$dir_path/$temp_plp'\n";

    # Generate phylogeny tree
    my $timeExe = time;
    if ( 1 ) {
        $timeExe = time;
        ($temp_new, $errcode) = Genotype_calctree::calc_tree( $faa, $tempdir, $temp_plp, $progs, $no_rerun, $dir_path, $debug);
        $debug && print STDERR "$subn: execution time of calc_tree is: ".(time-$timeExe)." sec\n";
        # Test to see if the PAUP alignment would give identical result wrt VBRC or rivm.nl
    } else {
            # In case there is alignment file from other source
            my $in = Bio::TreeIO->new(-file => "$dir_path/alignment_tree.txt", -format => 'nexus');
            my $out= Bio::TreeIO->new(-file => ">$tempdir/$temp_new", -format => 'newick');
            while( my $t = $in->next_tree ) { $out->write_tree($t); }
    }
    #$debug && print STDERR "$subn: \$deflines=\n".Dumper($deflines)."End of \$deflines\n";
    if (scalar(keys %$deflines)>0) {
        my $defline_fail = Genotype_calctree::restoreDefline($deflines, $seqExcluded, "$tempdir/$temp_new");
        (scalar(keys %$defline_fail)>0) && print STDERR "$subn: \$defline_fail=\n".Dumper($defline_fail)."End of \$defline_fail\n";
    }

    # check the phylogeny tree result
    if ($errcode) {
        print STDERR "$subn: \$errcode=$errcode\n";
        $genotype->{comment} = $errcode;
        return ($genotype, $temp_plp, $temp_new);
    } else {
    }
    #save the tree (newick format) file for output
    $debug && print STDERR "$subn: Saving $tempdir/$temp_new to $dir_path\n";
    system("cp $tempdir/$temp_new $dir_path");  # Save the Phylogeny tree
    print STDERR "$subn: \$taxon=$taxon Phylogeny tree: file='$dir_path/$temp_new'\n";


    # Run newickBI.pl to get branching index (BI) value and genotype
#    $cmd = "./phyloplace/newickBI.pl $temp_new \"AF011753\" ";
#    $cmd = "./phyloplace/newickBI.pl $tempdir/$temp_new \"$acc\" ";
#    $debug && print STDERR "Command to run newickBI.pl is:\n$cmd\n";
    my $msgs = "BI=\tgenotype=NA";
    if (0) { # Use the tree from VBRC.org, to see if newickBI generates same result. It did
        $temp_new = "AF289029-whole_genome-gaps_excluded.tre.new";
        system("cp $dir_path/$temp_new $tempdir");
    }
    $timeExe = time;
    $debug && print STDERR "$subn: Entering Genotype_newickBI::newickBI, \$seqName=$seqName\n";
    $msgs = Genotype_newickBI::newickBI("$tempdir/$temp_new", $seqName, $taxon, $debug);
    $debug && print STDERR "$subn: execution time of newickBI is: ".(time-$timeExe)." sec\n";
    chomp($msgs);

    # Change the format of genotype for ORF1 for NOV genomes
    if ($genotype->{'taxon'} =~ m/^NOV/ && $genotype->{'region'} eq 'ORF1' && $msgs =~ m/(I[.])[^P]/) {
        $old = $1;
        $msgs =~ s/I[.]/I.P/g;
        $debug && print STDERR "$subn: \$msgs changed: \$msgs='$msgs'\n";
    }

    print STDERR "$subn: sub newickBI returned: \$msgs='$msgs'\n";
    for my $msg (split /\n/, $msgs) {
        chomp($msg);
        $debug && print STDERR "$subn: \$msg='$msg'\n";
        $debug && print STDERR "$subn: before &convertGT: \$msg='$msg'\n";
        $msg = Genotype_Def::convertGT( $msg, $genotype->{taxon});
        $debug && print STDERR "$subn: after  &convertGT: \$msg='$msg'\n";

        $msg = Genotype_Def::convertSubtype( $msg, $genotype->{taxon});
        $debug && print STDERR "$subn: after  &convertSubtype: \$msg='$msg'\n";

        ($genotype->{status}, $msg, $genotype->{subtype}, $genotype->{comment}) = Genotype_util::parseBI( $msg, $taxon);
        $genotype->{Date} = Genotype_util::getTime();
        $genotype->{BI} = $msg;
#        if ($msg =~ /^(BI=([.0-9]+)\tgenotype=($Genotype_Def::fmtGenotype)\t(sub=([^\s]*)))\s*(.*)$/i) {
#            $debug && print STDERR "$subn: \$1='$1' \$4='$4'\n";
#            $genotype->{BI} = $1;
##            $genotype->{subtype} = $4;
#        }
        $debug && print STDERR "$subn: \$genotype=\n".Dumper($genotype)."End of \$genotype\n";
#        push @$genotypes, $genotype;

        $debug && print STDERR "$subn: checking Taxon: \$taxon='$taxon'\n";
        next if ($taxon !~ m/^NOV([I]*)/i);
        my $gg = $1;
        $debug && print STDERR "$subn: checking Taxon: \$taxon='$taxon', as \$gg='$gg'\n";
#        if ($genotype->{BI} =~ m/^BI=([.0-9]+)\tgenotype=(([0-9]\w+|\d-[A-Za-z]+|NA|.+))/i) {
        if ($genotype->{BI} =~ m/^BI=([.0-9]+)\tgenotype=($fmtGenotype)/i) {
            my $gt = $2;
            $debug && print STDERR "$subn: Taxon changed: \$taxon='$taxon', as \$gg='$gg' \$gt='$gt'\n";
            next if ($gt !~ m/^([I]+)[.]/i);
            $taxon = 'NOV' . $1 if ($gg ne $1);
            $genotype->{taxon} = $taxon;
            print STDERR "$subn: Taxon changed: \$taxon='$taxon', as \$gg='$gg' \$gt='$gt'\n";
        }

    }

    $debug && print STDERR "$subn: \n";
    for ($genotype) { $debug && print STDERR Genotype_Def::to_String( $_)."\n"; }

    $debug && print STDERR "$subn: leaving subroutine\n";
    return ($genotype, $temp_plp, $temp_new);
} # sub tree2genotype

=head1
sub summ_genotype_NOV looks through all genotypes, and determines the dominant genotypes among the analyses
  Prepends result to the input array
=cut

sub summ_genotype_NOV {
    my ($genotypes, $seq_length) = @_;
    my $debug = 0 || $debug_all || Genotype_util::getDebugAll();
    my $subn = "Genotype::summ_genotype_NOV";

    $debug && print STDERR "$subn: \$genotypes=\n".Dumper($genotypes);
    if ($#{$genotypes}!=3) {
        print STDERR "$subn: $#{$genotypes} genotypes in input, requires 4. Abort.\n";
        return [];
    }

    for my $i (1, 0) {
        $debug && print STDERR "$subn: genotype #$i\n";
        my $genotype = $genotypes->[$i];
        my $genotype2 = pop @$genotypes;
        $debug && print STDERR "$subn: genotype #$i='$genotype->{'BI'}'.\n";
        $debug && print STDERR "$subn: genotype #". ($i+2) ."='$genotype2->{'BI'}'.\n";
        if ($genotype->{'BI'} !~ m/genotype=II[.]P*4/) {
            $debug && print STDERR "$subn: genotype #$i='$genotype->{'BI'}' is not of II.4, skip.\n";
            next;
        }

        my $BI = '';
        my $BI2 = '';
        $BI = $1 if ($genotype->{'BI'} =~ m/sub=(.*)$/i);
        $BI2 = $1 if ($genotype2->{'BI'} =~ m/sub=(.*)$/i);
        $debug && print STDERR "$subn: \$BI='$BI' \$BI2='$BI2'.\n";
        if ($BI && $BI ne $BI2) {
            $genotype->{'BI'} =~ s/$BI/$BI2/;
            print STDERR "$subn: sub has been changed: genotype #$i='$genotype->{'BI'}'.\n";
        }
    }

    return $genotypes;
} # sub summ_genotype_NOV


=head1
sub summ_genotype_HCV looks through all genotypes, and determines the dominant genotypes among the analyses
  Prepends result to the input array
=cut

sub summ_genotype_HCV {
    my ($genotypes, $seq_length) = @_;
    my $debug = 0 || $debug_all || Genotype_util::getDebugAll();
    my $subn = "Genotype::summ_genotype_HCV";

    my $TH_BI = 0.711;
    $TH_BI = 0.4 if ($debug);
    $seq_length = 1;
    # take the required min. # of windows for a genotype to be considered as present as one less than square root of
    # total number of windows from the input sequence (not MSA), and the min. value is 1.
    #my $TH_min_num_windows = sqrt($seq_length/100 -2) -1 ;
    my $TH_min_num_windows = $seq_length/100 ;
    $TH_min_num_windows -= 2 if ($TH_min_num_windows>3);
    $TH_min_num_windows = sqrt($TH_min_num_windows) -1;
    $TH_min_num_windows = ($TH_min_num_windows<1) ? 1 : $TH_min_num_windows;
    $debug && print STDERR "$subn: \$seq_length=$seq_length \$TH_min_num_windows=$TH_min_num_windows\n";

    my $genotype_summary = [];
    my $sub_summary = {};
    my $comm = '';
    my $comm2 = '';
    my $genotype0 = Genotype_Def::newGenotype();
    $genotype0->{job_name}  = $genotypes->[$#{$genotypes}]->{job_name};
    $genotype0->{accession} = $genotypes->[$#{$genotypes}]->{accession};
    $genotype0->{long_name} = $genotypes->[$#{$genotypes}]->{long_name};
    $genotype0->{type} = 'genotype';
    $genotype0->{range1} = $genotypes->[$#{$genotypes}]->{range1};
    $genotype0->{range2} = $genotypes->[$#{$genotypes}]->{range2};
    $genotype0->{taxon}  = $genotypes->[$#{$genotypes}]->{taxon};
    $genotype0->{region} = 'NA';
    $genotype0->{Date} = Genotype_util::getTime();
    $genotype0->{status} = 'Failed';
    $genotype0->{comment} = $comm;
    $genotype0->{comment2} = $comm2;
    unshift @$genotypes, $genotype0;

    my $hasTrimmed = 0;
    for my $i (1 .. $#{$genotypes}) {
        next if ($genotypes->[$i]->{region} ne 'Trimmed');
        $hasTrimmed = 1;
        $debug && print STDERR "$subn: Found region 'Trimmed' \$i=$i \$genotype=\n".Dumper($genotypes->[$i])."\n";
    }

    my $maxBI = '0.000'; #-1;
    my $maxGT = '';
    my $allGT = '';
    for my $i (1 .. $#{$genotypes}) {
        my $genotype = $genotypes->[$i];
        # Skil 'All' if we have 'Trimmed'
        if ($hasTrimmed && $genotype->{region} eq 'All') {
            print STDERR "$subn: \$i=$i Skipping the 'All' result, in favored of result from 'Trimmed'\n";
            next;
        }
        $debug && print STDERR "$subn: \$i=$i \$genotype=\n".Dumper($genotype)."\n";

        # take cladinator output from region='All'
        if ($genotype->{region} eq 'All') {
          $genotype0->{cladinatorOutput} = $genotype->{cladinatorOutput};
        }

        my $BI = $genotype->{BI};
        $debug && print STDERR "$subn: \$i=$i region:$genotype->{region} range:$genotype->{range1}..$genotype->{range2} \$BI='$BI'\t";
#        if (!($BI =~ /^BI=([.0-9]+)\tgenotype=([0-9][A-Za-z]*|\d-[A-Za-z]+|.+)\tsub=(\w*)\t(.*)/
#                   || $BI =~ /^BI=([.0-9]+)\tgenotype=([0-9][A-Za-z]*|\d-[A-Za-z]+|.+)\tsub=(\w*)/)) {
        if (!($BI =~ /^BI=([.0-9]+)\tgenotype=($fmtGenotype|NA)\tsub=(\w*)\t(.*)/
                   || $BI =~ /^BI=([.0-9]+)\tgenotype=($fmtGenotype|NA)\tsub=(\w*)/)) {
            $debug && print STDERR ": Problem i=$i! BI=$BI\n";
            next;
        }

        $debug && print STDERR "\$1='$1' \$2='$2' \$3='$3'\n";
        ($maxBI, $maxGT) = ($1, $2) if ($maxBI<$1);
        # save the raw result in comment
        my $region = '';
          $region = $genotype->{region};
          $debug && print STDERR "$subn: \$i=$i \$genotype->{comment}='$genotype->{comment}'\n";
        $comm2 .= '|' if ($comm2);
        $comm2 .= "$region=$2:$1";
        if (exists($genotype->{comment}) && $genotype->{comment}) {
            $comm .= '|' if ($comm);
            $comm .= "$region=$genotype->{comment}";
            $debug && print STDERR "$subn: \$i=$i \$comm='$comm'\n";
        }
        if (0) {
            $genotype_summary->{$2}++ if ($1>=$TH_BI);
        } else {
            my $seen = 0;
            for my $j (0 .. $#{$genotype_summary}) {
              $debug && print STDERR "$subn: \$j=$j \$genotype_summary->[$j]->[0]=$genotype_summary->[$j]->[0] \$2=$2\n";
              my $oldGenotype = $genotype_summary->[$j]->[0];
if (0) {
              if (length($oldGenotype)>=length($2) && substr($oldGenotype, 0, length($2)) ne $2) {
                next;
              } elsif (length($oldGenotype)<length($2) && substr($2, 0, length($oldGenotype)) ne $oldGenotype) {
                next;
              }
              $genotype_summary->[$j]->[1] = $1 if ($1>$genotype_summary->[$j]->[1]);
              $seen = 1;
              last;
}
              if (length($oldGenotype)>=length($2)) {
                if (substr($oldGenotype, 0, length($2)) eq $2) {
                  #$genotype_summary->[$j]->[1] = $1 if ($1 > $genotype_summary->[$j]->[1]);
                  $seen = 1;
                  last;
                }
              } else {
                if (substr($2, 0, length($oldGenotype)) eq $oldGenotype) {
                  $genotype_summary->[$j]->[1] = $1;
                  $genotype_summary->[$j]->[0] = $2;
                  $seen = 1;
                  last;
                }
              }
            }
            if (!$seen) {
                push @$genotype_summary, [$2, $1];
            }
            $debug && print STDERR "$subn: \$genotype_summary=\n".Dumper($genotype_summary)."\n";
        }
        $sub_summary->{$3}++ if ($1>=$TH_BI);
    }
    $genotype0->{comment} = $comm;
    $genotype0->{comment2} = $comm2;
    $debug && print STDERR "$subn: \$genotype0=\n".Dumper($genotype0)."\n";
    $debug && print STDERR "$subn: \$genotype_summary=\n".Dumper($genotype_summary)."\n";
    $debug && print STDERR "$subn: \$genotype_summary=\n".Dumper($sub_summary)."\n";
    my $BI = "BI=0.000\tgenotype=";
    for my $j (0 .. $#{$genotype_summary}) {
        $k = $genotype_summary->[$j]->[0];
        $debug && print STDERR "$subn: \$k=$k\n";
        next if ($k eq 'NA');
        $genotype0->{status} = 'Success'; # as long as there is an assignment, it's a Success
        next if ($genotype_summary->[$j]->[1] < $TH_BI);
        $allGT .= '/' if ($allGT);
        $allGT .= $k;
        $debug && print STDERR "$subn: \$allGT='$allGT'\n";
    }
    if ($maxBI>=$TH_BI) {
        $BI = "BI=$maxBI\tgenotype=$allGT";
    } else {
        $BI = "BI=$maxBI\tgenotype=$maxGT";
    }
    $BI =~ s/genotype=\t/genotype=NA\t/i;
    $BI .= "\tsub=";
    for $k (sort keys %$sub_summary) {
        $debug && print STDERR "$subn: \$k=$k:$sub_summary->{$k} \$TH_min_num_windows=$TH_min_num_windows\n";
        next if ($k eq 'NA');
        next if ($sub_summary->{$k} < $TH_min_num_windows);
#        $BI .= ',' if ($BI =~ m/genotype=([0-9][A-Za-z]*|\d-[A-Za-z]+|.+)\s+sub=\w+$/);
        $BI .= '/' if ($BI =~ m/genotype=$fmtGenotype\s+sub=\w+$/);
        $BI .= $k;
#        $genotype0->{status} = 'Success';
        $debug && print STDERR "$subn: \$BI='$BI'\n";
    }
    $BI =~ s/genotype=\t/genotype=NA\t/i;
    $debug && print STDERR "$subn: \$BI='$BI'\n";
    $BI =~ s/sub=$/sub=NA/i;
    $genotype0->{BI} = $BI;
    $genotype0->{comment} = '' if ($genotype0->{status} !~ m/Fail/i);
    $debug && print STDERR "$subn: \$genotype0=\n".Dumper($genotype0)."\n";

    my $KEEP_analyses3 = 0;
    if ( !$KEEP_analyses3 ) {
        $genotypes = [ $genotypes->[0] ];
    }
    $debug && print STDERR "$subn: \$genotypes=\n".Dumper($genotypes)."\n";
    $debug && print STDERR "$subn: leaving subroutine\n";
    return $genotypes;
} # sub summ_genotype_HCV

=head2 msa_get_aln
Takes an alignment, and an id, return the gaps within the alignment with such id
=cut

sub msa_get_aln {
    my ($aln, $id) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'msa_get_aln';
    my $refaln = [ $aln->each_seq_with_id($id) ]; # 
    $debug && print "$subn: \$refaln=\n".Dumper($refaln)."End of \$refaln\n\n";
    if ($#{$refaln} < 0) {
        $debug && print "$subn: Couldn't find id=$id in alignment file.\n";
        return undef;
    }
    my $seq = $refaln->[0];
    $debug && print "$subn: \$refaln = ".$seq->display_id()."\n";
    $debug && print "$subn: \$refaln = ".$seq->start()."\n";
    $debug && print "$subn: \$refaln = ".$seq->end()."\n";
    $debug && print "$subn: \$refaln = ".$seq->alphabet()."\n";
    $debug && print "$subn: \$refaln = ".$seq->length."\n";
    $debug && print "$subn: \$refaln = ".$seq->seq."\n";

    my @gaps = (''); # we don't need element 0, but need it to suppress error from print
    $seq = $refaln->[0];
    $debug && print "$subn: \$seq='$seq'\n";

    return ($seq);
} # sub msa_get_aln


sub _getPostMappingFile {
  my ($file) = @_;

  my $postMapping = {};
  return $postMapping if ( !-e $file );
  $debug && print STDERR "_getPostMappingFile: Acquiring post mapping file:  $file\n";

  my $fh = new FileHandle;
  $fh->open( $file, '<' );
  while ( !$fh->eof ) {
    my $line = $fh->getline;
    chomp($line);
    next if ( util::Constants::EMPTY_LINE($line) );
    my ( $clade, $post_clade ) = split( /\t/, $line );
    $postMapping->{$clade} = $post_clade;
  }
  $fh->close;
  return $postMapping;
}


1;
