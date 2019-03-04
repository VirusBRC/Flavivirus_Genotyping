package Genotype_calctree;

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
use version; our $VERSION = qv('1.3.0'); # Jun 06, 2013
#use Getopt::Long;
use Annotate_Def;
use Genotype_util;

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

    my $newHCV_ICTV = 1; # Sep 30, 2013
    my $fmtGenotype = '\d[-._\dA-Za-z]*';
    $fmtGenotype = '([0-9]\w+|\d-[A-Za-z]+|NA|.+)' if (!$newHCV_ICTV);

############## Subroutines ################


=head1
sub run_MSA, runs MUSCLE or ClustalW,
 returns the file name of the alignment.
=cut

sub run_MSA {
    my ($tempdir, $dir_path, $progs, $ref_rmsa, $temp_faa, $temp_afa, $temp_plp, $deflines, $ALN_PROG) = @_;
    my $debug = 0 || $debug_all || Genotype_util::getDebugAll();
    my $subn = "run_MSA";

    $progs->{muscle}  = "muscle"  if (!exists($progs->{muscle}) || !$progs->{muscle});

    my $errcode = '';

    # Run MUSCLE with ref alignment and input sequence
    if ( 1 ) {
        if ( $ALN_PROG eq 'mafft' ) {
#            $ref_rmsa =~ s/muscle[.]aln/mafft.aln/;
        }
        print STDERR "$subn: ref MSA:'$ref_rmsa'\n";

        # Run MUSCLE with fasta output, which is then converted to phylip format via seqConverter.pl
        my $timeExe = time;
        my $run_result = '';
        # Copy rmsa to temdir
        my $fn = ($ref_rmsa =~ m/([^\\\/]+)$/i) ? $1 : '';
        $debug && print STDERR "$subn: \$fn='$fn'\n";
        system("cp $ref_rmsa $tempdir") if (!-e "$tempdir/$fn");
        my $temp_rmsa = ($ref_rmsa =~ m/([^\/]+)$/) ? $1 : "";
        $debug && print STDERR "$subn: \$temp_rmsa='$temp_rmsa'\n";

        # replace defline in phylip file with only accession, and only first 10 letters
        # save the mapping in order to restore the deflines in tree file 
        if (scalar(keys %$deflines)>0) {
            &cleanupDefline($deflines, "$tempdir/$temp_rmsa");
        }

        # Run MUSCLE or ClustalW or MAFFT
        my $cmd = "";
        if ( $ALN_PROG eq 'mafft' ) {
          if ( 0 ) {
            # mafft-profile
            $cmd = "$progs->{mafft} ";
            $cmd .= "$tempdir/$temp_rmsa ";
            $cmd .= "$tempdir/$temp_faa ";
            $cmd .= "> $tempdir/$temp_afa";
          } else {
            # mafft --add
            $cmd = "mafft --keeplength ";
#            $cmd .= "--maxiterate 1000 --localpair "; # according to Christian Zmasek
            $cmd .= "--add $tempdir/$temp_faa ";
            $cmd .= "$tempdir/$temp_rmsa ";
            $cmd .= "> $tempdir/$temp_afa";
          }
        } elsif ( $ALN_PROG eq 'muscle' ) {
            $cmd = "$progs->{muscle} ";
            $cmd .= "-profile -in1 $tempdir/$temp_rmsa -in2 $tempdir/$temp_faa ";
            $cmd .= "-out $tempdir/$temp_afa";
        } else {
            # clustalw -profile2=AB049153.input -profile1=AB049153.A_1.afa
            # -output=aln -align -gapopen=20 -gapext=1 -type=dna -outfile=AB049153.output
            $cmd = "$progs->{clustalw} ";
            $cmd .= "-profile1=$tempdir/$temp_rmsa -profile2=$tempdir/$temp_faa ";
            $cmd .= "-align -gapopen=20 -gapext=1 -type=dna ";
            $cmd .= "-output=fasta "; # -output=gcg OR gde OR pir OR phylip OR nexus OR fasta
            $cmd .= "-outfile=$tempdir/$temp_afa";
        }
        $cmd .= (!$debug && -e '/dev/null' || $ALN_PROG eq 'mafft') ? ' 2>/dev/null' : ' >&1'; # Throw away stdout and stderr if not debug
        print STDERR "$subn: \$cmd='$cmd'\n";

        $run_result = `$cmd`;

        $debug && print STDERR "$subn: \$run_result='$run_result'\n";
        $debug && print STDERR "$subn: execution time is: ".(time-$timeExe)." sec\n";
        if ( $ALN_PROG eq 'mafft'|| $ALN_PROG eq 'muscle'  || $ALN_PROG eq 'clustalw' ) {
          if ( 0 ) {
            $temp_afa = 'AF289029.fasta_aln'; # MSA downloaded from vbrc.org
            system("cp $dir_path/$temp_afa $tempdir"); # Use the alignment file from VBRC.org
          }
    
          if ( 0 ) {
            # re-run the alignment, trying to resolve the problem with Dengue genomes AY100465, AY732474, AY732477, & AY732478
            system("$progs->{muscle} -in $tempdir/$temp_afa -out $tempdir/${temp_afa}2");
            system("mv $tempdir/${temp_afa}2 $tempdir/$temp_afa");
          }
    
          $errcode = Genotype_util::check_file("$tempdir/$temp_afa", 3, $debug);
          if ($errcode) {
            $errcode = "muscle alignment $errcode";
            print STDERR "$subn: \$errcode=$errcode\n";
            $genotype->{comment} = $errcode;
            push @$genotypes, $genotype;
            return ($temp_plp, $errcode);
          }
          system("cp $tempdir/$temp_afa $dir_path"); # Save the alignment file in FASTA
          $debug && print STDERR "$subn: MSA in Fasta  : file='$dir_path/$temp_afa'\n";

          # Convert fasta alignment file to an extended phylip file
          if (0) {
            $cmd = "./seqConverter.pl -d$tempdir/$temp_afa -ope -v";
            $run_result = `$cmd`; # Run seqConverter.pl
            $debug && print STDERR "$subn: \$run_result='$run_result'\n";

          } else {
            # Test to see if the PAUP alignment would give identical result wrt VBRC or rivm.nl
#            `cp $dir_path/alignment.fasta $tempdir/$temp_afa`; # In case there is alignment file from other source
            Genotype_util::fasta2phylip("$tempdir/$temp_afa", "$tempdir/$temp_plp");
          }

          # Save fasta alignment file to a FastME 2.15 format file
          if ( 0) {
            my $fastme2file = Genotype_util::fasta2FastME("$tempdir/$temp_afa");
            $debug && print STDERR "$subn: \$fastme2file='$fastme2file'\n";
            $debug && print STDERR "$subn: \$dir_path='$dir_path'\n";
            if (-e $fastme2file) {
                print STDERR "$subn: copy FastME file: ". `cp $fastme2file $dir_path` ."\n";
            }
            #my $errcode = Genotype_util::check_file($fastme2file);
          }
        }
        my $msg = `cp $tempdir/$temp_plp $dir_path`; # Save the alignment file in PHYLIP
        ($msg) && print STDERR "$subn: $msg\n";
        $msg = "MSA(Fasta):";
        $msg .= (-e "$dir_path/$temp_afa") ? "$dir_path/$temp_afa" : "-N/A-";
        $msg .= "; MSA(phylip):";
        $msg .= (-e "$dir_path/$temp_plp") ? "$dir_path/$temp_plp" : "-N/A-";
        print STDERR "$subn: $msg\n";

    } else {
        # Run MUSCLE, with output in extended phylip format
        # However, the phylip output from MUSCLE can only have ids up to 10 letter-long, longer ids are truncated
        my $cmd = "$progs->{muscle} -profile -in1 $ref_rmsa -in2 $tempdir/$temp_faa -phyiout $tempdir/$temp_plp\n";
        system($cmd);
    }
    #save the alignment (phylip format) file for testing
    $debug && system("cp $tempdir/$temp_plp $dir_path");
    # Check the Phylip alignment file
    $errcode = Genotype_util::check_file("$tempdir/$temp_plp", 2, $debug);
    if ($errcode) {
            $errcode = "phylip alignment $errcode";
            print STDERR "$subn: \$errcode=$errcode\n";
#            $genotype->{comment} = $errcode;
#            push @$genotypes, $genotype;
#            next;
    }

    $debug && print STDERR "$subn: leaving subroutine\n";
    return ($temp_plp, $errcode);
} # sub run_MSA


=head1
sub run_dnadist, takes an input alignment in phylip format, runs either PhyML or dnadist/FastME,
 returns the file name of the phylogeny tree or undef.
=cut

sub run_dnadist {
    my ($tempdir, $temp_plp, $progs, $no_rerun, $dir_path, $debug) = @_;
    $debug = 0 || $debug_all || Genotype_util::getDebugAll();
    my $subn = "Genotype_calctree::run_dnadist";

    $progs->{dnadist} = "dnadist" if (!defined($progs) || !$progs->{dnadist});

    my $temp_new = '';
    $debug && print STDERR "$subn: \$temp_plp='$temp_plp'\n";
    my $errcode = '';

    # generate phylogeny tree from phylip-style alignment
    my $cmd;
    # Using dnadist/fastme to generate phylogeny tree
    # Note: FastME can only accommodate 10-character long ids, so be careful
    $temp_new = "$temp_plp"."_fastme_tree.new";
    if (-e "$tempdir/$temp_new" && $no_rerun) {
        system("cp $dir_path/$temp_new $tempdir");
        print STDERR "$subn: Found existing phylogeny tree file: $dir_path/$temp_new. FastME not run.\n";
    } else {
        my $pwd = `pwd`;
        chomp($pwd);
        $debug && print STDERR "$subn: \$tempdir='$tempdir'\n";
        chdir($tempdir) or die "$!";
        $debug && print STDERR "$subn: \$pwd=$pwd\n";
        my $cp_result = `cp $tempdir/$temp_plp $tempdir/infile`;
        ($cp_result) && print STDERR "$subn: cp result=$cp_result\n";
        $debug && system("cp $tempdir/$temp_plp $dir_path"); # Save the alignment file for debugging
        $errcode = Genotype_util::check_file("$tempdir/infile", 2, $debug);
        if ($errcode) {
            $errcode = "dnadist input $errcode";
            print "$subn: \$errcode=$errcode\n";
            chdir($pwd) or die "$!";
            return ($temp_new, $errcode);
        }

        # Run dnadist to get distance matrix from MSA, dnadist takes input from "infile" and outputs to "outfile"
        my $timeExe = time;
        `rm $tempdir/outfile` if (-e "$tempdir/outfile"); # Remove existing outfile, required by dnadist
        my $result;
        $cmd = "cd $tempdir; echo 'Y' | $progs->{dnadist}";
        $cmd .= (!$debug && -e '/dev/null') ? ' 2>/dev/null' : ' >&1'; # Throw away stdout and stderr if not debug
        $debug && print STDERR "$subn: starting to run: \$cmd='$cmd'\n";
        $run_result = `$cmd`;
        system("cp $tempdir/outfile $pwd/${temp_plp}_outfile"); # Save the distance matrix file for debugging
        $debug && print STDERR "$subn: execution time is: ".(time-$timeExe)." sec\n";
        $debug && print STDERR "$subn: \$run_result=$run_result\n";
        $errcode = Genotype_util::check_file("$tempdir/outfile", 0, $debug);
        if ($errcode) {
            $errcode = "dnadist output $errcode";
            print STDERR "$subn: \$errcode=$errcode\n";
            chdir($pwd) or die "$!";
            return ($temp_new, $errcode);
        }
        # When there is not overlap between 2 sequences, dnadist writes -1.00 as the distance,
        # which causes fastME to abort. To overcome this, all -1.00 will be changed to 0.99.
        # Based on testing, this change doesn't cause any unexpected side-effect. - 10/15/2013

        $errcode = `grep -P "[-]1[.]00" $tempdir/outfile`;
        $debug && print STDERR "$subn: \$errcode=$errcode\n";
#        if ($errcode) {
#            $debug && print STDERR "$subn: 'cat $tempdir/infile'\n".`cat $tempdir/infile`."End of $tempdir/infile\n";
#            $debug && print STDERR "$subn: `cat $tempdir/outfile`".`cat $tempdir/outfile`."End of $tempdir/outfile\n";
#            $errcode = "dnadist output has '-1.00' as distance";
#            print STDERR "$subn: \$errcode=$errcode for $temp_plp\n";
#            chdir($pwd) or die "$!";
#            return ($temp_new, $errcode);
#        }
        if ($errcode) {
            $debug && print STDERR "$subn: 'cat $tempdir/infile'\n".`cat $tempdir/infile`."End of $tempdir/infile\n";
            $debug && print STDERR "$subn: `cat $tempdir/outfile`".`cat $tempdir/outfile`."End of $tempdir/outfile\n";
            $errcode = "dnadist output has '-1.00' as distance";
            print STDERR "$subn: \$errcode=$errcode for $temp_plp\n";

            # Change distance=-1.00 to 0.99
            $run_result = `mv $tempdir/outfile $tempdir/outfile0`;
            $run_result = `sed "s/-1\.00/ 0.99/g" "$tempdir/outfile0" > "$tempdir/outfile" `;

            print STDERR "$subn: Changed -1.00 to 0.99 in $temp_plp, re-check for '-1.00'\n";
            $errcode = `grep -P "[-]1[.]00" $tempdir/outfile`;
            $debug && print STDERR "$subn: After changing -1.00 to 0.99 \$errcode=$errcode\n";
            if ($errcode) {
#                $debug && print STDERR "$subn: 'cat $tempdir/infile'\n".`cat $tempdir/infile`."End of $tempdir/infile\n";
                $debug && print STDERR "$subn: `cat $tempdir/outfile`".`cat $tempdir/outfile`."End of $tempdir/outfile\n";
                $errcode = "dnadist output has '-1.00' as distance";
                $debug && print STDERR "$subn: \$errcode=$errcode\n";
                chdir($pwd) or die "$!";
                return ('', $errcode);
            } else {
                $debug && print STDERR "$subn: `cat $tempdir/outfile`".`cat $tempdir/outfile`."End of $tempdir/outfile\n";
            }
        }
#        $debug && print STDERR "$subn: ". `cat $tempdir/outfile`;
        $debug && print STDERR "$subn: ". `head -n5 $tempdir/outfile`;

    }
    $debug_all && print STDERR "$subn: `cat $tempdir/outfile`".`cat $tempdir/outfile`."End of $tempdir/outfile\n";
    $debug && print STDERR "$subn: distance file: $tempdir/outfile\n";

    my $dnadist_output = 'outfile';
    $debug && print STDERR "$subn: leaving subroutine\n";
    return ($dnadist_output, $errcode);
} # sub run_dnadist


=head1
sub calc_tree, takes an input alignment in phylip format, runs either PhyML or dnadist/FastME,
 returns the file name of the phylogeny tree or undef.
=cut

sub calc_tree {
    my ($faa, $tempdir, $temp_plp, $progs, $no_rerun, $dir_path, $debug) = @_;
    $debug = 0 || $debug_all || Genotype_util::getDebugAll();
    my $subn = "Genotype_calctree::calc_tree";

    $progs->{fastme}  = "fastme_linux32" if (!defined($progs) || !$progs->{fastme});

    $debug && print STDERR "$subn: \$faa=\n".Dumper($faa)."End of \$faa\n";
    my $temp_new = '';
    $debug && print STDERR "$subn: \$temp_plp='$temp_plp'\n";
    my $errcode = '';

    # generate phylogeny tree from phylip-style alignment
    my $cmd = '';
    my $use_phyml = 0;
  if ($use_phyml) {
    $temp_new = "$temp_plp"."_phyml_tree.txt";
    $debug && print STDERR "$subn: \$temp_new=$temp_new\n";
    # Using PhyML to generate phylogeny tree from phylip-style alignment, tree file is .phylip_phyml_tree.txt
    if (!-e "$dir_path/$temp_new" || !$no_rerun) {
      if ( 1 ) { # the actual calculation
        $cmd = "phyml -i $tempdir/$temp_plp -d nt";
        print STDERR "$subn: Command for tree is: '".$cmd."'\n";
        my $result_string = Genotype_util::getTime();
        system($cmd);
        $result_string .= "\tUPDATE\tGENOTYPE\t". "$cmd\n". Genotype_util::getTime() . "\tCOMPLETE\tGENOTYPE\tCalculated whole sequence genotype.\n";
        $debug && system("cp $tempdir/$temp_new $dir_path");
      } else { # Skip the actual calculation
        system("cp $dir_path/$temp_new $tempdir"); # copy existing file to tempdir
        print STDERR "$subn: Didn't calculate the tree, using existing result: \n";
        my $result_string = Genotype_util::getTime();
        $result_string .= "\tUPDATE\tGENOTYPE\t". "$cmd\n". Genotype_util::getTime() . "\tCOMPLETE\tGENOTYPE\tCalculated whole sequence genotype.\n";
      }
    } else {
        print STDERR "$subn: Found existing phylogeny tree file: $dir_path/$temp_new. PhyML not run.\n";
    }

  } else {
    # Using dnadist/fastme to generate phylogeny tree
    # Note: FastME can only accommodate 10-character long ids, so be careful
    $temp_new = "$temp_plp"."_fastme_tree.new";
    if (-e "$tempdir/$temp_new" && $no_rerun) {
        system("cp $dir_path/$temp_new $tempdir");
        print STDERR "$subn: Found existing phylogeny tree file: $dir_path/$temp_new. FastME not run.\n";
    } else {
        my $pwd = `pwd`;
        chomp($pwd);
        $debug && print STDERR "$subn: \$tempdir='$tempdir'\n";
        $debug && print STDERR "$subn: \$dir_path='$dir_path'\n";
        chdir($tempdir) or die "$!";
        $debug && print STDERR "$subn: \$pwd=$pwd\n";

        my $cp_result = `cp $tempdir/$temp_plp $tempdir/infile`;
        ($cp_result) && print STDERR "$subn: cp result=$cp_result\n";
        $debug && system("cp $tempdir/$temp_plp $dir_path"); # Save the alignment file for debugging
        $errcode = Genotype_util::check_file("$tempdir/infile", 12, $debug);
        if ($errcode) {
            $errcode = "dnadist input $errcode";
            print "$subn: \$errcode=$errcode\n";
            chdir($pwd) or die "$!";
            return ($temp_new, $errcode);
        }

        my $dnadist_output = '';
        ($dnadist_output, $errcode) = Genotype_calctree::run_dnadist( $tempdir, $temp_plp, $progs, $no_rerun, $dir_path, $debug);
        $debug && system("cp $tempdir/$dnadist_output $dir_path/${temp_plp}_outfile"); # Save distance matrix from dnsdist
        $debug && print STDERR "$subn: distance file saved to '$dir_path/${temp_plp}_outfile'\n";

        if ($errcode) {
            $debug && print STDERR "$subn: 'cat $tempdir/infile'\n".`cat $tempdir/infile`."End of $tempdir/infile\n";
            $debug && print STDERR "$subn: `cat $tempdir/outfile`".`cat $tempdir/outfile`."End of $tempdir/outfile\n";
            print STDERR "$subn: \$errcode=$errcode for $temp_plp\n";
            chdir($pwd) or die "$!";
            return ($temp_new, $errcode);
        }
#        $debug && print STDERR "$subn: ". `cat $tempdir/outfile`;
        $debug && print STDERR "$subn: ". `head -n5 $tempdir/outfile`;

        # Run FastME for phylogeny tree
        my $timeExe = time;
if ( 1 ) {
        $cmd = "cd $tempdir; $progs->{fastme} -i outfile";
} else {
        $cmd = "cd $tempdir; $progs->{fastme} --NNI=O -i outfile";
}
        $cmd = $cmd . ' 2>/dev/null' if (!$debug && -e '/dev/null');
        $debug && print STDERR "$subn: starting to run: \$cmd='$cmd'\n";
        $run_result = `$cmd`; # Run FastME
        $debug && print STDERR "$subn: execution time is: ".(time-$timeExe)." sec\n";
        $debug && print STDERR "$subn: Command for FastME is: '".$cmd."'\n";
        $debug && print STDERR "$subn: \$run_result=$run_result\n";
        $debug && print STDERR "$subn: 'ls -l'=\n". `ls -l`;
        my $fastme_result = '';
        if (-e 'output.t') {
            $fastme_result = 'output.t';
        } elsif (-e 'outfile_fastme_tree.txt') {
            $fastme_result = 'outfile_fastme_tree.txt';
        } else {
            $debug && print STDERR "$subn: \$cmd='$cmd'\n";
            $debug && print STDERR "$subn: pwd='". `pwd` ."'\n";
            &mydie("ERROR: $subn: Couldn't find either output.t or outfile_fastme_tree.txt");
        }
#        $debug && print STDERR "$subn: ". `cat $tempdir/$fastme_result`;
        $debug && print STDERR "$subn: ". `head -n5 $tempdir/$fastme_result`;
        $errcode = Genotype_util::check_file("$tempdir/$fastme_result");
        if ($errcode) {
            $errcode = "FastME output $errcode";
            $debug && print STDERR "$subn: \$errcode=$errcode\n";
            chdir($pwd) or die "$!";
            return ($temp_new, $errcode);
        }
        system("mv $fastme_result $temp_new");
        $debug && print STDERR "$subn: $temp_new=" . `cat $temp_new`;
        $debug && system("cp $tempdir/$temp_new $dir_path"); # Save the phylogeny tree for debugging
        chdir($pwd) or die "$!"; # restore the original PWD
    }
  }

    # check the phylogeny tree result
    $errcode = Genotype_util::check_file("$tempdir/$temp_new");
    if ($errcode) {
        $errcode = "newick tree file $errcode";
        $debug && print STDERR "$subn: \$errcode=$errcode\n";
        return ($temp_new, $errcode);
    }

    $debug && print STDERR "$subn: leaving subroutine\n";
    return ($temp_new, $errcode);
} # sub calc_tree


=head1
sub restoreDefline, looks through a newick file, restors the ids
 from accession only back to genotype|accession
 returns number of changes
=cut
sub restoreDefline {
    my ($deflines, $seqExcluded, $tmpNewick) = @_;
    my $debug = 0 || $debug_all || Genotype_util::getDebugAll();
    my $subn = "Genotype_calctree::restoreDefline";
    my $deflineFail = {};
    $debug && print STDERR "$subn: \$tmpNewick=$tmpNewick\n";
    $debug && print STDERR "$subn: \$deflines=\n".Dumper($deflines)."End of \$deflines\n";

    if (!-e $tmpNewick) {croak("$subn: '$tmpNewick', couldn't find: $OS_ERROR")};
    open my $file, '<', "$tmpNewick"
       or croak("$subn: '$tmpNewick', couldn't open: $OS_ERROR");
    open my $ofile, '>', "$tmpNewick.1"
       or croak("$subn: '$tmpNewick.1', couldn't open: $OS_ERROR");
    my $line;

    while ($line = <$file>) {
      chomp($line);
      my $i = 0;
      if ($line =~ m/\S+/) {
        $debug && print STDERR "$subn: 0 \$line='$line'\n";
        my $changed = 0;
        my ($def0, $def1);
        foreach $def1 (sort keys %$deflines) {
            $i++;
            $def0 = $deflines->{$def1};
#            $debug && print STDERR "$subn: ID #$i: \$def1='$def1' \$def0='$def0'\n";
            if ($line =~ s/$def1/$def0/) {
                ($debug && $i%10==0) && print STDERR "$subn: 0 \$line='$line'\n";
                $changed = 1;
            } elsif (exists($seqExcluded->{$def1})) {
                #$deflineFail->{$def1} = 1;
                $debug && print STDERR "$subn: Line #$i: found in excluded: $def1 => $def0 in '$tmpNewick'\n";
                $changed = -1;
            } else {
                $deflineFail->{$def1} = 1;
#                print STDERR "$subn: Line #$i: didn't find Key: $def1 => $def0 in '$tmpNewick'\n";
            }
        }
        if ($debug) {
#            !$changed && print STDERR "$subn: UnChanged line='$line'\n";
            $changed && print STDERR "$subn:   Changed line='$line'\n";
        }
#        $debug && print STDERR "$subn: 1 \$line='$line'\n";
      }
      print $ofile "$line\n";
    }
    close $file
       or croak("$subn: '$tmpNewick', couldn't close: $OS_ERROR");
    close $ofile
       or croak("$subn: '$tmpNewick.1', couldn't close: $OS_ERROR");
    $debug && print STDERR `ls -l $tmpNewick $tmpNewick.1`;
    $debug && print STDERR `head -n5 $tmpNewick $tmpNewick.1` . "\n";
    system("mv $tmpNewick.1 $tmpNewick") if (!-z "$tmpNewick.1");
    ($debug) && print STDERR "$subn: \$deflineFail=\n".Dumper($deflineFail)."End of \$deflineFail\n";
#exit;
    return $deflineFail;
} # sub restoreDefline


=head1
sub mapDefline, looks through a fasta MSA file, gets the deflines,
 and replaces them with accession only,
 returns a hash of all deflines in new => old format.
=cut
sub mapDefline {
    my ($rmsafile) = @_;
    my $debug = 0 || $debug_all || Genotype_util::getDebugAll();
    my $subn = "Genotype_calctree::mapDefline";
    my $deflines = {};

    if (!-e $rmsafile) {croak("$subn: '$rmsafile', couldn't find: $OS_ERROR")};
    open my $file, '<', "$rmsafile"
       or croak("$subn: '$rmsafile', couldn't open: $OS_ERROR");
    my $i = 0;
    my $line;

    while ($line = <$file>) {
      chomp($line);
      if ($line =~ s/^>//) {
        $i++;
        #$debug && print STDERR "$subn: 0 \$line='$line'\n";
        my ($def0, $def1, $data);
        if ($line =~ m/^\s*(\S+)\s*/) {
            $def0 = $1;
            #$debug && print STDERR "$subn: Line #$i: \$def0='$def0'\n";
            $data = [ split(/[|]/, $def0) ];
            #$debug && print STDERR "$subn: Line #$i: \$data => '@$data'\n";
            $def1 = ($#{$data}>0) ? $data->[1]
                  : $data->[0];
            $def1 = substr($def1, 0, (length($def1)>10) ? 10 : length($def1));
            $def1 = '_'.$def1 if ($def1 !~ m/(^\d_|_$)/ && length($def1)<10);
            #$debug && print STDERR "$subn: Line #$i: $def1 => $def0\n";
            if (exists($deflines->{$def1})) {
                croak("$subn: Line #$i: Key exists: $def1 => $def0 in '$rmsafile'\n");
            }
            $deflines->{$def1} = $def0;
        }
        $def0 =~ s/[|]/\[\|\]/;
        #$debug && print STDERR "$subn: 1 \$def0='$def0 \$def1='$def1'\n";
        $line =~ s/$def0/$def1/;
        $line = '>' . $line;
        #$debug && print STDERR "$subn: 1 \$line='$line'\n";
      }
    }
    close $file
       or croak("$subn: '$rmsafile', couldn't close: $OS_ERROR");
    $debug && print STDERR "$subn: \$deflines=".scalar(keys %$deflines)."\n".Dumper($deflines)."End of \$deflines\n";
#exit;
    return $deflines;
} # sub mapDefline


=head1
sub cleanupDefline, looks through a fasta MSA file, gets the deflines,
 and replaces them with accession only,
 returns a hash of all deflines in new => old format.
=cut
sub cleanupDefline {
    my ($deflines, $rmsafile) = @_;
    my $debug = 0 || $debug_all || Genotype_util::getDebugAll();
    my $subn = "Genotype_calctree::cleanupDefline";

    if (!-e $rmsafile) {croak("$subn: '$rmsafile', couldn't find: $OS_ERROR")};
    open my $file, '<', "$rmsafile"
       or croak("$subn: '$rmsafile', couldn't open: $OS_ERROR");
    open my $ofile, '>', "$rmsafile.1"
       or croak("$subn: '$rmsafile.1', couldn't open: $OS_ERROR");
    my $i = 0;
    my $line;

    while ($line = <$file>) {
      chomp($line);
      if ($line =~ m/^>/) {
        $i++;
        $debug && print STDERR "$subn: 0 Line #$i: \$line='$line'\n";
        my ($def0, $def1, $data);
        foreach $def1 (sort keys %$deflines) {
            $def0 = $deflines->{$def1};
            $def0 =~ s/[|]/\[\|\]/g;
            if ($line =~ s/$def0/$def1/) {
                $debug && print STDERR "$subn:   Line #$i: \$def0='$def0 \$def1='$def1'\n";
                $debug && print STDERR "$subn: 1 Line #$i: \$line='$line'\n";
                last;
            }
        }
      }
      print $ofile "$line\n";
    }
    close $file
       or croak("$subn: '$rmsafile', couldn't close: $OS_ERROR");
    close $ofile
       or croak("$subn: '$rmsafile.1', couldn't close: $OS_ERROR");
    $debug && print STDERR `ls -l $rmsafile $rmsafile.1`;
    $debug && print STDERR `head -n5 $rmsafile $rmsafile.1` . "\n";
    system("mv $rmsafile.1 $rmsafile") if (!-z "$rmsafile.1");
    $debug && print STDERR "$subn: \$deflines=\n".Dumper($deflines)."End of \$deflines\n";
#exit;
    return $deflines;
} # sub cleanupDefline


1;
