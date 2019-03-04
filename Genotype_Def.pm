package Genotype_Def;

# This file analyzes a single HCV sequence for genotype and recombination
# Requires access to:
# 1. MUSCLE, to generate multiple sequence alignment with refseqs of all genotypes,
# 2. dnsdist, to generate distance matrix from phylip-formatted alignment (DNA only), and
# 3. FastME to generate phylogeny tree from distance matrix
#use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';
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
use version; our $VERSION = qv('1.3.3'); # Apr 17, 2017
#use Getopt::Long;
use Annotate_Def;

@ISA = qw(Exporter);
@EXPORT      = qw($RMSA $fmtGenotype);

my $debug_all = 0;

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

    # May create link genotype to actual dir of particular version
    my $GENOTYPE_HOME = undef;
    $GENOTYPE_HOME = (-l __FILE__) ? dirname(readlink(__FILE__))
                   : dirname(__FILE__);
    $GENOTYPE_HOME = abs_path($GENOTYPE_HOME) . '/';
    $debug_all && print STDERR "Genotyp_Def: \$GENOTYPE_HOME='$GENOTYPE_HOME'\n";

# Program locations, may leave blank if the programs are accessible from the prompt
# Default for fastme is fastme_linux32
my $PROGS = { # default is for front-end analysis servers
    clustalw => '', # Location of MUSCLE, eg /net/home/gsun/prog/muscle/mus37/muscle
    muscle   => '', # Location of MUSCLE, eg /net/home/gsun/prog/muscle/mus37/muscle
    dnadist  => '', # Location of dnadist, eg /net/home/gsun/prog/phylip/phylip-3.69/exe/dnadist
    fastme   => $GENOTYPE_HOME. "bin/fastme", # Location of fastme, eg /net/home/gsun/prog/fastme/FastME_2.07/fastme_linux32
};
my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
$debug_all && print STDERR "Genotyp_Def: \$username='$username'\n";
$debug_all && print STDERR "Genotyp_Def: \$username='".`hostname`."'\n";

if ($username eq 'dbadmin' || `hostname` =~ m/^brc-?(dev|prod)-?d/i) {
# Following is for the convenience of Northrop backend processing, or AWS
    $PROGS->{muscle} = '/home/dbadmin/loader/ext/muscle3.7/muscle',
    $PROGS->{dnadist}= '/home/dbadmin/loader/ext/phylip-3.69/exe/dnadist',
#    $PROGS->{fastme} = '/home/dbadmin/loader/ext/FastME_2.07/fastme_linux64',
    $PROGS->{mafft}    = "/usr/bin/mafft-profile";
    $PROGS->{taxit}    = "/usr/bin/taxit";
    $PROGS->{pplacer}  = "/home/dbadmin/loader/ext/pplacer-Linux-v1.1.alpha19/pplacer";
    $PROGS->{guppy}    = "/home/dbadmin/loader/ext/pplacer-Linux-v1.1.alpha19/guppy";

} elsif ($username eq 'tomcat' || `hostname` =~ m/^brc-?(dev|prod)-?ana/i) {
    # Following is for the convenience of front-end analysis, or ASW
    $PROGS->{clustalw} = '/usr/bin/clustalw';
    $PROGS->{muscle}   = '/usr/bin/muscle';
    $PROGS->{dnadist}  = '/usr/bin/dnadist';
    $PROGS->{mafft}    = "/usr/bin/mafft-profile";
    $PROGS->{taxit}    = "/usr/bin/taxit";
    # To install pplacer/guppy: unzip pplacer-linux-v1.1.alpha19.zip
    $PROGS->{pplacer}  = "/opt/tools/pplacer-Linux-v1.1.alpha19/pplacer";
    $PROGS->{guppy}    = "/opt/tools/pplacer-Linux-v1.1.alpha19/guppy";

} else {
    $PROGS->{clustalw} = (-x '/usr/bin/clustalw') ? '/usr/bin/clustalw'
                       : 'clustalw';
    $PROGS->{muscle}   = (-x '/home/dbadmin/loader/ext/muscle3.7/muscle') ? '/home/dbadmin/loader/ext/muscle3.7/muscle'
                       : (-x '/usr/bin/muscle') ? '/usr/bin/muscle'
                       : 'muscle';
    $PROGS->{dnadist}  = (-x '/home/dbadmin/loader/ext/phylip-3.69/exe/dnadist') ? '/home/dbadmin/loader/ext/phylip-3.69/exe/dnadist'
                       : (-x '/usr/bin/dnadist') ? '/usr/bin/dnadist'
                       : 'dnadist';
    $PROGS->{mafft}    = "/usr/bin/mafft-profile";
    $PROGS->{taxit}    = "/usr/bin/taxit";
    $PROGS->{pplacer}  = "/home/gsun/prog/pplacer-Linux-v1.1.alpha19/pplacer";
    $PROGS->{guppy}    = "/home/gsun/prog/pplacer-Linux-v1.1.alpha19/guppy";
}
    # To instal MAFFT: sudo rpm -Uvh mafft-7.307-gcc_fc6.x86_64.rpm
    $PROGS->{mafft}    = (-x '/usr/bin/mafft-profile') ? '/usr/bin/mafft-profile'
                       : 'mafft-profile';

    my $BLAST_DB = {
            'Flaviviridae' => $GENOTYPE_HOME.'rmsa/refseq_Flaviviridae85.faa',
#            'All'          => $GENOTYPE_HOME.'rmsa/refseq_Taxon_All100.faa',
            'All'          => $GENOTYPE_HOME.'rmsa/refseq_Taxon_All139.faa',
           };

    my $RMSA = [{ # The MSA of all refseqs for each species, 'shorthand'=>filename
            #
#            'HCV'         => $GENOTYPE_HOME. 'rmsa/HCV_ICTV_ref220_partial_mafft_Yun.aln', # ICTV from Yun, 4/2017
#            'HCV'         => $GENOTYPE_HOME. 'rmsa/HCV_ICTV_ref204_mafft_Yun.aln', # ICTV from Yun, 4/2017
            'HCV'         => $GENOTYPE_HOME. 'rmsa/HCV_ICTV_ref193_muscle.aln', # ICTV new standard, 10/2013
            #'HCV'         => $GENOTYPE_HOME. 'rmsa/7-15_HCV_3ref_MAFFTaln.fas',
            'DENGUE'      => $GENOTYPE_HOME. 'rmsa/refseq_Dengue_34.aln',
            'STLOUIS'     => $GENOTYPE_HOME. 'rmsa/refseq_stlouisenvephalitis_33.aln',
            'WESTNILE'    => $GENOTYPE_HOME. 'rmsa/refseq_westnile_16.aln',
            'JAPENCEPH'   => $GENOTYPE_HOME. 'rmsa/refseq_Japanenceph_13.aln',
    # note: seqs of tick borne encephalitis and louping and omsk are very similar. To make TKBE stand out,
    # the refseqs of louping and omsk are excluded in the blast database.
            'TKBENCEPH'   => $GENOTYPE_HOME. 'rmsa/refseq_Tickborneenceph_19.aln',
            'YELLOWFEVER' => $GENOTYPE_HOME. 'rmsa/refseq_Yellowfever_19.aln',
            'BOVDIARRHEA1'=> $GENOTYPE_HOME. 'rmsa/refseq_Bovineviraldiarrhea_24.aln',
            'MURRAY'      => $GENOTYPE_HOME. 'rmsa/refseq_MurrayValley_12.aln',
            'NOROVIRUS'   => '',
            'NORWALK'     => '',
            # combination of all refseqs for ORF1, including regular + those for II.4
# left out for Jan2016 release            'NOV'         => $GENOTYPE_HOME. 'rmsa/refseq_Norovirus_ORF1_104.aln',
# left out for Jan2016 release            'NOVI'        => $GENOTYPE_HOME. 'rmsa/refseq_Norovirus_ORF1_104.aln',
# left out for Jan2016 release            'NOVII'       => $GENOTYPE_HOME. 'rmsa/refseq_Norovirus_ORF1_104.aln',
            'ZIKA'        => $GENOTYPE_HOME. 'rmsa/refseq_ZikaVirus_33.aln',
          },

          { # These are used to find genotype based on ORF2  for Norovirus
            # combination of all refseqs for ORF2, including regular + those for II.4
# left out for Jan2016 release            'NOV'         => $GENOTYPE_HOME. 'rmsa/refseq_Norovirus_ORF2_107.aln',
# left out for Jan2016 release            'NOVI'        => $GENOTYPE_HOME. 'rmsa/refseq_Norovirus_ORF2_107.aln',
# left out for Jan2016 release            'NOVII'       => $GENOTYPE_HOME. 'rmsa/refseq_Norovirus_ORF2_107.aln',
          },

          { # These are used to find the subtupe for II.P4 of Norovirus
#           'HCV'         => $GENOTYPE_HOME. 'rmsa/7-15_HCV_3ref_MAFFTaln.fas',
            # following file is combination of all 165(-161-) refseq from DSmith + combination of accessions for
            # "Remaining provisionally assigned HCV subtypes" such as 1d_L39299_L38377. Three obsolete
            # accessions are excluded: D14189, D14200, D14203, in favor of newer version of same sequences
#            'HCV'         => $GENOTYPE_HOME. 'rmsa/HCV_ICTV_ref189_muscle.fasta', # ICTV new standard, 10/2013
##            'HCV'         => $GENOTYPE_HOME. 'rmsa/HCV_ICTV_ref193_muscle.aln', # ICTV new standard, 10/2013
#            'NOV'         => $GENOTYPE_HOME. 'rmsa/refseq_Norovirus_ORF1_21.fasta',
#            'NOVI'        => $GENOTYPE_HOME. 'rmsa/refseq_Norovirus_ORF1_21.fasta',
#            'NOVII'       => $GENOTYPE_HOME. 'rmsa/refseq_Norovirus_ORF1_21.fasta',
          },

          { # These are used to find the subtupe for II.4 of Norovirus
#            'NOV'         => $GENOTYPE_HOME. 'rmsa/refseq_Norovirus_ORF2_29.fasta',
#            'NOVI'        => $GENOTYPE_HOME. 'rmsa/refseq_Norovirus_ORF2_29.fasta',
#            'NOVII'       => $GENOTYPE_HOME. 'rmsa/refseq_Norovirus_ORF2_29.fasta',
          }];

    my $tax2Short = { # The taxid => 'shorthand'
            11103 => 'HCV',
            12637 => 'DENGUE',
            11080 => 'STLOUIS',
            11082 => 'WESTNILE',
            11072 => 'JAPENCEPH',
            11084 => 'TKBENCEPH',
            11089 => 'YELLOWFEVER',
            11099 => 'BOVDIARRHEA1',
            11079 => 'MURRAY',
            11983 => 'NOV',
            64320 => 'ZIKA',
          };

    my $taxon_ids = {
              idfile => $GENOTYPE_HOME. 'Genotype_refseq_def.txt',
              loaded => 0,
            };

    my $newHCV_ICTV = 1; # Sep 30, 2013
    my $fmtGenotype = '\d[-._\dA-Za-z/]*';
    $fmtGenotype = '[-._\dA-Za-z/]*'; # Took out '\d' at beginning of pattern for genotype='II.P4' of NOV
    $fmtGenotype = '([0-9]\w+|\d-[A-Za-z]+|NA|.+)' if (!$newHCV_ICTV);

############## Subroutines ################


sub getTime {
    return strftime "%Y%m%d%H%M%S", localtime;
} # sub getTime


sub getTaxon_Ids {
    return $taxon_ids;
} # sub getTaxon_Ids

sub getBLAST_DB {
    return $BLAST_DB;
} # sub getBLAST_DB

sub getGenotypeHome {
    return $GENOTYPE_HOME;
} # sub getGenotypeHome

sub getProgs {
    return $PROGS;
} # sub getProgs

sub getRMSA {
    return $RMSA;
} # sub getRMSA

sub getfmtGenotype {
    return $fmtGenotype;
} # sub getfmtGenotype


sub getShortSpecies {
    my ($tax) = @_;
    my $shortSpecies = ($tax && (exists $tax2Short->{$tax})) ? $tax2Short->{$tax} : $tax;
    return $shortSpecies;
} # sub getShortSpecies


=head1
sub get_taxon_allowed, 
 Returns the list of allowed taxons
=cut

sub get_taxon_allowed {
    my ($taxon) = @_;
    my $subname = "Genotype_Def::get_taxon_allowed";

    my $debug = 0 || $debug_all;
    my $allowed_taxons = [ sort keys %{$RMSA->[0]} ];

    if ( 1 ) {
        $allowed_taxons = [" taxonID: shorthand\n"];
        for my $t (sort {$tax2Short->{$a} cmp $tax2Short->{$b}} keys %$tax2Short) {
            push @$allowed_taxons, sprintf("%8d: $tax2Short->{$t}\n", $t);
        }
    }
    $debug && print STDERR "$subname: \$allowed_taxons=@$allowed_taxons\n";

    return $allowed_taxons;
} # get_taxon_allowed


=head1
sub check_taxon_allowed, takes a taxonm checks if such taxon is covered by the script.
 Returns 1 if covered, undef otherwise
=cut

sub check_taxon_allowed {
    my ($taxon) = @_;
    my $subname = "Genotype_Def::check_taxon_allowed";

    my $allowed = undef;
    if (!$taxon) {
        return $allowed;
    }
    my $debug = 0 || $debug_all;
    my $allowed_taxons = [ sort keys %{$RMSA->[0]} ];
    $debug && print STDERR "$subname: \$taxon=$taxon\n";
    $debug && print STDERR "$subname: \$allowed_taxons=@$allowed_taxons\n";
    # Return 1 if $taxon is covered
    for my $k (keys %{$RMSA->[0]}) {
        next if ($taxon !~ /^$k$/i);
        $allowed = 1;
        return $allowed;
    }
    $debug && print STDERR "$subname: leaving subroutine\n";
    return $allowed;
} # check_taxon_allowed


## Create a new record of genotype
sub newGenotype {
    my ($job_name, $acc, $long_name) = @_;
    my $debug = 0;
    my $subname = "newGenotype";

    $job_name = ($job_name) ? $job_name : '';
    $acc = ($acc) ? $acc : '';
    $long_name = ($long_name) ? $long_name : '';
    # a hash to hold the result
    my $genotype = {
                job_name => $job_name,
                accession => $acc,
                long_name => $long_name,
#                accession => '',
                type => 'genotype',
                region => 'NA',
                range1 => '000',
                range2 => '000',
                BI => "BI=0.000\tgenotype=NA\tsub=NA",
                subtype => 'sub=',
                taxon => '',
                status => 'Failed',
                MSA_phylip => '',
                rmsa => '',
#                comment => 'new genotype hash',
                comment => '',
                comment2 => '',
                Date => getTime(),
                cladinatorOutput => '',
                };

    $debug && print STDERR "$subname: \$genotype=\n".Dumper($genotype)."\n";
    return $genotype;
} # newGenotype


=head1
sub convertGT converts the genotype (1, 2 ...) returned from newickBI, to what the refseq was called by biologists
=cut

sub convertGT {
    my ($msg, $taxon) = @_;
    my $debug = 0 || $debug_all;
    my $subname = "Genotype_Def::convertGT";

    chomp($msg);
    $debug && print STDERR "$subname: \$msg='$msg' \$taxon=$taxon\n";

    my $taxon_labels = {
               ZIKA   => { # This species has to be all UPPERCASE
                           '1' => 'African',
                           '1a' => 'East_African',
                           '1b' => 'West_African',
                           '2' => 'Asian',
                            },
               DENGUE => {
                           '2-AsianAme' => '2_AsianAmerican',
                           '2-Cosmopol' => '2_Cosmopolitan',
                            },
               HCV => {
                           '1_HQ007' => '1_HQ537007',
                           '1_AJ228' => '1_AJ851228',
                           '1_KC195' => '1_KC248195',
                           '2_JF112' => '2_JF735112',
                           '2_JF116' => '2_JF735116',
                           '2_JF118' => '2_JF735118',
                           '2_JF117' => '2_JF735117',
                           '2_JF119' => '2_JF735119',
                           '2_JF110' => '2_JF735110',
                           '2_KC236' => '2_KC197236',
                           '2_KC237' => '2_KC197237',
                           '2_KC239' => '2_KC197239',
                           '3_JF124' => '3_JF735124',
                           '4_FJ854' => '4_FJ025854',
                           '4_JX964' => '4_JX227964',
                           '6_DQ891' => '6_DQ278891',
                           '6_JX558' => '6_JX183558',
                           '6_JX553' => '6_JX183553',
                           '6_JX554' => '6_JX183554',
                           '6_JX551' => '6_JX183551',
                           '6_JX552' => '6_JX183552',
                           '6_JX549' => '6_JX183549',
                           '6_JX550' => '6_JX183550',
                           '6_JX557' => '6_JX183557',
                            },
               JAPENCEPH => {
                           1 => 'GI',
                           2 => 'GII',
                           3 => 'GIII',
                           4 => 'GIV',
                            },
                TKBENCEPH => {
                           1 => '1 (European Subtype)',
                           2 => '2 (Far Eastern Subtype)',
                           3 => '3 (Siberian Subtype)',
                           4 => 'Baltic and Finish S-TBEV strains',
                             },
                YELLOWFEVER => {
                           1 => 'Angola Genotype',
                           2 => 'East/Central African Genotype',
                           3 => 'East Africa Genotype',
                           4 => 'West Africa Genotype II',
                           5 => 'West Africa Genotype I',
                             },
                MURRAY => {
                           1 => 'I',
                           2 => 'II',
                           3 => 'III',
                           4 => 'IV',
                             },
                       };

    $taxon = uc $taxon;
    return $msg if (!exists($taxon_labels->{$taxon}));
        $debug && print STDERR "$subname: before converting: '$msg'\n";
#    if ($msg =~ /^BI=[0-9.]+\s+(genotype=([0-9][a-zA-Za-z]*|\d-[a-zA-Za-z]+))/) {
    if ($msg =~ m/^BI=[0-9.]+\s+(genotype=($fmtGenotype))\s+sub=\w+$/ ) {
      if (exists($taxon_labels->{$taxon}->{$2})) {
        my @keys = keys(%{$taxon_labels->{$taxon}});
        $debug && print STDERR "$subname: \$1='$1' \$2='$2'\n";
        $debug && print STDERR "$subname: \@keys='@keys'\n";
        my $old = $1;
        foreach (sort @keys) {
          $debug && print STDERR "$subname: \$2=$2 \$_='$_'\n";
          next if ($2 ne $_);
          my $new = "genotype=" . $taxon_labels->{$taxon}->{$2};
          $debug && print STDERR "$subname: \$1='$1' \$2='$2'\n";
          $debug && print STDERR "$subname: \$old='$old' \$new='$new'\n";
          $debug && print STDERR "$subname: \$msg='$msg'\n";
          $msg =~ s/$old/$new/;
#          $_[0] = $msg;
          $debug && print STDERR "$subname: after conversion \$msg='$msg'\n";
          $debug && print STDERR "$subname: genotype changed to: '$new', from $old\n";
          last;
        }
      }
    } else {
        $debug && print STDERR "$subname: somehow the match didn't work: \$msg='$msg'\n";
    }

    $debug && print STDERR "$subname: leaving subroutine\n";
    return $msg;
} # sub convertGT


=head1
sub convertSubtype converts the genotype (1, 2 ...) returned from newickBI, to what the refseq was called by biologists
=cut

sub convertSubtype {
    my ($msg, $taxon) = @_;
    my $debug = 0 || $debug_all;
    my $subname = "Genotype_Def::convertSubtype";

    chomp($msg);
    $debug && print STDERR "$subname: \$msg='$msg' \$taxon=$taxon\n";

    my $subtype_table = {
               NOV => {
                           'Camb'   => 'Camberwell_1994',
                           'Brist'  => 'Bristol_1993',
                           '1995'   => 'US95_96',
                           '2003'   => 'Kaiso_2003',
                           '2002'   => 'Farmington_Hills_2002',
                           '2002C'  => 'Lanzou_2002', #2002CN
                           '2004'   => 'Hunter_2004',
                           '2006a'  => 'Yerseke_2006a',
                           '2006b'  => 'Den_Haag_2006b',
                           '2007J'  => 'Osaka_2007', #2007JP
                           '2007E'  => 'Apeldoorn_2007', #2007EU
                           '2009'   => 'New_Orleans_2009',
                           '2005'   => 'Asia_2003',
                           '2012'   => 'Sydney_2012',
                            },
                       };

    $taxon = uc $taxon;
    return $msg if (!exists($subtype_table->{$taxon}));
        $debug && print STDERR "$subname: before converting: '$msg'\n";
#    if ($msg =~ /^BI=[0-9.]+\s+(genotype=([0-9][a-zA-Za-z]*|\d-[a-zA-Za-z]+))/) {
    if ($msg =~ m/^BI=[0-9.]+\s+(genotype=$fmtGenotype)\s+sub=(\w+)$/ ) {
      if (exists($subtype_table->{$taxon}->{$2})) {
        my @keys = keys(%{$subtype_table->{$taxon}});
        $debug && print STDERR "$subname: \$1='$1' \$2='$2'\n";
        $debug && print STDERR "$subname: \@keys='@keys'\n";
        my $old = "sub=$2";
        foreach (sort @keys) {
          $debug && print STDERR "$subname: \$2=$2 \$_='$_'\n";
          next if ($2 ne $_);
          my $new = "sub=" . $subtype_table->{$taxon}->{$2};
          $debug && print STDERR "$subname: \$1='$1' \$2='$2'\n";
          $debug && print STDERR "$subname: \$old='$old' \$new='$new'\n";
          $debug && print STDERR "$subname: \$msg='$msg'\n";
          $msg =~ s/$old/$new/;
#          $_[0] = $msg;
          $debug && print STDERR "$subname: after conversion \$msg='$msg'\n";
          print STDERR "$subname: subtype changed to: '$new', from $old\n";
          last;
        }
      } else {
        $debug && print STDERR "$subname: Couldn't find matching subtype in datatable: \$2='$2'\n";
      }
    } else {
        $debug && print STDERR "$subname: \$1='$1' \$2='$2'\n";
        $debug && print STDERR "$subname: somehow the match didn't work: \$msg='$msg'\n";
    }

    $debug && print STDERR "$subname: leaving subroutine\n";
    return $msg;
} # sub convertSubtype


=head2 to_String
Turn a genotype object into a string
=cut

sub to_String {
    my ($genotype) = @_;

    my $debug = 0;
    my $subname = 'Genotype_Def::to_String';

#    $debug && print STDERR "$subname: \$genotype = \n".Dumper($genotype)."End of \$genotype\n\n";
    my $msg;
#    $msg .= "$genotype->{accession}";
    $msg .= "$genotype->{long_name}";
    $msg .= "\t$genotype->{type}";
    $msg .= "\t$genotype->{region}";
    $msg .= "\t$genotype->{range1}..$genotype->{range2}";
    $msg .= "\t$genotype->{BI}";
    $msg .= "\t$genotype->{taxon}";
#    $msg .= "\t$genotype->{Date}";
    $msg .= "\t$genotype->{status}";
    $msg .= "\t";
#    $msg .= "comment=";
    $msg .= "$genotype->{comment}";
    $msg .= "\t$genotype->{comment2}" if (exists($genotype->{comment2}) && $genotype->{comment2});
    $debug && print STDERR "$subname: \$msg='$msg'\n";

    return $msg;
} # sub to_String


=head2 parseGenotypeString
Turn a string into a genotype object, used to load the recomb data before drawing a graph
=cut

sub parseGenotypeString {
    my ($line, $ALGO) = @_;

    my $debug = 0;
    my $subname = 'Genotype_Def::parseGenotypeString';

    chomp $line;
    my $msg = [ split("\t", $line) ];
    $debug && print STDERR "$subname: \$msg=\n".Dumper($msg)."End of \$msg\n";
    if (!$msg || $#{$msg}<6 || $msg->[0] eq 'Job_name') {
        $debug && print STDERR "$subname: skipped line: '$line'\n";
        return undef;
    }

    my $genotype = {
                job_name => '',
                accession => '',
                long_name => '',
                type => '',
                range1 => '0',
                range2 => '0',
                BI => "BI=\tgenotype=NA",
                taxon => '',
                status => 'Failed',
                comment => '',
                comment2 => '',
                Date => '',
                algo => $ALGO,
                };
#    $msg .= "$genotype->{accession}";
    $genotype->{long_name} = $msg->[0];
    $genotype->{type} = $msg->[1];
    $genotype->{region} = $msg->[2];
    if ($msg->[3] =~ m/^\s*(\d+)\.\.\s*(\d+)$/i) {
        ($genotype->{range1}, $genotype->{range2}) = ($1, $2);
    }
    $genotype->{BI} = "$msg->[4]\t$msg->[5]\t$msg->[6]";
    $genotype->{taxon} = $msg->[7];
#    $genotype->{Date} = $msg->[7];
    $genotype->{status} = $msg->[8];
    $genotype->{comment} = $msg->[9];
=head2
    if ($msg->[2] =~ m/^\s*(\d+)\.\.\s*(\d+)$/i) {
        ($genotype->{range1}, $genotype->{range2}) = ($1, $2);
    }
    $genotype->{BI} = "$msg->[3]\t$msg->[4]";
    $genotype->{taxon} = $msg->[5];
    $genotype->{Date} = $msg->[6];
    $genotype->{status} = $msg->[7];
    $genotype->{comment} = $msg->[8];
=cut
    $debug && print STDERR "$subname: \$genotype=\n".Dumper($genotype)."End of \$genotype\n";

    return $genotype;
} # sub parseGenotypeString

1;
