package Genotype_util;

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
use Genotype_Def;

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

my $BLAST_DB = Genotype_Def::getBLAST_DB;

    my $taxon_ids = Genotype_Def::getTaxon_Ids;
#    my $taxon_ids = {
#                      idfile => 'Genotype_refseq_def.txt',
#                      loaded => 0,
#                    };

    my $RMSA = Genotype_Def::getRMSA;

#    my $newHCV_ICTV = 1; # Sep 30, 2013
#    my $fmtGenotype = '\d[-._\dA-Za-z]*';
#    $fmtGenotype = '([0-9]\w+|\d-[A-Za-z]+|NA|.+)' if (!$newHCV_ICTV);
    my $fmtGenotype = Genotype_Def::getfmtGenotype;

############## Subroutines ################

sub setDebugAll { my ($debug) = @_; $debug_all = $debug; } # sub setDebugAll

sub getDebugAll { return $debug_all; } # sub getDebugAll


=head2
  sub getPerlVersions finds out the version of the important BioPerl modules used
/perl -MBio::Root::Version -e 'print $Bio::Root::Version::VERSION,"\n"'/
=cut
sub getPerlVersions{
    my $debug = 0 || $debug_all;
    my $subname = "Genotype_util::getPerlVersions";
    my $msgs = [];
    my $cmds = [ ['perl', ''],
                 ['Bio::Root::Version', '$Bio::Root::Version::VERSION'],
                 ['Bio::SeqIO',   '$Bio::SeqIO::VERSION'],
                 ['Bio::AlignIO', '$Bio::AlignIO::VERSION'],
                 ['Bio::TreeIO',  '$Bio::TreeIO::VERSION'],
#                 ['Bio::Tools::Run',  'Bio::Tools::Run::VERSION'],
#                 ['Bio::Tools::Run::StandAloneBlast',  'Bio::Tools::Run::StandAloneBlast::VERSION'],
               ];
    for my $i (0 .. $#{$cmds}) {
        my $cmd;
        if ($cmds->[$i]->[0] eq 'perl') {
            $cmd = 'perl -v | head -n3';
        } else {
            $msg .= "\n";
            $cmd = "perl -M$cmds->[$i]->[0] -e 'print $cmds->[$i]->[1]'"
        }
        $cmd = ` $cmd `;
        $cmd =~ s/(^\n|\n$)//ig; # Remove the leading return
        push @$msgs, sprintf("%-14s'", "$cmds->[$i]->[0]="). "$cmd'";
    }

    for (@$msgs) {
        $debug && print STDERR "$subname: $_\n";
    }
    return $msgs;
} # sub getPerlVersion


sub getTime {
    return strftime "%Y%m%d%H%M%S", localtime;
} # sub getTime


=head1
sub list_accessions takes a file name that contains the list of accessions, read the info into an array
=cut
sub list_accessions {
    my ($list_fn, $dir_path ) = @_;
    my $debug = 0 || $debug_all;
    my $subname = "Genotype_util::list_accessions";

    my $list_file;
    open $list_file, '<', "$dir_path/$list_fn"
           or croak("$subname: Found $list_fn, but couldn't open");
    while (<$list_file>) {
        my ($number, $acc);
        chomp;
       print $_;
        if ($_ =~ /^\|\s+(\w+)\s+\|\s+\d+\s+\|/) { # '| AM910652  | 484894 | Hepatitis C virus'
            $number = $#{$accs} +1;
            $acc = $1;
            $debug && print STDERR "\$number=$number \$1=\'$1\'\n";
        } elsif ($_ =~ /^\s*([^\s.]+)[.]\d+\s*$/) { # 'FJ888392.1'
            $number = $#{$accs} +1;
            $acc = $1;
            $debug && print STDERR "\$number=$number \$1=\'$1\'\n";
        } elsif ($_ =~ /^\s*(\d+):\s+(\w+)\s*$/) { # '806: AY858050'
            $number = $1;
            $acc = $2;
            $debug && print STDERR "\$number=$number \$1=\'$1\' \$2=\'$2\'\n";
        } else {
           $debug && print STDERR "$subname: skipping line: '$_'\n" if ($_ || $_ ne "/n");
           next;
        }
        push @$accs, [$number, $acc];

    }
    close $list_file or croak "$0: Couldn't close $dir_path/$list_fn: $OS_ERROR";
    $debug && print STDERR "\n$subname: finished reading list file: $dir_path/$list_fn, found $#{$accs} gbk files.\n\n";

    return $accs;
} # sub list_accessions

=head1
sub get_faas, takes file name (or an IO::String object), check if it's a genbank file or a fasta file.
For genbank file, write each sequence to a separate fasta file and return the file names in an array
For fasta file, simply return the file name in an array
=cut

sub get_faas {
    my ($infile, $job_name, $dir_path, $taxon, $INCL_genotype) = @_;
    my $debug = 0 || Genotype_util::getDebugAll();
    my $subn = "Genotype_util::get_faas";

    $INCL_genotype = '' if (!$INCL_genotype);
    my $faas;

    $debug && print STDERR "$subn: \$infile=\n".Dumper($infile)."End of \$infile\n";
    if (!defined($infile)) {
        push @$faas, {job_name=>$job_name, accession=>$job_name, file=>$infile, long_name=>$job_name, taxon=>$taxon, informat=>''};
        return $faas;
    }

    # For a genbank file, find the taxon if the taxon if not defined
    if (!$taxon && $infile =~ m/[.](gb|gbk|genbank)/i) {
        $debug && print STDERR "$subn: \$infile=$infile \$dir_path='$dir_path'\n";
        $taxon = Genotype_util::taxonByGbk( $infile, $dir_path);
        $taxon = '' if (!defined($taxon));
        print STDERR "$subn: taxon obtained from genbank file $infile: \$taxon='$taxon'\n";
    }

    my $in;
if ( 0 ) {
    $job_name = $infile;
    $job_name =~ s/[.]\w{1,7}$//; # remove the file name extension
}
    # Open input file in genbank format
    if (UNIVERSAL::isa( $infile, 'IO::String')) {
        $debug && print STDERR "$subn: trying IO::String \$infile=$infile\n";
        $in = Bio::SeqIO->new( -fh   => $infile, -format => 'genbank',);
    } else {
        $debug && print STDERR "$subn: trying genbank \$infile=$infile\n";
        $in = Bio::SeqIO->new( -file   => "<$dir_path/$infile", -format => 'genbank',);
    }
    my $ct = 0;
    while (my $inseq = $in->next_seq()) {
        $debug && print STDERR "$subn: \$inseq=\n".Dumper($inseq)."End of \$inseq\n";
        $ct++;
#        my $acc = $inseq->accession_number;
        my $acc = sprintf("_seq%03d", $ct);
        my $long_name = $inseq->accession_number;
        # get the country of the sequence
        foreach my $feat ($inseq->get_SeqFeatures) {
            next if ($feat->primary_tag ne 'source');
            my $country = '';
            my $year = '';
            if ($feat->has_tag('country')) {
              $country = ($feat->get_tag_values('country'))[0];
            }
            if ($feat->has_tag('collection_date')) {
              $year = ($feat->get_tag_values('collection_date'))[0];
            }
            $debug && print STDERR "$subn: acc=$long_name\t$acc\tcountry=$country\tyear=$year\n";
        }
        $debug && print STDERR "$subn: \$acc=$acc \$INCL_genotype='$INCL_genotype'\n";
        if ( $INCL_genotype ) {
            my $gbk_genotype = 'na';
            my $feats = [ $inseq->get_SeqFeatures ];
            for (my $ct = 0; $ct<=$#{$feats}; $ct++) {
                my $feat = $feats->[$ct];
                $debug && printf STDERR ("$subn: \$ct=$ct \$feat=%-8s", $feat->primary_tag);
                $debug && printf STDERR (" \$loc=%s", $feat->location->to_FTstring);
                if ($feat->primary_tag ne 'source') {
                    $debug && print STDERR " is not source. Skip.\n";
                    next;
                }
                $debug && printf STDERR ("$subn: \$ct=$ct \$feat=%-8s", $feat->primary_tag);
                $debug && printf STDERR (" \$loc=%s", $feat->location->to_FTstring);
                if ($feat->has_tag('subtype')) {
                  $debug && print STDERR " has tag 'subtype'.\n";
                  for my $note ($feat->get_tag_values('subtype')) {
                    $debug && print STDERR "$subn: \$acc=$acc \$note='$note'\n";
                    next if ($note !~ m/\s*($fmtGenotype(\/$fmtGenotype)*)\s*/i);
                    $gbk_genotype = $1;
                    print STDERR "$subn: 'subtype' from genbank: \$acc=$acc \$gbk_genotype='$gbk_genotype'\n";
                  }
                } elsif ($feat->has_tag('genotype')) {
                  $debug && print STDERR " has tag 'genotype'.\n";
                  for my $note ($feat->get_tag_values('genotype')) {
                    $debug && print STDERR "$subn: \$acc=$acc \$note='$note'\n";
                    next if ($note !~ m/\s*($fmtGenotype(\/$fmtGenotype)*)\s*/i);
                    $gbk_genotype = $1;
                    print STDERR "$subn: 'genotype' from genbank: \$acc=$acc \$gbk_genotype='$gbk_genotype'\n";
                  }
                }
                 if ((!$gbk_genotype || $gbk_genotype =~ m/^na$/i) && $feat->has_tag('organism')) {
                  $debug && print STDERR " has tag 'organism'.\n";
                  for my $note ($feat->get_tag_values('organism')) {
                    $debug && print STDERR "$subn: \$acc=$acc \$note='$note'\n";
                    next if ($note !~ m/\s*(subtype|genotype):* ($fmtGenotype(\/$fmtGenotype)*)\s*/i);
                    $gbk_genotype = $2;
                    print STDERR "$subn: '$1' from genbank: \$acc=$acc \$gbk_genotype='$gbk_genotype'\n";
                  }
                }
                if ((!$gbk_genotype || $gbk_genotype =~ m/^na$/i) && $feat->has_tag('note')) {
                  $debug && print STDERR " has tag 'note'.\n";
                  for my $note ($feat->get_tag_values('note')) {
                    $debug && print STDERR "$subn: \$acc=$acc \$note='$note'\n";
                    if ($note =~ m/inter-genotypic [(]($fmtGenotype(\/$fmtGenotype)*)[)] recomb/i) {
                        $gbk_genotype = $1;
                        print STDERR "$subn:1 'note' from genbank: \$acc=$acc \$gbk_genotype='$gbk_genotype'\n";
                    } elsif ($note =~ m/\s*genotype:* ($fmtGenotype) isolate .* genotype:* ($fmtGenotype) isolate\s*/i) {
                        $gbk_genotype = "$1/$2";
                        print STDERR "$subn:2 'note' from genbank: \$acc=$acc \$gbk_genotype='$gbk_genotype'\n";
                    } elsif ($note =~ m/\s*(subtype|genotype):* ($fmtGenotype(\/$fmtGenotype)*)\s*/i) {
                        $gbk_genotype = $2;
                        print STDERR "$subn:3 '$1' from genbank: \$acc=$acc \$gbk_genotype='$gbk_genotype'\n";
                    } elsif ($note =~ m/\s*recombinant .*; ($fmtGenotype(\/$fmtGenotype)*)\s*/i) {
                        $gbk_genotype = $1;
                        print STDERR "$subn:4 'note' from genbank: \$acc=$acc \$gbk_genotype='$gbk_genotype'\n";
                    } elsif ($note =~ m/\s*type: ($fmtGenotype(\/$fmtGenotype)*)\s*/i) {
                        $gbk_genotype = $1;
                        print STDERR "$subn:5 'note' from genbank: \$acc=$acc \$gbk_genotype='$gbk_genotype'\n";
                    }
                  }
                }
                if ((!$gbk_genotype || $gbk_genotype =~ m/^na$/i) && $feat->has_tag('strain')) {
                  $debug && print STDERR " has tag 'strain'.\n";
                  for my $note ($feat->get_tag_values('strain')) {
                    $debug && print STDERR "$subn: \$acc=$acc \$note='$note'\n";
                    if ($note =~ m/\s*(HCV|M21)-*($fmtGenotype(\/$fmtGenotype)*)\s*/i) {
                        $gbk_genotype = $2;
                    } elsif ($note =~ m/\s*(RF\d+)\s*($fmtGenotype(\/$fmtGenotype)*)\s*/i) {
                        $gbk_genotype = $2;
#                    } elsif ($note =~ m/\s*(RF\d+)(\d+(\/\d+)*)\s*/i) {
                    } elsif ($note =~ m/\s*(type)\s*($fmtGenotype(\/$fmtGenotype)*)/i) { # D85516
                        $gbk_genotype = $2;
                    } else {
                        next;
                    }
                    print STDERR "$subn: 'strain' from genbank: \$acc=$acc \$gbk_genotype='$gbk_genotype'\n";
                  }
                }
                if ((!$gbk_genotype || $gbk_genotype =~ m/^na$/i) && $feat->has_tag('serotype')) {
                  $debug && print STDERR " has tag 'serotype'.\n";
                  for my $note ($feat->get_tag_values('serotype')) {
                    $debug && print STDERR "$subn: \$acc=$acc \$note='$note'\n";
                    next if ($note !~ m/\s*($fmtGenotype(\/$fmtGenotype)*)\s*/i);
                    $gbk_genotype = $1;
                    print STDERR "$subn: 'serotype' from genbank: \$acc=$acc \$gbk_genotype='$gbk_genotype'\n";
                  }
                }
                if ((!$gbk_genotype || $gbk_genotype =~ m/^na$/i) && $feat->has_tag('isolate')) {
                  $debug && print STDERR " has tag 'isolate'.\n";
                  for my $note ($feat->get_tag_values('isolate')) {
                    $debug && print STDERR "$subn: \$acc=$acc \$note='$note'\n";
                    next if ($note !~ m/\s*genotype:*\s*($fmtGenotype(\/$fmtGenotype)*)\s*/i);
                    $gbk_genotype = $1;
                    print STDERR "$subn: 'serotype' from genbank: \$acc=$acc \$gbk_genotype='$gbk_genotype'\n";
                  }
                }
                if (!$gbk_genotype || $gbk_genotype =~ m/^na$/i) {
                    print STDERR " has no tag 'genotype' in tags: subtype, genotype, note, organism, or strain. Skip.\n";
                    next;
                }
                last;
            }
            $long_name .= ":$gbk_genotype" if ($gbk_genotype);
        }
#        my $temp_faa = "$acc.faa";
        my $temp_faa = "${job_name}${acc}.faa";
        $debug && print STDERR "$subn: \$acc=$acc \$temp_faa=$temp_faa\n";
        my $out  = Bio::SeqIO->new( -file   => ">$dir_path/$temp_faa", -format => 'fasta',);
        my $outseq = Bio::PrimarySeq->new( -seq => $inseq->seq, -id => $acc,);
        $out->write_seq($outseq);
        $errcode = Genotype_util::check_file("$dir_path/$temp_faa");
        if ($errcode) {
            $errcode = "fasta input $errcode";
            print STDERR "$subn: \$errcode=$errcode\n";
        } else {
            push @$faas, {job_name=>$job_name, accession=>$acc, file=>$temp_faa, long_name=>$long_name, taxon=>$taxon, informat=>'genbank'};
        }
    }

    # Open input file in FASTA format
    if ($#{$faas}<0) {
    # See if the input file is a fasta file
        $debug && print STDERR "$subn: trying fasta \$infile=$infile\n";
        $in  = Bio::SeqIO->new( -file   => "<$dir_path/$infile", -format => 'fasta',);

        return $faas if (!$in);

        $job_name = $infile;
        $job_name =~ s/[.]\w{1,7}$//; # remove the file name extension
        $debug && print STDERR "$subn: \$job_name=$job_name \$infile=$infile\n";
        my $ct = 0;
        my $inseq = $in->next_seq();
        while ($inseq) {
          $ct++;
          my $acc = sprintf("_seq%03d", $ct);
          my $inseq2 = $in->next_seq();
          my $long_id = $inseq->primary_id;
          if ($long_id =~ m/[():;,\[\]]/i) {
            # Check the input file to make sure it suits dnadist, which doesn't allow "():;,[]"
            $errcode = "fasta defline contains '():;,[]', disallowed by dnadist";
            croak("$subn: \$errcode=$errcode, abort\n");
          }
          $debug && print STDERR "$subn: \$ct=$ct \$acc=$acc \$long_id=$long_id\n";
          my $temp_faa = "${job_name}${acc}.faa";
          $debug && print STDERR "$subn: \$ct=$ct \$temp_faa=$temp_faa\n";
          my $out  = Bio::SeqIO->new( -file   => ">$dir_path/$temp_faa", -format => 'fasta',);
          my $outseq = Bio::PrimarySeq->new( -seq => $inseq->seq, -id => $acc); # id = '_seq001'
#          my $outseq = Bio::PrimarySeq->new( -seq => $inseq->seq, -id => $long_id); # id = input defline
          $out->write_seq($outseq);
          $errcode = Genotype_util::check_file("$dir_path/$temp_faa", 3, 0);
          if ($errcode) {
            $errcode = "fasta input $errcode";
            print STDERR "$subn: \$errcode=$errcode\n";
          } else {
            push @$faas, {job_name=>$job_name, accession=>$acc, file=>$temp_faa, long_name=>$long_id, taxon=>$taxon, informat=>'fasta'};
          }
          $inseq = $inseq2;
      }

    }

    $debug && print STDERR "$subn: \$faas=\n".Dumper($faas)."End of \$faas\n";
    $debug && print STDERR "$subn: leaving subroutine\n";
    return $faas;
} # sub get_faas


=head1
sub writeSliceFasta, takes an MSA object, and the start/end values of a slice, then take such slice,
 and writes the slice to a FASTA file. Returns any errcode.
=cut

sub writeSliceFasta {
    my ($aln, $job_name, $acc, $seq_name, $start, $end, $tempdir) = @_;
    my $debug = 0 || $debug_all || Genotype_util::getDebugAll();
    my $subn = "Genotype_util::writeSliceFile";

    if (!$aln || !$aln->isa('Bio::SimpleAlign')) {
        print STDERR "$subn: Input \$aln is not of Bio::SimpleAlign, abort\n";
        exit(1);
    }
#   push @$faas, {job_name=>$job_name, accession=>$acc, file=>$temp_faa, long_name=>$long_id};
    $job_name = 'test' if (!$job_name);

    my $tally = { total => 0, empty => 0, longgap => 0};
    my $seqExcluded = {};
    my $recombs = [];
#    $debug && print STDERR "$subn: \$aln = \n".Dumper($aln)."End of \$aln\n\n";

    my $errcode = '';
    $debug && print STDERR "$subn: \$job_name=$job_name, \$acc=$acc, \$seq_name=$seq_name, \$start=$start \$end=$end\n";
    my $seg_afa = sprintf("%s_%d_%d.afa", "${job_name}_$acc", $start+1, $end+1);

    # Take a slice of the MSA and write to a new fasta file
    my $seq_out = Bio::SeqIO->new(
                          '-file'     => ">$tempdir/$seg_afa",
                          '-format' => 'fasta'
                          );
    my $sequences = {};
    my $target_empty = 0;
    my $RMSA_ct = 0; # check if there is any refseq in MSA
    foreach my $seq ($aln->each_seq) {
        $tally->{total}++;
        my $sequence = $seq->seq;
        $sequence = substr($sequence, $start, $end-$start);
        ($debug && $RMSA_ct<5) && print STDERR "$subn: \$start=$start \$seq->{display_id}=$seq->{display_id}\n";
#        $debug && print STDERR "$subn: \$start=$start id=$seq->{display_id} \$sequence=$sequence\n";
        my $GAP_PCT = 80;
        my $pct = ($sequence =~ tr/-Nn//) / length($sequence) *100;
        if ( $pct > $GAP_PCT ) {
            $tally->{longgap}++;
            $seqExcluded->{$seq->{display_id}} = 1;
            $debug && printf STDERR "$subn: %6.2f%% of id=%-10s is gap \@$start..$end, name=$job_name:$seq_name skip\n", $pct, $seq->{display_id};
            if ($seq->{display_id} eq $seq_name ) {
                $errcode = sprintf("%2d%% of $seq_name in $job_name is gap \@$start..$end", $pct);
                if ($pct>99.9) {
                    $errcode = "$seq_name in $job_name is empty \@$start..$end";
                    $debug && print STDERR "$subn: $errcode, abort\n";
                    $tally->{empty}++;
                    $target_empty = 1;
                    last;
                }
            }
            next;
        } elsif ($seq->{display_id} ne $seq_name) {
              $RMSA_ct++;
        }
        if ($seq->{display_id} ne $seq_name ) {
            $sequences->{$seq->display_id} += 1;
        }
        my $id = ($seq->{display_id} ne $seq_name ) ? ".$sequences->{$seq->{display_id}}" : '';
        $id = $seq->{display_id}. $id;
        ($debug && $RMSA_ct<5) && print STDERR "$subn: \$RMSA_ct=$RMSA_ct \$start=$start \$id=$id\n";
        my $seq1 = Bio::LocatableSeq->new( -seq => "$sequence",
                                           -id  => $id,
                                           -start => 0,
                                           -end => length($sequence),
                                          );
#        $debug && print STDERR "$subn: \$seq1 = \n".Dumper($seq1)."End of \$seq1\n\n";
        $seq_out->write_seq($seq1);
    }
    $seq_out = undef;

    $errcode = ($RMSA_ct==0) ? "ref MSA is empty \@$start"
             : ($RMSA_ct==1) ? "ref MSA has 1 seq \@$start"
             : ((keys %$sequences)<1) ? "MSA has ".(keys %$sequences)." seqs \@$start"
             : $errcode;
    $debug && print STDERR "$subn: errcode='$errcode'\n";
    $debug && print STDERR "$subn: MSA between $start..$end has $RMSA_ct good sequences\n";
    $debug && print STDERR "$subn: \$seqExcluded=\n".Dumper($seqExcluded)."End of \$seqExcluded\n\n";
    $debug && print STDERR "$subn: \$sequences=\n".Dumper($sequences)."End of \$sequences\n\n";

    $debug && print STDERR "$subn: cat $tempdir/$seg_afa=\n";
    $debug && print STDERR `cat $tempdir/$seg_afa`;
    $debug && print STDERR "$subn: name=$job_name:$start..$end total:$tally->{total} longgap:$tally->{longgap} empty:$tally->{empty}\n";

    return ($seg_afa, $seqExcluded, $errcode);
} # sub writeSliceFasta


=head1
sub taxonByGbk, takes a DNA sequence in fasta file, runs blast against the refseqs of flaviviridae,
 and deteremines the refseqs with highest similarity.
 Returns the species, if the similarity is high enough.
 Note: User specified taxon is preferred, as blast-based auto detection might miss, e.g. FV537325(Japenceph)
=cut

sub taxonByGbk {
    my ($infile, $dir_path) = @_;
    my $subname = "taxonByGbk";

    my $debug = 0 || $debug_all;
    my $taxon = '';

    $debug && print STDERR "$subname: \$infile=$infile \$dir_path=$dir_path\n";
    if ($infile !~ m/[.](gb|gbk|genbank)/i) {
        return $taxon;
    }
    my $in;
    if (UNIVERSAL::isa( $infile, 'IO::String')) {
        $debug && print STDERR "$subname: trying IO::String \$infile=$infile\n";
        $in = Bio::SeqIO->new( -fh   => $infile, -format => 'genbank',);
    } else {
        $debug && print STDERR "$subname: trying genbank \$infile=$infile\n";
        $in = Bio::SeqIO->new( -file   => "<$dir_path/$infile", -format => 'genbank',);
    }
    my $inseq = $in->next_seq();
    my $taxid = $inseq->species->ncbi_taxid;
#    $debug && print STDERR "$subname: \$taxid=$taxid \$exe_dir=$exe_dir\n";
    my $ti = Annotate_Def::getTaxonInfo( $taxid, $exe_dir);
    if ($ti) {
        $taxon = Genotype_Def::getShortSpecies( $ti->[1]);
    }
    $debug && print STDERR "$subname: \$infile=$infile \$taxid=$taxid \$ti='@$ti' \$taxon=$taxon\n";

    $debug && print STDERR "$subname: leaving subroutine with \$taxon='$taxon'\n";
    return $taxon;
} # sub taxonByGbk

=head1
sub determine_taxon, takes a DNA sequence in fasta file, runs blast against the refseqs of flaviviridae,
 and deteremines the refseqs with highest similarity.
 Returns the species, if the similarity is high enough.
 Note: User specified taxon is preferred, as blast-based auto detection might miss, e.g. FV537325(Japenceph)
=cut

sub determine_taxon {
    my ($infile, $dir_path, $tempdir) = @_;
    my $subname = "determine_taxon";

    my $debug = 0 || $debug_all;
    my $taxon = '';
    my $hit_refseq = [];

    my $in  = Bio::SeqIO->new( -file   => "<$dir_path/$infile", -format => 'fasta',);
    my $seq_obj = $in->next_seq();
    my $seq_length = length($seq_obj->seq());
    $debug && print STDERR "$subname: \$infile=$infile \$seq_length=$seq_length\n";
#    $debug && print STDERR "$subname: \$seq_obj=\n".Dumper($seq_obj)."End of \$seq_obj\n";

    my @params = (
                   -program  => 'blastn',
#                   -database => $BLAST_DB->{'Flaviviridae'},
                   -database => $BLAST_DB->{'All'},
#                   -database => 'refseq_flaviviridae54.faa',
                 );
    my $pid = $seq_obj->primary_id ? $seq_obj->primary_id : 'vipr_genotype';
    push @params, (-outfile  => "$tempdir/" . $pid . '_blast1.out') if ($debug);
    $debug && print STDERR "$subname: \@params='@params'\n";
    my $blast_obj = Bio::Tools::Run::StandAloneBlast->new(@params);
    my $report_obj = $blast_obj->blastall($seq_obj);
    my $result_obj = $report_obj->next_result;
    $debug && print STDERR "$subname: \$infile=$infile \$result_obj->num_hits=".$result_obj->num_hits."\n";

    my $ct_hit = -1;
    while( $hit = $result_obj->next_hit ) {
      $ct_hit++;
      my $ct_hsp = -1;
      while(my $hsp = $hit->next_hsp ) {
         push @$hit_refseq, $hit->name if ($hsp->percent_identity>=75 && $hsp->length('total')>$seq_length*0.80);
         $ct_hsp++;
         my $hsp_lent = $hsp->length('total');
         my $hsp_percid = $hsp->percent_identity;
         $debug && ($ct_hit+$ct_hsp<5) && print STDERR "$subname: hit=$ct_hit \thsp=$ct_hsp\t", $hit->name, "\t",
             sprintf("%%Length:%5d/$seq_length=%6.2f\t%%id=%6.2f\n", $hsp_lent, $hsp_lent/$seq_length*100, $hsp_percid);
         if ($hsp->percent_identity>=50 && $hsp->length('total')>$seq_length*0.50) {
         }
      }
    }
#    $debug && print STDERR "$subname: \$result_obj=\n".Dumper($result_obj)."End of \$result_obj\n";
    $debug && print STDERR "$subname: \$hit_refseq='@$hit_refseq'\n";

    if (!$taxon_ids->{loaded}) {
            if (!-e $taxon_ids->{idfile}) {croak("$subname: '$taxon_ids->{idfile}', couldn't find: $OS_ERROR")};
            open my $taxon_file, '<', "$taxon_ids->{idfile}"
               or croak("$subname: '$taxon_ids->{idfile}', couldn't open: $OS_ERROR");
            while (<$taxon_file>) {
                next if (/^#/); # Skip comments
                s/(^[|]*\s+|\s+$)//x; # Remove first '|' or trailing space
                my $words = [ split(/\s*\|\s*/x) ];
                next if (!$words->[5]); # Skip if no shorthand for species
                next if ($words->[0] eq 'accession'); # Skip header of MySQL output
#                next if (!$words->[0] || $words->[0] !~ /^NC_/x); # Skip if not refseq
                next if (!$words->[0]); # Skip if not refseq
                chomp($words->[5]);
                $debug && print STDERR "$subname: \$words='".join(', ',@$words)."'\n";
                $taxon_ids->{$words->[0]} = $words->[5];
            }
            close $taxon_file or croak "$subname: Couldn't close $taxon_ids->{idfile}: $OS_ERROR";
            if (scalar(keys %$taxon_ids)>1) {
                $taxon_ids->{loaded} = 1;
            }
            $debug && print STDERR "$subname: \$taxon_ids=\n".Dumper($taxon_ids)."End of \$taxon_ids\n";
    }

    my $species = {};
    for my $i (@$hit_refseq) {
        $debug && print STDERR "$subname: \$i=$i\n";
        next if (!$taxon_ids->{$i});
        next if (exists($species->{$taxon_ids->{$i}}));
        $species->{$taxon_ids->{$i}} = 1;
    }
    $debug && print STDERR "$subname: \$species=\n".Dumper($species)."End of \$species\n";

#    $taxon = $taxon_ids->{$hit_refseq->[0]} if (scalar(keys %$species)==1 && exists($taxon_ids->{$hit_refseq->[0]}));
    for my $k (keys %$species) {
        $taxon .= ' ' if ($taxon);
        $taxon .= $k;
    }
    $debug && print STDERR "$subname: \$seq_obj=".$seq_obj->primary_id." ".scalar(keys %$species)." \$taxon='$taxon'\n";
    $debug && print STDERR "$subname: leaving subroutine\n";
    return $taxon;
} # determine_taxon


=head1
sub determine_region, takes a DNA sequence in fasta file, runs blast against the refseqs of flaviviridae,
 and deteremines the refseqs with highest similarity.
 Returns the species, if the similarity is high enough.
 Note: User specified taxon is preferred, as blast-based auto detection might miss, e.g. FV537325(Japenceph)
=cut

sub determine_region {
    my ($infile, $taxon, $ref_fasta0, $dir_path, $tempdir) = @_;
    my $subname = "determine_region";

    my $debug = 0 || $debug_all;
    my $hit_refseq = [];

    my $region = 'NA';
#    if ($taxon =~ /^NOV/i) {
#        $region = ($ref_fasta eq $RMSA->{NOVI})  ? 'ORF1'
#                : ($ref_fasta eq $RMSA2->{NOVI}) ? 'ORF2'
#                :                                  'NA' ;
#    }
    my $in  = Bio::SeqIO->new( -file   => "<$dir_path/$infile", -format => 'fasta',);
    my $seq_obj = $in->next_seq();
    my $seq_length = length($seq_obj->seq());
    $debug && print STDERR "$subname: \$infile=$infile \$seq_length=$seq_length\n";
#    $debug && print STDERR "$subname: \$seq_obj=\n".Dumper($seq_obj)."End of \$seq_obj\n";

    # change the filename for the refseqs
    my $ref_fasta = $ref_fasta0;
    $ref_fasta =~ s/(21|29|104|107)/org/i;
    $ref_fasta =~ s/[.]aln/.fasta/i;
    if (!-f $ref_fasta) {
        print STDERR "$subname: Refseq file doesn't exist: \$ref_fasta='$ref_fasta'\n";
        return $region;
    }
    if ($taxon !~ /^NOV/i) {
        print STDERR "$subname: input taxon:$taxon, region is set to: \$region=$region\n";
        return $region;
    }
    $debug && print STDERR "$subname: \$ref_fasta='$ref_fasta'\n";
    my @params = (
                   -program  => 'blastn',
#                   -database => $Genotype_Def::BLAST_DB->{'Flaviviridae'},
                   -database => $ref_fasta,
                 );
    my $pid = $seq_obj->primary_id ? $seq_obj->primary_id : 'vipr_genotype';
    my $outfileName = $pid . '_blast_region.out';
    push @params, (-outfile  => "$tempdir/" . $pid . '_blast_region.out');
    $debug && print STDERR "$subname: \@params='@params'\n";
    my $blast_obj = Bio::Tools::Run::StandAloneBlast->new(@params);
    my $report_obj = $blast_obj->blastall($seq_obj);
    $debug && print STDERR "$subname: blast=\n". `head -n400 $tempdir/$outfileName`;
    my $result_obj = $report_obj->next_result;
    $debug && print STDERR "$subname: \$infile=$infile \$result_obj->num_hits=".$result_obj->num_hits."\n";
#    $debug && print STDERR "$subname: \$result_obj=\n".Dumper($result_obj)."End of \$result_obj\n";

    my $ct_hit = -1;
    my $length = 0;
    my $conserved = 0;
    HIT: while( $hit = $result_obj->next_hit ) {
      $ct_hit++;
      my $ct_hsp = -1;
      while(my $hsp = $hit->next_hsp ) {
         $debug && print STDERR "$subname: \$hsp=\n".Dumper($hsp)."End of \$hsp\n";
         $debug && print STDERR "$subname: total=".$hsp->{'QUERY_LENGTH'}." query=".$hsp->{'HIT_LENGTH'}."\n";
         $debug && print STDERR "$subname: \$conversed=".$hsp->length('conversed')."\n";

         $length = ($hsp->{'QUERY_LENGTH'}<$hsp->{'HIT_LENGTH'}) ? $hsp->{'QUERY_LENGTH'} : $hsp->{'HIT_LENGTH'};
         $debug && print STDERR "$subname: \$length=$length\n";
         $conserved = $hsp->length('conversed');
         $ct_hsp++;
         my $hsp_lent = $hsp->length('total');
         my $hsp_percid = $hsp->percent_identity;
         $debug && ($ct_hit+$ct_hsp<5) && print STDERR "$subname: hit=$ct_hit\thsp=$ct_hsp\t", $hit->name, "\t",
             sprintf("%%Length:%5d/$length=%5.1f\t%%id=%5.1f\n", $hsp_lent, $hsp_lent/$length*100, $hsp_percid);
         # This is the criteria to call the region
         # norovirus JN183159 needs the criteria to be 70%
         if ($hsp->percent_identity>=75 && $hsp->length('total')>$length*0.67) {
             my $name = $hit->name;
             push @$hit_refseq, $name;
#             $region = ($ref_fasta0 eq $Genotype_Def::RMSA->[0]->{NOVI})  ? 'ORF1'
#                     : ($ref_fasta0 eq $Genotype_Def::RMSA2->[1]->{NOVI}) ? 'ORF2'
#                     : ($ref_fasta0 eq $Genotype_Def::RMSA3->[2]->{NOVI}) ? 'ORF1'
#                     : ($ref_fasta0 eq $Genotype_Def::RMSA4->[3]->{NOVI}) ? 'ORF2'
#    my $RMSA = Genotype_Def::getRMSA;
             $region = ($ref_fasta0 eq $RMSA->[0]->{NOVI})  ? 'ORF1'
                     : ($ref_fasta0 eq $RMSA->[1]->{NOVI}) ? 'ORF2'
                     : ($ref_fasta0 eq $RMSA->[2]->{NOVI}) ? 'ORF1'
                     : ($ref_fasta0 eq $RMSA->[3]->{NOVI}) ? 'ORF2'
                     :                                   'NA';
             $debug && print STDERR "$subname: \$region='$region' \$ref_fasta0='$ref_fasta0'\n";
             last HIT;
         }
      }
    }
#    $debug && print STDERR "$subname: \$result_obj=\n".Dumper($result_obj)."End of \$result_obj\n";
    $debug && print STDERR "$subname: \$hit_refseq='@$hit_refseq'\n";

    $debug && print STDERR "$subname: \$region=$region\n";
    $debug && print STDERR "$subname: leaving subroutine\n";

    return $region;
} # determine_region


=head1
sub check_file, takes a filename, possibly including the path, tests if the file exists or have 0 length.
  Could add more tests.
=cut

sub check_file {
    my ($file, $n, $debug) = @_;

    my $subname = "Genotype_util::check_file";
    $debug && print STDERR "$subname: \$file=$file\n";
    $debug = 0 || $debug_all;
    $n = 3 if (!$n);

    my $errcode = '';
    if (!-e "$file") {
        $errcode = "Missing file $file";
        $debug && print STDERR "$subname: \$errcode=$errcode\n";
    } elsif (-z "$file") {
        $errcode = "File $file has 0 length";
        $debug && print STDERR "$subname: \$errcode=$errcode\n";
    } else {
        $debug && print STDERR "$subname: " . `ls -l $file`;
        $debug && print STDERR "$subname: " . `head -n$n $file`;
    }

    #$debug && print STDERR "$subname: leaving subroutine\n";
    return $errcode;
} # sub check_file


################################################################################
=head1
sub check_program, takes a program, tests if the program exists or can be executed.
  Could add more tests.
=cut
sub check_program {
    my ($program, $debug) = @_;
    my $subname = "Genotype_util::check_program";

    $debug = (0 || $debug_all) if (!defined($debug));
    my $errcode = '';
    my $result = '';

    my $cmd = "which $program";
    $result = ` $cmd `;
    chomp($result);
    #$debug && print STDERR "$subname: '$cmd'='$result'\n";
    my $status;
    if (!$result || $result =~ /which: no /) {
        $status = "NOT accessible";
        $errcode = "Program $status $program";
    } else {
        $status = "is accessible";
    }
    $debug && print STDERR "$subname: Program $status: '$result'\n";

    return $errcode;
} # sub check_program


################################################################################
=head1
sub check_programAll, takes hash of all programs, tests if each program exists or can be executed.
=cut
sub check_programAll {
    my ($progs, $debug) = @_;
    my $subn = "Genotype_util::check_programAll";

    $debug = (0 || $debug_all) if (!defined($debug));
    my $errcode = '';
    foreach my $k (keys %$progs) {
      if ($k =~ /cladinator/) {
        print STDERR "$subn: skips Program=cladinator";
        next;
      }
      my $err = Genotype_util::check_program( $progs->{$k});
      $errcode .= $err . "\n" if ($err);
    }
    if ($errcode) {
      print STDERR "$subn: \$errcode=\n$errcode";
      exit(1);
    }

} # sub check_programAll


################################################################################
=head1
sub fasta2phylip, takes an input alignment in fasta format, created a new file in phylip format
 returns the file name of new phylip file.
=cut
sub fasta2phylip {
    my ($fasta, $phylip) = @_;
        # Convert fasta MSA to phylip format
        my $in  = Bio::AlignIO->new('-file' => "$fasta", '-format' => 'fasta',);
        my $in1 = $in->next_aln;
#        $debug && print STDERR "$subname: \$in1=\n".Dumper($in1)."End of \$in1\n\n";
        my $out = Bio::AlignIO->new('-file' => ">$phylip", '-format' => 'phylip',);
        $out->write_aln($in1);
} # sub fasta2phylip


################################################################################
=head1
sub fasta2FastME, takes an input alignment in fasta format, created a new file in FastME 2.15 format
 returns the file name of new phylip file.
=cut
sub fasta2FastME {
    my ($fasta) = @_;
    my $subn = "Genotype_util::fasta2FastME";
    $debug = 0 || $debug_all || Genotype_util::getDebugAll();

    # Convert fasta MSA to FastME 2.15 format: 
    # Count-of-sequence  length-of-MSA
    # <ID>  <sequence>
    my $fastme2file = $fasta;
    $fastme2file =~ s/[.][^.]+$//;
    $fastme2file .= '.seq';
    $debug && print STDERR "$subn: \$fastme2file=$fastme2file\n";

    my $in  = Bio::AlignIO->new('-file' => "$fasta", '-format' => 'fasta',);
    my $aln = $in->next_aln;
    my $count = $aln->num_sequences;
    my $len = $aln->length;
    $debug && print STDERR "$subn: \$count=$count \$len=$len\n\n";
    open my $ofile, '>', "$fastme2file"
       or croak("$subname: '$fastme2file', couldn't open: $OS_ERROR");
    print $ofile " $count $len\n";

    foreach my $seq ($aln->each_seq) {
        $debug && print STDERR "$subn: \$seq=".ref($seq). " \$len=$len, length ".$seq->id()."\t".$seq->length()."\n";
        printf $ofile "%-12s %s\n", $seq->id(), $seq->seq();
    }
#    $debug && print STDERR "$subn: \$aln=\n".Dumper($aln)."End of \$aln\n\n";
    close $ofile
       or croak("$subname: '$fastme2file', couldn't close: $OS_ERROR");

    return $fastme2file;
} # sub fasta2FastME


################################################################################
=head1
sub parseBI parses the message returned from newickBI
=cut
sub parseBI {
    my ($msg, $taxon) = @_;
    my $debug = 0 || $debug_all;
    my $subname = "Genotype_util::parseBI";

    chomp($msg);
    $debug && print STDERR "$subname: \$msg='$msg'\n";
    my $status;
    my $comment = '';

#    $status = "Success" if ($msg =~ /^[0-9.]+\t[0-9][a-zA-Za-z]$/);
#    if ($msg =~ /^BI=[0-9.]+\s+genotype=([0-9][a-zA-Za-z]*|\d-[a-zA-Za-z]+)$/
#         || ( $taxon =~ /^Nov/i && $msg =~ /^BI=[0-9.]+\s+genotype=(\w+[.]\w+)$/)) {
#    if ($msg =~ /^BI=[0-9.]+\s+genotype=([0-9][a-zA-Za-z]*|\d-[a-zA-Za-z]+)/
#         || ( $taxon =~ /^NOV/i && $msg =~ /^BI=[0-9.]+\s+genotype=(\w+[.]\w+)/)) {
=head1
    if ($msg =~ /^BI=[0-9.]+\s+genotype=$fmtGenotype/
         || ( $taxon =~ /^NOV/i && $msg =~ /^BI=[0-9.]+\s+genotype=(\w+[.]\w+)/)) {
        $status = "Success";
        $debug && print STDERR "$subname: newickBI returned: '$msg'\n";
#        $comment = "Success: '$msg'";
        $comment = '';
    } elsif ($msg =~ /^\s*$/ || $msg =~ /^BI=\s+genotype=\s*/) {
        $status = "Failed";
        $debug && print STDERR "$subname: newickBI returned empty string: '$msg'\n";
#        $comment = "newickBI returned empty string: '$msg'";
        $comment = '';
    } elsif ($msg =~ /^BI=[0-9.]+\s+genotype=NA\s+(.*)/) {
        $status = "Failed";
        $debug && print STDERR "$subname: newickBI returned NA: '$msg'\n";
#        $comment = "newickBI returned NA: '$msg'";
        $comment = $1;
    } else {
        $status = "Failed";
        $debug && print STDERR "$subname: newickBI returned unrecognized string: '$msg'\n";
        $comment = "newickBI returned unrecognized string: '$msg'";
    }
=cut
    my $subtype = '';
    my @items = split(/\t/, $msg);
    $msg = $items[0]."\t".$items[1]."\t".$items[2];
    $subtype = $items[2];
    $comment = $items[3] if ($items[3]);
    $debug && print STDERR "$subname: \@items='@items'\n";
    $debug && print STDERR "$subname: \$msg='$msg' \$subtype='$subtype' \$comment='$comment'\n";
    if ($msg =~ /^BI=[0-9.]+\s+genotype=NA\s+(.*)/) {
        $status = "Failed";
        $debug && print STDERR "$subname: newickBI returned NA: '$msg'\n";
#        $comment = "newickBI returned NA: '$msg'";
#        $comment = $1;
    } elsif ($msg =~ /^BI=[0-9.]+\s+genotype=$fmtGenotype/
         || ( $taxon =~ /^NOV/i && $msg =~ /^BI=[0-9.]+\s+genotype=(\w+[.]\w+)/)) {
        $status = "Success";
        $debug && print STDERR "$subname: newickBI returned: '$msg'\n";
#        $comment = "Success: '$msg'";
        $comment = '';
    } elsif ($msg =~ /^\s*$/ || $msg =~ /^BI=\s+genotype=\s*/) {
        $status = "Failed";
        $debug && print STDERR "$subname: newickBI returned empty string: '$msg'\n";
        $comment = "newickBI result has no BI value: '$msg'";
#        $comment = '';
    } else {
        $status = "Failed";
        $debug && print STDERR "$subname: newickBI returned unrecognized string: '$msg'\n";
        $comment = "newickBI returned unrecognized string: '$msg'";
    }

    $debug && print STDERR "$subname: leaving subroutine. \$status='$status' \$comment='$comment'\n";
    return ($status, $msg, $subtype, $comment);
} # sub parseBI

################################################################################
=head2 list_dir_files
List the files with given pattern in a folder
=cut

sub list_dir_files {
    my ($list_fn, $ptn) = @_;

    my $debug = 0 || $debug_all;
    my $subname = 'Genotype::list_dir_files';

    $debug && print STDERR "$subname: \$list_fn=$list_fn\n";
    $debug && print STDERR "$subname: \$ptn=$ptn\n";
    my @files = ();
    my $accs = [];
    if (!-d $list_fn) {
        croak("$subname: Couldn't locate accession file/directory: $list_fn: $OS_ERROR");
    }

    # if input -l file is directory
    opendir(DIR, $list_fn)
           or croak("$subname: Couldn't open dir $list_fn: $OS_ERROR");
    @files = sort readdir(DIR)
           or croak("$subname: Couldn't read dir $list_fn: $OS_ERROR");
    closedir(DIR)
           or croak("$subname: Couldn't close dir $list_fn: $OS_ERROR");
    $debug && print STDERR "$subname: \@files='@files'\n";

    for (my $f = 0; $f<=$#files; $f++) {
            my ($number, $acc);
            my $file = $files[$f];
            chomp $file;
            if ($file !~ /$ptn/) { # Keep the gbk files
#            if ($file !~ /^\s*([^\s]+)(\.(gb|gbk|genbank))\s*$/) { # Keep the gbk files
               $debug && print STDERR "$subname: skipping file: '$file'\n";
               next;
            } else {
                $number = $#{$accs}+1;
                $acc = "$file";
            }
            push @$accs, [$number, $acc];
    }
    $debug && print STDERR "$subname: \@\$accs=".Dumper($accs)."\n";

    $debug && print STDERR "$subname: leaving subroutine\n";
    return $accs;
} # sub list_dir_files


1;
