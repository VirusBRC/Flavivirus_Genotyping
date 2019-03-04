#!/bin/sh

if [ -z "${1}" ]; # need a target
  then echo "No target is provided, please call this script as: ./pack_vipr_genotype.sh <target>";
  exit 1;

elif [ "${1}" != "${PWD##*/}" ]; # the target needs to be same as pwd
  then echo "The target has to be identical to current dir. Command: ./pack_vipr_genotype.sh <target>";
  exit 1;

else # step out of pwd, pack-up pwd, then return to pwd
  cd ..;

  # can't start with tar.gz, since that prevents more files to be added later.
  tar -cf ${1}.tar ${1}/vipr_genotype.pl  ${1}/vipr_genotype_parallel.pl  ${1}/vipr_genotype_recombgraph.pl # Production script
  tar -rf ${1}.tar ${1}/vipr_genotype_dev.pl  ${1}/vipr_genotype_test.pl  ${1}/vipr_genotype_testp.pl
  tar -rf ${1}.tar ${1}/Annotate_Def.pm  ${1}/Genotype_Def.pm  ${1}/Genotype_Draw.pm  ${1}/Genotype_util.pm
  tar -rf ${1}.tar ${1}/Genotype.pm  ${1}/Genotype_TPP.pm  ${1}/Genotype_calctree.pm  ${1}/Genotype_newickBI.pm  ${1}/Genotype_recomb.pm 
  tar -rf ${1}.tar ${1}/GBKUpdate/Configuration.pm  ${1}/GBKUpdate/Database.pm
  tar -rf ${1}.tar ${1}/Annotate_taxon_records.txt  ${1}/Genotype_refseq_def.txt  ${1}/out.txt  ${1}/err.txt
  tar -rf ${1}.tar ${1}/release_note.txt  ${1}/ViPR_Genotype_SOP.pdf
  tar -rf ${1}.tar ${1}/algo ${1}/bin ${1}/config ${1}/referenceData ${1}/util ${1}/XML
  tar -rf ${1}.tar ${1}/rmsa/HCV_ICTV_mapping.tsv
  tar -rf ${1}.tar ${1}/rmsa/refseq_Taxon_All139.faa.*
  tar -rf ${1}.tar ${1}/rmsa/refseq_Flaviviridae85.faa*  ${1}/rmsa/HCV_ICTV_*ref???_muscle.aln  ${1}/rmsa/refseq_Dengue_34.aln
  tar -rf ${1}.tar ${1}/rmsa/refseq_stlouisenvephalitis_33.aln  ${1}/rmsa/refseq_westnile_16.aln
  tar -rf ${1}.tar ${1}/rmsa/refseq_Japanenceph_13.aln  ${1}/rmsa/refseq_Tickborneenceph_19.aln
  tar -rf ${1}.tar ${1}/rmsa/refseq_Yellowfever_19.aln  ${1}/rmsa/refseq_Bovineviraldiarrhea_24.aln  ${1}/rmsa/refseq_MurrayValley_12.aln
  tar -rf ${1}.tar ${1}/rmsa/refseq_ZikaVirus_33.aln
  tar -rf ${1}.tar ${1}/rmsa/refseq_Norovirus_ORF?_*.aln  ${1}/rmsa/refseq_Norovirus_ORF?_org.fasta*
  tar -rf ${1}.tar ${1}/test/*.gb  ${1}/test/output/*_index.tsv  ${1}/test/output/usage.err
  # include FastMe V2.07. Necessary as phyml has a utility also named fastme
  tar -rf ${1}.tar ${1}/bin/fastme
  # include the packing shell script
  tar -rf ${1}.tar ${1}/pack_vipr_genotype.sh

  # zip the tar file
  if [ -e ${1}.tar.gz ]
  then
    echo "File=${1}.tar.gz exists, deleted."
    rm ${1}.tar.gz
  fi
  gzip ${1}.tar

  mv ${1}.tar.gz  ${1};
  cd ${1};
  ls -l ${1}.tar.gz;

  exit 0;
fi

