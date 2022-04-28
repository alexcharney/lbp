#!/bin/bash

## USAGE -------------------------------------------------------------------

set -e -o pipefail
if [ $# -lt 3 ]; then
    echo "About:   Annotate a VCF file"
    echo "Usage:   sh RUN_lbp_wgs_variant_annotation.sh --input <full path to input VCF> --output <output prefix with fullpath> --scratchdir <scratch directory>"
    exit 1
fi

## MODULES -------------------------------------------------------------------

module purge
module load vep/96
module load CPAN/5.28.1 perl5
module load mysql
module load bcftools java picard tabix

## ARGUMENTS -------------------------------------------------------------------

while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
	--input)
	    INPUT="$2"
	    shift
	    shift
	    ;;
	--scratchdir)
	    SCRATCH="$2"
	    shift
	    shift
	    ;;
    esac
done
if [[ ! -d ${SCRATCH} ]]
then mkdir ${SCRATCH}
fi
   
## SETUP -------------------------------------------------------------------

cachedir=/hpc/packages/minerva-centos7/vep/96/Cache
lofteedir=/hpc/packages/minerva-centos7/vep/96/loftee/
lofteeres=/sc/arion/projects/psychgen2/resources/loftee_resources
PERL5LIB=${PERL5LIB}:${lofteedir}
ref=/sc/arion/projects/H_PBG/REFERENCES/GRCh38/Ref.GRCh38.VariantCalling/GRCh38.primary_assembly.genome.fa
dbnsfp="/sc/arion/projects/psychgen2/resources/dbNSFP4.0a/dbNSFPv4.0a_custombuild.gz"
snpeff="/hpc/packages/minerva-common/snpeff/4.3/snpEff"
dbnsfp_fields="$(echo MPC_score,Ensembl_{geneid,transcriptid})"


## SITES ONLY -------------------------------------------------------------------

I=${INPUT}
O=${SCRATCH}/TMP_sites_only.vcf
if [[ ! -f ${O}.success ]]
then
    java -Xmx50g -jar ${PICARD} MakeSitesOnlyVcf INPUT=${I} OUTPUT=${O}
    STATUS=$?
    if [[ ${STATUS} -eq 0 ]]
    then
	touch ${O}.success
    else
	touch ${O}.fail
	exit 1
    fi
fi

I=${SCRATCH}/TMP_sites_only.vcf
O=${SCRATCH}/TMP_sites_only.vcf.gz
if [[ ! -f ${O}.success ]]
then
    cat ${I} | bgzip > ${O} && tabix ${O}
    STATUS=$?
    STATUS=$?
    if [[ ${STATUS} -eq 0 ]]
    then
	touch ${O}.success
    else
	touch ${O}.fail
	exit 1
    fi    
fi

I=${SCRATCH}/TMP_sites_only.vcf.gz
O=${SCRATCH}/TMP_sites_only.tsv
if [[ ! -f ${O}.success ]]
then
    echo 'CHROM POS REF ALT' | tr ' ' '\t' > ${O}
    bcftools query -f '%CHROM %POS %REF %ALT\n' ${I} | tr ' ' '\t' >> ${O}
    STATUS=$?
    if [[ ${STATUS} -eq 0 ]]
    then
	touch ${O}.success
    else
	touch ${O}.fail
	exit 1
    fi
fi

I=${SCRATCH}/TMP_sites_only.vcf.gz
O=${SCRATCH}/TMP_sites_only.noinfo.vcf.gz
if [[ ! -f ${O}.success ]]
then
    bcftools annotate -x INFO ${I} | bgzip > ${O} && tabix ${O}
    STATUS=$?
    if [[ ${STATUS} -eq 0 ]]
    then
	touch ${O}.success
    else
	touch ${O}.fail
	exit 1
    fi
fi

## LOFTEE -------------------------------------------------------------------

I=${SCRATCH}/TMP_sites_only.noinfo.vcf.gz
O=${SCRATCH}/TMP_sites_only.noinfo.loftee.vcf.gz
if [[ ! -f ${O}.success ]]
then
    vep -i ${I} --vcf -o ${O} --compress_output bgzip --cache --dir_cache ${cachedir} --offline --force_overwrite \
	--plugin LoF,loftee_path:${lofteedir},human_ancestor_fa:${lofteeres}/human_ancestor.fa.gz
    STATUS=$?
    if [[ ${STATUS} -eq 0 ]]
    then
	touch ${O}.success
    else
	touch ${O}.fail
	exit 1
    fi
fi
zgrep -v ^\# ${O} > ${O}.rinput

