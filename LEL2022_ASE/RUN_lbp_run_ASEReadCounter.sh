#!/bin/bash
 
ml bcftools tabix gatk/4.2.0.0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Parse command line arguments
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
	--iid-rna)
	    RID="$2"
	    shift
	    shift
	    ;;
	--iid-wgs)
	    WID="$2"
	    shift
	    shift
	    ;;
	--sid)
	    SID="$2"
	    shift
	    shift
	    ;;
	--bam-rna)
	    BAM="$2"
	    shift
	    shift
	    ;;	
	--results-dir)
            RDR="$2"
            shift
            shift
            ;;
    esac
done
DIR=/sc/arion/projects/psychgen/lbp/data/dna/wgs_GenotypeRefinement
REF=/sc/arion/projects/H_PBG/REFERENCES/GRCh38/Ref.GRCh38.VariantCalling/GRCh38.primary_assembly.genome.fa
ODR=/sc/arion/projects/psychgen/lbp/data/dna/wgs_ASEReadCounterInput/${SID}
if [[ ! -d ${ODR} ]]
then mkdir ${ODR}
fi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SUBSET VCF FOR INDIVIDUAL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

for i in {1,2,3,4,5,6,7,8,9,10,1{1..9},20,21,22,X,Y}
do
    I=${DIR}/lbp_wgs.chr${i}.vqsr_pipeline/strictAnnoFiltPass.vcf.gz
    O=${ODR}/chr${i}.vcf.gz
    if [[ ! -f ${O}.success ]]
    then
	bcftools view --samples ${WID} ${I} | bgzip > ${O} && tabix ${O}
	STATUS=$?
	if [[ ${STATUS} -eq 0 ]]
	then
	    touch  ${O}.success
	else
	    touch  ${O}.fail
	    exit 1
	fi
    fi
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# RENAME VCF IDs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

echo "${WID} ${SID}" | tr ' ' '\t' > ${ODR}/rename
for i in {1,2,3,4,5,6,7,8,9,10,1{1..9},20,21,22,X,Y}
do
    I=${ODR}/chr${i}.vcf.gz
    O=${ODR}/chr${i}.rename.vcf.gz
    if [[ ! -f ${O}.success ]]
    then
	bcftools view ${I} | bcftools reheader --samples ${ODR}/rename | bgzip > ${O} && tabix ${O}
	STATUS=$?
	if [[ ${STATUS} -eq 0 ]]
	then
	    touch  ${O}.success
	else
	    touch  ${O}.fail
	    exit 1
	fi
    fi
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CONCAT VCF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

I=${ODR}/flist
O=/sc/arion/projects/psychgen/lbp/data/dna/wgs_ASEReadCounterInput/${SID}.vcf.gz
echo ${ODR}/chr{1,2,3,4,5,6,7,8,9,10,1{1..9},20,21,22,X,Y}.rename.vcf.gz | tr ' ' '\n' > ${I}
if [[ ! -f ${O}.success ]]
then
    bcftools concat --file-list ${I} | bgzip > ${O} && tabix ${O}
    STATUS=$?
    if [[ ${STATUS} -eq 0 ]]
    then
	touch ${O}.success
    else
	touch ${O}.fail
	exit 1
    fi
fi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LIST HET SITES BY INDIVIDUAL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

I=/sc/arion/projects/psychgen/lbp/data/dna/wgs_ASEReadCounterInput/${SID}.vcf.gz
O=/sc/arion/projects/psychgen/lbp/data/dna/wgs_ASEReadCounterInput/${SID}.vcf.gz.hetsites
if [[ ! -f ${O}.success ]]
then
    bcftools view -i 'GT="het"' ${I} | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE\t%GT\t%AD]\n' > ${O}
    STATUS=$?
    if [[ ${STATUS} -eq 0 ]]
    then
	touch ${O}.success
    else
	touch ${O}.fail
	exit 1
    fi
fi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# RUN ASEReadCounter
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

VCF=/sc/arion/projects/psychgen/lbp/data/dna/wgs_ASEReadCounterInput/${SID}.vcf.gz
O=${RDR}/${SID}.output.table
if [[ ! -f ${O}.success ]]
then
    gatk --java-options "-Xmx100g" ASEReadCounter -R ${REF} -I ${BAM} -V ${VCF} -O ${O}
    STATUS=$?
    if [[ ${STATUS} -eq 0 ]]
    then
	touch ${O}.success
    else
	touch ${O}.fail
	exit 1
    fi
fi
