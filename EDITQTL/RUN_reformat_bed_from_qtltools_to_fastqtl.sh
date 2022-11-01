#!/bin/bash
 
module load tabix bedtools
PSC=/sc/arion/work/charna02/scripts/misc/subset_supplycol_space.pl

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Parse command line arguments
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
	--input)
	    BD1="$2"
	    shift
	    shift
	    ;;
	--output)
	    BD2="$2"
	    shift
	    shift
	    ;;
	--scratch)
	    SCR="$2"
	    shift
	    shift
	    ;;
    esac
done

zcat ${BD1} > ${SCR}.bed
echo "pid" > ${SCR}.bed.keepme
awk '{print $4}' ${SCR}.bed | sort | uniq -c | awk '$1==1 {print $2}' | awk '$1!="pid"' >> ${SCR}.bed.keepme
perl ${PSC} ${SCR}.bed.keepme ${SCR}.bed 3 > ${SCR}.bed2
mv ${SCR}.bed2 ${SCR}.bed
cut -f 1-4,7- ${SCR}.bed | sed s/pid/TargetID/g | bgzip > ${BD2}
tabix -p bed ${BD2}    
