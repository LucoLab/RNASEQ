#!/usr/bin/env bash

set -o errexit
set -o pipefail
set -o nounset

export LC_ALL=C

#
#################################################################
#
#date: March 3, 2017
#platform: Ubuntu 16.04
#author: Villemin Jean-Philippe
#team: Epigenetic Component of Alternative Splicing - IGH
#
# merge.splicing.fromRMATS.sh
# Usage : 
# merge.splicing.fromRMATS.sh
# Purpose :
# Run merging/fitering on RMATS output by specific config per event
# PASS EVENT NAME AS ARGUMENT !!
# Filter supposed reliable events.
#
#################################################################

#CODE_DIR=/home/jean-philippe.villemin/code/RNA-SEQ/src

EVENT=$1
SUBDIR=$2
READS=$3
PATH_TO_CONFIG=$4
BASE_PATH=$5
ID1=${SUBDIR}
ISCONTROL=$6
CODE_DIR=$7
PYTHON3=$8



echo ""
echo "=========>>START BASH SCRIPT: "

echo "EVENT : ${EVENT}"
echo "READS : ${READS}"
echo "PATH_TO_CONFIG : ${PATH_TO_CONFIG}"
echo "BASE_PATH : ${BASE_PATH}"
echo "SUBDIR : ${SUBDIR}"
echo "SAMPLE_NAME : ${ID1}"
echo "ISCONTROL : ${ISCONTROL}"
echo "CODE_DIR : ${CODE_DIR}"
echo "PYTHON3 : ${PYTHON3}"


ALT=""

if [ ${ISCONTROL} == "True" ] ; then
ALT="--control"
fi

echo "ALT : ${ALT}"


FINAL_OUT=${BASE_PATH}${SUBDIR}/${EVENT}/ 
echo "FINAL_OUT : ${FINAL_OUT}"


FILE_STAT_OUTPUT=${BASE_PATH}${SUBDIR}/${EVENT}/stats.txt

${PYTHON3} ${CODE_DIR}/src/mergeFinal.py -c ${PATH_TO_CONFIG} -r ${READS} -e ${EVENT} ${ALT}

echo "SANITIZE : LOOK FOR ${ID1} INSIDE DIR : ${BASE_PATH}${SUBDIR}/${EVENT}/"

${PYTHON3} ${CODE_DIR}/src/cleanfusion.py -s1 TMP -p ${FINAL_OUT} -o ${ID1} -ct ${ISCONTROL}


echo "REMOVE RAW BED"
rm ${BASE_PATH}${SUBDIR}/${EVENT}/*_TMP_*bed

if [ ${ISCONTROL} == "False" ]; then
	
	echo "SEPARE EVENTS BY DPSI"
	
	awk '{if ($4 > 0) print $0 }' ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.bed > ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.INC.${EVENT}.bed
	
	awk '{if ($4 < 0) print $0 }' ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.bed > ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.EXC.${EVENT}.bed
	
	echo "SORT EVENTS BY DPSI"
	
	
	sort -t$'\t' -nk4 ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.EXC.${EVENT}.bed > ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.sorted.EXC.${EVENT}.bed
	
	sort -t$'\t' -rgk4 ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.INC.${EVENT}.bed > ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.sorted.INC.${EVENT}.bed
	
	cat ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.sorted.INC.${EVENT}.bed ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.sorted.EXC.${EVENT}.bed >  ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.ALL.${EVENT}.bed
	sort -k1,1V -k2,2n ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.ALL.${EVENT}.bed > ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.sorted.ALL.${EVENT}.bed
	echo "REMOVE UNSORTED EVENTS"
	rm ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.ALL.${EVENT}.bed
	rm ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.INC.${EVENT}.bed
	rm ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.EXC.${EVENT}.bed

	echo "FIND EXCEL find ${BASE_PATH}${SUBDIR}/${EVENT} -name *clean.SPLICING.xlsx"
	EXCEL=$(find ${BASE_PATH}${SUBDIR}/${EVENT} -name "*clean.SPLICING.xlsx")
	echo ${EXCEL}
	
	echo "ADD CLEAN TAB TO EXCEL"
	readarray -t arr < <(find ${BASE_PATH}${SUBDIR}/${EVENT} -name "*sorted.INC*" -o -name "*sorted.EXC*")
	for FILE in "${arr[@]}";
	do
	echo ""
	echo "python3 ${CODE_DIR}/src/addToExcel.py --file ${FILE} --excel  ${EXCEL} "
	${PYTHON3} ${CODE_DIR}/src/addToExcel.py --file ${FILE} --excel  ${EXCEL} 
	
	done
	find ${BASE_PATH}${SUBDIR}/${EVENT} -name "*.bed" ! -name "*INC*.bed" ! -name "*EXC*.bed" ! -name "*ALL*.bed"
	echo  "CLEAN DIR"
	readarray -t filesToRemove < <(find ${BASE_PATH}${SUBDIR}/${EVENT} -name "*.bed"  ! -name "*INC*.bed" ! -name "*EXC*.bed" ! -name "*ALL*.bed")
	for FILE2REMOVE in "${filesToRemove[@]}";
	do
	echo ""
	echo  "REMOVE ${FILE2REMOVE}"
	rm ${FILE2REMOVE}
	done
	echo  "Finish MergeSplicingFromOnelist.sh"
	
fi

if [ ${ISCONTROL} == "True" ]; then
	
	echo "SEPARE EVENTS BY DPSI"
	
	awk '{print $0 }' ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.bed > ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.CONTROL.${EVENT}.bed
	
	
	echo "SORT EVENTS BY DPSI"
	
	sort -t$'\t' -rgk4 ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.CONTROL.${EVENT}.bed > ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.sorted.CONTROL.${EVENT}.bed
	
	
	echo "REMOVE UNSORTED EVENTS"
	
	rm ${BASE_PATH}${SUBDIR}/${EVENT}/${ID1}.CONTROL.${EVENT}.bed

	echo "FIND EXCEL find ${BASE_PATH}${SUBDIR}/${EVENT} -name *bad.SPLICING.xlsx"
	EXCEL=$(find ${BASE_PATH}${SUBDIR}/${EVENT} -name "*bad.SPLICING.xlsx")
	echo ${EXCEL}
	
	echo "ADD CLEAN TAB TO EXCEL"
	
	readarray -t arr < <(find ${BASE_PATH}${SUBDIR}/${EVENT} -name "*sorted.CONTROL*")
	for FILE in "${arr[@]}";
	do
	echo ""
	echo "python3 ${CODE_DIR}/src/addToExcel.py --file ${FILE} --excel  ${EXCEL} "
	${PYTHON3} ${CODE_DIR}/src/addToExcel.py --file ${FILE} --excel  ${EXCEL} 
	done
	
	echo  "CLEAN DIR"
	readarray -t filesToRemove < <(find ${BASE_PATH}${SUBDIR}/${EVENT} -name "*.bed"  ! -name "*INC*.bed" ! -name "*EXC*.bed" !  -name "*ALL*.bed" ! -name "*CONTROL*.bed")
	for FILE2REMOVE in "${filesToRemove[@]}";
	do
	echo ""
	echo  "REMOVE ${FILE2REMOVE}"
	rm ${FILE2REMOVE}
	done
	
fi
echo "=========>> END BASH SCRIPT: "


exit 0
