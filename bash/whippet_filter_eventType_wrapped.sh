#!/usr/bin/env bash

set -o errexit
set -o pipefail
set -o nounset
#set -o xtrace

# Set magic variables for current file & dir
__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
__file="${__dir}/$(basename "${BASH_SOURCE[0]}")"
__base="$(basename ${__file} .sh)"

#################################################################
#
#date: Fev 22, 2017
#platform: Ubuntu 16.04
#author: Villemin Jean-Philippe
#team: Epigenetic Component of Alternative Splicing - IGH
#
# whippet.sh
# Usage : 
# whippet.sh 
# Purpose :
# Filter supposed reliable events.
#
#################################################################


EVENT=$1
FILE=$2
PATH=$3


TYPE=""
echo ${PATH}${FILE}.diff

echo "Spliced "

#################################################################

if [ ${EVENT} == "CE" ] || [ ${EVENT} == "SE"  ] ; then

TYPE="CE"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="CE" && sqrt($8*$8) >= sqrt(0.1*0.1) && $9 >= 0.95) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.clean.${TYPE}.diff

fi

#################################################################

if [ ${EVENT} == "A3SS" ] || [ ${EVENT} == "AA"  ] ; then

TYPE="AA"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="AA" && sqrt($8*$8) >= sqrt(0.1*0.1) && $9 >= 0.95) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.clean.${TYPE}.diff

fi

#################################################################

if [ ${EVENT} == "A5SS" ] || [ ${EVENT} == "AD"  ] ; then

TYPE="AD"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="AD" && sqrt($8*$8) >= sqrt(0.1*0.1) && $9 >= 0.95) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.clean.${TYPE}.diff

fi

#################################################################

if [ ${EVENT} == "RI" ] || [ ${EVENT} == "IR"  ] ; then


TYPE="RI"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="RI" && sqrt($8*$8) >= sqrt(0.1*0.1) && $9 >= 0.95) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.clean.${TYPE}.diff

fi

#################################################################

if [ ${EVENT} == "TS"   ] ; then


TYPE="TS"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="TS" && sqrt($8*$8) >= sqrt(0.1*0.1) && $9 >= 0.95) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.clean.${TYPE}.diff

fi

#################################################################

if [ ${EVENT} == "TE"   ] ; then


TYPE="TE"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="TE" && sqrt($8*$8) >= sqrt(0.1*0.1) && $9 >= 0.95) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.clean.${TYPE}.diff

fi

#################################################################

if [ ${EVENT} == "AF"   ] ; then


TYPE="AF"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="AF" && sqrt($8*$8) >= sqrt(0.1*0.1) && $9 >= 0.95) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.clean.${TYPE}.diff

fi

#################################################################

if [ ${EVENT} == "AL"   ] ; then


TYPE="AL"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="AL" && sqrt($8*$8) >= sqrt(0.1*0.1) && $9 >= 0.95) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.clean.${TYPE}.diff

fi

#################################################################

/bin/sed -i $'1 i\\\ngene\tcoordinates\tstrand\tevent\tpsiA\tpsiB\tdpsi\tprobability\tcomplexity\tentropy' ${PATH}${FILE}.clean.${TYPE}.diff

# removed $9 < 0.95
echo "Control"

if [ ${EVENT} == "CE" ] || [ ${EVENT} == "SE"  ] ; then

TYPE="CE"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="CE" && $6 >= 0.90 && $7 >= 0.90 && sqrt($8*$8) <=  0.01 ) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.bad.${TYPE}.diff

fi

#################################################################

if [ ${EVENT} == "A3SS" ] || [ ${EVENT} == "AA"  ] ; then

TYPE="AA"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="AA" && $6 >= 0.90 && $7 >= 0.90  && sqrt($8*$8) <=  0.01) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.bad.${TYPE}.diff

fi

#################################################################

if [ ${EVENT} == "A5SS" ] || [ ${EVENT} == "AD"  ] ; then

TYPE="AD"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="AD" && $6 >= 0.90 && $7 >= 0.90  && sqrt($8*$8) <=  0.01) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.bad.${TYPE}.diff

fi

#################################################################

if [ ${EVENT} == "RI" ] || [ ${EVENT} == "IR"  ] ; then


TYPE="RI"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="RI" && $6 >= 0.90 && $7 >= 0.90 && sqrt($8*$8) <=  0.01 ) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.bad.${TYPE}.diff

fi

#################################################################

if [ ${EVENT} == "TE"  ] ; then


TYPE="TE"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="TE" && $6 >= 0.90 && $7 >= 0.90 && sqrt($8*$8) <=  0.01 ) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.bad.${TYPE}.diff

fi

#################################################################

if [ ${EVENT} == "TS"  ] ; then


TYPE="TS"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="TS" && $6 >= 0.90 && $7 >= 0.90 && sqrt($8*$8) <=  0.01 ) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.bad.${TYPE}.diff

fi

#################################################################

if [ ${EVENT} == "AF"  ] ; then


TYPE="AF"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="AF" && $6 >= 0.90 && $7 >= 0.90 && sqrt($8*$8) <=  0.01 ) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.bad.${TYPE}.diff

fi

#################################################################

if [ ${EVENT} == "AL"  ] ; then


TYPE="AL"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="AL" && $6 >= 0.90 && $7 >= 0.90 && sqrt($8*$8) <=  0.01 ) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH}${FILE}.diff  > ${PATH}${FILE}.bad.${TYPE}.diff

fi

#################################################################

/bin/sed -i $'1 i\\\ngene\tcoordinates\tstrand\tevent\tpsiA\tpsiB\tdpsi\tprobability\tcomplexity\tentropy' ${PATH}${FILE}.bad.${TYPE}.diff


exit 0

#################################################################
#################################################################
#################################################################
