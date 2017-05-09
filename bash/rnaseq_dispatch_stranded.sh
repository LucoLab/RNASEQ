set -ue
SAMTOOLS=$(which samtools)
#################################################################
#
#date: apr 04, 2017
#platform: Ubuntu 16.04
#author: Villemin Jean-Philippe
#team: Epigenetic Component of Alternative Splicing - IGH
#
# rnaseq_dispatch_stranded.sh
# Usage : 
# rnaseq_dispatch_stranded.sh pathToBam nameOutput DirOutput
#
# For stranded RNASEQ, library type ISR(Salmon),fr-secondstrand(RMATS) separe reads per strand
#
#Source : 
#
#https://www.biostars.org/p/92935/
#
# Bamtools (other solution)
#https://www.biostars.org/p/88582/
#
#################################################################

# Get the bam file from the command line
DATA=$1
NAME=$2
DIR=$3
# Forward strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
$SAMTOOLS view -b -f 128 -F 16 ${DATA} > ${DIR}/${NAME}_fwd1.bam
$SAMTOOLS index ${DIR}/${NAME}_fwd1.bam

$SAMTOOLS view -b -f 80 ${DATA} > ${DIR}/${NAME}_fwd2.bam
$SAMTOOLS index ${DIR}/${NAME}_fwd2.bam

$SAMTOOLS merge -f ${DIR}/${NAME}_pos.bam ${DIR}/${NAME}_fwd1.bam ${DIR}/${NAME}_fwd2.bam
$SAMTOOLS index ${DIR}/${NAME}_pos.bam
$SAMTOOLS flagstat ${DIR}/${NAME}_pos.bam > ${DIR}/${NAME}_pos.txt

rm ${DIR}/${NAME}_fwd1.bam
rm ${DIR}/${NAME}_fwd2.bam
rm ${DIR}/${NAME}_fwd1.bam.bai
rm ${DIR}/${NAME}_fwd2.bam.bai

# Reverse strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#
$SAMTOOLS view -b -f 144 ${DATA} > ${DIR}/${NAME}_rev1.bam
$SAMTOOLS index ${DIR}/${NAME}_rev1.bam

$SAMTOOLS view -b -f 64 -F 16 ${DATA} > ${DIR}/${NAME}_rev2.bam
$SAMTOOLS index ${DIR}/${NAME}_rev2.bam

#
# Combine alignments that originate on the reverse strand.
#
$SAMTOOLS merge -f ${DIR}/${NAME}_neg.bam ${DIR}/${NAME}_rev1.bam ${DIR}/${NAME}_rev2.bam
$SAMTOOLS index ${DIR}/${NAME}_neg.bam
$SAMTOOLS flagstat ${DIR}/${NAME}_neg.bam > ${DIR}/${NAME}_neg.txt

rm ${DIR}/${NAME}_rev2.bam
rm ${DIR}/${NAME}_rev1.bam
rm ${DIR}/${NAME}_rev1.bam.bai
rm ${DIR}/${NAME}_rev2.bam.bai
#
#If keeping full pairs is essential then the filtering steps need to be performed to operate on mates like so:
#
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if their mate maps to the forward strand
#
#$SAMTOOLS view -b -f 128 -F 16 $DATA > fwd1.bam
#$SAMTOOLS index fwd1.bam

#$SAMTOOLS view -b -f 64 -F 32 $DATA > fwd2.bam
#$SAMTOOLS index fwd2.bam
#and similarly:

# 1. alignments of the second in pair if it maps to the reverse strand
# 2. alignments of the first in pair if their mates map to the reverse strand
#
#$SAMTOOLS view -b -f 144 $DATA > rev1.bam
#$SAMTOOLS index rev1.bam

#$SAMTOOLS view -b -f 96 $DATA > rev2.bam
#$SAMTOOLS index rev2.bam

