#!/bin/bash

## Useful shortcuts
export refdir=/data/REFS

## Change --sjdbOverhang to length of your sequence data /2 minus 1

echo "\n\n I TOLD YOU NOT TO RUN THIS ONE NOW! \n\n (unless you're in the future and trying to run this for real, in which case you need to edit this script and remove the # characters from the command)"

STAR 	\
	--runThreadN 8 \
        --limitGenomeGenerateRAM 321563573 \
	--genomeSAindexNbases 12 \
	--runMode genomeGenerate \
	--genomeDir  $refdir \
	--genomeFastaFiles $refdir/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa \
	--sjdbGTFfile $refdir/Arabidopsis_thaliana.TAIR10.53.gtf \
	--sjdbOverhang 49
