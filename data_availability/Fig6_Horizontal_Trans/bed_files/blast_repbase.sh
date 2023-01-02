#!/bin/bash

sixteen="4 15 28 39 57 63 2 31 56"
for NUM in $sixteen; do
	(
		blastn -query DN_mariner_sample$NUM.fasta -subject RepeatMasker.lib -outfmt 6 -max_hsps 1 -out HT_DN_mariner_$NUM.outfmt6
) 
done
