#!/bin/bash

sixteen="2 31"
for NUM in $sixteen; do
	(
		bedtools getfasta -fi cleaned6_sample$NUM.fasta -bed marinersample$NUM.bed -fo DN_mariner_sample$NUM.fasta
) 
done

