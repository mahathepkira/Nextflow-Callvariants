#!/bin/bash
FILE="all_gvcf.list"

for x in $(cat ${FILE});
do
	tabix -p vcf $x
done
