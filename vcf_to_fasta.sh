mkdir fasta_vcf

mkdir fasta_FINAL

for K in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32;
do

	~/Programs/gatk-4.2.0.0/gatk \
	FastaAlternateReferenceMaker \
	-R /Users/diego/Desktop/covid_IMT2/REF/sars_MN908947.fa \
	-V vcf/COVID_IMT.$K.filterDP.selected.vcf \
	-O fasta_vcf/COVID_IMT_1.$K.fasta


	cp fasta_vcf/COVID_IMT_1.$K.fasta fasta_FINAL/COVID_IMT_1.$K.fasta

	python ~/Dropbox/posdoc-imt/COVID_APODI/adjust_fasta.py -f fasta_FINAL/COVID_IMT_1.$K.fasta -d DP/COVID_IMT_1.$K.DEPTH.txt

	cd fasta_FINAL/

	awk '/^>/ {gsub(/.fa(sta)?$/,"",FILENAME);printf(">%s\n",FILENAME);next;} {print}' COVID_IMT_1.$K.fasta > COVID_IMT_1.$K.FINAL.fasta

	cd ..

	rm fasta_FINAL/COVID_IMT_1.$K.fasta fasta_FINAL/COVID_IMT_1.$K.fasta.fai low_cov_table.txt

done


