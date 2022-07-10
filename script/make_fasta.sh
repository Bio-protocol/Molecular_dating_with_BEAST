## Generate fasta for a gene in a sample (for parallelization)
## run as: cat samplelist.txt | parallel -j 10 bash make_fasta.sh {BAM_DIR} {ref}
## Kaichi Huang 2022 May

BAM_DIR=$1
ref=$2
gene_name=$3
sample=$4 # passed from pipe

bam="$BAM_DIR/$sample.bam"
program/gatk-4.0.8.1/gatk HaplotypeCaller -R $ref -I $bam --intervals cache/tmp.$gene_name.bed -ERC BP_RESOLUTION -O cache/tmp.$gene_name.$sample.g.vcf
cat cache/tmp.$gene_name.$sample.g.vcf | perl script/gvcf2fasta_nogaps.pl > cache/tmp.$gene_name.$sample.fa
