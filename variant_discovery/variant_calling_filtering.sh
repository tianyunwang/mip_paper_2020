mipgen_dir=/path/to/mipgen
work_dir=/path/to/work/directory

# load software modules
module load python/2.7.3 numpy/1.8.1 scipy/0.14.0 bwa/0.7.3 pear/0.9.5 samtools/0.1.7 freebayes/1.0.2-6-g3ce827d

# MIPgen analysis
python $mipgen_dir/mipgen_fq_cutter_pe.py $work_dir/read1.fastq.gz $work_dir/read2.fastq.gz -i $work_dir/index.fastq.gz -tb /path/to/barcode_file.txt -o $work_dir/sample.barcoded
pear -f $work_dir/sample.barcoded.r1.indexed.fq -r $work_dir/sample.barcoded.r2.indexed.fq -o $work_dir/sample.barcoded
python $mipgen_dir/mipgen_fq_cutter_se.py $work_dir/sample.barcoded.assembled.fastq -tb /path/to/barcode_file.txt -m 0,5 -o $work_dir/sample.barcoded
bwa mem -t 6 $work_dir/hg19.fa $work_dir/sample.barcoded.indexed.fq > $work_dir/sample.barcoded.indexed.sam
samtools view -bS $work_dir/sample.barcoded.indexed.sam | samtools sort - $work_dir/sample.barcoded.indexed.sort
samtools view -h $work_dir/sample.barcoded.indexed.sort.bam | python $mipgen_dir/mipgen_smmip_collapser.py 5 $work_dir/sample.barcoded.indexed.sort.collapse -m /path/to/mip_design_file.txt -f 1 -c -r -b /path/to/barcode_file.txt -s
samtools view -bS -T $work_dir/hg19.fa $work_dir/sample.barcoded.indexed.sort.collapse.all_reads.unique.sam | samtools sort - $work_dir/sample.barcoded.indexed.sort.collapse.all_reads.unique.sam.sort
samtools index $work_dir/sample.barcoded.indexed.sort.collapse.all_reads.unique.sam.sort.bam

# variant calling
freebayes -f $work_dir/hg19.fa $work_dir/sample.barcoded.indexed.sort.collapse.all_reads.unique.sam.sort.bam > unfiltered.vcf

# variant filtering
vcffilter -g "DP > 8" -f "QUAL > 20" unfiltered.vcf | vcffixup - > DP_Qual_filtered.vcf
vcffilter -f "AC > 0 & AC < 4" -s DP_Qual_filtered.vcf > DP_Qual_freq_filtered.vcf

echo "analysis commands have terminated"
