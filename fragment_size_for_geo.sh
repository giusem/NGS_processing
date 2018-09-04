

infolder="$1"
current_dir=$(pwd)

#fastq subsampling

echo "fastq subsampling"

mkdir $current_dir/fastq_subsamples

for i in $infolder/*.fastq.gz; do base=$(basename $i); \
zcat $i | head -n 1000000 > fastq_subsamples/${base%.fastq.gz}_sub.fastq; \
done

for i in fastq_subsamples/*_sub.fastq; do gzip $i; done

#mapping

echo "mapping"

mkdir $current_dir/bowtie2_mapping

for i in $current_dir/fastq_subsamples/*R1_sub.fastq.gz; \
do base=$(basename $i); \
bam=${base%_R1_sub.fastq.gz}_sub; logname=${base%_R1_sub.fastq.gz}_sub_bowtie2log.txt; \
bowtie2 -x /data/repository/organisms/GRCh37_ensembl/BowtieIndex/genome   \
-1 $i -2 ${i%R1_sub.fastq.gz}R2_sub.fastq.gz \
--end-to-end -p 32 2> $current_dir/bowtie2_mapping/$logname | \
samtools view -b -q 10 - | samtools sort - $current_dir/bowtie2_mapping/$bam; \
done

for i in $current_dir/bowtie2_mapping/*.bam; do samtools index $i; done

#calculate fragment size

echo "calculate fragment size"

mkdir $current_dir/PE_fragment_size

> $current_dir/PE_fragment_size/summary_inner_distance.txt

echo -e "Sample\tMean\tStdev" >> $current_dir/PE_fragment_size/summary_inner_distance.txt


for i in $current_dir/bowtie2_mapping/*.bam; \
do base=$(basename $i); \
inner_distance.py -i $i -o $current_dir/PE_fragment_size/${base%.bam}_temp -k 250000 \
-r /data/repository/organisms/GRCh37_ensembl/ensembl/release-75/genes.bed  | \
tee $current_dir/PE_fragment_size/${base%.bam}_inner_distance.txt; \
rm $current_dir/PE_fragment_size/${base%.bam}_temp*; \
mean=$(sed -n '2p' $current_dir/PE_fragment_size/${base%.bam}_inner_distance.txt); \
mean=$(printf "%.0f\n" $mean)
sd=$(sed -n '6p' $current_dir/PE_fragment_size/${base%.bam}_inner_distance.txt); \
sd=$(printf "%.0f\n" $sd); \
sample=${base%_sub.bam}; \
echo -e "$sample\t$mean\t$sd" >> $current_dir/PE_fragment_size/summary_inner_distance.txt; \
done

