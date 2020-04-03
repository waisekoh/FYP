
#to trim the first 10 base pairs of the reads
awk -v s=10 -v e=0 '{if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; }' MT1_1.fastq > Trim_MT1_1.fastq
awk -v s=10 -v e=0 '{if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; }' MT1_2.fastq > Trim_MT1_2.fastq
awk -v s=10 -v e=0 '{if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; }' MT2_1.fastq > Trim_MT2_1.fastq
awk -v s=10 -v e=0 '{if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; }' MT2_2.fastq > Trim_MT2_2.fastq
awk -v s=10 -v e=0 '{if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; }' MT3_1.fastq > Trim_MT3_1.fastq
awk -v s=10 -v e=0 '{if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; }' MT3_2.fastq > Trim_MT3_2.fastq
awk -v s=10 -v e=0 '{if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; }' WT1_1.fastq > Trim_WT1_1.fastq
awk -v s=10 -v e=0 '{if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; }' WT1_2.fastq > Trim_WT1_2.fastq
awk -v s=10 -v e=0 '{if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; }' WT2_1.fastq > Trim_WT2_1.fastq
awk -v s=10 -v e=0 '{if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; }' WT2_2.fastq > Trim_WT2_2.fastq
awk -v s=10 -v e=0 '{if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; }' WT3_1.fastq > Trim_WT3_1.fastq
awk -v s=10 -v e=0 '{if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; }' WT3_2.fastq > Trim_WT3_2.fastq

#build HISAT2 genome index
#need to download the genome.fa and snp.fa from HISAT2 website to the HISAT2_HOME folder

#hisat2-build /mnt/gtklab01/waise/Test/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --snp 
#pwd is HISAT2_HOME
/mnt/gtklab01/waise/hisat2-2.2.0/hisat2-build reference/genome.fa --snp reference/genome.snp mm10_snp
#OR
/mnt/gtklab01/waise/hisat2-2.2.0/hisat2-build -p 16 /mnt/gtklab01/waise/Test/Mus_musculus/Ensembl/NCBIM37/Sequence/WholeGenomeFasta/genome.fa genome

/mnt/gtklab01/waise/hisat2-2.2.0/hisat2-build -p 16 /mnt/gtklab01/waise/Test/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa genome

#stay in the same directory as the previous, where all the mm10_snp index files are at
/mnt/gtklab01/waise/hisat2-2.2.0/hisat2 -q --summary-file summary.txt -x genome -1 /mnt/gtklab01/waise/Test/Fastq/Trim_MT1_1.fastq -2 /mnt/gtklab01/waise/Test/Fastq/Trim_MT1_2.fastq -S MT1.sam 

#HISAT2 pipe
for file in MT1 MT2 MT3 WT1 WT2 WT3
do
	/mnt/gtklab01/waise/hisat2-2.2.0/hisat2 -p 6 -q --summary-file summary.txt -x /mnt/cbis/home/waise/redo/HISAT/genome -1 /mnt/gtklab01/waise/Test/new/fastq/${file}_1.fastq.gz -2 /mnt/gtklab01/waise/Test/new/fastq/${file}_2.fastq.gz -S /mnt/cbis/home/waise/redo/HISAT/sam/${file}.sam
	samtools view -bS /mnt/cbis/home/waise/redo/HISAT/sam/${file}.sam > /mnt/cbis/home/waise/redo/HISAT/bam/${file}.bam
	rm /mnt/cbis/home/waise/redo/HISAT/sam/${file}.sam
	samtools sort -o /mnt/cbis/home/waise/hnRNP/HISAT2/bam/${file}.sorted.bam /mnt/cbis/home/waise/hnRNP/HISAT2/bam/${file}.bam
	samtools index /mnt/cbis/home/waise/hnRNP/HISAT2/bam/${file}.sorted.bam
	rm /mnt/cbis/home/waise/hnRNP/HISAT2/bam/${file}.bam
	/mnt/gtklab01/waise/Test/regtools/build/./regtools junctions extract -s 0 -o /mnt/cbis/home/waise/hnRNP/HISAT2/bam/${file}.bed /mnt/cbis/home/waise/hnRNP/HISAT2/bam/${file}.sorted.bam
done

#STAR built already
#STAR pipe
# need to call STAR from its source folder

for file in MT1 MT2 MT3 WT1 WT2 WT3
do
	/mnt/gtklab01/waise/STAR-2.7.3a/source/./STAR --genomeDir /mnt/gtklab01/waise/Test/Mus_musculus/UCSC/mm10/Sequence/STARIndex --runThreadN 10 --readFilesIn <(gunzip -c /mnt/gtklab01/waise/Test/new/fastq/${file}_1.fastq.gz /mnt/gtklab01/waise/Test/new/fastq/${file}_2.fastq.gz) --outFileNamePrefix /mnt/cbis/home/waise/hnRNP/STAR/bam/${file}/ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes All --outSAMstrandField intronMotif
	mv /mnt/cbis/home/waise/hnRNP/STAR/bam/${file}/Aligned.sortedByCoord.out.bam /mnt/cbis/home/waise/hnRNP/STAR/bam/${file}.sorted.bam
	#rm /mnt/cbis/home/waise/hnRNP/STAR/bam/${file}/*
	#rmdir /mnt/cbis/home/waise/hnRNP/STAR/bam/${file}
	samtools index /mnt/cbis/home/waise/hnRNP/STAR/bam/${file}.sorted.bam
	/mnt/gtklab01/waise/Test/regtools/build/./regtools junctions extract -s 0 -o /mnt/cbis/home/waise/hnRNP/STAR/bam/${file}.bed /mnt/cbis/home/waise/hnRNP/STAR/bam/${file}.sorted.bam
done

#build GSNAP genome index
gmap_build -D MM10USCS /mnt/gtklab01/waise/Test/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa
#GSNAP pipe
for file in MT1 MT2 MT3 WT1 WT2 WT3
do
	gsnap -d MM10USCS -t 6 -A sam -N 1 <(gunzip -c /mnt/gtklab01/waise/Test/hnRNP/fastq/${file}_1.fastq.gz /mnt/gtklab01/waise/Test/hnRNP/fastq/${file}_2.fastq.gz) > /mnt/cbis/home/waise/hnRNP/GSNAP/sam/${file}.sam
	samtools view -bS /mnt/cbis/home/waise/hnRNP/GSNAP/sam/${file}.sam > /mnt/cbis/home/waise/hnRNP/GSNAP/bam/${file}.bam
	rm /mnt/cbis/home/waise/hnRNP/GSNAP/sam/${file}.sam
	samtools sort -o /mnt/cbis/home/waise/hnRNP/GSNAP/bam/${file}.sorted.bam /mnt/cbis/home/waise/hnRNP/GSNAP/bam/${file}.bam
	rm /mnt/cbis/home/waise/hnRNP/GSNAP/bam/${file}.bam
	samtools index /mnt/cbis/home/waise/hnRNP/GSNAP/bam/${file}.sorted.bam
	/mnt/gtklab01/waise/Test/regtools/build/./regtools junctions extract -s 0 -o /mnt/cbis/home/waise/hnRNP/GSNAP/bam/${file}.bed /mnt/cbis/home/waise/hnRNP/GSNAP/bam/${file}.sorted.bam
done