for file in MT1 MT2 MT3 WT1 WT2 WT3
do
	stringtie -p 8 -l ${file} -o /mnt/cbis/home/waise/hnRNP/GSNAP/Stringtie/denovo/${file}/transcripts.gtf /mnt/gtklab01/waise/Test/hnRNP/GSNAP/bam/${file}.sorted.bam
done

cd /mnt/cbis/home/waise/hnRNP/GSNAP/Stringtie/denovo/
ls -1 *T*/transcripts.gtf > assembly_GTF_list.txt
stringtie --merge -p 8 -o stringtie_merged.gtf assembly_GTF_list.txt
gffcompare -r /mnt/gtklab01/waise/Test/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o gffcompare stringtie_merged.gtf

for file in MT1 MT2 MT3 WT1 WT2 WT3
do
	stringtie -p 8 -G /mnt/cbis/home/waise/hnRNP/GSNAP/Stringtie/denovo/stringtie_merged.gtf -e -B -o /mnt/cbis/home/waise/hnRNP/GSNAP/Stringtie/denovo_merged/${file}/transcripts.gtf /mnt/gtklab01/waise/Test/hnRNP/GSNAP/bam/${file}.sorted.bam
done

cd /mnt/cbis/home/waise/hnRNP/GSNAP/Stringtie/
printf "\"ids\",\"type\",\"path\"\n\"MT1\",\"MT\",\"/mnt/cbis/home/waise/hnRNP/GSNAP/Stringtie/denovo_merged/MT1\"\n\"MT2\",\"MT\",\"/mnt/cbis/home/waise/hnRNP/GSNAP/Stringtie/denovo_merged/MT2\"\n\"MT3\",\"MT\",\"/mnt/cbis/home/waise/hnRNP/GSNAP/Stringtie/denovo_merged/MT3\"\n\"WT1\",\"WT\",\"/mnt/cbis/home/waise/hnRNP/GSNAP/Stringtie/denovo_merged/WT1\"\n\"WT2\",\"WT\",\"/mnt/cbis/home/waise/hnRNP/GSNAP/Stringtie/denovo_merged/WT2\"\n\"WT3\",\"WT\",\"/mnt/cbis/home/waise/hnRNP/GSNAP/Stringtie/denovo_merged/WT3\"\n" > MT_vs_WT.csv

#   
#  35119 reference transcripts loaded.
#  843 duplicate reference transcripts discarded.
#  34791 query transfrags loaded.