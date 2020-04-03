for file in MT1 MT2 MT3 WT1 WT2 WT3
do
	stringtie -p 8 -l ${file} -o /mnt/cbis/home/waise/hnRNP/STAR/Stringtie/denovo/${file}/transcripts.gtf /mnt/gtklab01/waise/Test/new/STAR/bam/${file}.sorted.bam
done

cd /mnt/cbis/home/waise/hnRNP/STAR/Stringtie/denovo/
ls -1 *T*/transcripts.gtf > assembly_GTF_list.txt
stringtie --merge -p 8 -o stringtie_merged.gtf assembly_GTF_list.txt
gffcompare -r /mnt/gtklab01/waise/Test/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o gffcompare stringtie_merged.gtf

for file in MT1 MT2 MT3 WT1 WT2 WT3
do
	stringtie -p 8 -G /mnt/cbis/home/waise/hnRNP/STAR/Stringtie/denovo/stringtie_merged.gtf -e -B -o /mnt/cbis/home/waise/hnRNP/STAR/Stringtie/denovo_merged/${file}/transcripts.gtf /mnt/gtklab01/waise/Test/new/STAR/bam/${file}.sorted.bam
done

cd /mnt/cbis/home/waise/hnRNP/STAR/Stringtie/
printf "\"ids\",\"type\",\"path\"\n\"MT1\",\"MT\",\"/mnt/cbis/home/waise/hnRNP/STAR/Stringtie/denovo_merged/MT1\"\n\"MT2\",\"MT\",\"/mnt/cbis/home/waise/hnRNP/STAR/Stringtie/denovo_merged/MT2\"\n\"MT3\",\"MT\",\"/mnt/cbis/home/waise/hnRNP/STAR/Stringtie/denovo_merged/MT3\"\n\"WT1\",\"WT\",\"/mnt/cbis/home/waise/hnRNP/STAR/Stringtie/denovo_merged/WT1\"\n\"WT2\",\"WT\",\"/mnt/cbis/home/waise/hnRNP/STAR/Stringtie/denovo_merged/WT2\"\n\"WT3\",\"WT\",\"/mnt/cbis/home/waise/hnRNP/STAR/Stringtie/denovo_merged/WT3\"\n" > MT_vs_WT.csv

#  35119 reference transcripts loaded.
#  843 duplicate reference transcripts discarded.
#  26153 query transfrags loaded.

cat stringtie_merged.gtf | awk '{print $1, $4, $5, $9, $10}' > temp_stringtie.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5}' temp_stringtie.bed > stringtie1.bed
grep -v '#' stringtie1.bed > stringtie.bed
#choose only the cols you want
cat genes.gtf | awk '{print $1, $4, $5, $9, $10}' > temp_genes.bed
#make tab delimited
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5}' temp_genes.bed > genes.bed
bedtools intersect -a genes.bed -b stringtie.bed -wa -wb > annotate.txt
awk 'BEGIN{OFS="\t"} {print $10, $5}' annotate.txt > genename.txt
# first sed removes all "", second sed removes all ;
cat genename.txt| sed 's/"//g' | sed 's/;//g' > final.txt
#removes duplicates
awk '!seen[$0]++' final.txt > last.txt
echo -e "id\tgene_name" | cat - last.txt > end.txt
awk 'NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{print $0, a[$1,$2]}' end.txt MT_vs_WT
awk 'FNR==NR{A[$1]=$NF;next} FNR!=NR && FNR>1{Q=$1;$0=(Q in A)?$0 FS A[$1]:$0 FS;print;next} {print}' FS="\t" denovo/end.txt MT_vs_WT_gene_results_sig.tsv > output.txt







