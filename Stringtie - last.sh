

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







