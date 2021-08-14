#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 120G
#SBATCH -c 1
#SBATCH -t 72:0:0


awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t""Rep1peakFile-"NR}' HeLa_ac4c_rep1_peaks.narrowPeak > HeLa_ac4c_rep1_peaks.narrowPeak_id
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t""Rep2peakFile-"NR}' HeLa_ac4c_rep2_peaks.narrowPeak > HeLa_ac4c_rep2_peaks.narrowPeak_id
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t""Rep1peakFile-"NR}' NAT10_rep1_peaks.narrowPeak > NAT10_rep1_peaks.narrowPeak_id
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t""Rep2peakFile-"NR}' NAT10_rep2_peaks.narrowPeak > NAT10_rep2_peaks.narrowPeak_id

cat HeLa_ac4c_rep1_peaks.narrowPeak_id HeLa_ac4c_rep2_peaks.narrowPeak_id | sort -k1,1 -k2,2n | mergeBed -i stdin -o collapse -c 8,11 >> HeLa_ac4c_All_peaks.narrowPeak.bed
cat NAT10_rep1_peaks.narrowPeak_id NAT10_rep2_peaks.narrowPeak_id | sort -k1,1 -k2,2n | mergeBed -i stdin -o collapse -c 8,11 >> NAT10_All_peaks.narrowPeak.bed

for file in /home/anksi/faststorage/Modification/Human/Arr_new_new/new/RIPseq/Peak_calling/Combined/1levelfile/*All_peaks.narrowPeak.bed; do
  intersectBed -a $file -b /home/anksi/faststorage/reference/UCSC/human/GRCh37/gencode.v36lift37.annotation.gtf -wa -wb >> ${file/.bed/.gencode_new.bed}
done

#filter
for file in /home/anksi/faststorage/Modification/Human/Arr_new_new/new/RIPseq/Peak_calling/Combined/1levelfile/*gencode_new.bed; do
  awk '{for (i = 1; i <= NF; i++) {if ($i ~ /gene_name/) {printf "%s ", $(i+1)}}print """\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}' $file | sort | uniq | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1"\t"$6}' | sed -e 's/"//g' -e 's/;//g' -e 's/ /\t/' >> ${file/.gencode_new.bed/.gencode_new_filter.bed}
done

#for HeLa and NAT10 Combined
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""NAT10"}' NAT10_All_peaks.narrowPeak.gencode_new_filter.bed >> NAT10_All_peaks.narrowPeak.gencode_new_name.bed
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""HeLa"}' HeLa_ac4c_All_peaks.narrowPeak.gencode_new_filter.bed >> HeLa_ac4c_All_peaks.narrowPeak.gencode_new_name.bed
cat NAT10_All_peaks.narrowPeak.gencode_new_name.bed HeLa_ac4c_All_peaks.narrowPeak.gencode_new_name.bed >> NAT10_HeLa_com.bed
awk '{print $1"_"$2"_"$3"\t"$4"\t"$5"\t"$6"\t"$7}' NAT10_HeLa_com.bed >> NAT10_HeLa_com_ID.bed

intersectBed -a NAT10_HeLa_com.bed -b /home/anksi/faststorage/reference/UCSC/human/hg19/3UTR_Exon_hg19.bed -wa -wb >> NAT10_HeLa_com_3UTR.bed
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""3UTR"}' NAT10_HeLa_com_3UTR.bed >> NAT10_HeLa_com_3UTR_new.bed
intersectBed -a NAT10_HeLa_com.bed -b /home/anksi/faststorage/reference/UCSC/human/hg19/5UTR_Exon_hg19.bed -wa -wb >> NAT10_HeLa_com_5UTR.bed
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""5UTR"}' NAT10_HeLa_com_5UTR.bed >> NAT10_HeLa_com_5UTR_new.bed
intersectBed -a NAT10_HeLa_com.bed -b /home/anksi/faststorage/reference/UCSC/human/hg19/Exon_hg19.bed -wa -wb >> NAT10_HeLa_com_exon.bed
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""exon"}' NAT10_HeLa_com_exon.bed >> NAT10_HeLa_com_exon_new.bed
intersectBed -a NAT10_HeLa_com.bed -b /home/anksi/faststorage/reference/UCSC/human/hg19/Intron_hg19.bed -wa -wb >> NAT10_HeLa_com_intron.bed
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""intron"}' NAT10_HeLa_com_intron.bed >> NAT10_HeLa_com_intron_new.bed

cat NAT10_HeLa_com*new.bed >> NAT10_HeLa_com_new.bed
cat NAT10_HeLa_com_new.bed | sort | uniq >> NAT10_HeLa_com_new_SU.bed
awk '{print $1"_"$2"_"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' NAT10_HeLa_com_new_SU.bed >> NAT10_HeLa_com_new_SU_ID.bed

# do GSEA analysis here

for file in /home/anksi/faststorage/Modification/Human/Arr_new_new/new/RIPseq/Peak_calling/Combined/1levelfile/*gencode_new_filter.bed; do
  intersectBed -a $file -b /home/anksi/faststorage/reference/UCSC/human/hg19/3UTR_Exon_hg19.bed -wa -wb >> ${file/.bed/.3UTR.bed}
done

for file in /home/anksi/faststorage/Modification/Human/Arr_new_new/new/RIPseq/Peak_calling/Combined/1levelfile/*gencode_new_filter.bed; do
  intersectBed -a $file -b /home/anksi/faststorage/reference/UCSC/human/hg19/5UTR_Exon_hg19.bed -wa -wb >> ${file/.bed/.5UTR.bed}
done

for file in /home/anksi/faststorage/Modification/Human/Arr_new_new/new/RIPseq/Peak_calling/Combined/1levelfile/*gencode_new_filter.bed; do
  intersectBed -a $file -b /home/anksi/faststorage/reference/UCSC/human/hg19/Exon_hg19.bed -wa -wb >> ${file/.bed/.exon.bed}
done

for file in /home/anksi/faststorage/Modification/Human/Arr_new_new/new/RIPseq/Peak_calling/Combined/1levelfile/*gencode_new_filter.bed; do
  intersectBed -a $file -b /home/anksi/faststorage/reference/UCSC/human/hg19/Intron_hg19.bed -wa -wb >> ${file/.bed/.intron.bed}
done
