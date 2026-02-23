#!/bin/bash

#slurm options
#SBATCH -p intel-sc3			#选择多个分区用逗号隔开
#SBATCH -q normal				#Qos只能选一个，否则会报错
#SBATCH -J Oli_CFP			#作业名称LDSC
#SBATCH -c 1 					#申请8个CPU核心
#SBATCH --mem 30024 					#申请100G内存
#SBATCH -o /storage/yangjianLab/jiayuran/Coloc/result/brainmeta_4_ct/log/Oli_coloc_format_process-%j.log				#%j表示实际运行时的作业号


for num in {1..22}; do
    input_file="/storage/yangjianLab/jiayuran/bulk_seq/Brainmeta/four_ct/eQTL_mapping/Oli/BrainMeta_Oli_chr${num}.cis_qtl_pairs.${num}.parquet.txt"
    output_file="/storage/yangjianLab/jiayuran/Coloc/result/brainmeta_4_ct/data/Oli/Oli_chr${num}.txt"
    
    #change2query.pl是依据reference，补齐每条eqtl构建query格式所需的信息。其中，如果snp的等位基因信息有错误(即等位基因列不为纯大写字母构成)，则舍弃这个snp
    perl /storage/yangjianLab/jiayuran/Coloc/result/brainmeta_4_ct/script/tensorqtl2query_change2query.pl /storage/yangjianLab/jiayuran/Data/Gene_annotation_hg37/hg19_ENSGGene_LncRNA_annotation.bed /storage/yangjianLab/jiayuran/bulk_seq/Brainmeta/genotype/Brainmeta.chr1-22.bim "$input_file" "$output_file"
    #将具有冗余query格式所需信息的文件，转化成query格式
    python /storage/yangjianLab/jiayuran/Coloc/result/brainmeta_4_ct/script/tensorqtl2query_covert2query.py "$output_file" "$output_file"
done
