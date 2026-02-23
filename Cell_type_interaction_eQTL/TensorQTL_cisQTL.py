import pandas as pd
import sys
import numpy as np
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')
 


def tensorEQTL_process(plink_path_genotypes,
                       expression_bed_phenotypes,
                       cell_proportion_covariates,
                       covarate,
                       cohort_name,
                       cell_name):
    
    plink_prefix_path = plink_path_genotypes   		
    expression_bed = expression_bed_phenotypes 		
    covariates_cell_proportion = pd.read_csv(cell_proportion_covariates,
								sep='\t',
								index_col=0)
    
    prefix = f"{cohort_name}_{cell_name}"
    
    
    # load phenotypes and covariates
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
    print("phenotype: ")
    print(phenotype_df)
    print("phenotype_pos: ")
    print(phenotype_pos_df)
    #pc = list(phenotype_df.columns)
    covariates_df = pd.read_csv(covarate,
								sep='\t',
								index_col=0)
    interaction = covariates_cell_proportion.loc[:, [cell_name]] 
    #ci=list(covariates_df.index)
    print("interaction cell proportion:")
    print(interaction)
    
    
    # PLINK reader for genotypes
    pr = genotypeio.PlinkReader(plink_prefix_path)
    genotype_df = pr.load_genotypes()
    print("genotype:")
    print(genotype_df)
    
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    print("genotype variant:")
    print(variant_df)
    
    
    for num in range(22):
        num1=num+1
        str1=str(num1)
        prefix = f"{cohort_name}_{cell_name}"
        prefix = prefix + "_chr" + str1
        cis.map_nominal(genotype_df, variant_df, 
					phenotype_df.loc[phenotype_pos_df['chr']==str1], 
					phenotype_pos_df.loc[phenotype_pos_df['chr']==str1],
                    prefix,
					maf_threshold_interaction=0.05,   
					covariates_df=covariates_df,
                    interaction_df=interaction,    
                    run_eigenmt=True
					)  
        
        filename=f'{prefix}.cis_qtl_pairs.'+str1
        fileaddress=f'/storage/eQTL_mapping/{cell_name}/'+filename +'.parquet'    
        pairs_df = pd.read_parquet(fileaddress)
        pairs_df.to_csv(fileaddress+".txt",sep="\t",index=0)



if __name__ == "__main__":
    tensorEQTL_process("Brainmeta.chr1-22",
                       "Brainmeta.RINT.sort.bed.gz",
                       "BayesPrism_estimate_brainmeta_RINT.txt",
                       "merge_covar.txt",
                       "BrainMeta",  
                       "Ex")  #Ex	Oli	In	Ast

	


                 
                

