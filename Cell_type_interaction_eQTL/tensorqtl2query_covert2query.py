import pandas as pd
import argparse



def convert_to_query(eqtl_file,output):
    eqtl = pd.read_csv(eqtl_file,sep="\t")
    print(eqtl)
    
    eqtl.rename(columns={'phenotype_id': 'Probe'}, inplace=True)
    eqtl.rename(columns={'variant_id': 'SNP'}, inplace=True)
    eqtl.rename(columns={'b_gi': 'beta'}, inplace=True)
    eqtl.rename(columns={'pval_gi': 'p'}, inplace=True)
    eqtl.rename(columns={'af': 'Freq'}, inplace=True)
    eqtl.rename(columns={'b_gi_se': 'se'}, inplace=True)
    eqtl.rename(columns={'start_distance': 'tss_distance'}, inplace=True)
    
    # 指定要保留的列
    columns_to_keep = ['SNP', 'Chr', 'BP','tss_distance' ,'ma_samples','ma_count','A1', 'A2', 'Freq', 'Probe', 'Probe_bp', 'beta', 'se', 'p']
    # 使用loc方法仅保留指定的列
    eqtl = eqtl.loc[:, columns_to_keep]
    eqtl['varbeta'] = eqtl['se'] ** 2
    new_order = ['SNP', 'Chr', 'BP','tss_distance' ,'ma_samples','ma_count','A1', 'A2', 'Freq', 'Probe', 'Probe_bp', 'beta', 'se','varbeta', 'p']
    eqtl = eqtl[new_order]
    

    print(eqtl)
    eqtl.to_csv(output,index=False)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Keep specific columns of eQTL file')
    
    parser.add_argument('eqtl_file', help='Input eQTL file')
    
    parser.add_argument('output', help='Output file')
    
    args = parser.parse_args()
    
    convert_to_query(args.eqtl_file, args.output)
