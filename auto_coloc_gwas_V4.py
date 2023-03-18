"""
Used to perform coloc between 2 GWAS traits
"""


##
import pandas as pd
import math
import numpy as np
import sys
import argparse
import warnings
warnings.filterwarnings('ignore')
from itertools import repeat
import sumstats
import h5py
import random
from hashlib import md5
import pyranges as pr
import os.path
from tqdm import tqdm


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Colocalisation test between SNPs from 2 GWAS traits ' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g1', '--gwas1', help='HDF5 Summary statistics from GWAS trait 1',required='True')
    parser.add_argument('-g2', '--gwas2', help='HDF5 Summary statistics from GWAS trait 2',required='True')
    parser.add_argument('-l', '--list', help='text file of SNPs to test - one SNP per a line' , default="list.txt" )
    parser.add_argument('-o', '--output', help='File name of the results', default='results.txt')
    parser.add_argument('-kb', '--kb_window', help='The window size around the SNP, DEFAULT is 100000 ', default=100)
    parser.add_argument('-c', '--cytoband', help='file that contains cytoband information with headers chr,start,end,cytoband for the current reference genome - default is cytoband_hg19.txt', default='cytoband_hg19.txt')
    results = parser.parse_args(args)
    return (results.gwas1 , results.gwas2 , results.list  , results.output , results.kb_window , results.cytoband)

def coloc1(
    trait1_lnbfs,
    trait2_lnbfs,
    prior1: float = 1e-4,
    prior2: float = 1e-4,
    prior12: float = 1e-5 ):
    """
    COLOC function take from Anthony Aylward
    https://github.com/anthony-aylward/coloc
    """
    log_numerators = (
        0,
        math.log(prior1) + sumstats.log_sum(trait1_lnbfs),
        math.log(prior2) + sumstats.log_sum(trait2_lnbfs),
        math.log(prior1) + math.log(prior2) + sumstats.log_sum(
            trait1_lnbf + trait2_lnbf
            for i, trait1_lnbf in enumerate(trait1_lnbfs)
            for j, trait2_lnbf in enumerate(trait2_lnbfs)
            if i != j
        ),
        math.log(prior12) + sumstats.log_sum(
            trait1_lnbf + trait2_lnbf
            for trait1_lnbf, trait2_lnbf in zip(trait1_lnbfs, trait2_lnbfs)
        )
    )
    yield from  (
        math.exp(log_numerator - sumstats.log_sum(log_numerators))
        for log_numerator in log_numerators
    )
    
def get_PP4(
    title: str,
    pp0: float,
    pp1: float,
    pp2: float,
    pp3: float,
    pp4: float):
    return pp4


def get_PP(df, snp):
    '''
    This is to do the PP for the coloc
    '''
    PP4 = get_PP4("snp", *coloc1((df.gwas_1_lbf.to_numpy()),(df.gwas_2_lbf.to_numpy() )))
    PP4_rounded = round(PP4,2)
    raw = [[ snp , PP4_rounded]]
    df2 = pd.DataFrame(data=raw, columns=["SNP", "PP4"])
    return df2



def get_snp_window( kb_window , SNP, SNP_file , gwas1 , gwas2, cyto):
    """ 
    extracts the given window and return the df
    """
    small = SNP_file[SNP_file["SNP"] == SNP]
    small = small.reset_index()
    position = small.at[0, 'Start']
    upper = int(position) + kb_window
    lower = int(position) - kb_window
    chrom = small.Chromosome.astype(int)[0]
    data = { 'Chromosome' : [chrom] , 'Start' : [lower] , 'End' : [upper]}
    gr3 = pr.from_dict(data)
    cyto_temp = cyto.join(gr3)
    cyto_temp = cyto_temp.df
    cytoband_list = cyto_temp.cytoband.to_list()
    gwas1_df = pd.DataFrame()
    gwas2_df = pd.DataFrame()
    try:
        for i in cytoband_list : 
            temp1 = pd.read_hdf(gwas1, i)
            temp2 = pd.read_hdf(gwas2, i)
            gwas1_df = gwas1_df.append(temp1)
            gwas2_df = gwas2_df.append(temp2)
    except:
        return
    
    if gwas1_df.isin([SNP]).any().any() :
        gwas1_clean = gwas1_df[(gwas1_df["Start"] >= lower ) & (gwas1_df["Start"] <= upper )]
        gwas1_clean = gwas1_clean[["SNP", "lbf", "Z" , "ref" , "alt", "test_allele"]].copy()
        gwas1_clean.columns=["SNP","gwas_1_lbf", "gwas_1_Z" , "gwas_1_ref" , "gwas_1_alt", "gwas_1_test_allele"]
       
    if gwas2_df.isin([SNP]).any().any() : 
        gwas2_clean = gwas2_df[(gwas2_df["Start"] >= lower ) & (gwas2_df["Start"] <= upper )]
        gwas2_clean = gwas2_clean[["SNP", "lbf", "Z" , "ref" , "alt", "test_allele", "se"]].copy()
        gwas2_clean.columns=["SNP","gwas_2_lbf", "gwas_2_Z" , "gwas_2_ref" , "gwas_2_alt", "gwas_2_test_allele", "gwas_2_se"]
        
    
    ## deal with allele flips -> if not drop snp
    merged = gwas1_clean.merge(gwas2_clean, on="SNP")
    bad_snps = merged[merged["gwas_1_test_allele"]!= merged["gwas_2_test_allele"]  ]
    good_snps = merged[~merged['SNP'].isin(bad_snps["SNP"])]
    bad_snps2 = bad_snps[(bad_snps["gwas_1_test_allele"] == bad_snps["gwas_2_ref"]) | (bad_snps["gwas_1_test_allele"] == bad_snps["gwas_2_alt"])]
    bad_snps2["gwas_2_Z"] = bad_snps2["gwas_2_Z"] * -1 
    #recal the LBF
    W = 0.2
    bad_snps2["V"] = bad_snps2["gwas_2_se"] ** 2 
    bad_snps2["R"] =  W **2 / (W**2 + bad_snps2["V"])
    bad_snps2["gwas_2_lbf"] =  0.5 * (np.log(1 - bad_snps2["R"]) + (bad_snps2["R"] * bad_snps2["gwas_2_Z"]**2))   
    bad_snps2 = bad_snps2.drop(["V", "R"], 1 )
    final = pd.concat([good_snps, bad_snps2])
    final = final[["SNP", "gwas_1_lbf" , "gwas_2_lbf" ]].copy()    
    return final


def error_snps(SNP):
    """
    send back a naughty SNP
    """
    raw = [[ SNP , np.nan]]
    df2 = pd.DataFrame(data=raw, columns=["SNP", "PP4"])
    return df2

def process_snps(kb_window , SNP, SNP_file , gwas1 , gwas2, cyto):
    """
    Process each SNP 
    """
    try:
        merged = get_snp_window(kb_window , SNP, SNP_file , gwas1 , gwas2, cyto)
        results = get_PP(merged, SNP)
        return results
    except: 
       results = error_snps(SNP)
       return results    

def main(gwas1 , gwas2 , list_file ,  output , kb_window, cyto):
    kb_window = round((int(kb_window) * 1000) / 2, 1 )
    cyto_df = pd.read_csv(cyto, sep="\t")
    gr = pr.PyRanges(cyto_df)
    SNP_list_raw = pd.read_csv(list_file , sep="\t", header=None)
    SNP_list_raw.columns = ["Chromosome" , "Start" , "SNP"]
    SNP_list_raw["End"] = SNP_list_raw.Start.astype(int) + 1 
    SNP_list_raw = SNP_list_raw[["Chromosome" , "Start", "End" , "SNP"]]
    gr2 = pr.PyRanges(SNP_list_raw)
    SNP_list = gr2.join(gr)
    SNP_list = SNP_list.df
    flat_list = SNP_list.SNP.tolist()
    results = pd.DataFrame()
    if os.path.isfile(gwas1) and os.path.isfile(gwas2) :
        for snp in tqdm(flat_list):
            temp = process_snps(kb_window ,snp, SNP_list , gwas1 , gwas2 , gr) 
            results = results.append(temp)
        results.to_csv(output, sep="\t", index=None)
    else:
        print ("Error - check that the GWAS files exist or the path is correct")
       
            
if __name__ == '__main__':
    gwas1 , gwas2 , list_file  , output , kb_window, cyto = check_arg(sys.argv[1:])
    main(gwas1 , gwas2 , list_file  , output , kb_window, cyto)
