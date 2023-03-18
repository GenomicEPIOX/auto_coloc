import h5py
import warnings
from tables import NaturalNameWarning
warnings.filterwarnings('ignore', category=NaturalNameWarning)
import numpy as np
import pandas as pd
import pyranges as pr
import sys
import argparse


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Colocalisation test between SNPs from 2 GWAS traits ' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='Summary statistics from GWAS trait ',required='True')
    parser.add_argument('-o', '--output', help='File name of the results', default='results.txt')
    parser.add_argument('-b', '--beta', help='Name of the beta column in GWAS1 if its not named as beta', default='beta')
    parser.add_argument('-se', '--se', help='Name of the se column in GWAS1 if its not named as se', default='se')
    parser.add_argument('-var', '--variant_col', help='Name of the rsid column in gwas1 - default is SNP ', default='SNP')
    parser.add_argument('-c', '--cytoband', help='file that contains cytoband information with headers chr,start,end,cytoband for the current reference genome - default is cytoband_hg19.txt', default='cytoband_hg19.txt')
    parser.add_argument('-chr', '--chrom_col', help='Name of bp column ', default='chr')
    parser.add_argument('-p', '--position_col', help='Name of bp column ', default='pos')
    parser.add_argument('-r', '--ref_col', help='Name of the reference allele column ', default='ref')
    parser.add_argument('-a', '--alt_col', help='Name of the alt_col allele column ', default='alt')
    parser.add_argument('-t', '--test_col', help='Name of the test allele column ', default='minor_allele')
    
    results = parser.parse_args(args)
    return (results.input   , results.output , results.beta , results.se , results.variant_col , results.cytoband, results.chrom_col , results.position_col , results.ref_col , results.alt_col , results.test_col)


def store_hd5f(df,name):
    for i in df.cytoband.unique() :
        temp = df[df["cytoband"] == i]
        if not temp.empty : 
            key_name = str(i)
            temp.to_hdf(name, key_name, format='t', complevel=4)
        del temp


def main(input_file  , output , beta ,  se  , variant_col,  cyto, chrom , pos , ref , alt , test):
    cyto_df = pd.read_csv(cyto, sep="\t")
    gr = pr.PyRanges(cyto_df)
    for gm_chunk in pd.read_csv(input_file, sep="\t" , chunksize=2000000, low_memory=False):
        temp = gm_chunk[[variant_col , beta , se , chrom , pos , ref , alt , test ]].copy()
        temp.columns = ["SNP", "beta", "se", "Chromosome" , "Start", "ref" , "alt" , "test_allele"]
        temp["End"] = temp.Start + 1 
        W = 0.2 
        temp["Z"] = temp["beta"] / temp["se"]
        temp["V"] = temp["se"] ** 2 
        temp["R"] =  W **2 / (W**2 + temp["V"])
        temp["lbf"] =  0.5 * (np.log(1 - temp["R"]) + (temp["R"] * temp["Z"]**2))
        temp_again = temp[["Chromosome", "Start", "End", "SNP", "beta", "se", "Z" , "lbf", "ref" , "alt" , "test_allele"]].copy()
        gr2 = pr.PyRanges(temp_again)
        results = gr2.join(gr)
        results = results.df
        results = results[["Chromosome" , "Start" , "SNP" , "beta" , "se", "cytoband", "Z" , "lbf", "ref" , "alt" , "test_allele"]].copy()
        store_hd5f(results,output)
         
if __name__ == '__main__':

    input_file  , output , beta ,  se  , variant_col,  cyto , chrom_col , position_col , ref_col , alt_col , test_col = check_arg(sys.argv[1:])

    main(input_file  , output , beta ,  se  , variant_col,  cyto , chrom_col , position_col , ref_col , alt_col , test_col)