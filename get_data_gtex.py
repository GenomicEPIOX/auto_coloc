import os 

import pandas as pd 
import numpy as np
import sys
import argparse

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='extract gtex eqtl information' ,formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-c', '--chrom', help='Current chromsome', required='True')
    parser.add_argument('-e', '--ensembl', help='Ensemble Gene ID',required='True')
    parser.add_argument('-r', '--rsid', help='SNP rsid',required='True')
    parser.add_argument('-p', '--plink', help='full path for plink file')
    parser.add_argument('-t', '--tissue', help='full path to the gtex file (allpairs)',required='True')
    parser.add_argument('-n', '--tissue_name', help='Tissue name for graph and file output',required='True')
    parser.add_argument('-ref', '--reference', help='the path of the directory to Map file for gtex make sure it ends with /',required='True')
    parser.add_argument('-d', '--dataset', help='AGs master dataset')
    parser.add_argument('-g', '--gene', help='gene name for graph',required='True')
    
    results = parser.parse_args(args)
    return (results.chrom, results.ensembl , results.rsid , results.plink , results.tissue , results.tissue_name ,results.reference, results.dataset, results.gene )




##[[0, "#E0E0E0"], [0.1, "#9E9E9E"], [0.6, "yellow"],[1, "rgb(227,26,28)"]]

def SetColor(x):
    if(float(x) < 0.1):
        return "#E0E0E0"
    elif(float(x) >= 0.1) & (float(x) <= 0.4):
        return "#9E9E9E"
    elif(float(x) >= 0.4) & (float(x) <= 0.8):
        return "yellow"
    elif(float(x) > 0.8):
        return "rgb(227,26,28)"

#combined_allDatasets_noDup_16062020.txt

def automate(ENSEMBLE, CHR, RSID, GENE_NAME , PLINK, TISSUE, TISSUE_NAME, PATH_REF, DATASET) : 
    file_name_lung = RSID + ".find."+GENE_NAME+ "_" + TISSUE_NAME + ".txt"

    command1 = "zgrep " + str(ENSEMBLE)  + " " + TISSUE + " > " + file_name_lung
    os.system(command1)
    header = ['gene_id', 'variant_id', 'tss_distance', 'ma_samples', 'ma_count','maf', 'pval_nominal', 'slope', 'slope_se']

    lung_find = pd.read_csv(file_name_lung, sep="\t", header=None)
    lung_find.columns = header
    file_name = PATH_REF + "ref_GTEX_hg38_chr" + str(CHR) + ".txt"
    chrom_ref = pd.read_csv(file_name, sep="\t")
    lung_find = lung_find.merge(chrom_ref,on="variant_id") 
    lung_find.to_csv(file_name_lung, sep="\t", index=None)

#     print ("zgrep complete")
#     lung_list = pd.DataFrame()
#     for gm_chunk in pd.read_csv(DATASET, sep="\t" , chunksize=100000, low_memory=False):
#         new = gm_chunk[["UKFH_SNP", "META_ZLC", "META_ref_new", "META_alt_new", "META_swap", "META_strand" ]]
#         lung_temp = lung_find.merge(new, left_on="rs_id_dbSNP151_GRCh38p7", right_on="UKFH_SNP")
#         lung_list= lung_list.append(lung_temp)

#     print ("merge with master GWAS file complete")
#     lung_list.rename(columns={'rs_id_dbSNP151_GRCh38p7':'SNP'}, inplace=True)
#     lung_Good_snps = lung_list[lung_list["alt"]==lung_list["META_alt_new"]]
#     lung_needs_flipping = lung_list[lung_list["ref"]==lung_list["META_alt_new"]]
#     lung_clean1 = (lung_list[~lung_list["SNP"].isin(lung_Good_snps.SNP)])
#     lung_clean2 = (lung_clean1[~lung_clean1["SNP"].isin(lung_needs_flipping.SNP)])
#     lung_clean2[['new_ref', "new_effect"]] =  lung_clean2.apply(lambda row : pd.Series(fix_strand(row)), axis=1)
#     lung_Good_snps2 = lung_clean2[lung_clean2["alt"]==lung_clean2["new_effect"]]
#     lung_needs_flipping2 = lung_clean2[lung_clean2["ref"]==lung_clean2["new_effect"]]
#     lung_Good_snps2 = lung_Good_snps2.drop(["new_ref", "new_effect"], axis=1)
#     lung_Good = pd.concat([lung_Good_snps,lung_Good_snps2])
#     lung_Good["meta_Z"] = lung_Good.META_ZLC
#     lung_needs_flipping2 = lung_needs_flipping2.drop(["new_ref", "new_effect"], axis=1)
#     lung_flip = pd.concat([lung_needs_flipping,lung_needs_flipping2])
#     lung_flip["meta_Z"] = lung_flip.META_ZLC * -1 
#     lung_complete = pd.concat([lung_Good,lung_flip])
#     name_t = TISSUE_NAME + "_T"
#     lung_complete[name_t] = round(lung_complete.slope / lung_complete.slope_se, 3 )
#     lung_small_comp = lung_complete[["SNP", "meta_Z", name_t]].copy()
#     lung_SNP_list = lung_complete[["SNP"]].copy()
#     d = {'SNP': [RSID] } 
#     temp_df = pd.DataFrame(data=d)
#     lung_SNP_list = lung_SNP_list.append(temp_df, ignore_index=True)
#     lung_SNP_list.to_csv("for_plink_r2.txt", sep="\t",index=None)
#     print ("flipping complete")
#     lung_plink1 = "plink --bfile "+ PLINK + " --extract for_plink_r2.txt --make-bed --out temp"
#     lung_plink2 = "plink --bfile temp --ld-snp "+ RSID + " --ld-window 10000 --ld-window-kb 10000 --ld-window-r2 0 --out "+RSID+" --r2 " 
#     os.system(lung_plink1)
#     os.system(lung_plink2)


#     ld_file = RSID + ".ld"
#     current_snp = pd.read_csv(ld_file, delimiter=r"\s+")
#     current_snp = current_snp[["SNP_B", "R2"]].copy()
#     snp_col = RSID + "_R2"
#     current_snp.columns=["SNP", snp_col]
#     print ("Plink R2 complete")
#     lung_current_snp = current_snp.merge(lung_small_comp, on="SNP") 
#     lung_current_snp.sort_values([snp_col,"meta_Z",name_t] , ascending=False)
#     lung_current_snp["scale"] = lung_current_snp[snp_col].round(4).astype(float)
#     lung_current_snp["reference"] = "SNP: " + lung_current_snp.SNP + " <br> R2 with "+RSID+" : " + lung_current_snp.scale.astype(str)
#     lung_current_snp = lung_current_snp.sort_values("scale", ascending=True)
#     #lung_current_snp = lung_current_snp[lung_current_snp["SNP"] !=RSID]

    
#     ## plotting lung
    
#     fig = go.Figure()
#     fig.update_layout(plot_bgcolor='white' ) 
#     title_plot = RSID + " - " + GENE_NAME 
#     fig.update_layout(
#         title={
#             'text': title_plot,
#             'y':0.9,
#             'x':0.5,
#             'xanchor': 'center',
#             'yanchor': 'top' , 'font_color':"black"})
#     fig.add_trace(go.Scatter( name="R-square", mode="markers"  , x = lung_current_snp["meta_Z"].tolist(), y = lung_current_snp[name_t].tolist(), text=lung_current_snp["reference"].tolist(),
#     showlegend=False , marker_size=5, marker_color=list(map(SetColor, lung_current_snp['scale']))) )  
#     #fig.update_layout(yaxis=dict(range=[-6,6]))
#     #fig.update_layout(xaxis=dict(range=[-6,6]))
#     y_name = TISSUE_NAME + "GTEX t-stat"
#     fig.update_xaxes(ticks="outside", title="FH-meta Z-score" , color='black', tickwidth=1, tickcolor='black', ticklen=5, showline=True, linewidth=2, linecolor='Black')
#     fig.update_yaxes(ticks="outside", title=y_name , color='black',tickwidth=1, tickcolor='black', ticklen=5, showline=True, linewidth=2, linecolor='Black')
#     lung_png = GENE_NAME + TISSUE_NAME +  ".png"
#     lung_svg = GENE_NAME + TISSUE_NAME +  ".svg"
#     lung_html = GENE_NAME + TISSUE_NAME + ".html"
#     fig.write_image(lung_png,width=956, height=770, scale=10)
#     fig.write_image(lung_svg,width=956, height=770, scale=10)
#     plotly.offline.plot(fig, filename = lung_html, auto_open=False)
#     #
#     print ("plot complete")

    
    
def fix_strand(row):
    ''' clean data
    1. IF strand flip, change each allele 
    2. IF allele flip, flip alleles and times Z by -1 
    '''
    base_T = "A"
    base_A = "T"
    base_C = "G"
    base_G = "C"
    
    new_effect = ""
    new_ref = ""
    
    if row["META_strand"] == True:
        if row["META_ref_new"] == "T":
            new_ref = base_T
        elif row["META_ref_new"] == "A":
            new_ref = base_A
        elif row["META_ref_new"] == "C":
            new_ref = base_C    
        elif row["META_ref_new"] == "G":
            new_ref = base_G
            
        if row["META_alt_new"] == "T":
            new_effect = base_T  
        elif row["META_alt_new"] == "A":
            new_effect = base_A  
        elif row["META_alt_new"] == "C":
            new_effect = base_C  
        elif row["META_alt_new"] == "G":
            new_effect = base_G  
    
    if row["META_strand"] == False:
        new_ref = row["META_ref_new"]
        new_effect = row["META_alt_new"]
    
    return new_ref, new_effect



def main(chrom, ensembl , rsid , plink , tissue , tissue_name ,reference, dataset , gene ):
    automate(ensembl, chrom, rsid, gene, plink, tissue, tissue_name, reference, dataset)



if __name__ == '__main__':

    chrom, ensembl , rsid , plink , tissue , tissue_name ,reference, dataset, gene  = check_arg(sys.argv[1:])
    main(chrom, ensembl , rsid , plink , tissue , tissue_name ,reference, dataset, gene )

