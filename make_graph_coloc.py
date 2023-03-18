import numpy as np
import math
import sumstats
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
import sys
import argparse
import os 
import plotly.graph_objects as go
import plotly

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Make Z-Z plots from GWAS summary statistics and eqtls file ' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--chrom', help='Current chromsome', required='True')
    parser.add_argument('-e', '--ensembl', help='Ensemble Gene ID',required='True')
    parser.add_argument('-r', '--rsid', help='SNP rsid',required='True')
    parser.add_argument('-p', '--plink', help='full path for plink file')
    parser.add_argument('-t', '--tissue', help='full path to the gtex file (allpairs)',required='True')
    parser.add_argument('-n', '--tissue_name', help='Tissue name for graph and file output',required='True')
    parser.add_argument('-ref', '--reference', help='the path of the directory to Map file for gtex make sure it ends with /',required='True')
    parser.add_argument('-d', '--dataset', help='GWAS summary stats file' , required='True')
    parser.add_argument('-g', '--gene', help='gene name for graph',required='True')
    parser.add_argument('-rsid_col', '--rsid_col', help='RSID column from Summary statistics',default="SNP")
    parser.add_argument('-ref_col', '--ref_col', help='reference allele column from Summary statistics',default= "REF")
    parser.add_argument('-alt_col', '--alt_col', help='alternative allele column from Summary statistics',default='ALT')
    parser.add_argument('-test_al_col', '--test_al_col', help='Test allele column if Summary statistics has both ref,alt and test - leave blank if not require',default=0)
    parser.add_argument('-beta', '--beta_col', help='Beta column from Summary statistics', default="BETA")
    parser.add_argument('-OR', '--odds_ratio', help='If Odds ratios are given instead of beta in the Summary statistics file',default=0)
    parser.add_argument('-SE', '--SE_col', help='Standard error from the Summary statistics',default="SE")
    results = parser.parse_args(args)
    return (results.chrom, results.ensembl , results.rsid , results.plink , results.tissue , results.tissue_name ,results.reference, results.dataset, results.gene , results.rsid_col ,results.ref_col , results.alt_col, results.test_al_col, results.beta_col, results.odds_ratio , results.SE_col   )
def make_bayes_factors(df):
    """
    returns a df with just the SNP and the bayes factor for each SNP
    INPUT DF contains SNP, beta, se
    """ 
    W = 0.2 
    df["Z"] = df["beta"] / df["se"]
    df["V"] = df["se"] ** 2 
    df["R"] =  W **2 / (W**2 + df["V"])
    df["lbf"] =  0.5 * (np.log(1 - df["R"]) + (df["R"] * df["Z"]**2))
    df_small = df[["SNP", "lbf" , "Z"]].copy()
    return df_small 
def get_PP4(
    title: str,
    pp0: float,
    pp1: float,
    pp2: float,
    pp3: float,
    pp4: float):
    """part of calling the below function"""
    return pp4
def coloc1(trait1_lnbfs, trait2_lnbfs, prior1: float = 1e-4, prior2: float = 1e-4, prior12: float = 1e-5 ):
    """
    COLOC function taken(stolen) from Anthony Aylward
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
def get_PP(df, snp):
    '''
    This is to do the PP for the coloc
    '''
    PP4 = get_PP4("snp", *coloc1((df.GWAS_1_lbf.to_numpy()),(df.GWAS_2_lbf.to_numpy() )))
    PP4_rounded = round(PP4,4)
    raw = [[ snp , PP4_rounded]]
    df2 = pd.DataFrame(data=raw, columns=["SNP", "PP4"])
    return df2
def combine_df(df1 , df2):
    """
    merged the 2 smaller DFs together, just incase SNPs are missing form one, it will push out the position each snp in the np array
    """
    df1.columns = ["SNP", "GWAS_1_lbf", "GWAS_1_Z"]
    df2.columns = ["SNP", "GWAS_2_lbf", "GWAS_2_Z"]
    merged = df1.merge(df2,on="SNP")
    return merged
def run_plink(DF, RSID, PLINK, tissue_name): 
    lung_SNP_list = DF[["rsid"]].copy()
    lung_SNP_list.to_csv(RSID + "_" + tissue_name + "_extractrs.txt", sep="\t",index=None)
    lung_plink1 = "plink --memory 32000 --bfile " + PLINK + " --extract " + RSID + "_" + tissue_name + "_extractrs.txt --make-bed --out " + RSID +"_" + tissue_name
    lung_plink2 = "plink --memory 32000 --bfile "+ RSID +"_" + tissue_name + " --ld-snp "+ RSID + " --ld-window 10000 --ld-window-kb 10000 --ld-window-r2 0 --out " +RSID+" --r2 "
    os.system(lung_plink1)
    os.system(lung_plink2)
    temp = add_rs(DF,RSID)
    return temp
def add_rs(df, RSID):
    ld_file = RSID + ".ld"
    current_snp = pd.read_csv(ld_file, delimiter=r"\s+")
    current_snp = current_snp[["SNP_B", "R2"]].copy()
    current_snp.columns=["rsid", "scale"]
    temp = df.merge(current_snp , on="rsid")
    return temp
def process_gtex_file(file): 
    temp = pd.read_csv(file, sep="\t")
    #gene_id	variant_id	tss_distance	ma_samples	ma_count	maf	pval_nominal	slope	slope_se	chr	variant_pos	ref	alt	num_alt_per_site	rs_id_dbSNP151_GRCh38p7
    temp.columns = ['gene_id', 'variant_id', 'tss_distance', 'ma_samples', 'ma_count','maf', 'pval_nominal', 'slope', 'slope_se', "chr",	"variant_pos"	,"ref",	"gtex_alt"	,"num_alt_per_site","rsid" ]
    temp["gtex_Z"] = round((temp.slope / temp.slope_se),4) 
    temp = temp.drop(['tss_distance', 'variant_id', 'ma_samples', 'ma_count','maf', 'pval_nominal' ,"ref","chr",	"variant_pos" ,"num_alt_per_site" ],1)
    return temp
def automate(ENSEMBLE, CHR, RSID, GENE_NAME , PLINK, TISSUE, TISSUE_NAME, PATH_REF, DATASET) : 
    '''
    This is a hackme function as I have not written a function to process the 20GB GTEx files
    ENSEMBLE is the ID for the gene of interest (not the isoform - I guess it would work)
    GTEX_FILE_NAME is the *.all_pairs.txt.gz for a given tissue of interest
    GENE_NAME and TISSUE_NAME are used to name the grep output file
    '''
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
def first_merge(data_df ,gtex_df, rsid_col ,ref_col , alt_col, test_al_col, beta_col, odds_ratio , SE_col):
    temp2 = gtex_df.merge(data_df, left_on="rsid", right_on=rsid_col)
    ##convert to beta if OR is set
    if odds_ratio != 0 :
        beta_col = np.log(odds_ratio)
    ##to deal with new plinks output where you have ref alt and test allele
    if test_al_col != 0: 
        if test_al_col == ref_col : 
            ref_col = alt_col
            alt_col = test_al_col
    smaller1 = temp2[["rsid", "gtex_alt", "gtex_Z" , ref_col, alt_col, beta_col, SE_col, 'slope', 'slope_se']].copy()
    smaller1.columns = ["rsid", "gtex_alt", "gtex_Z" , "ref", "alt", "beta", "se", 'slope', 'slope_se']
    smaller1["Z_score"] = round((smaller1.beta / smaller1.se),4) 
    return smaller1
def flip(df): 
    needs_flipping = df[df["gtex_alt"] != df["alt"]]
    cob_clean = df[~df['rsid'].isin(needs_flipping["rsid"])]
    needs_flipping["Z_score"] = needs_flipping.Z_score * -1 
    cob_clean["Z_score"] = cob_clean.Z_score
    combined = pd.concat([cob_clean, needs_flipping])
    return combined
def SetColor(x):
    if(float(x) < 0.01):
        return "#E0E0E0"
    elif(float(x) >= 0.01) & (float(x) <= 0.4):
        return "green"
    elif(float(x) >= 0.4) & (float(x) <= 0.8):
        return "yellow"
    elif(float(x) > 0.8) &  (float(x) < 0.9999):
        return "rgb(227,26,28)"
    elif(float(x) > 0.9999):
        return "purple"
def make_graph(df, tissue , gene, rsid): 
    final1 = df.sort_values("scale",ascending=True)
    final1['marker_size'] = (final1['scale'] + 1 ) * 9
    combined2 = final1[final1["scale"] != 1 ]
    combined3 = final1[final1["scale"] == 1 ]
    Y_name = "<b> "+ gene + " " + tissue + " expression t-score </b>"
    title_plot = rsid + " in " + gene
    fig = go.Figure()
    #fig.update_layout(yaxis=dict(range=[-6,6]))
    #fig.update_layout(xaxis=dict(range=[-6,6]))
    fig.update_layout(plot_bgcolor='white' ) 
    fig.update_layout(title={'text': title_plot,'x':0.5 , 'xanchor': 'center','yanchor': 'top' , 'font_color':"black"})
    fig.update_xaxes(ticks="outside", title="<b> GWAS trait </b>"  ,zeroline = True, zerolinewidth = 3, zerolinecolor="black"  ,showgrid = False,  gridwidth = 1,gridcolor = '#bdbdbd', color='black', tickwidth=5, tickcolor='black', ticklen=10, showline=True, linewidth=5, linecolor='Black')
    fig.update_yaxes(ticks="outside", title=Y_name , zeroline = True, zerolinewidth = 3,zerolinecolor="black", showgrid = False, gridwidth = 1,gridcolor = '#bdbdbd',color='black',tickwidth=5, tickcolor='black', ticklen=10, showline=True, linewidth=5, linecolor='Black')
    fig.update_layout(font=dict(size=20) ) 
    fig.add_trace(go.Scatter( name="R-square", mode="markers" , text=final1["rsid"] , x = combined2["Z_score"].tolist(), y = combined2["gtex_Z"].tolist(), 
    showlegend=False , marker_size=combined2['marker_size'],opacity=1, marker_line_width=0.8,line=dict(color="black"), marker_color=list(map(SetColor, combined2['scale']))))
    fig.add_trace(go.Scatter( name="SNP", text=combined3["rsid"], mode="markers" ,marker_symbol="asterisk" , x = combined3["Z_score"].tolist(), y = combined3["gtex_Z"].tolist(),
    showlegend=False , marker_size=20,opacity=1, marker_line_width=3,line=dict(color="black"), marker_color=list(map(SetColor, combined3['scale']))))   
    output_name = rsid + "_" + gene + "_" + tissue + ".html"
    #fig.write_image(output_name,width=1000, height=800, scale=10)
    plotly.offline.plot(fig, filename = output_name, auto_open=False)
def make_df_coloc(df): 
    gwas_file = df[["rsid", "beta", "se"]].copy()
    gwas_file.columns=["SNP", "beta", "se"]
    gtex_file = df[["rsid", "slope" , "slope_se"]].copy()
    gtex_file.columns=["SNP", "beta", "se"]
    return gwas_file , gtex_file

def main(chrom, ensembl , rsid , plink , tissue , tissue_name ,reference, dataset, gene, rsid_col ,ref_col , alt_col, test_al_col, beta_col, odds_ratio , SE_col):
    print ("reading in summary stats")
    data = pd.read_csv(dataset , sep="\t") 
    print ("Done")
    print ("extracting gene from GTEX file")
    automate(ensembl, chrom, rsid, gene , plink, tissue, tissue_name, reference, dataset)
    print ("Done")
    print ("Processing and merging datasets")
    file_name_gtex = rsid + ".find."+gene+ "_" + tissue_name + ".txt"
    gene_gtex_file = process_gtex_file(file_name_gtex)
    clean = first_merge(data,gene_gtex_file,rsid_col ,ref_col , alt_col, test_al_col, beta_col, odds_ratio , SE_col)
    print (clean.head(2))
    print ("Done")
    print ("running plink to get r2 for snps")
    clean_temp = run_plink(clean, rsid  , plink, tissue_name)
    print ("Done")
    print ("flipping alleles")
    clean1 = flip(clean_temp)
    print ("Done")
    print ("plotting")
    make_graph(clean1, tissue_name , gene , rsid )
    print ("Done")
    print ("doing coloc for PP4")
    gwas_file , gtex_file = make_df_coloc(clean1)
    moo = make_bayes_factors(gwas_file)
    moo2 = make_bayes_factors(gtex_file)
    combined = combine_df(moo,moo2)
    combined = combined.dropna()
    
    print (get_PP(combined, rsid))
    print ("Done")
if __name__ == '__main__':
    print ( """ZxZ plots from GTEX data and coloc """ )
    chrom, ensembl , rsid , plink , tissue , tissue_name ,reference, dataset, gene,  rsid_col ,ref_col , alt_col, test_al_col, beta_col, odds_ratio , SE_col  = check_arg(sys.argv[1:])
    main(chrom, ensembl , rsid , plink , tissue , tissue_name ,reference, dataset, gene , rsid_col ,ref_col , alt_col, test_al_col, beta_col, odds_ratio , SE_col)