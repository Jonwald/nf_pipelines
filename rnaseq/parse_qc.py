## parse almac RNAseq and get summary files qc files and generate a summary table
import os
import gzip
import pandas as pd
import datetime
import csv
import sys
from argparse import ArgumentParser

def parse_command(command):
    """parses inputs"""
    parser = ArgumentParser(description='parse almac RNAseq and get summary files qc files and generate a summary table')
    parser.add_argument('-f', '--folder', help='folder with all qc files', required=True)
    args = parser.parse_args()
    return args

def parse_bamqc(bamqc):
    with open(bamqc) as f:
        d = {}
        for line in f:
            if '|' not in line:
                continue
            else:
                line = line.strip(" ").strip("\n").split("|\t")
                d[line[0]] = line[1]
    return d

def parse_rrna(rrna):
    with open(rrna) as f:
        rrna_count = f.readlines()[2].split(" ")[0]
    return int(rrna_count)
        
def parse_fastqc(fastqc1, fastqc2):
    with open(fastqc1) as f:
        for line in f:
            if line.startswith('%GC'):
                gc1 = line.strip("%GC\t").strip("\n")
            else:
                continue
    with open(fastqc2) as f:
        for line in f:
            if line.startswith('%GC'):
                gc2 = line.strip("%GC\t").strip("\n")
            else:
                continue
    return float(gc1), float(gc2)

def parse_duplication(duplication):
    with open(duplication) as f:
        metrics = f.readlines()[7].split("\t")
        dups = round(float(metrics[8]) * 100, 2)
        opt_dups = round(float(metrics[7]) / float(metrics[2]) * 100, 2)
    return dups, opt_dups
            

def main():
    input_command = parse_command(sys.argv[1:])
    qc_fol = input_command.folder
    ## get timestamp
    timestamp=datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    ## define output file header
    header  = ['File',
               'Number_of_input_read_pairs',
               'Uniquely_mapped_read_pairs_number',
               'Uniquely_mapped_reads_perc',
               'Number_of_read_pairs_mapped_to_too_many_loci',
               'Perc_of_reads_mapped_to_too_many_loci',
               'Mapping_rate',
               'Perc_GC_Read1',
               'Perc_GC_Read2',
               'Ribosomal_RNA_contamination',
               'Duplication_rate',
               'Optical_Duplication_rate',
               'HK_coverage']
    
    ## list files and get sample names
    
    all_files = os.listdir(qc_fol)
    snames = set([x.split("_Log")[0] for x in all_files if x.endswith("_Log.final.out")])
    
    ## read housekeeping genes
    with open(qc_fol + "hk_gene_ids.txt") as f:
            hk_genes = [x.split("\t")[0] for x in f.readlines() if x.startswith("ENSG")]
    
    ## read cov file and get median expression of HK genes
    #tpm_file = "ADX20114_B1-5_gene_RSEM_TPM.tsv"
    cov_file = [x for x in all_files if x.endswith("COV.tsv")][0]
    df =pd.read_csv(qc_fol + cov_file, sep="\t")
    df = df[df["GENE"].isin(hk_genes)]
    med_exp = df.iloc[:, 2:].median(axis=0)
    
    ## parse QC files and compile into list of lists
    outlines = []
    for sample in list(snames):
        bamqc = qc_fol + sample + "_Log.final.out"
        rrna = qc_fol + sample + "_rrna_stats.txt"
        fastqc1 = qc_fol + sample + "_fastqc" + "/fastqc_data.txt"
        fastqc2 = qc_fol + sample.replace("_R1_", "_R2_") + "/fastqc_data.txt"
        duplication = qc_fol + sample + "_marked_dup_metrics.txt"

        bamqc_dict = parse_bamqc(bamqc)
        rrna_count = parse_rrna(rrna)
        gc1, gc2 = parse_fastqc(fastqc1, fastqc2)
        dups, opt_dups = parse_duplication(duplication)
        hk_cov = med_exp.loc[sample] 
        
        outline = [sample, bamqc_dict['Number of input reads '],
        bamqc_dict['Uniquely mapped reads number '],
        bamqc_dict['Uniquely mapped reads % '].strip("%"),
        bamqc_dict['Number of reads mapped to too many loci '],
        bamqc_dict['% of reads mapped to too many loci '].strip("%"),
        round(float(bamqc_dict['Uniquely mapped reads % '].strip("%")) + float(bamqc_dict['% of reads mapped to too many loci '].strip("%")),2),
        gc1, gc2, rrna_count / int(bamqc_dict['Number of input reads ']) , dups, opt_dups, hk_cov]
        
        outlines.append(outline)
    outname = qc_fol + timestamp + "_qc_metrics.txt"
    ## write output file
    with open(outname, "w") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header)
        for i in outlines:
            writer.writerow(i)

main()