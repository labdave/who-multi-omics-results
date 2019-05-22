#!/usr/bin/Rscript                                                                                                                                    
                                                                                                                                                      
# Script Name: meta_analysis_CLI.R                                                                                                          
#                                                                                                                                                     
# Author: DaveLab Bioinformatics Team                                                                                                                 
#                                                                                                                                                     
# Created on: 05/21/2019                                                                                                                              
#                                                                                                                                                     
# Description: Analysis of vcf and expr count matrices to produce results for variant filtering, expression normalization, cell of origin calls, genetic subgroups, survival analysis, translocations,  etc.                                    
#                                                                                                                                                     
# Syntax: Rscript --vanilla meta_analysis_CLI.R --exprfile <raw_count_matrix> --mutfile <recoded_vcf_file> --sampleid <sample_id> --outputdir <outputdir/sample.id>
#                                                                                                                                                     
# Usage: Rscript --vanilla meta_analysis_CLI.R <TBD>
                                                                                                                                                      
# load require library
library("survival")
library("ggplot2")
library("jsonlite")
library("optparse")


#make a list of arguments
options_list <- list(
  make_option(c("-e", "--exprfile"), dest = "expr_file", type = "character", help = "raw read count file from STAR", default = F),
  make_option(c("-m", "--mutfile"), dest = "mut_file", type = "character", help = "recoded annovar vcf", default = F),
  make_option(c("-s", "--sampleid"), dest = "sample_id", type = "character", help = "a sample name/id"),
  make_option(c("-o", "--outputdir"), dest = "out_dir", type = "character", help = "output dir name"),
  make_option("--recodemutlevel", dest = "recode_mut_level", type = "character", help = "threshold used for recoded AC for calling mutation", default = 0.3),
  make_option("--recodewtlevel", dest = "recode_wt_level", type = "character", help = "threshold used for recoded AC for calling WT", default = -0.4),
  make_option("--mycexprthreshold", dest = "myc_expr_threshold", type = "character", help = "threshold used for defining high-MYC expression", default = 0.5), 
  make_option("--bcl2exprthreshold", dest = "bcl2_expr_threshold", type = "character", help = "threshold used for defining high-BCL2 expression", default = 0.5), 
  make_option(c("-r", "--refexpr"), dest = "ref_expr", type = "character", help = "a reference R object with 1kdlbcl expr data", default = "data/1kdlbcl_counts_fpkm_norm_1MB_gene_panel.RData"),
  make_option(c("-g", "--refcellmodel"), dest = "ref_cell_model", type = "character", help = "a reference R object with the genomic risk model from Cell", default = "data/Genomic_risk_model_Ref_data.RData")
)

#parse the arguments
opt_parser <- OptionParser(option_list = options_list)
opt <- parse_args(opt_parser)


############## Input
sample.id <- opt$sample_id
variants.file=opt$mut_file
tna.expr.file=opt$expr_file

recode.mut.level=opt$recode_mut_level
recode.wt.level=(opt$recode_wt_level)

ref.expr.rdata.file=opt$ref_expr
ref.genomic.risk.model.file=opt$ref_cell_model

myc.cutoff=opt$myc_expr_threshold
bcl2.cutoff=opt$bcl2_expr_threshold

##################
output.prefix <- opt$out_dir


################## Run filtering variants script
if (variants.file!=F) {
  source("variant_filtering.R")
}



##########Run expression analysis script
if (tna.expr.file!=F) {
  source("abc_gcb_classifier.R")
}



###########Run genomic risk model (Cell paper) & survival analysis 
if (variants.file!=F & tna.expr.file!=F) {
  source("survival_analysis_Cell_paper.R")
}




