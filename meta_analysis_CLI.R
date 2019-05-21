#!/usr/bin/Rscript                                                                                                                                    
                                                                                                                                                      
# Script Name: meta_analysis_CLI.R                                                                                                          
#                                                                                                                                                     
# Author: DaveLab Bioinformatics Team                                                                                                                 
#                                                                                                                                                     
# Created on: 05/21/2019                                                                                                                              
#                                                                                                                                                     
# Description: TBD                                    
#                                                                                                                                                     
# Syntax: Rscript --vanilla meta_analysis_CLI.R <TBD>
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
  make_option(c("-r", "--refexpr"), dest = "ref_expr", type = "character", help = "a reference R object with 1kdlbcl expr data", default = F),
  make_option(c("-g", "--refcellmodel"), dest = "ref_cell_model", type = "character", help = "a reference R object with the genomic risk model from Cell", default = F),
  make_option(c("-o", "--outputdir"), dest = "out_dir", type = "character", help = "output dir name")
)

#parse the arguments
opt_parser <- OptionParser(option_list = options_list)
opt <- parse_args(opt_parser)


############## Input
sample.id <- opt$sample_id
variants.file=opt$mut_file
tna.expr.file=opt$expr_file


mut.flag=if(variants.file!=F) {TRUE} else {FALSE}
expr.flag=T
cnv.flag=F
translocation.flag=F
survival.flag=mut.flag & expr.flag

recode.mut.level=0.3
recode.wt.level=(-0.4)

ref.expr.rdata.file=opt$ref_expr
ref.genomic.risk.model.file=opt$ref_cell_model
myc.cutoff=0.5
bcl2.cutoff=0.5

##################

result.dir <- opt$out_dir
output.prefix=paste0(result.dir,"/",sample.id)


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




