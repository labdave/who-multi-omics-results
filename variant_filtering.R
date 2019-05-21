


###Preprocess variants - from recoded vcf

  read.table(variants.file,sep="\t",header=T,check.names=F,na.strings=c(".","","NA"),comment.char="",quote="",as.is=T) -> variants.all
  rownames(variants.all) <- paste0(variants.all$CHROM,"-",variants.all$POS,"-",variants.all$REF,"-",variants.all$ALT)
  
  variant.annot <- variants.all[,c(1:138)]
  variants <- as.matrix(variants.all[,c(139),drop=F]); colnames(variants) <- tolower(colnames(variants))
  variants.depth <- as.matrix(variants.all[,140]); colnames(variants.depth) <- tolower(colnames(variants.depth))
  
  pop.cols <- grep("gnomAD|ExAC|popfreq_max",colnames(variant.annot),value=T)
  variant.annot[,pop.cols] -> pop.freq.data
  pop.freq.data[is.na(pop.freq.data)] <- 0
  pop.freq.max.all <- apply(pop.freq.data,1,max)
  
  cbind(variant.annot,pop.freq.max.all) -> variant.annot
  
###############
  
  
  
#### Convert genotypes to 0=WT, 1=MUT and NA
  variants[variants>=recode.mut.level] <- 1
  variants[variants<recode.mut.level & variants>recode.wt.level] <- NA
  variants[variants<=recode.wt.level] <- 0
  
  filter.variants <- variants==1 & !is.na(variants)
  variants <- variants[filter.variants]
  variant.annot <- variant.annot[filter.variants,]

##############
  
  
  
####Filter variants
  
  quality.variants.moderate <- !((as.numeric(variant.annot$QD)<2 & !is.na(variant.annot$QD)) | 
                                   (as.numeric(variant.annot$FS)>100 & !is.na(variant.annot$FS)) | 
                                   (as.numeric(variant.annot$SOR)>3 & !is.na(variant.annot$SOR)) | 
                                   (as.numeric(variant.annot$MQ)<40 & !is.na(variant.annot$MQ)) |
                                   (as.numeric(variant.annot$MQRankSum)<(-2.5) & !is.na(variant.annot$MQRankSum)) | 
                                   (as.numeric(variant.annot$ReadPosRankSum)<(-8) & !is.na(variant.annot$ReadPosRankSum)))
  
  rare.variants <-  (as.numeric(variant.annot$pop.freq.max.all)<=0.001)
  not.syn <- variant.annot[,'ExonicFunc.refGene']!="synonymous_SNV" & !is.na(variant.annot[,'ExonicFunc.refGene'])
  remove.fs <- !grepl("frameshift",variant.annot[,'ExonicFunc.refGene']) & !is.na(variant.annot[,'ExonicFunc.refGene'])
  
  not.superdups <- is.na(variant.annot$genomicSuperDups)
  not.dbsnp <- is.na(variant.annot$avsnp150)
  
  filter.variants <- rare.variants & not.syn & not.superdups & quality.variants.moderate & remove.fs

  variant.annot.filtered <- variant.annot[filter.variants,]
  variants.filtered <- variants[filter.variants]

  write.table(cbind(variant.annot.filtered,variants.filtered), paste0(output.prefix,"_filtered_variants.txt"), sep="\t", col.names=NA, row.names=T)
################
  
  
####Convert to json
  
  variant.annot.filtered.sel <- variant.annot.filtered[, c("CHROM","POS","REF","ALT","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","pop.freq.max.all","CADD_phred","cosmic87_coding")]
  colnames(variant.annot.filtered.sel) <- c("CHROM","POS","REF","ALT","Gene","Mutation Type","Protein Change","Population frequency","CADD_phred","cosmic87_coding")
  
  write.table(variant.annot.filtered.sel, paste0(output.prefix,"_filtered_variants_subset.txt"), sep="\t", col.names=NA, row.names=T)
  
  write_json(variant.annot.filtered.sel, paste0(output.prefix,"_filtered_variants.json"))
  
  write.table(cbind(variant.annot,variants,filter.variants,rare.variants,not.syn,not.superdups,quality.variants.moderate),paste0(output.prefix,"_all_variants_with_filters.txt"),sep="\t",col.names=NA,row.names=T)
  
#########################
  
  
  
  