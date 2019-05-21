
###############Load reference data set
load(ref.expr.rdata.file)



###############Input parameters and read file
tna.expr <- read.table(tna.expr.file,sep="\t",header=T,row.names=1,as.is=T)



#################Select same genes as ref data (gene panel)
tna.expr <- tna.expr[rownames(ref.1kdbcl),]


##################Compute FPKM
library.size <- sum(tna.expr[,3])

tna.expr$fpkm <- (tna.expr[,3]/library.size)*10^6
tna.expr$fpkm <- (tna.expr$fpkm/ref.1kdbcl.annot$gene.length)*10^3

##################Log normalize
tna.expr$fpkm.log <- log2(tna.expr$fpkm+1)

##################Z-normalize
tna.expr$fpkm.log.z <- (tna.expr$fpkm.log - mean.sd.1kdlbcl[,"mean"])/mean.sd.1kdlbcl[,"sd"]

##################Select subtype genes
tna.expr[rownames(subtype.genes),] -> tna.expr.subset

#################Signature scores
abc.score <- mean(tna.expr.subset$fpkm.log.z[subtype.genes[,"subtype"]=="ABC"])
gcb.score <- mean(tna.expr.subset$fpkm.log.z[subtype.genes[,"subtype"]=="GCB"])
sig.score <- abc.score - gcb.score

cbind(tna.expr.subset,"subtype"=factor(subtype.genes[,"subtype"],levels=c("GCB","ABC")),"genes"=subtype.genes[,"gene_name"]) -> tna.expr.subset
tna.expr.subset$genes <- factor(tna.expr.subset$genes,levels=as.vector(tna.expr.subset[order(tna.expr.subset$subtype,tna.expr.subset$fpkm.log.z),]$genes))

dlbcl.subtype <- "Unclassified"
if (abc.score>0 & sig.score>0.25) {
  dlbcl.subtype="ABC DLBCL"
} else if (gcb.score>0 & sig.score<(-0.25)) {
  dlbcl.subtype="GCB DLBCL"
}


################Output files
pdf(paste0(output.prefix,"_COO_genes.pdf"),width=6,height=6)
  my.colors <- c("#7fbf7b","#af8dc3")
  p <- ggplot(tna.expr.subset,aes(genes,fpkm.log.z,fill=subtype)) + scale_fill_manual(values=my.colors) + geom_bar(stat="identity")  + coord_flip() + ggtitle(paste0(sample.id," (Num Reads=",library.size,")\n",dlbcl.subtype)) + geom_hline(yintercept = c(gcb.score,abc.score), size=1.5, linetype="dashed",color=my.colors)
  print(p)
dev.off()


write.table(tna.expr.subset[,c("gene","fpkm.log.z","subtype")], paste0(output.prefix,"_expr_subset.txt"),sep="\t", col.names=NA, row.names=T)
write_json(tna.expr.subset[,c("gene","fpkm.log.z","subtype")], paste0(output.prefix,"_expr_subset_COO.json"))
write_json(data.frame("ABC_score"=abc.score,"GCB_score"=gcb.score, "Cell_of_origin_classification"=dlbcl.subtype), paste0(output.prefix,"_COO_classification.json"))

########################

#############MYC & BCL2 high expression
tna.expr[tna.expr$gene %in% c("MYC","BCL2"),] -> tna.expr.myc.bcl2

pdf(paste0(output.prefix,"_MYC_BCL2_expr.pdf"),width=6,height=6)
  p <- ggplot(tna.expr.myc.bcl2,aes(gene,fpkm.log.z,fill=gene)) + geom_bar(stat="identity")  + coord_flip() + ggtitle(sample.id) 
  print(p)
dev.off()


