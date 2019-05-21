

##############Read reference genomic risk model
load(ref.genomic.risk.model.file)

###################


##############Format TNA data into input features

all.input.features <- unique(unlist(strsplit(rownames(genomic.risk.model)," & ")))
tna.input.data <- matrix(data=0,nrow=length(all.input.features),ncol=1,dimnames=list(all.input.features,sample.id))

common.genes <- intersect(unique(variant.annot.filtered$Gene.refGene),all.input.features)

tna.input.data[common.genes,] <- 1

if (dlbcl.subtype=="ABC DLBCL"){
  tna.input.data["ABC.775",] <- 1
} else if (dlbcl.subtype=="GCB DLBCL") {
  tna.input.data["GCB.775",] <- 1
}

if (tna.expr.myc.bcl2[tna.expr.myc.bcl2$gene=="MYC","fpkm.log.z"]>=myc.cutoff){
  tna.input.data["MYC.high.775",] <- 1
} 

if (tna.expr.myc.bcl2[tna.expr.myc.bcl2$gene=="BCL2","fpkm.log.z"]>=bcl2.cutoff){
  tna.input.data["BCL2.high.775",] <- 1
} 

write.table(tna.input.data,paste0(output.prefix,"_input_features.txt"),sep="\t",col.names=NA,row.names=T)

##############################



###########TNA combo data
tna.combo.data <- matrix(data=0,nrow=nrow(genomic.risk.model),ncol=1,dimnames=list(rownames(genomic.risk.model),sample.id))

for (sel.combo.feature in rownames(genomic.risk.model)) {
  
  ind.features <- unlist(strsplit(sel.combo.feature," & "))
  tna.combo.data[sel.combo.feature,] <- prod(tna.input.data[ind.features,])  
}

write.table(tna.combo.data,paste0(output.prefix,"_input_combos.txt"),sep="\t",col.names=NA,row.names=T)

######################



#########TNA risk score

tna.risk.score <- sum(tna.combo.data*genomic.risk.model)

closest.samples <- names(sort(abs(pred.ref[,1]-tna.risk.score))[1:50])
risk.df <- rbind(data.frame("Risk.score"=pred.ref[,1],"Samples"="All"),data.frame("Risk.score"=pred.ref[closest.samples,1],"Samples"="Similar_risk"))

pdf(paste0(output.prefix,"_distribution_hist_similar_patients.pdf"))
  p <- ggplot(risk.df, aes(x=Risk.score,fill=Samples)) + geom_density(alpha=0.4, adjust=1) + scale_fill_manual(values=c("darkgrey","red")) + geom_vline(xintercept = tna.risk.score, col="red", size=1) + ggtitle(paste0(sample.id,": risk score=", round(tna.risk.score,2)))
  print(p)
dev.off()

###################




###############Survival analysis - function
survival_plot=function(surv.obj,plot.feature,name.plot.feature,col=rev(seq_along(levels(plot.feature))),plot.p=T) {
  par(mar=c(5,5,2,4))
  #plot.feature <- factor(plot.feature)
  #col=rev(seq_along(levels(plot.feature)))
  plot(survfit(surv.obj~plot.feature),col=col,lwd=5,main=name.plot.feature,xlab="Overall Survival (years)",ylab="Survival probability",las=1,cex.axis=1.25,cex.lab=1.5,cex.main=2)
  #title(cex=1.5)
  sdf <- survdiff(surv.obj~plot.feature); 
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  if (plot.p==T) {
    if (p.val<0.0001) {
      legend("topright",expression(P<10^-4),cex=1.75,bty="n");
    } else {
      legend("topright",paste0("P=",prettyNum(p.val,digits=3)),cex=1.75,bty="n");
    }
  }		
  legend("bottomleft", legend=paste0(levels(factor(plot.feature))," (N=",table(plot.feature),")"),col=col,pch="-",lwd=2,cex=1.25)
}
###########


#########Survival plot 
surv.obj <- Surv(ref.input.survival[,'Overall_Survival'],1-ref.input.survival[,'Censored'])

plot.feature <- rep(0,nrow(ref.input.survival)); names(plot.feature) <- rownames(ref.input.survival)
plot.feature[closest.samples] <- 1

pdf(paste0(output.prefix,"_survival_similar_risk.pdf"))
  survival_plot(surv.obj,plot.feature,paste0(sample.id,": genomic risk score = ",round(tna.risk.score,2)),col=c("grey","red"))
dev.off()  

################



