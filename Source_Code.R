setwd("~/swirl")
library(limma)
targets <- readTargets("SwirlSample.txt")
RG <- read.maimages(targets, source="spot")
RG

setwd("~/swirl")
RG$genes <- readGAL("fish.gal")
RG$printer <- getLayout(RG$genes)

imageplot(log2(RG$Rb[,1]), RG$printer, low="white", high="red")
imageplot(log2(RG$Gb[,1]), RG$printer, low="white", high="green")

MA_raw <- normalizeWithinArrays(RG, method="none")
imageplot(MA_raw$M[,1], RG$printer, zlim=c(-3,3))

plotMD(MA_raw, column = 1) ## column = 1 is the default
abline(0,0,col="blue")
plotPrintTipLoess(MA_raw, array = 1)  ## array = 1 is the default

library(ggokabeito)
boxplot(MA_raw$M[,1]~MA_raw$genes$Block,xlab="Print-Tip",ylab="M",col=rep(palette_okabe_ito(1:8),each=2),main="Raw M vs Print-Tip for swirl.1")
boxplot(MA_raw$M[,2]~MA_raw$genes$Block,xlab="Print-Tip",ylab="M",col=rep(palette_okabe_ito(1:8),each=2),main="Raw M vs Print-Tip for swirl.2")
boxplot(MA_raw$M[,3]~MA_raw$genes$Block,xlab="Print-Tip",ylab="M",col=rep(palette_okabe_ito(1:8),each=2),main="Raw M vs Print-Tip for swirl.3")
boxplot(MA_raw$M[,4]~MA_raw$genes$Block,xlab="Print-Tip",ylab="M",col=rep(palette_okabe_ito(1:8),each=2),main="Raw M vs Print-Tip for swirl.4")

MA_print_tip <- normalizeWithinArrays(RG,method="printtiploess") 
# ^^^ method="printtiploess" is the default
plotMD(MA_print_tip)
abline(0,0,col="blue")
plotPrintTipLoess(MA_print_tip)
boxplot(MA_print_tip$M[,1]~MA_print_tip$genes$Block,xlab="Print-Tip",ylab="M",col=rep(palette_okabe_ito(1:8),each=2),main="Normlized M (printtiploess) vs Print-Tip for swirl.1")
boxplot(MA_print_tip$M[,2]~MA_print_tip$genes$Block,xlab="Print-Tip",ylab="M",col=rep(palette_okabe_ito(1:8),each=2),main="Normlized M (printtiploess) vs Print-Tip for swirl.2")
boxplot(MA_print_tip$M[,3]~MA_print_tip$genes$Block,xlab="Print-Tip",ylab="M",col=rep(palette_okabe_ito(1:8),each=2),main="Normlized M (printtiploess) vs Print-Tip for swirl.3")
boxplot(MA_print_tip$M[,4]~MA_print_tip$genes$Block,xlab="Print-Tip",ylab="M",col=rep(palette_okabe_ito(1:8),each=2),main="Normlized M (printtiploess) vs Print-Tip for swirl.4")

MA_print_tip_scale <- normalizeBetweenArrays(MA_print_tip,method="scale")
plotMD(MA_print_tip_scale)
abline(0,0,col="blue")
plotPrintTipLoess(MA_print_tip_scale)
boxplot(MA_print_tip_scale$M[,1]~MA_print_tip_scale$genes$Block,xlab="Print-Tip",ylab="M",col=rep(palette_okabe_ito(1:8),each=2),main="Normlized M (printtiploess + scale) vs Print-Tip for swirl.1")
boxplot(MA_print_tip_scale$M[,2]~MA_print_tip_scale$genes$Block,xlab="Print-Tip",ylab="M",col=rep(palette_okabe_ito(1:8),each=2),main="Normlized M (printtiploess + scale) vs Print-Tip for swirl.2")
boxplot(MA_print_tip_scale$M[,3]~MA_print_tip_scale$genes$Block,xlab="Print-Tip",ylab="M",col=rep(palette_okabe_ito(1:8),each=2),main="Normalized M (printtiploess + scale)vs Print-Tip for swirl.3")
boxplot(MA_print_tip_scale$M[,4]~MA_print_tip_scale$genes$Block,xlab="Print-Tip",ylab="M",col=rep(palette_okabe_ito(1:8),each=2),main="Normalized M (printtiploess + scale) vs Print-Tip for swirl.4")

design <- modelMatrix(targets, ref="wild type")
design
fit_print_tip_scale <- lmFit(MA_print_tip_scale,design)

ebayes_print_tip_scale <- eBayes(fit_print_tip_scale) ## the current function is eBayes
qqt(ebayes_print_tip_scale$t,
    df=ebayes_print_tip_scale$df.prior+ebayes_print_tip_scale$df.residual,pch=16,cex=0.2)
plotMD(ebayes_print_tip_scale)
abline(0,0,col="blue")
top30 <- order(ebayes_print_tip_scale$lods,decreasing=TRUE)[1:30]
text(ebayes_print_tip_scale$Amean[top30],ebayes_print_tip_scale$coef[top30],labels=ebayes_print_tip_scale$genes[top30,"Name"],cex=0.8,col="blue")
volcanoplot(ebayes_print_tip_scale,ylab="Surprisal (-logP)")

ordinary.t <- fit_print_tip_scale$coef / fit_print_tip_scale$stdev.unscaled / fit_print_tip_scale$sigma
head(ordinary.t)
library(tidyverse)
tibble(moderated_t=unname(ebayes_print_tip_scale$t),ordinary_t=unname(ordinary.t)) |> mutate(shrink=(abs(moderated_t) < abs(ordinary_t))) |> ggplot(aes(x=ordinary_t,y=moderated_t,color=shrink)) + geom_point(alpha=0.5) + labs(title="Effect of Limma Moderation on Gene-Wise t") +  scale_color_okabe_ito()

options(digits=3)
topTable(ebayes_print_tip_scale,number=30)
tt <- topTable(ebayes_print_tip_scale,n=Inf)
head(tt[tt$adj.P.Val <= 0.1,])
dim(tt[tt$adj.P.Val <= 0.1,])
dt <- decideTests(ebayes_print_tip_scale,p.value = 0.1)
summary(dt)
plotMD(ebayes_print_tip_scale,status=dt)

library(FDRestimation)
hist(ebayes_print_tip_scale$p.value)

ebayes.fit.fdr <- p.fdr(pvalues=ebayes_print_tip_scale$p.value,threshold = 0.1)
plot(ebayes.fit.fdr, main="Benjamini-Hochberg FDRs for Swirl")
plot(ebayes.fit.fdr,xlim = c(0,1500),ylim=c(0,0.2),main="Blowup for Benjamini-Hochberg FDRs for Swirl")

sum(ebayes.fit.fdr$fdrs <= 0.1)
sum(ebayes.fit.fdr$`Results Matrix`$`BH FDRs`<= 0.1)
sum(ebayes.fit.fdr$`Results Matrix`$`Adjusted p-values`<= 0.1)
sum(ebayes.fit.fdr$`Reject Vector` == "Reject.H0")


ebayes.fit.fdr.set.pi0 <- p.fdr(pvalues=ebayes_print_tip_scale$p.value,
                                set.pi0 = get.pi0(ebayes_print_tip_scale$p.value),threshold = 0.1)
plot(ebayes.fit.fdr.set.pi0,main="Benjamini-Hochberg FDRs for Swirl with pi0 estimation")
plot(ebayes.fit.fdr.set.pi0,xlim = c(0,1500),ylim=c(0,0.2),main="Blowup for Benjamini-Hochberg FDRs for Swirl with pi0 estimation")


sum(ebayes.fit.fdr.set.pi0$fdrs <= 0.1)
sum(ebayes.fit.fdr.set.pi0$`Results Matrix`$`BH FDRs`<= 0.1)
sum(ebayes.fit.fdr.set.pi0$`Results Matrix`$`Adjusted p-values`<= 0.1)
sum(ebayes.fit.fdr.set.pi0$`Reject Vector` == "Reject.H0")
