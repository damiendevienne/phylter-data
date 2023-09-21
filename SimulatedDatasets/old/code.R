require(phylter)
require(ape)
require(ggplot2)
source("../functions.R")

############################################
############# SIMULATIONS FOR PAPER ########
############################################
###
### SD1
###
trees.sd1<-simtrees(25, 20, Nb.cell.outlier=10, out.type="topology", brlen.sd=0.15)
write.tree(trees.sd1, file="SD1.tre")
pdf("trees.sd1.pdf")
par(mfrow=c(5,5))
par(mar=c(0,0,0,0))
for (i in 1:25) plot(trees.sd1[[i]])
dev.off()

phy.sd1<-phylter(trees.sd1)
pdf("mat.sd1.pdf")
par(mfrow=c(5,5))
par(mar=c(1,1,1,1))
for (i in 1:25) image(phy.sd1$Initial$matrices[[i]], axes=FALSE)
dev.off()

pdf("RV.sd1.pdf")
heatmap(phy.sd1$Initial$RV, scale="none", Rowv=NA, Colv=NA)
dev.off()

pdf("weights.sd1.pdf", height=3)
plot(phy.sd1$Initial$weights, type="h", lwd=10, frame=FALSE, axes=F, ylim=c(0, max(phy.sd1$Initial$weights)))
axis(2)
dev.off()

pdf("compromise.sd1.pdf")
heatmap(phy.sd1$Initial$compromise, scale="none", Rowv=NA, Colv=NA)
dev.off()


gn.of.interest<-1
gn.no.interest<-2
bigmat<-do.call(rbind,phy.sd1$Initial$PartialF)
allmat<-lapply(phy.sd1$Initial$PartialF, function(x,y) cbind(x,y),y=phy.sd1$Initial$F)
ALLMAT<-allmat
interestmat<-allmat[[gn.of.interest]]
allmat<-allmat[-gn.of.interest]
startdraw<-function() {
	plot(bigmat, type="n", frame=FALSE, axes=FALSE)
	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ebebeb", border=NA)
	abline(h=seq(-1,1,by=0.2),v=seq(-1,1,by=0.2),col="white",lwd=3)
	abline(h=seq(-1,1,by=0.1),v=seq(-1,1,by=0.1),col="white",lwd=1.5)
	axis(1,lwd=0, lwd.ticks=3, col.ticks="grey", col.axis="grey59", cex.axis=0.75)
	axis(2,lwd=0, lwd.ticks=3, col.ticks="grey", col.axis="grey59", cex.axis=0.85)
}

pdf("proj1_1.sd1.pdf", height=10)
startdraw()
for (i in 1:20) {
	points(phy.sd1$Initial$F[i,1],phy.sd1$Initial$F[i,2], pch=21, bg="#357db3",col="black", cex=3)
	text(phy.sd1$Initial$F[i,1],phy.sd1$Initial$F[i,2], label=rownames(phy.sd1$Initial$F)[i], col="white", cex=0.8)
}
dev.off()

pdf("proj1_2.sd1.pdf", height=10)
startdraw()
apply(allmat[[4]],1,function(y) segments(y[1],y[2],y[3],y[4], col="black"))
points(allmat[[4]][,1:2], col="black", pch=19, cex=0.5)
for (i in 1:20) {
	points(phy.sd1$Initial$F[i,1],phy.sd1$Initial$F[i,2], pch=21, bg="#357db3",col="black", cex=3)
	text(phy.sd1$Initial$F[i,1],phy.sd1$Initial$F[i,2], label=rownames(phy.sd1$Initial$F)[i], col="white", cex=0.8)
}
dev.off()

pdf("proj1_3.sd1.pdf", height=10)
startdraw()
apply(allmat[[gn.no.interest]],1,function(y) segments(y[1],y[2],y[3],y[4], col="black"))
points(allmat[[gn.no.interest]][,1:2], col="black", pch=19, cex=0.5)

apply(interestmat,1,function(y) segments(y[1],y[2],y[3],y[4], col="#d11b26"))
points(interestmat[,1:2], col="#d11b26", pch=19, cex=0.5)

for (i in 1:20) {
	points(phy.sd1$Initial$F[i,1],phy.sd1$Initial$F[i,2], pch=21, bg="#357db3",col="black", cex=3)
	text(phy.sd1$Initial$F[i,1],phy.sd1$Initial$F[i,2], label=rownames(phy.sd1$Initial$F)[i], col="white", cex=0.8)
}
dev.off()


pdf("proj1_4.sd1.pdf", height=10)
startdraw()
lapply(allmat, function(x) apply(x,1,function(y) segments(y[1],y[2],y[3],y[4])))
apply(interestmat,1,function(y) segments(y[1],y[2],y[3],y[4], col="red"))
points(do.call(rbind, allmat)[,1:2], col="black", pch=19, cex=0.5)
points(interestmat[,1:2], col="red", pch=19, cex=0.5)
for (i in 1:100) {
	points(phy.sd1$Initial$F[i,1],phy.sd1$Initial$F[i,2], pch=21, bg="#357db3",col="black", cex=3)
	text(phy.sd1$Initial$F[i,1],phy.sd1$Initial$F[i,2], label=rownames(phy.sd1$Initial$F)[i], col="white", cex=0.8)
}
dev.off()


pdf("proj1_4_extra.sd1.pdf", height=10)
startdraw()
lapply(ALLMAT[-5], function(x) apply(x,1,function(y) segments(y[1],y[2],y[3],y[4])))
apply(ALLMAT[[5]][-5,],1,function(y) segments(y[1],y[2],y[3],y[4], col="black"))
points(do.call(rbind, ALLMAT[-5])[,1:2], col="black", pch=19, cex=0.5)
points(ALLMAT[[5]][-5,1:2], col="black", pch=19, cex=0.5)
segments(ALLMAT[[5]][5,][1],ALLMAT[[5]][5,][2],ALLMAT[[5]][5,][3],ALLMAT[[5]][5,][4], col="red")
points(ALLMAT[[5]][5,1],ALLMAT[[5]][5,2], col="red", pch=19, cex=0.5)
for (i in 1:20) {
	points(phy.sd1$Initial$F[i,1],phy.sd1$Initial$F[i,2], pch=21, bg="#357db3",col="black", cex=3)
	text(phy.sd1$Initial$F[i,1],phy.sd1$Initial$F[i,2], label=rownames(phy.sd1$Initial$F)[i], col="white", cex=0.8)
}
dev.off()


pdf("WR.sd1.pdf")
heatmap(phy.sd1$Initial$WR, scale="none", Rowv=NA, Colv=NA)
dev.off()


###
### SD2
###
NBT<-20
sd2<-replicate(NBT, simtrees(100,40, Nb.cell.outlier=10, out.type="topology", brlen.sd=0.15),simplify=FALSE)
trees.sd2<-lapply(sd2, function(x) x$trees)
outl.sd2<-lapply(sd2, function(x) x$outl[,2:1])

system("rm SD2/*")
for (i in 1:NBT) {
	write.tree(trees.sd2[[i]], file=paste("SD2/SD2_",i,".tre",sep=""))
}
#on each dataset we run phyler
PHY<-lapply(trees.sd2, phylter, k2=10)
OUTPHY<-lapply(PHY, function(x) x$Final$Outliers)
OUTTS<-lapply(trees.sd2, runtreeshrink.default.getoutl)
phy.nbout<-unlist(lapply(PHY, function(x) nrow(x$Final$Outliers)))
ts.nbout<-unlist(lapply(OUTTS, nrow))
#compute precision = TP/(TP+FP)
TruePositive<-function(o1,ref) { #o1 dn o2 are two mat of outliers
	if (!is.null(o1)) {
		O1<-apply(o1,1,paste, collapse="-")
		REF<-apply(ref,1,paste, collapse="-")
		TP<-sum(is.element(O1, REF))
		FP<-sum(!is.element(O1, REF))
		FN<-length(setdiff(REF,O1))

	}
	else {
		TP<-0
		FP<-0
		FN<-nrow(ref)
	}
	precision<-TP/(TP+FP)
	recall<-TP/(TP+FN)
	return(c(precision,recall))
}

PR.PHY<-t(sapply(1:length(OUTPHY), function(x, o,r) TruePositive(o[[x]],r[[x]]),o=OUTPHY, r=outl.sd2))
PR.TS<-t(sapply(1:length(OUTTS), function(x, o,r) TruePositive(o[[x]],r[[x]]),o=OUTTS, r=outl.sd2))
DF<-rbind(data.frame(x=PR.PHY[,1],what="precision",method="PhylteR"),data.frame(x=PR.TS[,1],what="precision",method="TreeShrink"), data.frame(x=PR.PHY[,2],what="recall",method="PhylteR"),data.frame(x=PR.TS[,2],what="recall",method="TreeShrink"))
plotsd2<-ggplot(DF, aes(x=what,y=x,fill=method)) + geom_boxplot() + xlab("") + ylab("mean") + ggtitle("Simulated Dataset 2 (SD2)")

DF.nbout<-rbind(data.frame(nbout=phy.nbout, method="PhylteR", dataset="SD2"), data.frame(nbout=ts.nbout, method="TreeShrink", dataset="SD2"))


###
### SD3
###
NBT<-20
sd3<-replicate(NBT, simtrees(100,40, Nb.cell.outlier=10, out.type="brlength", brlen.sd=0.15, bl.mult=20),simplify=FALSE)
trees.sd3<-lapply(sd3, function(x) x$trees)
outl.sd3<-lapply(sd3, function(x) x$outl[,2:1])

system("rm SD3/*")
for (i in 1:NBT) {
	write.tree(trees.sd3[[i]], file=paste("SD3/SD3_",i,".tre",sep=""))
}
#on each dataset we run phyler
PHY3<-lapply(trees.sd3, phylter, k2=10)
OUTPHY3<-lapply(PHY3, function(x) x$Final$Outliers)
OUTTS3<-lapply(trees.sd3, runtreeshrink.default.getoutl)
phy.nbout3<-unlist(lapply(PHY3, function(x) nrow(x$Final$Outliers)))
ts.nbout3<-unlist(lapply(OUTTS3, nrow))


PR.PHY3<-t(sapply(1:length(OUTPHY3), function(x, o,r) TruePositive(o[[x]],r[[x]]),o=OUTPHY3, r=outl.sd3))
PR.TS3<-t(sapply(1:length(OUTTS3), function(x, o,r) TruePositive(o[[x]],r[[x]]),o=OUTTS3, r=outl.sd3))
DF3<-rbind(data.frame(x=PR.PHY3[,1],what="precision",method="PhylteR"),data.frame(x=PR.TS3[,1],what="precision",method="TreeShrink"), data.frame(x=PR.PHY3[,2],what="recall",method="PhylteR"),data.frame(x=PR.TS3[,2],what="recall",method="TreeShrink"))
plotsd3<-ggplot(DF3, aes(x=what,y=x,fill=method)) + geom_boxplot() + xlab("") + ylab("mean") + ggtitle("Simulated Dataset 3 (SD3)")

DF3.nbout<-rbind(data.frame(nbout=phy.nbout3, method="PhylteR", dataset="SD3"), data.frame(nbout=ts.nbout3, method="TreeShrink", dataset="SD3"))

###
### SD4
###
NBT<-20
sd4<-replicate(NBT, simtrees(100,40, Nb.cell.outlier=20, out.type="both", brlen.sd=0.15, bl.mult=20),simplify=FALSE)
trees.sd4<-lapply(sd4, function(x) x$trees)
outl.sd4<-lapply(sd4, function(x) x$outl[,2:1])

system("rm SD4/*")
for (i in 1:NBT) {
	write.tree(trees.sd4[[i]], file=paste("SD4/SD4_",i,".tre",sep=""))
}
#on each dataset we run phyler
PHY4<-lapply(trees.sd4, phylter, k2=10)
OUTPHY4<-lapply(PHY4, function(x) x$Final$Outliers)
OUTTS4<-lapply(trees.sd4, runtreeshrink.default.getoutl)
phy.nbout4<-unlist(lapply(PHY4, function(x) nrow(x$Final$Outliers)))
ts.nbout4<-unlist(lapply(OUTTS4, nrow))


PR.PHY4<-t(sapply(1:length(OUTPHY4), function(x, o,r) TruePositive(o[[x]],r[[x]]),o=OUTPHY4, r=outl.sd4))
PR.TS4<-t(sapply(1:length(OUTTS4), function(x, o,r) TruePositive(o[[x]],r[[x]]),o=OUTTS4, r=outl.sd4))
DF4<-rbind(data.frame(x=PR.PHY4[,1],what="precision",method="PhylteR"),data.frame(x=PR.TS4[,1],what="precision",method="TreeShrink"), data.frame(x=PR.PHY4[,2],what="recall",method="PhylteR"),data.frame(x=PR.TS4[,2],what="recall",method="TreeShrink"))
plotsd4<-ggplot(DF4, aes(x=what,y=x,fill=method)) + geom_boxplot() + xlab("") + ylab("mean") + ggtitle("Simulated Dataset 4 (SD4)")

DF4.nbout<-rbind(data.frame(nbout=phy.nbout4, method="PhylteR", dataset="SD4"), data.frame(nbout=ts.nbout4, method="TreeShrink", dataset="SD4"))


#THE TWO PLOTS FROM THERE !
DFALL<-rbind(cbind(DF,dataset="SD2"), cbind(DF3,dataset="SD3"), cbind(DF4,dataset="SD4"))
p<-ggplot(DFALL, aes(x=what,y=x,fill=method)) + geom_boxplot() + xlab("") + ylab("mean") + ylab("mean") + facet_wrap(~dataset) +ylim(0,1)
ggsave("simulations-plots.pdf", p, height=5)


DFALL.nbout<-rbind(DF.nbout, DF3.nbout, DF4.nbout)
nbout.theoretical<-data.frame(dataset=c("SD2","SD3","SD4"), z = c(10, 10,20))
ggplot(DFALL.nbout, aes(x=method,y=nbout,fill=nbout)) + geom_boxplot() + xlab("") + ylab("mean") + ylab("mean") + facet_wrap(~dataset) + geom_hline(data=nbout.theoretical, aes(yintercept=nbout), linetype=3)



