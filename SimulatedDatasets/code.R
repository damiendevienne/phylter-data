##LOAD PACKAGES
require(ggplot2)
##SOURCE FUNCTIONS
source("../functions.R")

## This scripts allows performing the simulations
## for the article presenting PhylteR. 

## 1-TRANSFORM CARNIVORA TREE TO SIMPHY-COMPATIBLE TREE
spu<-preparetree2()
## 2 - RUN SIMPHY
## increase sp (nb of species per pop) to increase the level of ILS
## cs is the seed for each simulation. Needs to change each time otherwise all simulations are identical for a given set of parameters. 
## we generate 500 gene trees and we sample 100 of them afterwards, keeping only cases with 0 or 1 HGT only.

##########################################################################
### SIMULATION 1: comparison of TSH and PHY                            ###
### in terms of correctly detecting misplaced species (result of HGT)  ###
### we vary the number of                                              ###
##########################################################################

RESULTS<-matrix(NA, nrow=1, ncol=13)
hg<-100
hgt<-1e-8
sp<-c(10, 100000, 200000, 500000)
maxnb<-c(1,10,53)
rep<-1:100
method<-c("phy", "tsh")
cpt<-0
for (spi in sp) {
	for (ma in maxnb) {
			for (re in rep) {
				cpt<-cpt+1
				print(cpt)
				print(nrow(RESULTS))
				TREES<-runSimPhy(rl=500,sp=spi, s=write.tree(spu), lt=hgt, cs=sample.int(1e10,1), hg=hg)
				formatted_results<-analyse2(TREES, ma, 100)
				res1<-phylter(formatted_results$GENETREES, k2=10000)$Final$Outliers
				res2<-runtreeshrink.default.getoutl(formatted_results$GENETREES)
				nbout1<-nrow(res1)
				nbout2<-nrow(res2)
				nboutreal<-nrow(formatted_results$FORMATTED_OUTLIERS)
				tpfpfnprecrec1<-TruePositive(res1, formatted_results$FORMATTED_OUTLIERS)
				tpfpfnprecrec2<-TruePositive(res2, formatted_results$FORMATTED_OUTLIERS)
				RES<-rbind(c(re,"phy", spi, hg, hgt, ma, nbout1, tpfpfnprecrec1, nboutreal),c(re,"tsh", spi, hg, hgt, ma, nbout2, tpfpfnprecrec2, nboutreal))
				RESULTS<-rbind(RESULTS, RES)
#				print(RESULTS)
			}
	}
}
RESULTS<-RESULTS[-1] #remove first row made of NA.

## THESE SIMULATIONS WHERE DONE SIMULTANEOUSLY ON 4 MACHINES, WITH VALUE FOR SP
## SO THAT WE NOW NEED TO READ THE FOUR TABLES AND COMBINE THEM. THIS IS 
## WHAT IS DONE HERE
RES1<-read.table("RESULTS-NEWSIMUL-SP-10.tab")
RES2<-read.table("RESULTS-NEWSIMUL-SP-1e+05.tab")
RES3<-read.table("RESULTS-NEWSIMUL-SP-2e+05.tab")
RES4<-read.table("RESULTS-NEWSIMUL-SP-5e+05.tab")
# RESULTS<-rbind(RES1,RES2,RES3,RES4)

RESDF<-as.data.frame(RESULTS)
colnames(RESDF)<-c("rep","Method","sp","hg","hgt","maxnb","nbout","TP","FP","FN","Precision","Recall","realnbout")
RESDF[,c(1,3:13)] <- lapply(RESDF[,c(1,3:13)], as.numeric)


RESDFsave<-RESDF
RESDF2<-RESDF


RESDF2$TP<-(RESDF$realnbout*RESDF$nbout)/5300
RESDF2$FP<-RESDF$nbout-RESDF2$TP
RESDF2$FN<-RESDF$realnbout-RESDF2$TP
RESDF2$Precision<-RESDF2$TP/(RESDF2$TP+RESDF2$FP)
RESDF2$Recall<-RESDF2$TP/(RESDF2$TP+RESDF2$FN)
RESDF2$Method<-"random"

RESDFfinal<-rbind(RESDF,RESDF2)

###REFORMAT DATA and COMPUTE PRECISION AND RECALL UNDER RANDOM SELECTION OF OUTLIERS

RESDFfinalprecision<-RESDFfinal[,-12]
names(RESDFfinalprecision)[11]<-"value"
RESDFfinalprecision$what<-"precision"
RESDFfinalrecall<-RESDFfinal[,-11]
names(RESDFfinalrecall)[11]<-"value"
RESDFfinalrecall$what<-"recall"
RESDFfinalOK<-rbind(RESDFfinalrecall, RESDFfinalprecision)
RESDFfinalOK$ILS<-RESDFfinalOK$sp
RESDFfinalOK$ILS[RESDFfinalOK$ILS==10]<-"NO-ILS"
RESDFfinalOK$ILS[RESDFfinalOK$ILS==1e5]<-"LOW-ILS"
RESDFfinalOK$ILS[RESDFfinalOK$ILS==2e5]<-"MODERATE-ILS"
RESDFfinalOK$ILS[RESDFfinalOK$ILS==5e5]<-"HIGH-ILS"

RESDFfinalOK$ILS<-factor(RESDFfinalOK$ILS, levels=c("NO-ILS","LOW-ILS","MODERATE-ILS","HIGH-ILS"), ordered=TRUE)
RESDFfinalOK$Method<-factor(RESDFfinalOK$Method, levels=c("phy","tsh","random"), ordered=TRUE)

pdf("results-final-nbtransfmax-ALL-allils.pdf")
for (i in c(1,10,53)) {
	ggplot(RESDFfinalOK, aes(x=what,y=value,fill=Method)) + geom_boxplot(data=subset(RESDFfinalOK,maxnb==i)) + facet_wrap(~ILS) + scale_fill_manual(values = c("#f8766d", "#00bfc4","lightgrey"), labels=c('PhylteR','TreeShrink','random')) + xlab("") + ylab("mean")
}
dev.off()

#############################################################################################################################
## SIMULATION 2: SAME AS BEFORE BUT SIMPLY COMPUTE TREES AND GET TOPOLOGICAL DISTANCES BETWEEN GENE TREES AND SPECIES TREE ##
## We could have done this directly in the previous script, but we didn't.                                                 ##
#############################################################################################################################
spu<-preparetree2()
spu00<-spu
spu0<-spu
spu00$tip.label<-paste(spu00$tip.label,"_0_0", sep="")	
spu0$tip.label<-paste(spu0$tip.label,"_0", sep="")	
TOPODISTANCES<-NULL
hg<-100
hgt<-1e-8
sp<-c(10, 100000, 200000, 500000)
maxnb<-c(53)
rep<-1:100
cpt<-0
for (spi in sp) {
	for (ma in maxnb) {
			for (re in rep) {
				cpt<-cpt+1
				cat("\n\n\n0000000000000000000000000000000000000000000000000000000000000000000\n")
				cat(cpt)
				cat("\n\n\n0000000000000000000000000000000000000000000000000000000000000000000\n")
				TREES<-runSimPhy(rl=500,sp=spi, s=write.tree(spu), lt=hgt, cs=sample.int(1e10,1), hg=hg)
				formatted_results<-analyse2(TREES, ma, 100)
				topodist.gntr<-unlist(lapply(formatted_results$GENETREES, function(x,y) dist.topo(unroot(x),unroot(y)), y=spu00))
				RES<-c(re, spi, hg, hgt, ma, mean(topodist.gntr), sd(topodist.gntr))
				TOPODISTANCES<-rbind(TOPODISTANCES, RES)
#				print(RESULTS)
			}
	}
}

RESTOPO<-as.data.frame(TOPODISTANCES)
colnames(RESTOPO)<-c("rep","sp","hg","hgt","maxnb","meantopodist","sdtopodist")

RESTOPO$ILS<-RESTOPO$sp
RESTOPO$ILS[RESTOPO$ILS==10]<-"NO-ILS"
RESTOPO$ILS[RESTOPO$ILS==100000]<-"LOW-ILS"
RESTOPO$ILS[RESTOPO$ILS==200000]<-"MODERATE-ILS"
RESTOPO$ILS[RESTOPO$ILS==500000]<-"HIGH-ILS"
RESTOPO$ILS<-factor(RESTOPO$ILS, levels=c("NO-ILS","LOW-ILS","MODERATE-ILS","HIGH-ILS"), ordered=TRUE)
ggplot(RESTOPO, aes(x=ILS, y=meantopodist)) + geom_boxplot()






