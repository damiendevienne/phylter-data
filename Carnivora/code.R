
############################################
######### PHY vs TS COMP. FOR PAPER ########
############################################
require(ape)
require(phylter)
require(phangorn)
require(ggplot2)
source("../functions.R")

######### 1. get data computed on ALL sequences by T. Tricou
## 

data<-read.table("Outlier_Duplication.txt", header=TRUE)
##recompute GC3 BETTER and length because original values were incorrect
GC3_2<-apply(data,1,function(x,y) sum(as.numeric(x[2:65])*y)/sum(as.numeric(x[2:65])),y=sapply(colnames(data)[2:65], function(x) strsplit(x,"")[[1]][3] %in% c("C","G")))
data<-cbind(data,GC3_2)
##count number of nodes to root in species tree
sptr<-read.nexus("../phylter-data/Carnivora/eLife/Supermatrix_14307_genes_53_spp_UFBS_TESTNEW_Codon.tre")
nodes2rootperspp<-data.frame(spp=sptr$tip.label, nodes2root=unlist(sapply(1:Ntip(sptr),function(x,tr) length(Ancestors(tr,x)),tr=sptr)))
n2r<-nodes2rootperspp$nodes2root[match(data$ID, nodes2rootperspp$spp)]
data<-cbind(data,n2r)
data<-cbind(data, norm_score_dup=data$score_dup/data$n2r)
#compute number of ZNF-labeled genes in dataset
znf<-rep("non-ZNF",nrow(data))
znf[grep("ZNF",data$FAM)]<-"ZNF"
data$znf<-znf


######### 2. Get (for Phylter) and compute (for treeshrink) list of outliers. 

trees<-readCARNIVORAzero()

getPhyRes<-function(inittrees, outfile="file1_v2.txt") {
	PHY<-read.table(outfile)
	K<-paste(PHY[,2],PHY[,1],"_final_align_NT.aln",sep="")
	PHY[,1]<-gsub(".treefile","",PHY[,1])
	k<-unique(PHY[,1])
	trees2<-inittrees
	trees2remove<-NULL
	for (j in k) {
		where<-which(names(inittrees)==j)
		tips2remove<-PHY[PHY[,1]==j,2]
		if(sum(duplicated(tips2remove))>0) cat ("ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo")
		if (length(tips2remove)<Ntip(inittrees[[where]])) trees2[[where]]<-drop.tip(inittrees[[where]], tips2remove)
		else {
			trees2remove<-c(trees2remove, where)
		}
	}
	TREES2REMOVE<-c(trees2remove, which(lapply(trees2, Ntip)<=2))
	return(list(tr=trees2, tr2rm=TREES2REMOVE, outl=PHY, k=K))
}
getPhyRes2<-function(file) {
	phyres<-read.table(file)
	phyres<-cbind(phyres,k=paste(phyres[,2],phyres[,1],"_final_align_NT.aln",sep=""))
	phyres
}
getTSres<-function(inittrees, TS) {
	tsres<-NULL
	for (i in 1:length(TS)) {
		out<-setdiff(inittrees[[i]]$tip.label, TS[[i]]$tip.label)
		if (length(out)>0)	{
			tsres<-rbind(tsres, cbind(names(inittrees)[i], out))
		}
	}
	tsres<-as.data.frame(tsres)
	tsres$k<-paste(tsres[,2],tsres[,1],"_final_align_NT.aln",sep="")
	tsres
}

#### 2.1. get Phylter Results. 

##Latest version ofPhylter was run elsewhere and result files were imported (for ex. file1_v2.txt). look at the header of the files to find parameters
phyres1<-getPhyRes(trees, "file1_v2.txt")
phyres4<-getPhyRes(trees, "file_k1.55_v2.txt") #gives the closest nb of outliers as compared to treeshrink with default values

#### 2.2. run TeeShrink and get results

TS0.012<-runtreeshrink(trees, q=0.012) #gives 7032 outliers. using 0.013 gives 7491. We keep 0.012
TS0.05<-runtreeshrink(trees, q=0.05)
TS1<-TS0.012
TS4<-TS0.05
tsres1<-getTSres(trees, TS0.012)
tsres4<-getTSres(trees, TS0.05)


######### 2BIS. Compare sets of outliers, 
######### with proportionnal venn diagrams
require(eulerr)

tsphy1<-cbind(PhylteR=is.element(unique(c(phyres1$k,tsres1$k)), phyres1$k), TreeShrink=is.element(unique(c(phyres1$k,tsres1$k)), tsres1$k))
tsphy4<-cbind(PhylteR=is.element(unique(c(phyres4$k,tsres4$k)), phyres4$k), TreeShrink=is.element(unique(c(phyres4$k,tsres4$k)), tsres4$k))

pdf("euler-plots-like-venns.pdf")
par(mfrow=c(1,2))
plot(euler(tsphy1, shape="ellipse"), quantities=TRUE, fills=c("#f8766d","#00bfc4"), edges=FALSE)
plot(euler(tsphy4, shape="ellipse"), quantities=TRUE, fills=c("#f8766d","#00bfc4"), edges=FALSE)
dev.off()
##This plot is modified afterwards on Inkscape.



######### 3. Compute and plot distributions showing phylter-outliers, 
#########    treeshrink-outliers and random-outliers 

#### 3.1. GC3

gedatagc3<-function(tsres, phyres) {
	GC3_treeshrink<-data[match(tsres$k,data$paste),]$GC3_2
	GC3_outlier<-data[match(phyres$k,data$paste),]$GC3_2
	# Distribution of GC3 of random sequences
	GC3_random<-data[sample(1:nrow(data),length(phyres$k)),]$GC3_2
	DFGC3<-data.frame(GC3=c(GC3_outlier, GC3_random, GC3_treeshrink), type=c(rep(c("PhylteR outliers","Random outliers"),each=length(GC3_outlier)),rep("TreeShrink outliers",length(GC3_treeshrink))))
	DFGC3$type<-factor(DFGC3$type, levels=c("PhylteR outliers","TreeShrink outliers","Random outliers"))
#	p<-ggplot(DFGC3, aes(x=GC3, fill=type)) + geom_density(alpha=.3)
	DFGC3
}
D1<-cbind(gedatagc3(tsres1,phyres1), dataset="small")
D4<-cbind(gedatagc3(tsres4,phyres4), dataset="large")
D14<-rbind(D1,D4)
D14$dataset<-factor(D14$dataset, levels=c("small","large"))
gc3plot<-ggplot(D14, aes(x=GC3, fill=type)) + geom_density(alpha=.6) + facet_wrap(~dataset) + theme(legend.title = element_blank(), legend.text=element_text(size=7)) + scale_fill_manual(values=c("#f8766d","#00bfc4", "lightgrey")) + theme(legend.position="none")
ggsave("gc3plot.png", gc3plot, height=5, dpi=600)

#### 3.2. sequence length

getdatalen<-function(tsres, phyres) {
	length_outlier<-data[match(phyres$k,data$paste),]$length
	length_treeshrink<-data[match(tsres$k,data$paste),]$length
	length_random<-data[sample(1:nrow(data),length(phyres$k)),]$length
	DFlength<-data.frame(length=c(length_outlier, length_random, length_treeshrink), type=c(rep(c("PhylteR outliers","Random outliers"),each=length(length_outlier)),rep("TreeShrink outliers",length(length_treeshrink))))
	DFlength$type<-factor(DFlength$type, levels=c("PhylteR outliers","TreeShrink outliers","Random outliers"))
	DFlength
}
D1<-cbind(getdatalen(tsres1,phyres1), dataset="small")
D4<-cbind(getdatalen(tsres4,phyres4), dataset="large")
D14<-rbind(D1,D4)
D14$dataset<-factor(D14$dataset, levels=c("small","large"))
lenplot<-ggplot(D14, aes(y=length, fill=type)) + geom_boxplot(alpha=.6) + facet_wrap(~dataset) + ylim(1,15000) + xlab("") + ylab("sequence length (bp)") + scale_fill_manual(values=c("#f8766d","#00bfc4", "lightgrey")) + theme(axis.text.x = element_text(color="white"), axis.ticks = element_blank(), legend.position="none") 
ggsave("lenplot.png", lenplot, height=5, dpi=600)


#### 3.3. Normalized duplication score


getdatanormdup<-function(tsres, phyres) {
	normdup_treeshrink<-data[match(tsres$k,data$paste),]$norm_score_dup
	normdup_outlier<-data[match(phyres$k,data$paste),]$norm_score_dup
	normdup_random<-data[sample(1:nrow(data),length(phyres$k)),]$norm_score_dup
	DFnormdup<-data.frame(normdup=c(normdup_outlier, normdup_random, normdup_treeshrink), type=c(rep(c("PhylteR outliers","Random outliers"),each=length(normdup_outlier)),rep("TreeShrink outliers",length(normdup_treeshrink))))
	DFnormdup$type<-factor(DFnormdup$type, levels=c("PhylteR outliers","TreeShrink outliers","Random outliers"))
	DFnormdup
#	ggplot(DFnormdup, aes(x=normdup, fill=type)) + geom_density(alpha=.3)
	#ggplot(DFnormdup, aes(x=normdup, fill=type)) + geom_boxplot(alpha=.3)
}
D1<-cbind(getdatanormdup(tsres1,phyres1), dataset="small")
D4<-cbind(getdatanormdup(tsres4,phyres4), dataset="large")
D14<-rbind(D1,D4)
D14$dataset<-factor(D14$dataset, levels=c("small","large"))
dupplot<-ggplot(D14, aes(x=normdup, fill=type)) + geom_density(alpha=.6) +  facet_wrap(~dataset) + theme(legend.title = element_blank(), legend.text=element_text(size=7)) + xlab("Duplication score") + scale_fill_manual(values=c("#f8766d","#00bfc4", "lightgrey")) + theme(legend.position="none")
ggsave("dupplot.png", dupplot, height=5, dpi=600)




#### 3.4. ZNF genes
getdataznf<-function(tsres, phyres) {
	df1<-data.frame(as.data.frame(table(data$znf)),type="ALL")
	df2<-data.frame(as.data.frame(table(data$znf[match(phyres$k,data$paste)])),type="PhylteR outliers")
	df3<-data.frame(as.data.frame(table(data$znf[match(tsres$k,data$paste)])),type="TreeShrink outliers")
	df4<-data.frame(as.data.frame(table(data$znf[sample(1:nrow(data),length(phyres$k))])),type="Random outliers")
	prop1<-df1$Freq[df1$Var1=="ZNF"]/sum(df1$Freq)
	prop2<-df2$Freq[df2$Var1=="ZNF"]/sum(df2$Freq)
	prop3<-df3$Freq[df3$Var1=="ZNF"]/sum(df3$Freq)
	prop4<-df4$Freq[df4$Var1=="ZNF"]/sum(df4$Freq)
	pdfznf<-data.frame(prop=c(prop2, prop3,prop4), type=c("PhylteR outliers","TreeShrink outliers","Random outliers"))
	pdfznf$type<-factor(pdfznf$type, levels=c("PhylteR outliers","TreeShrink outliers","Random outliers"))
	DFznf<-rbind(df1,df2, df3,df4)
#	ggplot(DFznf, aes(x=type,y=Freq, fill=Var1)) + geom_bar(position = "fill",stat = "identity") + scale_y_continuous(labels = scales::percent_format())
	#ggplot(DFznf[DFznf$Var1=="ZNF",], aes(x=type,y=Freq)) + geom_bar(position = "fill",stat = "identity")

	print(length(grep("ZNF",unique(data$FAM[match(phyres$k,data$paste)]))))
	print(length(grep("ZNF",unique(data$FAM[match(tsres$k,data$paste)]))))
	print(length(grep("ZNF",unique(data$FAM[sample(data$paste, length(phyres$k))]))))

	print(length(grep("ZNF",unique(data$FAM))))


	return(list(pdfznf=pdfznf, DFznf=DFznf))
}
D1<-cbind(getdataznf(tsres1,phyres1)$pdfznf, dataset="small")
D4<-cbind(getdataznf(tsres4,phyres4)$pdfznf, dataset="large")
D14<-rbind(D1,D4)
D14$dataset<-factor(D14$dataset, levels=c("small","large"))
znfplot<-ggplot(D14, aes(x=type,y=prop, fill=type)) + geom_bar(stat="identity",alpha=.6, colour="black") + scale_fill_manual(values=c("#f8766d","#00bfc4", "lightgrey")) + facet_wrap(~dataset) + ylim(0,0.08) + xlab("") + theme(axis.text.x = element_text(color="white"), axis.ticks = element_blank(), legend.position="none") + ylab("Proportion of outliers being ZNF")


###PLOT ALL GRAPHS IN ONE! 

require("ggpubr")
THEPLOT<-ggarrange(lenplot, dupplot, znfplot,labels = c("A", "B", "C"),ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom") + theme(legend.title = element_blank())
ggsave("THEBIGPLOT.pdf",THEPLOT,width=6, height=10, dpi=600)


######### 4. Gene Concordance Factor (gCF analysis) 

gCFinit<-computeGCF(trees, sptr)
gCFphysmall<-computeGCF(phyres1$tr[-phyres1$tr2rm], sptr)
gCFphylarge<-computeGCF(phyres4$tr[-phyres4$tr2rm], sptr)
gCFtssmall<-computeGCF(TS1,sptr)
gCFtslarge<-computeGCF(TS4,sptr)

##prepare dataset for ggplot
smallphy<-data.frame(Method="PhylteR",delta=sort(gCFphysmall$gCF-gCFinit$gCF), x=1:length(gCFinit$gCF), paired=1:length(gCFinit$gCF), size=gCFinit$gCF[order(gCFphysmall$gCF-gCFinit$gCF)])
smallts<-data.frame(Method="TreeShrink",delta=(gCFtssmall$gCF-gCFinit$gCF)[order(gCFphysmall$gCF-gCFinit$gCF)], x=1:length(gCFinit$gCF), paired=1:length(gCFinit$gCF), size=gCFinit$gCF[order(gCFphysmall$gCF-gCFinit$gCF)])
small<-rbind(smallphy,smallts)
small<-cbind(small, dataset="small")
small$shape<-cut(small$size, c(0,70,100), labels=c("0%-70%","71%-100%"))
#ggplot(small,aes(x=x,y=delta, color=Method)) + geom_line(aes(group = paired), color="grey") + geom_point(aes(shape=shape), size=3) + scale_shape_manual(values = c(21, 16))


largephy<-data.frame(Method="PhylteR",delta=sort(gCFphylarge$gCF-gCFinit$gCF), x=1:length(gCFinit$gCF), paired=1:length(gCFinit$gCF), size=gCFinit$gCF[order(gCFphylarge$gCF-gCFinit$gCF)])
largets<-data.frame(Method="TreeShrink",delta=(gCFtslarge$gCF-gCFinit$gCF)[order(gCFphylarge$gCF-gCFinit$gCF)], x=1:length(gCFinit$gCF), paired=1:length(gCFinit$gCF), size=gCFinit$gCF[order(gCFphylarge$gCF-gCFinit$gCF)])
large<-rbind(largephy,largets)
large<-cbind(large, dataset="large")
large$shape<-cut(large$size, c(0,70,100), labels=c("0%-70%","71%-100%"))
ggplot(large,aes(x=x,y=delta, color=Method)) + geom_line(aes(group = paired), color="grey") + geom_point(aes(shape=shape), size=3) + scale_shape_manual(values = c(21, 16))

#ggplot(large,aes(x=x,y=delta, color=Method)) + geom_line(aes(group = paired), color="grey") + geom_point(size=3)

SL<-rbind(small, large)
SL$dataset<-factor(SL$dataset, levels=c("small","large"))
#names(SL)[which(names(SL)=="shape")]<-"Initial gCF"

gcfplot<-ggplot(SL,aes(x=x,y=delta, color=Method)) + geom_line(aes(group = paired), color="grey") + geom_point(aes(shape=shape), size=3) + scale_shape_manual(values = c(21, 16)) + facet_wrap(~dataset) + xlab("branches of the species tree") + ylab(expression(Delta*gCF~(gCF[init]-gCF))) + theme(legend.position=c(0.12,0.8)) + labs(shape="Initial gCF") + guides(color = guide_legend(order = 1))

ggsave("gcfplot.png",gcfplot, dpi=600)

