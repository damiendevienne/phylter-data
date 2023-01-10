require(ape)
require(phylter)

######### FONCTIONS
#run treeshrink with varrious q values and return pruned trees
runtreeshrink<-function(trees, q=0.05) {
	filename<-paste("fortreeshrink",as.character(runif(1)),".tre", sep="")
	foldname<-paste(gsub("\\.tre","",filename),"_treeshrink", sep="")
	write.tree(trees, file=filename)
	system(paste("run_treeshrink.py -t ",filename," -q ",q," -m per-species", sep=""))
	cleaned<-read.tree(paste(foldname,"/output.tre", sep=""))
	names(cleaned)<-names(trees)
	system(paste("rm -r ",foldname,"/",sep="")) #clean everything
	system(paste("rm ",filename, sep="")) #clean everything
	return(cleaned)
}

#run treeshrink default parameters and return outliers
runtreeshrink.default.getoutl<-function(trees) {
	filename<-paste("fortreeshrink",as.character(runif(1)),".tre", sep="")
	foldname<-paste(gsub("\\.tre","",filename),"_treeshrink", sep="")
	write.tree(trees, file=filename)
	system(paste("run_treeshrink.py -t ",filename, sep=""))
	cleaned<-read.tree(paste(foldname,"/output.tre", sep=""))
	names(cleaned)<-names(trees)
	system(paste("rm -r ",foldname,"/",sep="")) #clean everything
	system(paste("rm ",filename, sep="")) #clean everything
	L1<-lapply(trees, function(x) x$tip.label)
	L2<-lapply(cleaned, function(x) x$tip.label)
	OUT<-do.call(rbind, sapply(1:length(L1), function(x,a,b) if (length(setdiff(a[[x]],b[[x]]))>0) cbind(x, setdiff(a[[x]],b[[x]])),a=L1,b=L2))
	return(OUT)
}

# compute gene Concordance Factor for one tree (sptree) and multiple gene trees (trees)
computeGCF<-function(trees, sptree) {
	trees<-lapply(trees, unroot)
	sptree<-unroot(sptree)
	speciestreefilename<-paste("speciestree",as.character(runif(1)),".tre", sep="")
	write.tree(sptree, file="species.tre")
	write.tree(trees, file="genes.tre")
	system("iqtree2 -t species.tre --gcf genes.tre --prefix concord")
	GCF<-read.table("concord.cf.stat", header=TRUE, fill=TRUE)
	system("rm concord.*")
	return(GCF)
}

#read the carnivora dataset, naming the gebne trees accordng to the file names
readCARNIVORAzero<-function() {
	#list all trees
	system("ls data/trees/* > listoftreesCarnivorazero")
	fi<-readLines("listoftreesCarnivorazero")
	genelist<-gsub("data/trees/","",fi)
	genelist<-gsub(".treefile","",genelist)
	T<-list()
	for (i in 1:length(fi)) {
		print(i)
		T[[i]]<-read.tree(fi[i])
	}	
	class(T)<-"multiPhylo"
	genes<-genelist
	trees<-T
	names(trees)<-genes
	return(trees)
}
