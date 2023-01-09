
#read the carnivora dataset, naming the gebne trees accordng to the file names
readCARNIVORAzero<-function() {
	#list all trees
	system("ls Carnivora/trees/* > listoftreesCarnivorazero")
	fi<-readLines("listoftreesCarnivorazero")
	genelist<-gsub("Carnivora/trees/","",fi)
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


