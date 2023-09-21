require(ape)
require(phylter)
require(phytools)
require(phangorn)
require(scales)
require(TreeDist)

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
	system(paste("run_treeshrink.py -t ",filename, " -m per-species > trash", sep=""))
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

#read the carnivora dataset, naming the gene trees according to the file names
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

#run Simphy for simulations
runSimPhy<-function(sp,s,rs=1,rl=1000,rg=1,sg=7.62125,su=2.2e-9,ld=0,lb=0,lt=0,lg=0,cs=223041,om=1,ot=0,od=1,op=1,oc=1, lk=0, hg=100) {
	# Explanation of each argument of Simphy: 
	# -rs 1 //Number of replicates
	# -rl f:1000 //1000 locus per replicate
	# -rg 1 //1 gene tree per locus tree
	# -s (Manis_javanica:455847749,(((Paradoxurus_hermaphroditus:169205608,((Cryptoprocta_ferox:116840692,(Suricata_suricatta:43236738,(Helogale_parvula:36575073,Mungos_mungo:36575073):6661665):73603954):22011556,((Hyaena_hyaena:44525423,Crocuta_Crocuta:44525423):3478061,(Proteles_cristatus:28735578,Proteles_septentrionalis:28735578):19267906):90848764):30353360):12602776,((Neofelis_nebulosa:101733296,((Panthera_leo:89276016,Panthera_onca:89276016):8868364,(Panthera_pardus:94378790,Panthera_tigris:94378790):3765590):3588916):7905315,(Prionailurus_bengalensis:105207893,(Felis_catus:102419208,((Lynx_pardinus:91939989,Lynx_canadensis:91939989):9121251,(Acinonyx_jubatus:95307196,Puma_concolor:95307196):5754044):1357968):2788685):4430718):72169773):76105110,(((Canis_familiaris:61953836,Lycaon_pictus:61953836):10642036,(Otocyon_megalotis:69873007,Vulpes_vulpes:69873007):2722865):150117840,((Ailuropoda_melanoleuca:120620140,(Ursus_maritimus:95750616,(Ursus_thibetanus:86573938,(Ursus_americanus:84893799,Ursus_arctos:84893799):1680139):9176678):24869524):70794484,(((Phoca_vitulina:118427298,(Neomonachus_schauinslandi:108616671,(Leptonychotes_weddellii:106358010,Mirounga_angustirostris:106358010):2258661):9810627):20265597,(Odobenus_rosmarus:112506973,(Callorhinus_ursinus:96502348,(Arctocephalus_gazella:93159391,(Eumetopias_jubatus:91288721,Zalophus_californianus:91288721):1870670):3342957):16004625):26185922):46219417,(Spilogale_gracilis:155549714,(Ailurus_fulgens:146040473,((Potos_flavus:110913439,(Nasua_narica:92036988,(Bassariscus_sumichrasti:51703254,Procyon_lotor:51703254):40333734):18876451):25926107,(Taxidea_taxus:84318741,(Mellivora_capensis:72547374,(Gulo_gulo:68269081,((Mustela_putorius:35781771,Neovison_vison:35781771):21439962,(Pteronura_brasiliensis:47674385,(Enhydra_lutris:34080846,Lutra_lutra:34080846):13593539):9547348):11047348):4278293):11771367):52520805):9200927):9509241):29362598):6502312):31299088):35199782):197934255); //Fixed species tree
	# -sg f:7.62125 //Generation time
	# -sp f:1000000 //Population size
	# -su f:2.2e-9 //Tree-wide substitution rate
	# -ld f:0.0 //Loss rate
	# -lb f:0.0 //Duplication rate
	# -lt f:0.0 //HGT rate
	# -lg f:0.0 //Gene conversion rate
	# -cs 223041 //Seed for the random number generator, in order to make the experiment repetible.
	# -o ILSonly_1M
	# -om 1 //Tree mapping output
	# -ot 0 //The species and locus tree branches are written in number of generations
	# -od 1 //Database
	# -op 1 //Output with the general sampled options (describes the simulation run)
	# -oc 1 //Activates the backup of the original command line and configuration file (we recommend to always activate this option)
	# -lk 0 //transfer recipient does not depend on distance to donor of transfer

	#generate a unique folder output name
	o<-paste("simphysimul_",paste(sample(c(LETTERS,letters, 1:10),10, replace=TRUE), collapse=""),sep="")
	print(o)
	genefiles<-paste(o,"_g_trees_files", sep="")
	rl<-paste("f:",rl,sep="")
	sg<-paste("f:",sg,sep="")
	sp<-paste("f:",sp,sep="")
	su<-paste("f:",su,sep="")
	ld<-paste("f:",ld,sep="")
	lb<-paste("f:",lb,sep="")
	lt<-paste("f:",lt,sep="")
	lg<-paste("f:",lg,sep="")
	hg<-paste("f:",hg,sep="")
	s<-paste("\'",s,"\'", sep="")

	#prepare call
	call<-paste("simphy_lnx64 -rs", rs,"-rl",rl ,"-rg", rg,"-s", s,"-sg", sg,"-sp", sp,"-su",su ,"-ld",ld ,"-lb", lb,"-lt", lt,"-lg", lg,"-cs", cs,"-o", o,"-om", om,"-ot", ot,"-od", od,"-op", op,"-oc",oc, "-lk",lk, "-hg",hg, "-ol 1", sep=" ")
	print(call)
	#run simphy
	system(call)
	#GET RESULTS
	#list_trees
	system(paste("ls ",o,"/1/g_trees* > ",genefiles, sep=""))
	print(o)
	allg<-readLines(genefiles)
	allmat<-gsub(".trees",".mapsl",gsub("g_trees","",allg))
	T<-list()
	M<-list()
	for (i in 1:length(allg)) {
		T[[i]]<-read.tree(allg[i])
		if (file.exists(allmat[i])) M[[i]]<-read.table(allmat[i], header=TRUE)
		else M[[i]]<-NULL
	}
	class(T)<-"multiPhylo"

	#GET LOCUS TREES AS WELL (useful for geting transfers)
	locustrees<-read.tree(paste(o,"/1/l_trees.trees", sep=""))
	sptree<-read.tree(paste(o,"/1/s_tree.trees", sep=""))
	##clean everything
	system(paste("rm -r",o,sep=" "))
	return(list(gntr=T,ltr=locustrees, sptr=sptree, maps=M))
}

# Analyse locustrees and genetrees output by Simphy to produce
# a list of "true" outliers
analyse2<-function(alltrees, MAXNB, NBTR) {
	sptr<-alltrees$sptr
	#get trees with 0 or 1 transfer(s)
	nbtransf<-unlist(lapply(alltrees$ltr, function(x) length(grep("Rtransf",x$tip.label))))
	trees2keep<-nbtransf<=1
	LOCUTREES<-alltrees$ltr[trees2keep]
	GENETREES<-alltrees$gntr[trees2keep]
	MATRICES<-alltrees$maps[trees2keep]
	nbtransf2<-nbtransf[trees2keep]

	GetGroupOfMovedSpecies<-function(lctr, sptr, mtr) {
		moved<-lctr$tip.label[grep("Rtransf",lctr$tip.label)]
		if (length(moved)>0) {
			moved<-gsub("_0","",moved)
			moved<-gsub("Rtransf-","",moved)
			if (is.element(moved, sptr$tip.label)) MOVEDspecies<-moved
			else {
				nodemoved<-as.numeric(moved)
				M<-mtr$Lt_node[mtr$St_node==nodemoved]
				nodereallymoved<-M[setdiff(1:2,grep("RTransf",mtr$Lt_node[mtr$St_node==nodemoved]))]
				spnb<-c(lctr$tip.label, lctr$node.label)
				MOVEDspecies<-lctr$tip.label[Descendants(lctr, match(nodereallymoved,spnb), type="tips")[[1]]]
			}
			MOVEDspecies<-gsub("_0","",MOVEDspecies)
			## here we look at the topological distance between lctr and sptr trees. 
			## we remove any transfered species that does not change the topology of
			## the locustree
			lctrtmp<-drop.tip(lctr,grep("Rtransf",lctr$tip.label))
			lctrtmp$tip.label<-gsub("_0","",lctrtmp$tip.label)
			sizeofmastnbout<-Ntip(sptr)-Ntip(mast(lctrtmp, sptr))
			if (sizeofmastnbout==0) {
				MOVEDspecies<-NULL
			}
		}
		else MOVEDspecies<-NULL
		return(MOVEDspecies)
		print(MOVEDspecies)
	}
	TRUE_OUTLIERS<-sapply(1:length(LOCUTREES), function(x,l,s,m) GetGroupOfMovedSpecies(l[[x]],s,m[[x]]),l=LOCUTREES,s=alltrees$sptr, m=MATRICES)
	TRUE_OUTLIERS[unlist(lapply(TRUE_OUTLIERS, is.null))]<-"f"

	##KEEP TRANSFERS INVOLVING NO MORE THAN MAXNB SPECIES
	nspeciesintransfer<-unlist(lapply(TRUE_OUTLIERS, function(x) length(x[x!="f"])))
	# print(nspeciesintransfer)
	TRUE_OUTLIERS<-TRUE_OUTLIERS[nspeciesintransfer<=MAXNB]
	GENETREES<-GENETREES[nspeciesintransfer<=MAXNB]
	LOCUTREES<-LOCUTREES[nspeciesintransfer<=MAXNB]

	##SELECT THE DESIRED NUMBER OF TREES
	stopifnot(length(TRUE_OUTLIERS) >= NBTR)
	TRUE_OUTLIERS<-TRUE_OUTLIERS[1:NBTR]
	GENETREES<-GENETREES[1:NBTR]
	LOCUTREES<-LOCUTREES[1:NBTR]

	FORMATTED_OUTLIERS<-do.call(rbind, sapply(seq_along(TRUE_OUTLIERS), function(x,y) if(y[[x]][1]!="f") {cbind(as.character(x),y[[x]])},y=TRUE_OUTLIERS))	
	return(list(GENETREES=GENETREES, LOCUTREES=LOCUTREES, FORMATTED_OUTLIERS=FORMATTED_OUTLIERS))
}


## Compute precision and recall from phylter, treeshrink and true list of outliers
TruePositive<-function(o1,ref) { #o1 and o2 are two mat of outliers
	o1[,2]<-gsub("_0_0","",o1[,2])
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
	print(paste(TP, FP, FN, sep=" - "))
	precision<-TP/(TP+FP)
	if (is.na(precision)) precision<-0
	recall<-TP/(TP+FN)
	return(c(TP,FP,FN,precision,recall))
}

## Transform Carnivora species to ultrametric tree in generation time usable by simphy
preparetree<-function(sptree="../../Carnivora/Supermatrix_14463_genes_53_spp_UFBS_TESTNEW_Codon.tre") {
	if (is.numeric(sptree)) {
		sp<-rphylo(sptree, 1,0)
	}
	else {
		sp<-read.nexus(sptree)
	}
	spu<-chronos(sp) #see doc to see what it does
	class(spu)
	#then we change branch length so that the root to tip is 74 million years (age recovered from timetree)
	# multiplied by the mean nb of generation per year (7.62125)

	root2tipdist<-max(cophenetic(spu))/2
	multiplier<-74e6/root2tipdist
	spu$edge.length<-round(spu$edge.length*multiplier*7.62125)
	root<-Ntip(spu)+1
	print(root)
	alldist<-dist.nodes(spu)
	tip2root<-alldist[root,1:(root-1)]
	maxtip2root<-max(tip2root)
	diff2add<-maxtip2root-tip2root #what needs to be added to terminal branches to get EXACT ultrametricity
	spu$edge.length[match(1:(root-1),spu$edge[,2])]<-spu$edge.length[match(1:(root-1),spu$edge[,2])]+diff2add
	return(spu)
}

preparetree2<-function(sptree="../../Carnivora/Supermatrix_14463_genes_53_spp_UFBS_TESTNEW_Codon.tre") {
	if (is.numeric(sptree)) {
		sp<-rphylo(sptree, 1,0)
	}
	else {
		sp<-read.nexus(sptree)
	}
	spu<-chronos(sp) #see doc to see what it does
	class(spu)
	#then we change branch length so that the root to tip is 74 million years (age recovered from timetree)
	# multiplied by the mean nb of generation per year (7.62125)

	root2tipdist<-max(cophenetic(spu))/2
	root2tipgens<-74e6/8.315 # time/gen_time=gens
	multiplier<-root2tipgens/root2tipdist
	spu$edge.length<-round(spu$edge.length*multiplier)
	root<-Ntip(spu)+1
	print(root)
	alldist<-dist.nodes(spu)
	tip2root<-alldist[root,1:(root-1)]
	maxtip2root<-max(tip2root)
	diff2add<-maxtip2root-tip2root #what needs to be added to terminal branches to get EXACT ultrametricity
	spu$edge.length[match(1:(root-1),spu$edge[,2])]<-spu$edge.length[match(1:(root-1),spu$edge[,2])]+diff2add
	return(spu)
}


