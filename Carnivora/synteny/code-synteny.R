library(ggplot2)

####  Compute pval of synteny outlier scores for ALL pairs of species

sps<-c("Ailuropoda_melanoleuca", "Callorhinus_ursinus", "Canis_familiaris", "Eumetopias_jubatus", "Felis_catus", "Mustela_putorius", "Neomonachus_schauinslandi", "Odobenus_rosmarus", "Panthera_pardus", "Panthera_tigris", "Puma_concolor", "Suricata_suricatta", "Ursus_maritimus", "Zalophus_californianus")



# this function computes all pairwise simple and double breaks and the associated pvalues
# for all possible pairs of species (in both direction)
ComputePairwisePvalue<-function(outfile) {
    RESU<-NULL
    NAM<-NULL
    for (i in 1:length(sps)) {
        for (j in sps[-i]) {
            sp1<-sps[i]
            sp2<-j
            nam<-paste(sp1,sp2,sep=" - ")
            NAM<-c(NAM, nam)
            print(nam)
            run<-paste("python3 test_phylter_synteny_3.py genomes",outfile,sp1,sp2, sep=" ")
            resu<-system(run, intern=TRUE)
            RESU<-rbind(RESU, as.numeric(strsplit(resu[10],"\t")[[1]]))
        }
    }
    colnames(RESU)<-c("markers","nb.markers.phyltered","breaks","dbl.breaks","breaks.remain", "dbl.breaks.remain","pval.breaks","pval.dbl.breaks")
    rownames(RESU)<-NAM
    RESU<-as.data.frame(RESU)
    RESU$breaks.phyltered<-RESU$breaks-RESU$breaks.remain
    RESU$dbl.breaks.phyltered<-RESU$dbl.breaks-RESU$dbl.breaks.remain
    return(RESU)
}

# this function prepares data for plotting a heatmap containing the pvalues associated to the number of outliers associated to single or double breaks.
GetDataFrameForHeatmapAndPlot<-function(RESU, what="breaks") {
    sp1<-unlist(lapply(strsplit(rownames(RESU)," - "), function(x) x[1]))
    sp2<-unlist(lapply(strsplit(rownames(RESU)," - "), function(x) x[2]))
    if (what=="breaks") DF<-data.frame(sp1=sp1, sp2=sp2,pvalue=RESU$pval.breaks)
    if (what=="dblbreaks") DF<-data.frame(sp1=sp2, sp2=sp1,pvalue=RESU$pval.dbl.breaks)

    p<-ggplot(DF, aes(x = sp1, y = sp2, fill = pvalue)) + geom_tile(color="black") + geom_text(aes(label = pvalue), color = "white", size = 4) + coord_fixed() +  theme(axis.text.x = element_text(angle = 90))
    p
}

R1.phy<-ComputePairwisePvalue("../data/phylter-results/file1_v2.txt") #small
R2.phy<-ComputePairwisePvalue("../data/phylter-results/file_k1.55_v2.txt") #large
R1.ts<-ComputePairwisePvalue("../data/treeshrink-results/treeshrink0.012.txt") #small
R2.ts<-ComputePairwisePvalue("../data/treeshrink-results/treeshrink0.05.txt") #large

p1<-GetDataFrameForHeatmapAndPlot(R1.phy,"breaks")
p2<-GetDataFrameForHeatmapAndPlot(R1.ts,"breaks")
p3<-GetDataFrameForHeatmapAndPlot(R1.phy,"dblbreaks")
p4<-GetDataFrameForHeatmapAndPlot(R1.ts,"dblbreaks")
p5<-GetDataFrameForHeatmapAndPlot(R2.phy,"breaks")
p6<-GetDataFrameForHeatmapAndPlot(R2.ts,"breaks")
p7<-GetDataFrameForHeatmapAndPlot(R2.phy,"dblbreaks")
p8<-GetDataFrameForHeatmapAndPlot(R2.ts,"dblbreaks")

require("ggpubr")
#THEPLOT<-ggarrange(p1, p2, p3,p4, p5, p6, p7, p8, labels=c("A","B","C","D","E","F","G","H"),ncol = 2, nrow = 4, common.legend = TRUE, legend = "bottom") + theme(legend.title = element_blank())
THEPLOT.PHY<-ggarrange(p1, p3, p5,p7, labels=c("A","B","C","D"), font.label = list(size = 30, color = "black"), ncol = 2, nrow = 2, common.legend = TRUE, legend = "right") + theme(legend.title = element_blank())
THEPLOT.TS<-ggarrange(p2, p4, p6,p8, labels=c("A","B","C","D"), font.label = list(size = 30, color = "black"), ncol = 2, nrow = 2, common.legend = TRUE, legend = "right") + theme(legend.title = element_blank())
ggsave("HEATMAPRESULTSYNTENY.PHYLTER.pdf",THEPLOT.PHY,width=16, height=17, dpi=600)
ggsave("HEATMAPRESULTSYNTENY.TREESHRINK.pdf",THEPLOT.TS,width=16, height=17, dpi=600)




# This function produces the circular plot presented in the manuscript where 
# orthologous genes of all genes in one scaffold of one species are connected with links
# and those that are outliers according to phylter are in red (others in blue)
require(circlize)
require(scales)

plotit<-function(sp1="Suricata_suricatta", sp2="Zalophus_californianus", i, outfile="../data/phylter-results/file1_v2.txt") {
    run<-paste("python3 test_phylter_synteny_3.py genomes",outfile,sp1,sp2, "1", sep=" ") #"1" is for 1 replicate.
    resu<-system(run, intern=TRUE)
    fi<-paste(sp1,"-",sp2,".markers", sep="")
    ok<-read.table(fi)

    scaff<-unique(ok$V7)[i]
    sub<-ok[ok$V7==scaff,]
    scaffs2<-unique(sub$V2)
    group1<-which(!is.na(match(ok$V2,c(scaff,scaffs2))))
    group2<-which(!is.na(match(ok$V7,c(scaff,scaffs2))))
    sub0<-ok[unique(c(group1, group2)),]
    sub0<-sub0[sub0$V7!=scaff,]
    DF<-data.frame(sectors=c(rep(sub$V2,2),rep(sub$V7,2), rep(sub0$V2,2)), x=c(sub$V3,sub$V4,sub$V8,sub$V9, sub0$V3, sub0$V4), species=c(rep(sub$V1,2),rep(sub$V6,2), rep(sub0$V1,2)))
    circos.par(cell.padding = c(0.02, 0, 0.02, 0))
    circos.initialize(DF$sectors, x=DF$x)
    circos.track(DF$sectors,ylim=c(0,0.5), track.height = 0.05, bg.col = c(rep("#bebebe",length(unique(DF$sectors))-1),"#3d6282"),bg.border=NA,
    panel.fun = function(x, y) {
           circos.text(CELL_META$xcenter, 
               CELL_META$cell.ylim[2] + mm_y(2), 
               CELL_META$sector.index, cex=0.5)
    #       circos.axis(labels.cex = 0.4, major.at=CELL_META$xlim[1]:CELL_META$xlim[2])

            }
        )
    apply(sub, 1, function(x) circos.link(x[2],as.integer(x[3:4]), x[7],as.integer(x[8:9]), col=alpha(ifelse(x[11],"#D81B60","#1E88E5"),0.5), lwd=3))

    out<-sub[which(sub$V11=="True"),c(2,3,4)]
    out$mid<-out$V3+(out$V4-out$V3)/2
    #add stars to show outliers
    apply(out, 1, function(x) circos.text(as.integer(x[4]), 0.5 + mm_y(0.5), labels="*", x[1], cex=0.8, col="#999999"))
}


pdf("circplot.pdf")
par(mfrow=c(2,2))
for (i in c(1,2,3,5)) plotit(i=i)
dev.off()
