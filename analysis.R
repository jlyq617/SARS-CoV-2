######PCA#####

suppressPackageStartupMessages(library(ade4))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(sva))

mapfile='./calu3.design.txt'
taxafile='./VSN-normalized.txt'
#taxafile='../04.normalization/Calu3.txt'


dat = read.csv(taxafile, 
               sep="\t",
               row.names = 1,
               check.names = F,
               comment.char = "",
               stringsAsFactors = F)


data1 = read.table(mapfile,header = T,stringsAsFactors = F,sep='\t',comment.char = "",check.names = F)

dat = dat[,data1[,1]]
dat[is.na(dat)]=0

dat<-apply(dat)
dat=apply(dat,2,function(x){x/sum(x)})

dat<-ComBat(dat=dat,batch = data1$group)


data = dat

data1[,2]=data1[,3]
colnames(data1)[2] ="Group"
data1$Group<-factor(data1$Group)


incw = 1.5
inlength = .1
inlift = "Axis1"
inright = "Axis2"
newpalette<-colorRampPalette(brewer.pal(9,"Set1"))(12)
newpalette<-c("#E41A1C", "#FF8200", "#FCC801", "#BCD401", "#58B53C", "#00835D",
              "#59C3E1","#004896","#920882","#E05695", "#EFAEB9","#999999")


data<-t(data)

data.dudi <- dudi.pca(data, center=TRUE, scale=T, scan=F, nf=10)
data2 <- data.dudi$li
incw = as.double(incw)
inlength = as.double(inlength)
classified_c = as.character(unique(data1[,2]))


#Adonis test
adonis1<-adonis(t(dat) ~ data1[,2],permutations = 999,method = "euclidean")


phenotype <- data1[,2]
f=classified_c
Type <- factor(phenotype,levels=f)
m = data.dudi$li
n = data.dudi$c1

lift_data = m[as.character(inlift)]
row.names(lift_data) = gsub(pattern = "[.]",replacement = "-",row.names(lift_data))
right_data = m[as.character(inright)]
row.names(right_data) = gsub(pattern = "[.]",replacement = "-",row.names(right_data))
data.dudi$li = cbind(lift_data,right_data)
num1 = substring(as.character(inlift),5,6)
num2 = substring(as.character(inright),5,6)
num1_data = n[paste("CS",num1,sep = '')]
num2_data = n[paste("CS",num2,sep = '')]
data.dudi$c1 = cbind(num1_data,num2_data)

right_data_box= cbind(data1[,2],right_data)
colnames(right_data_box)[1] = "Group" 
lift_data_box = cbind(data1[,2],lift_data)
colnames(lift_data_box)[1] = "Group" 

x1 <- min(data.dudi$li[,1]) - 0.3
y1 <- min(data.dudi$li[,2]) - 0.3
x2 <- max(data.dudi$li[,1]) + 0.3
y2 <- max(data.dudi$li[,2]) + 0.3
bb <- head(data.dudi$c1[order(sqrt((data.dudi$c1[,1])^2+(data.dudi$c1[,2])^2),decreasing=T),],n=7L)[1:7,]
rownames(bb) <- gsub("^X", "", rownames(bb))
rownames(bb) <- gsub("\\S+o__", "o__", rownames(bb))
cutoff <- (x2-0.3) / abs(bb[1,1]) * inlength
d2 <- data.frame(X=bb[1:dim(bb)[1],1]*cutoff, Y=bb[1:dim(bb)[1],2]*cutoff, LAB=rownames(bb)[1:dim(bb)[1]])
#d2[[3]]<- gsub('.*f__','f__',as.character(d2[[3]]))
d2[[3]] <- gsub('.*o__','o__',as.character(d2[[3]]))
d2[[3]] <- gsub('.*c__','c__',as.character(d2[[3]]))
d2[[3]] <- gsub('.*p__','p__',as.character(d2[[3]]))

d <- data.dudi$li
eig <- ( data.dudi$eig / sum( data.dudi$eig ) )

#track <- read.table("2.beta_div/plots/pca/2017.Jun11.Throat.Stool.16S.otu.Tonsil.prof.track", head=T, sep="	")
#points <- read.table("2.beta_div/plots/pca/2017.Jun11.Throat.Stool.16S.otu.Tonsil.prof.points", head=T, sep="	")

ggdata <- data.frame(d)

label_title <- paste("Calu-3 cell line","\n","P value = ",
                     round(adonis1[[1]][6][[1]][1],4)," R^2 = ",round(adonis1[[1]][5][[1]][1],4),sep="")

p<-ggplot(ggdata) +
  xlab("") +
  ylab("") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(aes(x=d[,1], y=d[,2], color=Type), size=3,alpha=0.75) +
  scale_color_manual(values=newpalette)+
  scale_fill_manual(values=newpalette)+
  guides(color=guide_legend(colnames(data1)[2]),
         fill=guide_legend(colnames(data1)[2]),
         shape=guide_legend(colnames(data1)[2]) ) +
  theme_classic()+
  theme(
    panel.background = element_rect(linetype=1,size=1,color="black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.background = element_blank(),
    legend.position=c(0.26,0.15),
    legend.title=element_text(
      #family="Arial",
      colour="black",
      size=6
    ),
    legend.text=element_text(
      #family="Arial",
      colour="black",
      size=10
    ))+ 
  xlim(min(d[,1], 0)*incw, max(d[,1])*incw)+
  ylim(min(d[,2], 0)*incw, max(d[,2])*incw)+
  geom_text(aes(x=-22,y=1,label=paste("PC",num1," (",round(eig[as.numeric(num1)]*100,2),"%)",sep="")))+
  geom_text(aes(x=-1,y=-15,label=paste("PC",num2," (",round(eig[as.numeric(num2)]*100,2),"%)",sep="")),angle=90)+
  guides(color=guide_legend(nrow=3,title=label_title,
                            title.theme=element_text(size=12,face = "bold"),
  ))



pdf('./Calu3.pdf',width = as.double(7), height = as.double(5))
plot(p)
dev.off()

#GO enrichmenti#####
dat<-read.csv('./allproteins_analysis.txt',stringsAsFactors = F,
              check.names = F,sep='\t',header = T,skip = 11)
dat<-dat[order(dat$`upload_1 (FDR)`),]
plotdata<-dat[1:15,]
plotdata<-plotdata[c(1,2,6,7,13),]
plotdata$`GO cellular component complete`<-factor(plotdata[,1],levels=plotdata[,1])
plotdata$term<-c("RNP complex","Membrane","Nucleus","Nucleoplasm","Spliceosome")
plotdata$term<-factor(plotdata$term,levels=rev(plotdata$term))

library(ggplot2)

pdf('./GO-CP.pdf',width = 4,height = 2.5)
ggplot(plotdata,aes(x=term,y=-log10(`upload_1 (FDR)`)))+
  geom_bar(stat="identity",fill="#9A7CA3",width=0.6)+
  theme_bw()+
  coord_flip()+
  xlab("")+
  ylab("-log10(FDR adjusted p-value)")+
  theme(
    panel.grid.major  = element_line(size=0.5),
    panel.grid.minor  = element_line(size=0.5)
    
  )
dev.off()

#Diff analysis######
library(fdrtool)
mapfile='../05.PCA/calu3.design.txt'
taxafile='../05.PCA/VSN-normalized-combat.txt'

dat = read.csv(taxafile, 
               sep="\t",
               row.names = 1,
               check.names = F,
               stringsAsFactors = F)

data1 = read.table(mapfile,header = T,stringsAsFactors = F,sep='\t',comment.char = "",check.names = F)

dat = dat[,data1[,1]]
dat = t(dat)

difftaxa<-function(x,y){
  result = matrix(rep("",ncol(x)*5),ncol(x),5)
  colnames(result)<-c('protein',substitute(x),substitute(y),'pvalue','Enrich')
  for(i in c(1:ncol(x))){
    t<-wilcox.test(x[,i],y[,i],paired=F)
    
    result[i,1]<-colnames(x)[i]
    result[i,2]<-as.double(mean(x[,i],na.rm = T))
    result[i,3]<-as.double(mean(y[,i],na.rm = T))
    result[i,4]<-t$p.value
    abundance<-result[i,2]
    type2<-colnames(result)[2]
    if(as.double(result[i,3])>as.double(abundance)){
      abundance=result[i,3]
      type2<-colnames(result)[3]
    }
    result[i,5]=type2
    print(i)
  }
  return(result)
}

NRC1<-dat[data1[data1$treatment=="NRC1",]$prep,]
NRC2<-dat[data1[data1$treatment=="NRC2",]$prep,]
NRC3<-dat[data1[data1$treatment=="NRC3",]$prep,]
NRC4<-dat[data1[data1$treatment=="NRC4",]$prep,]
NRC5<-dat[data1[data1$treatment=="NRC5",]$prep,]
NRC6<-dat[data1[data1$treatment=="NRC6",]$prep,]
NRC7<-dat[data1[data1$treatment=="NRC7",]$prep,]
NRC8<-dat[data1[data1$treatment=="NRC8",]$prep,]
NRC9<-dat[data1[data1$treatment=="NRC9",]$prep,]
NRC10<-dat[data1[data1$treatment=="NRC10",]$prep,]
NRC11<-dat[data1[data1$treatment=="NRC11",]$prep,]
NRC12<-dat[data1[data1$treatment=="NRC12",]$prep,]


result1<-difftaxa(NRC1,NRC12)

result1<-data.frame(result1,check.names = F,stringsAsFactors = F)

result1$pvalue<-as.numeric(result1$pvalue)
result1$qvalue<-fdrtool(result1$pvalue,statistic = c("pvalue"),plot=F)$qval

result1<-result1[order(result1$qvalue),]


filename='NRC12.txt'

write.table(result1,filename,quote=F,sep='\t')

#WGCNA#####
library(WGCNA)
library(dplyr)

#dat<-read.table('~/Downloads/Projects/COVID-19-3Cellline/Calu3/05.PCA/VSN-normalized.combat.txt',
#                 stringsAsFactors = F,check.names = F,sep='\t',header=T,row.names = 1)

dat<-read.table('~/Downloads/Projects/COVID-19-3Cellline/Calu3/05.PCA/Calu3-combat.txt',
                stringsAsFactors = F,check.names = F,sep='\t',header=T,row.names = 1)

de<-read.table('~/Downloads/Projects/COVID-19-3Cellline/Calu3/05.PCA/calu3.design.txt',
               stringsAsFactors = F,check.names = F,sep='\t',header=T )



#prepare 
datExpr<-dat%>%
  .[,de$prep]%>%
  t()%>%
  data.frame(stringsAsFactors = F,check.names = F)

SubGeneNames=colnames(datExpr)

#Choosing a soft-threshold to fit a scale-free topology to the network
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft<-pickSoftThreshold(datExpr,
                       dataIsExpr = TRUE,
                       powerVector = powers,
                       corFnc = cor,
                       corOptions = list(use = 'p'),
                       networkType = "unsigned")

# Plot the results
pdf("./softpower.pdf",width = 9, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

#Generating adjacency and TOM similarity matrices based on the selected softpower
softPower = sft$powerEstimate
softPower
#softPower = 7

net=blockwiseModules(datExpr,power=softPower,
                     maxBlockSize = 5000,
                     TOMType = "signed", 
                     minModuleSize = 20,
                     reassignThreshold = 0, 
                     mergeCutHeight = 0.25,
                     numericLabels = TRUE, 
                     pamRespectsDendro = FALSE,
                     saveTOMs = TRUE,
                     saveTOMFileBase = "TOM",
                     networkType = "signed hybrid",
                     verbose = 3)

table(net$colors)

# Convert labels to colors for plotting
dynamicColors  = labels2colors(net$colors)
table(dynamicColors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], dynamicColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#remove grey
restGenes=(dynamicColors!="grey")
diss1=1-TOMsimilarityFromExpr(datExpr[,restGenes], power = sft$powerEstimate)
hier1=hclust(as.dist(diss1), method="average" )
plotDendroAndColors(hier1, dynamicColors[restGenes],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


module.order <- unlist(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))
m<-datExpr[,module.order]

#plot heatmap to show pattern between samples
heatmap(t(m),
        col =  blueWhiteRed(50),
        Rowv=NA,
        Colv=NA,
        labRow=NA,
        scale="row",
        RowSideColors=dynamicColors[module.order])

#similarity plot
nGenes<-ncol(datExpr)
nSamples<-nrow(datExpr)

datExpr_tree<-hclust(dist(datExpr),method="average")
plot(datExpr_tree, main = "Sample clustering")

geneTree=net$dendrograms[[1]]

colnames(diss1) =rownames(diss1) =SubGeneNames[restGenes]

diag(diss1) = NA
TOMplot(diss1, hier1, as.character(dynamicColors[restGenes]))

#design table 
design<-model.matrix(~0+de$treatment)
rownames(design)=de$prep
colnames(design)=gsub("de\\$treatment","",colnames(design))
design<-design[,unique(de$treatment)]
dynamicColors  = labels2colors(net$colors)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, dynamicColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design,method = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
#sizeGrWindow(10,6)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

pdf("moudle.correlation.pdf",width = 10,height = 6)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors =  blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#extract gene list from different modules
module_colors= unique(dynamicColors)
for (color in module_colors){
  module=SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

NM=as.data.frame(cor(datExpr,MEs,use="p"))


save.image(file="WGCNA.RData")

#
ADJ1=abs(cor(datExpr,use="p"))^6
colorh1=dynamicColors
Alldegress1=intramodularConnectivity(ADJ1,colorh1)
write.table(Alldegress1,"alldegree.txt",quote=F,sep='\t')

#
