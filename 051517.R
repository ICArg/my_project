pepe
pr1=names(gene.symbs)[jj] #probes
pr1
View(asd)
genes2look=c(GeneList)
jj=which(gene.symbs%in%genes2look)
pr1=names(gene.symbs)[jj] #probes
gns=gene.symbs[jj]        #genes
all.dat.genes = all.dat[pr1,]
all.dat.genes.df = data.frame(all.dat.genes)
all.dat.genes.df["GeneNames"] = gene.symbs[jj]
all.dat.ddply = ddply(all.dat.genes.df, .(GeneNames), colwise(mean))
rownames(all.dat.ddply) = all.dat.ddply$GeneNames
all.dat.ddply$GeneNames = factor(all.dat.ddply$GeneNames, levels = GeneOrder)
all.dat2 = all.dat.ddply[order(all.dat.ddply$GeneNames),]
tall.dat2 = setNames(data.frame(t(all.dat2[,-1])), all.dat2[,1])
cor.dat = cor(tall.dat2)
cor.datm = melt(cor.dat)
all.dat.matrix = as.matrix(all.dat.ddply)
# Get appropiate gene lists
Genes = read.csv("GeneList.csv", stringsAsFactors = FALSE, header = TRUE, na.strings = "")
GlycList = as.character(c(Genes$GlyGenes))
GlycList = na.omit(GlycList)
ImmList = as.character(c(Genes$ImmGenes))
ImmList = na.omit(ImmList)
ConList = as.character(c(Genes$ConGenes))
ConList = na.omit(ConList)
GeneList = c(GlycList, ImmList, ConList)
all.dat.ddply$GeneNames = factor(all.dat.ddply$GeneNames, levels = GeneList)
all.dat2 = all.dat.ddply[order(all.dat.ddply$GeneNames),]
tall.dat2 = setNames(data.frame(t(all.dat2[,-1])), all.dat2[,1])
cor.dat = cor(tall.dat2)
cor.datm = melt(cor.dat)
all.dat.matrix = as.matrix(all.dat.ddply)
View(cor.dat)
# Correlation plot of genes ordered by manual list - using ggplot. input is a melted correlation matrix.
ggplot_corr = function(x) {
setwd("C:/Users/coheni/Documents/")
GGPlotObject = ggplot(x, aes(x = Var1, y = Var2, fill = value)) +
geom_tile() +
coord_equal() +
xlab("") +
ylab("") +
scale_fill_gradient2(low = "red", high = "blue", mid = "white",
midpoint = 0, limit = c(-1,1), space = "Lab",
name="Pearson\nCorrelation") +
theme(text = element_text(size = 20),
axis.text = element_text(size = 30),
#legend.text = element_text(size=20),
axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
scale_y_discrete(limits = rev(levels(x$Var1)))
GGDrawObject = ggdraw(switch_axis_position(GGPlotObject, axis = 'x'))
return(GGDrawObject)
}
# Heatmap of patient expression data using ggplot. Input is a melted patient matrix.
ggplot_heatmap = function(x, y) {
setwd("C:/Users/coheni/Documents/")
GGPlotObject = ggplot(x, aes(x = Var2, y = Var1, fill = value)) +
geom_tile() +
xlab("") +
ylab("") +
scale_fill_gradient2(low = "red", high = "blue", mid = "white",
midpoint = 0, limit = c(-5,5), space = "Lab",
name="Pearson\nCorrelation") +
theme(text = element_text(size = 20),
axis.text = element_text(size = 30),
#legend.text = element_text(size=20),
axis.text.x = element_blank()) +
scale_x_discrete(limits = y) +
scale_y_discrete(limits = rev(GeneList))
return(GGPlotObject)
}
ggplot_corr(cor.datm)
setwd("C:/Users/coheni/Documents/TCGA")
# Get appropiate gene lists
Genes = read.csv("GeneList.csv", stringsAsFactors = FALSE, header = TRUE, na.strings = "")
GlycList = as.character(c(Genes$GlyGenes))
GlycList = na.omit(GlycList)
ImmList = as.character(c(Genes$ImmGenes))
ImmList = na.omit(ImmList)
ConList = as.character(c(Genes$ConGenes))
ConList = na.omit(ConList)
GeneList = c(GlycList, ImmList, ConList)
genes2look=c(GeneList)
jj=which(gene.symbs%in%genes2look)
pr1=names(gene.symbs)[jj] #probes
gns=gene.symbs[jj]        #genes
gns
pr1
gns["CD8A"]
which(gns%in%"CD8A")
which(gns%in%"CD8B")
pepe = which(gns%in%"CD8B")
pepe
pepe[[1]]
pepe[1]
all.dat.genes = all.dat[pr1,]
all.dat.genes.df = data.frame(all.dat.genes)
all.dat.genes.df["GeneNames"] = gene.symbs[jj]
all.dat.ddply = ddply(all.dat.genes.df, .(GeneNames), colwise(mean))
rownames(all.dat.ddply) = all.dat.ddply$GeneNames
all.dat.ddply$GeneNames = factor(all.dat.ddply$GeneNames, levels = GeneList)
all.dat2 = all.dat.ddply[order(all.dat.ddply$GeneNames),]
tall.dat2 = setNames(data.frame(t(all.dat2[,-1])), all.dat2[,1])
cor.dat = cor(tall.dat2)
cor.datm = melt(cor.dat)
all.dat.matrix = as.matrix(all.dat.ddply)
View(cor.dat)
GeneList
ggplot_corr(cor.datm)
all.dat.genes.df
all.dat.genes.df$GeneNames
length(all.dat.genes.df$GeneNames)
length(all.dat2$GeneNames)
all.dat2
all.dat.genes.df$GeneNames
all.dat.genes = all.datERn[pr1,]
all.dat.genes.df = data.frame(all.dat.genes)
all.dat.genes.df["GeneNames"] = gene.symbs[jj]
all.dat.ddply = ddply(all.dat.genes.df, .(GeneNames), colwise(mean))
rownames(all.dat.ddply) = all.dat.ddply$GeneNames
all.dat.ddply$GeneNames = factor(all.dat.ddply$GeneNames, levels = GeneList)
all.dat2 = all.dat.ddply[order(all.dat.ddply$GeneNames),]
tall.dat2 = setNames(data.frame(t(all.dat2[,-1])), all.dat2[,1])
cor.dat = cor(tall.dat2)
cor.datm = melt(cor.dat)
all.dat.matrix = as.matrix(all.dat.ddply)
ggplot_corr(cor.datm)
all.dat.genes = all.datERp[pr1,]
all.dat.genes.df = data.frame(all.dat.genes)
all.dat.genes.df["GeneNames"] = gene.symbs[jj]
all.dat.ddply = ddply(all.dat.genes.df, .(GeneNames), colwise(mean))
rownames(all.dat.ddply) = all.dat.ddply$GeneNames
all.dat.ddply$GeneNames = factor(all.dat.ddply$GeneNames, levels = GeneList)
all.dat2 = all.dat.ddply[order(all.dat.ddply$GeneNames),]
tall.dat2 = setNames(data.frame(t(all.dat2[,-1])), all.dat2[,1])
cor.dat = cor(tall.dat2)
cor.datm = melt(cor.dat)
all.dat.matrix = as.matrix(all.dat.ddply)
ggplot_corr(cor.datm)
all.dat.genes = all.dat[pr1,]
all.dat.genes = all.dat[pr1,]
all.dat.genes.df = data.frame(all.dat.genes)
all.dat.genes.df["GeneNames"] = gene.symbs[jj]
all.dat.ddply = ddply(all.dat.genes.df, .(GeneNames), colwise(mean))
rownames(all.dat.ddply) = all.dat.ddply$GeneNames
all.dat.ddply$GeneNames = factor(all.dat.ddply$GeneNames, levels = GeneList)
all.dat2 = all.dat.ddply[order(all.dat.ddply$GeneNames),]
tall.dat2 = setNames(data.frame(t(all.dat2[,-1])), all.dat2[,1])
cor.dat = cor(tall.dat2)
cor.datm = melt(cor.dat)
all.dat.matrix = as.matrix(all.dat.ddply)
####
GlySig = c(GlycList)
ImmSig = c(ImmList)
BRCAOrderList = get_order(all.dat.matrix, GlySig, ImmSig) # Get pt list ordered by Gly or Imm Sig score (mean)
BRCAOrderGly = BRCAOrderList[[1]] # List of all pts ordered by Gly
BRCAOrderImm = BRCAOrderList[[2]]
BRCAOrderList.Gly = get_thirds(BRCAOrderList)
BRCAOrderList.Imm = get_thirds(BRCAOrderList)
BRCAGlyHigh = BRCAOrderList.Gly[[1]] # List of pts in the "Gly-High"
BRCAGlyMid = BRCAOrderList.Gly[[2]]
BRCAGlyLow = BRCAOrderList.Gly[[3]]
BRCAImmHigh = BRCAOrderList.Imm[[1]]
BRCAImmMid = BRCAOrderList.Imm[[2]]
BRCAImmLow = BRCAOrderList.Imm[[3]]
# This FUN reads tables. Just use CAPS for 'x' and no-caps for 'y'.
read_table = function(x, y) {
setwd("C:/Users/coheni/Documents/TCGA")
dire = paste("TCGA-", x, "/Test/", y, "/tcga/data_RNA_Seq_v2_mRNA_median_Zscores.txt", sep = "")
Table = read.table(dire, stringsAsFactors = FALSE, header = TRUE)
return(Table)
}
# This FUN reads tables. Just use CAPS for 'x' and no-caps for 'y'.
read_table_exp = function(x, y) {
setwd("C:/Users/coheni/Documents/TCGA")
dire = paste("TCGA-", x, "/Test/", y, "/tcga/data_RNA_Seq_v2_expression_median.txt", sep = "")
Table = read.table(dire, stringsAsFactors = FALSE, header = TRUE)
return(Table)
}
# Makes all kinds of tables, matrices, and correlations.
# Input must be the original table from the read_table/original file.
make_tables = function(x) {
DisDF = as.data.frame(x)
rownames(DisDF) = DisDF$Hugo_Symbol
DisDF = na.omit(DisDF)
tDisDF = setNames(data.frame(t(DisDF[,-1])), DisDF[,1])
DisGOI = DisDF[c(GlycList, ImmList, ConList),]
DisGOI = na.omit(DisGOI)
DisGOI$Entrez_Gene_Id = NULL
tDisGOI = setNames(data.frame(t(DisGOI[,-1])), DisGOI[,1])
DisMatrix = data.matrix(tDisGOI)
tDisMatrix = data.matrix(DisGOI[,-1])
DisMatrixm = melt(DisMatrix)
tDisMatrixm = melt(tDisMatrix)
DisCorr = cor(tDisGOI)
DisCorrm = melt(DisCorr)
FinalTables = list(DisMatrix, tDisMatrix, DisMatrixm, tDisMatrixm, DisCorr, DisCorrm)
return(FinalTables)
}
# Get order of patients according to Gly or Imm. Input here is the DisMatrix, as well as Gly or Imm signature genes
get_order = function(x, y, z) {
DisGlyMatrix = x[,c(y)]
DisImmMatrix = x[,c(z)]
DisGlyMeans = rowMeans(DisGlyMatrix)
DisImmMeans = rowMeans(DisImmMatrix)
DisRowMeans = data.frame(cbind(DisGlyMeans, DisImmMeans))
DisGlyOrderNumber = rev(order(DisRowMeans$DisGlyMeans))
DisImmOrderNumber = rev(order(DisRowMeans$DisImmMeans))
DisGlyOrder = rownames(DisRowMeans[DisGlyOrderNumber,])
DisImmOrder = rownames(DisRowMeans[DisImmOrderNumber,])
OrderList = list(DisGlyOrder, DisImmOrder)
return(OrderList)
}
# Input a patient order (or any list of things I guess), and divide it into thirds.
get_thirds = function(x) {
High = x[1:round(length(x)/3)]
Mid = x[(1+round(length(x)/3)):(2*round(length(x)/3))]
Low = x[(2*round(length(x)/3)):length(x)]
ThirdsList = list(High, Mid, Low)
return(ThirdsList)
}
#FUNMatrix = head(BRCARowMeansDF)
#rownames(FUNMatrix2) = rownames(FUNMatrix)
#FUNMatrix = FUNMatrix[order(FUNMatrix$BRCAGlyMeans),]
####
GlySig = c(GlycList)
ImmSig = c(ImmList)
BRCAOrderList = get_order(all.dat.matrix, GlySig, ImmSig) # Get pt list ordered by Gly or Imm Sig score (mean)
BRCAOrderGly = BRCAOrderList[[1]] # List of all pts ordered by Gly
BRCAOrderImm = BRCAOrderList[[2]]
BRCAOrderList.Gly = get_thirds(BRCAOrderList)
BRCAOrderList.Imm = get_thirds(BRCAOrderList)
BRCAGlyHigh = BRCAOrderList.Gly[[1]] # List of pts in the "Gly-High"
BRCAGlyMid = BRCAOrderList.Gly[[2]]
BRCAGlyLow = BRCAOrderList.Gly[[3]]
BRCAImmHigh = BRCAOrderList.Imm[[1]]
BRCAImmMid = BRCAOrderList.Imm[[2]]
BRCAImmLow = BRCAOrderList.Imm[[3]]
GlySig
all.dat.matrix["HK1"]
all.dat.matrix
all.dat.matrix["HK1",]
all.dat.matrix
####
BRCATable = read_table("GBM", "gbm") # Make original table from file (GEM)
BRCATableList = make_tables(BRCATable) # Make a list of 6 tables/DFs/Matrices - see function for details
BRCAMatrix = BRCATableList[[1]]
tBRCAMatrix = BRCATableList[[2]]
BRCAMatrix.m = BRCATableList[[3]]
tBRCAMatrix.m = BRCATableList[[4]]
BRCACorr = BRCATableList[[5]]
BRCACorr.m = BRCATableList[[6]]
BRCAMatrix
all.dat.genes = all.dat[pr1,]
all.dat.genes.df = data.frame(all.dat.genes)
all.dat.genes.df["GeneNames"] = gene.symbs[jj]
all.dat.ddply = ddply(all.dat.genes.df, .(GeneNames), colwise(mean))
rownames(all.dat.ddply) = all.dat.ddply$GeneNames
all.dat.ddply$GeneNames = factor(all.dat.ddply$GeneNames, levels = GeneList)
all.dat2 = all.dat.ddply[order(all.dat.ddply$GeneNames),]
tall.dat2 = setNames(data.frame(t(all.dat2[,-1])), all.dat2[,1])
cor.dat = cor(tall.dat2)
cor.datm = melt(cor.dat)
all.dat.matrix = as.matrix(all.dat.ddply)
t.all.dat.matrix = setNames(data.frame(t(all.dat.matrix[,-1])), all.dat.matrix[,1])
t.all.dat.matrix
####
GlySig = c(GlycList)
ImmSig = c(ImmList)
BRCAOrderList = get_order(all.dat.matrix, GlySig, ImmSig) # Get pt list ordered by Gly or Imm Sig score (mean)
asd = t.all.dat.matrix[1:10,1:10]
asd
asd["ACTB"]
asd["ALDOA"]
asd["CD8A"]
asd[,"ALDOA"]
asd[,c("ACTB", "ALDOA")]
asd = t.all.dat.matrix[1:10,1:10]
asd2 = asd[,c("ACTB", "ALDOA")]
asd2
asd2 = asd[,c("ACTB", "ALDOA", "CD8A")]
asd2
asd3 = asd[,c("ACTB", "ALDOA", "CD8A")]
asd3
GlySig
all.dat.genes.df$GeneNames
all.dat.ddply$GSM50119
all.dat.ddply$GeneNames
asd2 = GlySig[-1]
asd2
ImmSig
which(ImmSig%in%"CD8B")
ImmSig = ImmSig[-5]
ImmSig
all.dat.genes = all.dat[pr1,]
all.dat.genes.df = data.frame(all.dat.genes)
all.dat.genes.df["GeneNames"] = gene.symbs[jj]
all.dat.ddply = ddply(all.dat.genes.df, .(GeneNames), colwise(mean))
rownames(all.dat.ddply) = all.dat.ddply$GeneNames
all.dat.ddply$GeneNames = factor(all.dat.ddply$GeneNames, levels = GeneList)
all.dat2 = all.dat.ddply[order(all.dat.ddply$GeneNames),]
tall.dat2 = setNames(data.frame(t(all.dat2[,-1])), all.dat2[,1])
cor.dat = cor(tall.dat2)
cor.datm = melt(cor.dat)
all.dat.matrix = as.matrix(all.dat.ddply)
t.all.dat.matrix = setNames(data.frame(t(all.dat.matrix[,-1])), all.dat.matrix[,1])
BRCAOrderList = get_order(all.dat.matrix, GlySig, ImmSig) # Get pt list ordered by Gly or Imm Sig score (mean)
BRCAOrderList = get_order(t.all.dat.matrix, GlySig, ImmSig) # Get pt list ordered by Gly or Imm Sig score (mean)
asd = t.all.dat.matrix[1:10,1:10]
asd
asd2 = BRCAMatrix[1:10,1:10]
asd2
str(asd)
str(asd2)
asd3 = all.dat.matrix[1:10,1:10]
asd3
asd4 = all.dat.ddply[1:10,1:10]
asd4
all.dat2
all.dat
asd = all.dat[1:10,1:10]
asd
asd2 = all.dat2[1:10,1:10]
asd2
asd = tall.dat2[1:10,1:10]
asd
BRCAOrderList = get_order(tall.dat2, GlySig, ImmSig) # Get pt list ordered by Gly or Imm Sig score (mean)
BRCAOrderList = get_order(tall.dat2, GlySig, ImmSig) # Get pt list ordered by Gly or Imm Sig score (mean)
BRCAOrderGly = BRCAOrderList[[1]] # List of all pts ordered by Gly
BRCAOrderImm = BRCAOrderList[[2]]
BRCAOrderList.Gly = get_thirds(BRCAOrderList)
BRCAOrderList.Imm = get_thirds(BRCAOrderList)
BRCAGlyHigh = BRCAOrderList.Gly[[1]] # List of pts in the "Gly-High"
BRCAGlyMid = BRCAOrderList.Gly[[2]]
BRCAGlyLow = BRCAOrderList.Gly[[3]]
BRCAImmHigh = BRCAOrderList.Imm[[1]]
BRCAImmMid = BRCAOrderList.Imm[[2]]
BRCAImmLow = BRCAOrderList.Imm[[3]]
ggplot_heatmap(BRCAMatrix, BRCAOrderList.Gly)
asd = tBRCAMatrix[1:10,1:10]
asd
str(tBRCAMatrix)
asd2 = all.dat2[1:10,1:10]
asd2
all.dat2$GeneNames = NULL
asd2 = all.dat2[1:10,1:10]
asd2
str(asd2)
asd3 = melt(asd2)
asd3
View(asd2)
View(asd)
asd2 = as.numeric(all.dat2[1:10,1:10])
asd3
asd2
str(asd2)
asd3 = matrix(asd2)
str(asd3)
asd4 = melt(asd3)
asd4
rownames(asd4 = rownames(asd2))
rownames(asd4) = rownames(asd2))
rownames(asd4) = rownames(asd2)
colnames(asd4) = colnames(asd2)
asd4
asd5 = melt(asd4)
gse = getGEO('GSE20624')
edata = exprs(gse[[1]])
pdata = pData(gse[[1]])
plot(density(log2(edata+1)))
plot(density(edata))
library(GEOquery)
gse = getGEO('GSE20624')
edata = exprs(gse[[1]])
pdata = pData(gse[[1]])
plot(density(log2(edata+1)))
plot(density(edata))
plot(density(na.omit(edata)))
plot(density(asd))
asd = na.omit(edata)
plot(density(asd))
View(pdata)
gse = getGEO('GSE3143')
edata = exprs(gse[[1]])
pdata = pData(gse[[1]])
plot(density(log2(edata+1)))
plot(density(edata))
asd = na.omit(edata)
plot(density(asd))
View(pdata)
gse = getGEO('GSE20624')
edata = exprs(gse[[1]])
pdata = pData(gse[[1]])
plot(density(log2(edata+1)))
plot(density(edata))
asd = na.omit(edata)
plot(density(asd))
View(pdata)
gse
summary(pdata$characteristics_ch2.1)
gse = getGEO('GSE22133')
edata = exprs(gse[[1]])
pdata = pData(gse[[1]])
plot(density(log2(edata+1)))
plot(density(edata))
asd = na.omit(edata)
plot(density(asd))
gse
View(pdata)
gse = getGEO('GSE10886')
edata = exprs(gse[[1]])
pdata = pData(gse[[1]])
plot(density(log2(edata+1)))
plot(density(edata))
asd = na.omit(edata)
plot(density(asd))
View(pdata)
gse
edata = exprs(gse[[3]])
pdata = pData(gse[[3]])
View(pdata)
gse = getGEO('GSE18229')
edata = exprs(gse[[1]])
pdata = pData(gse[[1]])
plot(density(log2(edata+1)))
plot(density(edata))
asd = na.omit(edata)
plot(density(asd))
gse
gse = getGEO('GSE87411')
edata = exprs(gse[[1]])
pdata = pData(gse[[1]])
plot(density(log2(edata+1)))
plot(density(edata))
asd = na.omit(edata)
plot(density(asd))
View(asd)
View(edata)
gse
View(pdata)
gse = getGEO('GSE80999')
edata = exprs(gse[[1]])
pdata = pData(gse[[1]])
plot(density(log2(edata+1)))
plot(density(edata))
asd = na.omit(edata)
plot(density(asd))
gse
View(pdata)
summary(pdata$characteristics_ch1.1)
View(asd)
View(pdata)
summary(pdata$characteristics_ch1.2)
gse = getGEO('GSE81000')
edata = exprs(gse[[1]])
pdata = pData(gse[[1]])
plot(density(log2(edata+1)))
plot(density(edata))
asd = na.omit(edata)
plot(density(asd))
View(pdata)
setwd("C:/Users/coheni/Documents")
getwd()
setwd("C:/Users/coheni/Documents/Biostats")
dd = read.delim("051517/ProstatePartCEvents.txt")
dd = read.delim(unz("class2.zip", "ProstatePartCEvents.txt"))
dd2 = read.delim("http://cbio.mskcc.org/~socci/RCourse2017/ProstatePartCEvents.txt")
View(dd)
View(dd2)
data()
dir("data")
dir("051517")
getwd()
data()
data("MSIscores")
dir("051517")
data("051517/MSIscores")
MSIscores = read.delim("051517/MSIscores.csv")
View(MSIscores)
MSIscores = read.csv("051517/MSIscores.csv")
ls()
dir()
dir("051517")
patient = read.csv("051517/patientManifest.csv")
View(patient)
dim(patient)
dat = data.frame(Sample.ID = patient$Sample.ID, Patient.ID = patient$Patient.ID)
View(dat)
dir("051517")
library(gdata)
msi = read.xls("051517/MSIscores.xlsx")
install.packages("xlsx")
data("051517/MSIscores")
msi = read.csv("051517/MSIscores.csv")
View(msi)
dim(msi)
library(dplyr)
mm = match(dat$Sample.ID, MSIscores$DMP_ASSAY_ID)
mm
dat$MSI.score = MSIscores$Score[mm]
View(dat)
xx = cbind(dat$Sample.ID, MSIscores[mm])
xx = cbind(dat$Sample.ID, MSIscores$DMP_ASSAY_ID[mm])
sum(dat$MSI.score, na.rm = T)
is.na(dat$MSI.score)
count(is.na(dat$MSI.score))
mean(dat$MSI.score)
mean(dat$MSI.score, na.rm = T)
mm2 = grep(dat$Sample.ID, MSIscores$DMP_ASSAY_ID)
mm2
savehistory("051517.R")
