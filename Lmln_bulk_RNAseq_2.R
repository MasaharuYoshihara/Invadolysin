# Lmln bulk RNA seq TAE 051

# KO_1.  KO     869 30W
# KO_2   KO     871 30W
# KO_3   KO     805 42W
# cont_1 hetero 870 30W
# cont_2 hetero 872 30W
# cont_3 WT     807 42W

# https://olvtools.com/documents/edger

library(edgeR)

data <- read.csv("TAE_051_quantile_RPKM_normalized.csv", sep=",", row.names=1)
#data_gene <- data[,1]
KO1 <- data[,14]
KO2 <- data[,17]
KO3 <- data[,20]
HET1 <- data[,23]
HET2 <- data[,26]
WT <- data[,29]
counts <- data.frame(KO1, KO2, KO3, HET1, HET2, WT)
group <- factor(c("KO", "KO", "KO", "HET", "HET", "WT"))
y <- DGEList(counts=counts, group=group)

nrow(y) # total gene number 49585
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
nrow(y) # DEG number 5212

y <- calcNormFactors(y)

design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=2)
topTags(qlf)
# Coefficient:  groupKO 
# data_gene      logFC    logCPM        F
# 9427        Gm10377   4.058147  4.695639 15.29712
# 31157          Gpx5 -15.655565 13.002601 12.16222
# 45222        Spink8 -12.120479 10.023118 11.69660
# 45228        Spint5  -7.061732  7.283712 11.46809
# 2606          Adam7 -13.697444 11.027143 11.46252
# 33660          Lcn5 -14.017915 11.471603 11.42337
# 1925  9230113P08Rik -10.539915  8.157241 11.33076
# 28457        Gm4846 -10.271735  9.503482 11.26664
# 45102       Spag11b -12.710326 10.044916 11.25122
# 6876         Defb28 -10.287709  7.981875 11.22368
# PValue       FDR
# 9427  0.01319263 0.4544459
# 31157 0.01996470 0.4544459
# 45222 0.02138001 0.4544459
# 45228 0.02212781 0.4544459
# 2606  0.02214651 0.4544459
# 33660 0.02227856 0.4544459
# 1925  0.02259550 0.4544459
# 28457 0.02281873 0.4544459
# 45102 0.02287287 0.4544459
# 6876  0.02297007 0.4544459

result <- as.data.frame(topTags(qlf, n=nrow(y)))
result[result$FDR<0.05,]
# [1] data_gene logFC     logCPM    F        
# [5] PValue    FDR      
# <0 rows> (or 0-length row.names)

plotMDS(y, gene.selection="common", ylim = c(-0.1,0.05), xlim = c(-2,6))


