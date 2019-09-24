setwd("/home/javier/Documents/Repositorios/Doctorado-Platelmintos/Peptidos_antimicrobianos/RNAseq_Fast_calculations/Egranulosus")
# Analisis de expresión diferencial entre los distintos estadíos de H. microstoma
# Consultas
#   - ¿Es combeniente usar otras aproximaciones para la dispersión (tendencias/tags)?  
#   - ¿Normalizar por GC%?

# 1) Cargar edgeR
library(edgeR)

# 2) Cargar los datos a R:
raw_counts <- read.delim("ERP004459.counts_per_run.tsv", row.names="gene_id")

#gene_id	ERR337909	ERR337922	ERR337912	ERR337925	ERR337961	ERR337973	ERR337937	ERR337949

# 3) Grupo Experimental: Metacestode (M), Non-activated protoscolex (P)
Groups <- factor(c("M","M","P","P","P","P","P","P"))

# 4) Producción del objeto a utilizar por edgeR
edgeR_Input <- DGEList(counts=raw_counts, group=Groups)

# 5) Filtrado de la Librería:
# filterByExpr.default <- function(y, design=NULL, group=NULL, lib.size=NULL, min.count=10, min.total.count=15, ...)
keep<- filterByExpr(edgeR_Input, group=Groups)
edgeR_Input <- edgeR_Input[keep, , keep.lib.sizes=FALSE]





# Re adaptar script: 
##################################################################################
##################################################################################
##################################################################################
##################################################################################

# 5) Filtrado de la Librería:
# filterByExpr.default <- function(y, design=NULL, group=NULL, lib.size=NULL, min.count=10, min.total.count=15, ...)
keep<- filterByExpr(edgeR_Input, group=Groups)
edgeR_Input <- edgeR_Input[keep, , keep.lib.sizes=FALSE]

# 6) Normalizado
edgeR_Input<-calcNormFactors(edgeR_Input)

# 8) Diseño de matriz experimental
design_matrix <- model.matrix(~0+Groups, data=edgeR_Input$samples) # No termino de entender como funciona...
colnames(design_matrix) <- levels(edgeR_Input$samples$group)

# 7) Calculos de expresión diferencial
edgeR_Input <- estimateDisp(edgeR_Input,design_matrix)

# 8) Ajuste de la distribución binomial
fit <- glmQLFit(edgeR_Input,design_matrix)

# 9) Expresión diferencial:
# Adulto vs Cycsticeroid
CvsA <- makeContrasts(C-A, levels=design_matrix)
CvsA_Test <- glmQLFTest(fit, contrast=CvsA)
topTags(CvsA_Test,10000)-> CyVsAd_2019_09_10
# write.table(CyVsAd_2019_09_10, file="Results_CyVsAd_2019-09-10.tab", sep = "\t", col.names = NA)

# Oncosphera vs Adulto
OvsA <- makeContrasts(O-A, levels=design_matrix)
OvsA_Test <- glmQLFTest(fit, contrast=OvsA)
topTags(OvsA_Test,10000)-> OnVsAd_2019_09_10
# write.table(OnVsAd_2019_09_10, file="Results_OnVsAd_2019-09-10.tab", sep = "\t", col.names = NA)

# Cycsticeroid vs Oncosphera
CvsO <- makeContrasts(O-C, levels=design_matrix)
CvsO_Test <- glmQLFTest(fit, contrast=CvsO)
topTags(CvsO_Test,10000)-> CyVsOn_2019_09_10
# write.table(CyVsOn_2019_09_10, file="Results_CyVsOn_2019-09-10.tab", sep = "\t", col.names = NA)



# plotMDS: Multidimensional scaling plot 
col.cond <- c(rep("Gold4",4),rep("blue",4),rep("green4",4)) # Colores de las condiciones (Se toman las columnas en orden)
plotMDS(edgeR_Input, col= col.cond, cex=5, pch = ".") # Gràfico en si
title("Adultos Vs Cisticercoides Vs Oncosferas") # Tíulo
legend(x=1.8,y=0.4,c("Adult","Cycsticeroid","Oncosphera"),cex=0.8,col=c("Gold4","blue","green4"),pch=".") # Leyenda (no he logrado que quede bien)

# Grafica de dispersión (Tengo que preguntar como se interpretan)
plotBCV(edgeR_Input)

# Resumen e resultaos
summary(decideTests(CvsA_Test, p.value=0.001))
summary(decideTests(OvsA_Test, p.value=0.001))
summary(decideTests(CvsO_Test, p.value=0.001))

test_results<-decideTests(CvsO_Test, p.value=0.00001)

# Plot log-fold change against log-counts per million, with DE genes highlighted:
plotMD(CvsA_Test)
abline(h=c(-1, 1), col="blue")

plotMD(OvsA_Test)
abline(h=c(-1, 1), col="blue")

plotMD(CvsO_Test)
abline(h=c(-1, 1), col="blue")
