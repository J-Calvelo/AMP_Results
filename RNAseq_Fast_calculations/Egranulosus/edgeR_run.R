setwd("/home/javier/Documents/Repositorios/Doctorado-Platelmintos/Peptidos_antimicrobianos/RNAseq_Fast_calculations/Egranulosus")
# Analisis de expresión diferencial entre los distintos estadíos de H. microstoma
# Consultas
#   - ¿Es combeniente usar otras aproximaciones para la dispersión (tendencias/tags)?  
#   - ¿Normalizar por GC%?

# 1) Cargar edgeR
library(edgeR)

# 2) Cargar los datos a R:
raw_counts <- read.delim("ERP004459.work_set.csv", row.names="gene_id")

#gene_id	ERR337909	ERR337922	ERR337912	ERR337925	ERR337961	ERR337973	ERR337937	ERR337949

# 3) Grupo Experimental: Metacestode (M), Non-activated protoscolex (P)
Groups <- factor(c("M","M","P","P","P","P"))

# 4) Producción del objeto a utilizar por edgeR
edgeR_Input <- DGEList(counts=raw_counts, group=Groups)

# 5) Filtrado de la Librería:
# filterByExpr.default <- function(y, design=NULL, group=NULL, lib.size=NULL, min.count=10, min.total.count=15, ...)
keep<- filterByExpr(edgeR_Input, group=Groups)
edgeR_Input <- edgeR_Input[keep, , keep.lib.sizes=FALSE]

# 6) Normalizado
edgeR_Input<-calcNormFactors(edgeR_Input)

# 7) Diseño de matriz experimental
design_matrix <- model.matrix(~0+Groups, data=edgeR_Input$samples)
colnames(design_matrix) <- levels(edgeR_Input$samples$group)

# 8) Calculos de expresión diferencial
edgeR_Input <- estimateDisp(edgeR_Input,design_matrix)

# 9) Ajuste de la distribución binomial
fit <- glmQLFit(edgeR_Input,design_matrix)

# 9) Expresión diferencial:
PvsM <- makeContrasts(P-M, levels=design_matrix)
PvsM_Test <- glmQLFTest(fit, contrast=PvsM)
topTags(PvsM_Test,10000)-> PvsM_Test_2019_09_13
write.table(PvsM_Test_2019_09_13, file="Results_PvsM_2019-09-13.tab", sep = "\t", col.names = NA)
  

summary(decideTests(PvsM_Test))
plotMD(PvsM_Test)
abline(h=c(-1, 1), col="blue")

