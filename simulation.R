library("limma")
library("edgeR")
library("DESeq2")
library("splines")
library("Biobase")
library("BiocParallel")
files <- c("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/cheung_eset.RData",
           "http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData",
           "http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
for (file in files) if (!file.exists(basename(file))) download.file(file, basename(file))
for (file in files) load(basename(file))
bottomly.eset <- bottomly.eset[ rowMeans(exprs(bottomly.eset)) >= 1,]
dim(bottomly.eset)
cheung.eset <- cheung.eset[ rowMeans(as.matrix(exprs(cheung.eset))) >= 1,]
dim(cheung.eset)
montpick.eset <- montpick.eset[ rowMeans(exprs(montpick.eset)) >= 1,]
dim(montpick.eset)

runDESeq2 <- function(e, alpha=1e-4) {
  dds <- DESeqDataSetFromMatrix(as.matrix(exprs(e)), pData(e), ~1)
  m <- ncol(e)
  dds$condition <- factor(rep(1:2,each=m/2))
  design(dds) <- ~ condition
  dds <- DESeq(dds, quiet=TRUE)
  res <- results(dds)
  sum(res$pvalue < alpha, na.rm=TRUE)
}

runEdgeR <- function(e, alpha=1e-4) {
  m <- ncol(e)
  condition <- factor(rep(1:2,each=m/2))
  design <- model.matrix(~ condition)
  dgel <- DGEList(as.matrix(exprs(e)))
  dgel <- calcNormFactors(dgel)
  dgel <- estimateGLMCommonDisp(dgel,design=design)
  dgel <- estimateGLMTrendedDisp(dgel,design=design)
  dgel <- estimateGLMTagwiseDisp(dgel,design=design)
  fit <- glmFit(dgel,design=design)
  lrt <- glmLRT(fit,coef=2)
  pvals <- lrt$table$PValue
  sum(pvals < alpha, na.rm=TRUE)
}

runLimmaVoom <- function(e, alpha=1e-4) {
  m <- ncol(e)
  condition <- factor(rep(1:2,each=m/2))
  design <- model.matrix(~ condition)
  dgel <- DGEList(as.matrix(exprs(e)))
  dgel <- calcNormFactors(dgel)
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  pvals <- tt$P.Value
  sum(pvals < alpha, na.rm=TRUE)
}

bottomly1 <- bottomly.eset[,bottomly.eset$strain == levels(bottomly.eset$strain)[1]]
bottomly2 <- bottomly.eset[,bottomly.eset$strain == levels(bottomly.eset$strain)[2]]
cheung1 <- cheung.eset[,cheung.eset$gender == levels(cheung.eset$gender)[1]]
cheung2 <- cheung.eset[,cheung.eset$gender == levels(cheung.eset$gender)[2]]
mont1 <- montpick.eset[,montpick.eset$population == levels(montpick.eset$population)[1]]
mont2 <- montpick.eset[,montpick.eset$population == levels(montpick.eset$population)[2]]

# subset how many total samples
m <- 6

# how many subsets
# divide by two: symmetric partitions count once
chooseTwoSets <- function(n,m) choose(n,m/2) * choose(n-m/2,m/2) / 2
dim(bottomly1)
chooseTwoSets(ncol(bottomly1), m)
dim(cheung1)
chooseTwoSets(ncol(cheung1), m)
dim(mont1)
chooseTwoSets(ncol(mont1), m)

register(MulticoreParam(20))

# master seed
set.seed(1)
nsim <- 200
seeds <- round(runif(nsim,1,1e6))

runAlgos <- function(i, e) {
  set.seed(seeds[i])
  sub <- sample(ncol(e), m, FALSE)
  nz <- rowSums(exprs(e)[,sub]) > 0
  d <- runDESeq2(e[nz,sub])
  ed <- runEdgeR(e[nz,sub])
  v <- runLimmaVoom(e[nz,sub])
  data.frame(DESeq2=d,edgeR=ed,limma.voom=v,nonzero=sum(nz),subset=paste(sub, collapse=","),stringsAsFactors=FALSE)
}

# run within-group comparisons for group 1 then for group 2
res.b1 <- bplapply(1:nsim, runAlgos, e=bottomly1)
res.b1 <- do.call(rbind, res.b1)
res.b2 <- bplapply(1:nsim, runAlgos, e=bottomly2)
res.b2 <- do.call(rbind, res.b2)

res.c1 <- bplapply(1:nsim, runAlgos, e=cheung1)
res.c1 <- do.call(rbind, res.c1)
res.c2 <- bplapply(1:nsim, runAlgos, e=cheung2)
res.c2 <- do.call(rbind, res.c2)

res.m1 <- bplapply(1:nsim, runAlgos, e=mont1)
res.m1 <- do.call(rbind, res.m1)
res.m2 <- bplapply(1:nsim, runAlgos, e=mont2)
res.m2 <- do.call(rbind, res.m2)

summary(res.b1) # bottomly
summary(res.b2) # bottomly
summary(res.c1) # cheung
summary(res.c2) # cheung
summary(res.m1) # montgomery
summary(res.m2) # montgomery

(session.info <- sessionInfo())

save(res.b1, res.b2, res.c1, res.c2, res.m1, res.m2, session.info, file="results.rda")

