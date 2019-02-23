library(xlsx)

data = read.csv("/Users/charlesmurphy/Desktop/Research/1216_bishoy/results/DNA/Mutations/somatic.csv",header=T,check.names=F)
samples = colnames(data)[12:ncol(data)]

vaf = c()
samples2 = c()

for (s in samples) {
  for (i in 1:nrow(data)) {
    if (data[i,s]!="-") {
      tmp = strsplit(as.character(data[i,s])," ")[[1]]
      vaf <- c(vaf,as.numeric(tmp[1]))
      samples2 <- c(samples2,strsplit(s,"_")[[1]][1])
    }
  }
}

mutation_data = data.frame(sample=samples2,VAF=vaf)
initial_dna_conc <- read.xlsx("~/Desktop/Research/1216_bishoy/docs/initial_dna_amount_REPLIg.xlsx","Sheet1",check.names=F)
row.names(initial_dna_conc) <- initial_dna_conc[,"Sample name"]
common_samples <- intersect(mutation_data[,1],initial_dna_conc[,"Sample name"])
initial_dna_conc <- initial_dna_conc[initial_dna_conc[,"Sample name"] %in% common_samples,]
mutation_data <- mutation_data[mutation_data[,"sample"] %in% common_samples,]

coverage <- read.csv("../../DNA/GATKCoverage_CDH1_coding/GATKCoverage.sample_summary.csv",header=T,check.names=F,row.names=1)
colnames(coverage) <- unlist(lapply(colnames(coverage),function(x) return(strsplit(x,"_")[[1]][1])))

d1 <- c()
d2 <- c()
d3 <- c()
for (i in common_samples) {
  tmp <- mutation_data_low[mutation_data_low[,"sample"]==i,]
  tmp_low <- tmp[tmp$VAF<=0.02,]

  d1 <- c(d1,nrow(tmp_low))
  d2 <- c(d2,initial_dna_conc[i,"pre-WGA"])
  d3 <- c(d3,coverage["mean",i])
}

rr <- data.frame(mut=d1,dna=d2,sample=common_samples,coverage=d3)

pdf("initialgDNA-vs-numMutations.pdf")
plot(rr$dna,rr$mut,xlab="Initial gDNA (ng)",ylab="Number mutations (VAF <= 5%)")
dev.off()

summary(glm(mut~log2(dna)+coverage,rr))
