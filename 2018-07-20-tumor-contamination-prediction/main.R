colremove <- c("chr","start","end")

#load the data

coverage.genome <- as.matrix(read.csv("./data/BedtoolsCount.counts.csv",header=T,check.names=F))
coverage.trp53 <- read.csv("./data/BedtoolsCount.Trp53.counts.csv",header=T,check.names=F)

coverage.genome <- coverage.genome[,setdiff(colnames(coverage.genome),colremove)]
coverage.trp53 <- coverage.trp53[,setdiff(colnames(coverage.trp53),colremove)]

coverage.genome <- apply(coverage.genome,2,as.numeric)
coverage.trp53 <- apply(coverage.trp53,2,as.numeric)

livers <- c("HL-WES-39", "HL-WES-40", "HL-WES-41")
tumors <- setdiff(colnames(coverage.genome),livers)


# filter the data so rows with values 0 and remove two 95% quantile of data
# to avoid outlier data

data <- data.frame(
  tumor=coverage.genome[,"HL-WES-28"],
  liver=apply(coverage.genome[,livers],1,mean)
)

upper.cutoff <- quantile(apply(coverage.genome,1,mean),0.99)[[1]]
data <- data[(data[,1]>0) & (data[,2]>0) & (data[,1]<upper.cutoff) & (data[,2]<upper.cutoff),]

# poisson regression

r <- glm(tumor~log(liver),data=data,poisson(link="log"))
trp53_predicted_coverage <- predict(r,data.frame(liver=mean(as.numeric(coverage.trp53[1,livers]))),type="response")[[1]]
#1 - (coverage.trp53[1,tumor] / trp53_predicted_coverage)

x <- seq(min(data[,2]),max(data[,2]),10)

png("./plots/tumor-vs-liver.png",width=800,height=400)
par(mfrow=c(1,2),mar=c(5,5,5,5))
plot(
  data[,2],
  data[,1],
  xlab="Liver read count",
  ylab="Tumor read count",
  pch=".",
  cex.lab=2.0,
  cex.main=2.0,
  cex.axis=2.0
)
lines(x,predict(r,data.frame(liver=x),type="response"),col="red")
plot(
  data[,2],
  r$residuals,
  xlab="Liver read count",
  ylab="Residuals",
  pch=".",
  ylim=c(-5,5),
  cex.lab=2.0,
  cex.main=2.0,
  cex.axis=2.0
)
dev.off()


# gaussian regression

data2 <- log(data+1)

r <- glm(tumor~liver,data=data2,gaussian)
trp53_predicted_coverage <- predict(r,data.frame(liver=log(mean(as.numeric(coverage.trp53[1,livers]))+1)),type="response")[[1]]
#1 - (log(coverage.trp53[1,tumor]+1) / trp53_predicted_coverage)


x <- seq(min(data2[,2]),max(data2[,2]),0.1)

png("./plots/tumor-vs-liver-gaussian.png",width=800,height=400)
par(mfrow=c(1,2),mar=c(5,5,5,5))
plot(
  data2[,2],
  data2[,1],
  xlab="Liver read count",
  ylab="Tumor read count",
  pch=".",
  cex.lab=2.0,
  cex.main=2.0,
  cex.axis=2.0
)
lines(x,predict(r,data.frame(liver=x),type="response"),col="red")
plot(
  data2[,2],
  r$residuals,
  xlab="Liver read count",
  ylab="Residuals",
  pch=".",
  ylim=c(-5,5),
  cex.lab=2.0,
  cex.main=2.0,
  cex.axis=2.0
)
dev.off()
