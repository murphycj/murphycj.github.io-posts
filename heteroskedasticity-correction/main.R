library(dglm)
library(caret)
library(ggplot2)

parameratize_dglm_method <- function(data) {

  #compare non-nested models of increasing powers

  r1 <- dglm(log2~gc,~gc,data=data,family=gaussian,phistart=1)
  r2 <- dglm(log2~gc,~I(gc^2),data=data,family=gaussian,phistart=1)
  r3 <- dglm(log2~gc,~I(gc^3),data=data,family=gaussian,phistart=1)
  r4 <- dglm(log2~gc,~I(gc^4),data=data,family=gaussian,phistart=1)
  r5 <- dglm(log2~gc,~I(gc^5),data=data,family=gaussian,phistart=1)
  r6 <- dglm(log2~gc,~I(gc^6),data=data,family=gaussian,phistart=1)
  r7 <- dglm(log2~gc,~I(gc^7),data=data,family=gaussian,phistart=1)
  r8 <- dglm(log2~gc,~I(gc^8),data=data,family=gaussian,phistart=1)

  rmses <- c(
    sqrt(sum((r1$dispersion.fit$residuals^2)/length(r1$dispersion.fit$residuals))),
    sqrt(sum((r2$dispersion.fit$residuals^2)/length(r2$dispersion.fit$residuals))),
    sqrt(sum((r3$dispersion.fit$residuals^2)/length(r3$dispersion.fit$residuals))),
    sqrt(sum((r4$dispersion.fit$residuals^2)/length(r4$dispersion.fit$residuals))),
    sqrt(sum((r5$dispersion.fit$residuals^2)/length(r5$dispersion.fit$residuals))),
    sqrt(sum((r6$dispersion.fit$residuals^2)/length(r6$dispersion.fit$residuals))),
    sqrt(sum((r7$dispersion.fit$residuals^2)/length(r7$dispersion.fit$residuals))),
    sqrt(sum((r8$dispersion.fit$residuals^2)/length(r8$dispersion.fit$residuals)))
  )

  aics <- c(
    2 * (3+3) + r1$m2loglik,
    2 * (3+3) + r2$m2loglik,
    2 * (3+3) + r3$m2loglik,
    2 * (3+3) + r4$m2loglik,
    2 * (3+3) + r5$m2loglik,
    2 * (3+3) + r6$m2loglik,
    2 * (3+3) + r7$m2loglik,
    2 * (3+3) + r8$m2loglik
  )

  png("./plots/sequencing_bias_correction/dglm-models.png",width=1200,height=600)
  par(mfrow=c(1,2),mar=c(5,5,5,5))
  plot(
    1:8,
    aics,
    xlab=expression(paste("n in ",log(sigma^2) == beta[0] + beta[1] * x^n)),
    ylab="AIC",
    main="Comparison of non-nested models"
  )
  lines(1:8,aics)
  plot(
    1:8,
    rmses,
    xlab=expression(paste("n in ",log(sigma^2) == beta[0] + beta[1] * x^n)),
    ylab="RMSE",
    main="Comparison of non-nested models"
  )
  lines(1:8,rmses)
  dev.off()


  # compare nested models of increasing powers
  # n = 5 does not converge

  r1 <- dglm(log2~gc,~gc,data=data,family=gaussian)
  r2 <- dglm(log2~gc,~gc+I(gc^2),data=data,family=gaussian)
  r3 <- dglm(log2~gc,~gc+I(gc^2)+I(gc^3),data=data,family=gaussian)
  r4 <- dglm(log2~gc,~gc+I(gc^2)+I(gc^3)+I(gc^4),data=data,family=gaussian)
  r5 <- dglm(log2~gc,~gc+I(gc^2)+I(gc^3)+I(gc^4)+I(gc^5),data=data,family=gaussian)
  r6 <- dglm(log2~gc,~gc+I(gc^2)+I(gc^3)+I(gc^4)+I(gc^5)+I(gc^6),data=data,family=gaussian)
  r7 <- dglm(log2~gc,~gc+I(gc^2)+I(gc^3)+I(gc^4)+I(gc^5)+I(gc^6)+I(gc^7),data=data,family=gaussian)
  r8 <- dglm(log2~gc,~gc+I(gc^2)+I(gc^3)+I(gc^4)+I(gc^5)+I(gc^6)+I(gc^7)+I(gc^8),data=data,family=gaussian)

  rmses <- c(
    sqrt(sum((r1$dispersion.fit$residuals^2)/length(r1$dispersion.fit$residuals))),
    sqrt(sum((r2$dispersion.fit$residuals^2)/length(r2$dispersion.fit$residuals))),
    sqrt(sum((r3$dispersion.fit$residuals^2)/length(r3$dispersion.fit$residuals))),
    sqrt(sum((r4$dispersion.fit$residuals^2)/length(r4$dispersion.fit$residuals))),
    sqrt(sum((r5$dispersion.fit$residuals^2)/length(r5$dispersion.fit$residuals))),
    sqrt(sum((r6$dispersion.fit$residuals^2)/length(r6$dispersion.fit$residuals))),
    sqrt(sum((r7$dispersion.fit$residuals^2)/length(r7$dispersion.fit$residuals))),
    sqrt(sum((r8$dispersion.fit$residuals^2)/length(r8$dispersion.fit$residuals)))
  )

  aics <- c(
    2 * (3+3) + r1$m2loglik,
    2 * (3+3) + r2$m2loglik,
    2 * (3+3) + r3$m2loglik,
    2 * (3+3) + r4$m2loglik,
    2 * (3+3) + r5$m2loglik,
    2 * (3+3) + r6$m2loglik,
    2 * (3+3) + r7$m2loglik,
    2 * (3+3) + r8$m2loglik
  )

  png("./plots/sequencing_bias_correction/dglm-models-nested.png",width=1400,height=700)
  par(mfrow=c(1,2),mar=c(5,5,5,5))
  plot(
    1:8,
    aics,
    xlab=expression(paste("n in ",log(sigma^2) == beta[0] + beta[1] * x^1 + ... + beta[1] * x^n)),
    ylab="AIC",
    main="Comparison of nested models"
  )
  lines(
    1:8,
    aics
  )
  plot(
    1:8,
    rmses,
    xlab=expression(paste("n in ",log(sigma^2) == beta[0] + beta[1] * x^1 + ... + beta[1] * x^n)),
    ylab="RMSE",
    main="Comparison of nested models"
  )
  lines(
    1:8,
    rmses
  )
  dev.off()

}

parameratize_loess_method <- function(window_data) {

  # 10-fold cross validation over loess parameters
  d1 <- c()
  d2 <- c()
  d3 <- c()

  for (sp in seq(0.1,1,0.05)) {
    for (deg in c(1,2)) {
      ss <- c()
      for (train_ids in createDataPartition(1:nrow(window_data),10,p=0.7)) {
        train <- window_data[train_ids,]
        test <- window_data[setdiff(1:nrow(window_data),train_ids),]
        r <- loess(window_sd~gc,span=sp,degree=deg,data=window_data)
        window_sd_pred <- predict(r,test)
        ss <- c(ss,sqrt(sum((window_sd_pred-test$window_sd)^2)/nrow(test)))
      }

      d1 <- c(d1,sp)
      d2 <- c(d2,deg)
      d3 <- c(d3,mean(ss))
    }
  }

  result <- data.frame(span=d1,degree=d2,ss=d3)
  result$degree <- factor(result$degree)

  p <- ggplot(result,aes(y=ss,x=span,col=degree,group=degree)) + geom_line() + xlab("Span") + ylab("Average RMSE (test data)")
  png("./plots/sequencing_bias_correction/loess-train.png",width=700,height=700)
  print(p)
  dev.off()

  write.csv(result,"./data/sequencing_bias_correction/loess-cross-validation.csv")

}

cnr <- read.table("./data/sequencing_bias_correction/HL-WES-18.cnr",sep="\t",header=T)
liver <- read.table("./data/sequencing_bias_correction/livers.cnn",sep="\t",header=T)

row.names(liver) <- paste(liver$chromosome,liver$start,liver$end,sep="-")
row.names(cnr) <- paste(cnr$chromosome, cnr$start, cnr$end,sep="-")
common <- intersect(row.names(liver),row.names(cnr))

cnr <- cnr[common,]
liver <- liver[common,]

data <- data.frame(gc=1-liver$gc,log2=cnr$log2)

# filter the extreme data and
# subsample the data since DGLM method does not seem to like too much data

data  <- data[(data$log2 > -5) & (data$log2 < 5),]
data <- data[sample(1:nrow(data),10000),]


png("./plots/sequencing_bias_correction/bias-data.png",width=700,height=700)
smoothScatter(
  data$gc,
  data$log2,
  ylim=c(-5,5),
  xlim=c(0,1),
  main="Biased data with heteroskedasticity",
  ylab="log2(tumor/liver)",
  xlab="GC content"
)
dev.off()

#rolling variance estimation

sds <- c()
gc_range <- seq(0.25,0.73,0.005)

for (window in gc_range) {
  tmp <- data[(data$gc>=window) & (data$gc<(window+0.005)),]
  sds <- c(sds,sd(tmp$log2)^2)
}

window_data <- data.frame(gc=gc_range,window_sd=sds)
print(window_data)

#parameratize the model and make some plots

#parameratize_dglm_method(data)
#parameratize_loess_method(window_data)

#create the best DGLM model

r1 <- dglm(log2~gc,~gc+I(gc^2)+I(gc^3),data=data,family=gaussian)
data$dlgm_predicted <- predict(r1$dispersion,data,type="response")

# create the best loess model

lr <- loess(window_sd~gc,data=window_data,span=0.1,degree=2)
data$lr_predicted <- predict(lr,data)
data$lr_predicted[is.na(data$lr_predicted)] <- median(data$lr_predicted,na.rm=T)

#correct the data by multiplying it by a ratio of standard deviations
# numerator: minimum estimated standard deviation from the data
# denominator: predicted standard deviations at each point

data$y_loess_corrected <- data$log2*sqrt(min(data$lr_predicted,na.rm=T))/sqrt(data$lr_predicted)
data$y_dglm_corrected <- data$log2*sqrt(min(data$dlgm_predicted,na.rm=T))/sqrt(data$dlgm_predicted)

png("./plots/sequencing_bias_correction/variance.png",width=700,height=700)
par(mar=c(6,6,3,3))
plot(
  gc_range,
  sds,
  main="Variance estimation",
  ylab="log2(tumor/liver) variance",
  xlab="GC content",
  pch=16,
  xlim=c(0,1),
  ylim=c(0,5)
)
dev.off()


png("./plots/sequencing_bias_correction/variance-estimate.png",width=700,height=700)
par(mar=c(6,6,3,3))
plot(
  gc_range,
  sds,
  main="Variance estimation",
  ylab="log2(tumor/liver) variance",
  xlab="GC content",
  pch=16,
  xlim=c(0,1),
  ylim=c(0,5)
)
lines(
  gc_range,
  predict(r1$dispersion.fit,data.frame(gc=gc_range),type="response"),
  col="red",
  pch=16,
  lwd=3
)
lines(
  gc_range,
  predict(lr,window_data),
  col="green",
  pch=16,
  lwd=3
)
legend("topleft",c("Window estimate","DGLM estimate","LOESS estimate"),col=c("black","red","green"),pch=16)
dev.off()


png("./plots/sequencing_bias_correction/corrected-data.png",width=1400,height=700)
par(mfrow=c(1,2))
smoothScatter(
  data$gc,
  data$y_loess_corrected,
  ylim=c(-5,5),
  xlim=c(0,1),
  main="LOESS corrected data",
  ylab="log2(tumor/liver)",
  xlab="GC content"
)
smoothScatter(
  data$gc,
  data$y_dglm_corrected,
  ylim=c(-5,5),
  xlim=c(0,1),
  main="DGLM corrected data",
  ylab="log2(tumor/liver)",
  xlab="GC content"
)
dev.off()
