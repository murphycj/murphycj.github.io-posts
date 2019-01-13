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

  png("./plots/dglm-models.png",width=1200,height=600)
  par(mfrow=c(1,2),mar=c(5,5,5,5))
  plot(
    1:8,
    aics,
    xlab=expression(paste("n in ",log(sigma^2) == beta[0] + beta[1] * x^n)),
    ylab="AIC",
    main="Comparison of non-nested models",
    cex.lab=2.0,
    cex.main=2.0,
    cex.axis=2.0
  )
  lines(1:8,aics)
  plot(
    1:8,
    rmses,
    xlab=expression(paste("n in ",log(sigma^2) == beta[0] + beta[1] * x^n)),
    ylab="RMSE",
    main="Comparison of non-nested models",
    cex.lab=2.0,
    cex.main=2.0,
    cex.axis=2.0
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

  png("./plots/dglm-models-nested.png",width=1400,height=700)
  par(mfrow=c(1,2),mar=c(5,5,5,5))
  plot(
    1:8,
    aics,
    xlab=expression(paste("n in ",log(sigma^2) == beta[0] + beta[1] * x^1 + ... + beta[1] * x^n)),
    ylab="AIC",
    main="Comparison of nested models",
    cex.lab=2.0,
    cex.main=2.0,
    cex.axis=2.0
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
    main="Comparison of nested models",
    cex.lab=2.0,
    cex.main=2.0,
    cex.axis=2.0
  )
  lines(
    1:8,
    rmses
  )
  dev.off()

}

parameratize_loess_method <- function(window_data) {

  # 10-fold cross validation over loess parameters
  spans = seq(0.1,1,0.01)
  degrees = c(1,2)
  
  i <- 1
  n = length(spans) * length(degrees) * 2
  result <- data.frame(span=rep(0,n),degree=rep(0,n),rmse=rep(0,n),rmse_type=rep(0,n))

  for (sp in spans) {
    for (deg in degrees) {
      train_rmse <- c()
      test_rmse <- c()
      
      for (train_ids in createDataPartition(1:nrow(window_data),10,p=0.7)) {
        train <- window_data[train_ids,]
        test <- window_data[setdiff(1:nrow(window_data),train_ids),]
        
        #train
        r <- loess(window_sd~gc,span=sp,degree=deg,data=train,control=loess.control(surface="direct"))
        window_sd_pred <- predict(r,train)
        train_rmse <- c(train_rmse,sqrt(sum((window_sd_pred-train$window_sd)^2)/nrow(train)))
        
        #test
        #r <- loess(window_sd~gc,span=sp,degree=deg,data=window_data)
        window_sd_pred <- predict(r,test)
        test_rmse <- c(test_rmse,sqrt(sum((window_sd_pred-test$window_sd)^2)/nrow(test)))
        
      }
      
      result[i,"span"] <- sp
      result[i,"degree"] <- deg
      result[i,"rmse"] <- mean(train_rmse)
      result[i,"rmse_type"] <- "train"
      
      i <- i + 1
      
      result[i,"span"] <- sp
      result[i,"degree"] <- deg
      result[i,"rmse"] <- mean(test_rmse)
      result[i,"rmse_type"] <- "test"
      
      i <- i + 1
    }
  }
  
  #result$span <- factor(result$span,levels=spans)
  result$degree <- factor(result$degree,levels=degrees)
  
  png("./plots/loess-train.png",width=1400,height=700)
  print(cowplot::plot_grid(
    ggplot2::ggplot(result[result$rmse_type=="train",],aes(y=rmse,x=span,col=degree,group=degree)) + geom_line() + xlab("Span") + ylab("Average RMSE (train)") + theme(text = element_text(size=20)),
    ggplot2::ggplot(result[result$rmse_type=="test",],aes(y=rmse,x=span,col=degree,group=degree)) + geom_line() + xlab("Span") + ylab("Average RMSE (test)") + theme(text = element_text(size=20)),
    nrow=1
  ))
  dev.off()

  write.csv(result,"./data/loess-cross-validation.csv")

}

#load data

cnr <- read.table("./data/HL-WES-18.cnr",sep="\t",header=T)
liver <- read.table("./data/livers.cnn",sep="\t",header=T)

#create row names and filter what is common between tumor (cnr) and control (liver)

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


png("./plots/bias-data.png",width=700,height=700)
par(mar=c(5,5,5,5))
smoothScatter(
  data$gc,
  data$log2,
  ylim=c(-5,5),
  xlim=c(0,1),
  main="Biased data with heteroskedasticity",
  ylab="log2(tumor/liver)",
  xlab="GC content",
  cex.lab=2.0,
  cex.main=2.0,
  cex.axis=2.0
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

parameratize_dglm_method(data)
parameratize_loess_method(window_data)

#create the best DGLM model

r1 <- dglm(log2~gc,~gc+I(gc^2)+I(gc^3),data=data,family=gaussian)
data$dlgm_predicted <- predict(r1$dispersion,data,type="response")

# create the best loess model

lr <- loess(window_sd~gc,data=window_data,span=0.42,degree=2)
data$lr_predicted <- predict(lr,data)
data$lr_predicted[is.na(data$lr_predicted)] <- median(data$lr_predicted,na.rm=T)

#correct the data by multiplying it by a ratio of standard deviations
# numerator: minimum estimated standard deviation from the data
# denominator: predicted standard deviations at each point

data$y_loess_corrected <- data$log2*sqrt(min(data$lr_predicted,na.rm=T))/sqrt(data$lr_predicted)
data$y_dglm_corrected <- data$log2*sqrt(min(data$dlgm_predicted,na.rm=T))/sqrt(data$dlgm_predicted)

png("./plots/variance.png",width=700,height=700)
par(mar=c(6,6,3,3))
plot(
  gc_range,
  sds,
  main="Variance estimation",
  ylab="log2(tumor/liver) variance",
  xlab="GC content",
  pch=16,
  xlim=c(0,1),
  ylim=c(0,5),
  cex.lab=2.0,
  cex.main=2.0,
  cex.axis=2.0
)
dev.off()


png("./plots/variance-estimate.png",width=700,height=700)
par(mar=c(6,6,3,3))
plot(
  gc_range,
  sds,
  main="Variance estimation",
  ylab="log2(tumor/liver) variance",
  xlab="GC content",
  pch=16,
  xlim=c(0,1),
  ylim=c(0,5),
  cex.lab=2.0,
  cex.main=2.0,
  cex.axis=2.0
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


png("./plots/corrected-data.png",width=1400,height=700)
par(mfrow=c(1,2),mar=c(5,5,5,5))
smoothScatter(
  data$gc,
  data$y_loess_corrected,
  ylim=c(-5,5),
  xlim=c(0,1),
  main="LOESS corrected data",
  ylab="log2(tumor/liver)",
  xlab="GC content",
  cex.lab=2.0,
  cex.main=2.0,
  cex.axis=2.0
)
smoothScatter(
  data$gc,
  data$y_dglm_corrected,
  ylim=c(-5,5),
  xlim=c(0,1),
  main="DGLM corrected data",
  ylab="log2(tumor/liver)",
  xlab="GC content",
  cex.lab=2.0,
  cex.main=2.0,
  cex.axis=2.0
)
dev.off()
