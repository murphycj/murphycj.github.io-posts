library(ggplot2)

data <- read.csv("by_initPop_sampling.csv")
data[data$InitialPop==278,"InitialPop"] <- "1 ng"
data[data$InitialPop==1392,"InitialPop"] <- "5 ng"
data[data$InitialPop==2785,"InitialPop"] <- "10 ng"
data[data$InitialPop==13927,"InitialPop"] <- "50 ng"
data[data$InitialPop==27855,"InitialPop"] <- "100 ng"
data[data$InitialPop==65710,"InitialPop"] <- "200 ng"
data$InitialPop <- factor(data$InitialPop,levels=c("1 ng","5 ng","10 ng","50 ng","100 ng","200 ng"))
data$Sampling <- factor(data$Sampling)
data$PCREfficiency <- factor(data$PCREfficiency)
data$PolymeraseError <- factor(data$PolymeraseError)

pdf("by_starting_amount.pdf",width=10,height=6)
print(ggplot(data=data,aes(factor(InitialPop),Max)) + geom_boxplot() + ylab("MVAF") + xlab("Starting gDNA amount"))
dev.off()

pdf("by_starting_amount_bySampling.pdf",width=10,height=6)
print(ggplot(data=data,aes(factor(InitialPop),Max,fill=Sampling)) + geom_boxplot() + ylab("MVAF") + xlab("Starting gDNA amount"))
dev.off()

pdf("by_starting_amount_PCREfficiency.pdf",width=10,height=6)
print(ggplot(data=data,aes(factor(InitialPop),Max,fill=PCREfficiency)) + geom_boxplot() + ylab("MVAF") + xlab("Starting gDNA amount"))
dev.off()

pdf("by_starting_amount_PolymeraseError.pdf",width=10,height=6)
print(ggplot(data=data,aes(factor(InitialPop),Max,fill=PolymeraseError)) + geom_boxplot() + ylab("MVAF") + xlab("Starting gDNA amount"))
dev.off()
