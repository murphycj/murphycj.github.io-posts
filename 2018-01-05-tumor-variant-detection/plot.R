library(ggplot2)
library(cowplot)

data <- read.csv("results.csv",header=T)
data$Tn <- factor(data$Tn, levels=sort(unique(data$Tn)))
data$Tr <- factor(data$Tr, levels=sort(unique(data$Tr)))
data$Cn <- factor(data$Cn, levels=sort(unique(data$Cn)))
data$Cc <- factor(data$Cc, levels=sort(unique(data$Cc)))
data$Cr <- factor(data$Cr, levels=sort(unique(data$Cr)))


pdf("Tumor-totalMutations.pdf",width=14,height=7)
print(cowplot::plot_grid(
  ggplot(data,aes(y=total_mutations, x=Tn)) + geom_boxplot() + ylab("Total mutations") + theme(text = element_text(size=20)),
  ggplot(data,aes(y=total_mutations, x=Tr)) + geom_boxplot() + ylab("Total mutations") + theme(text = element_text(size=20)),
  nrow=1
))
dev.off()

pdf("Normal-totalMutations.pdf",width=12,height=12)
print(cowplot::plot_grid(
  ggplot(data,aes(y=total_mutations, x=Cn)) + geom_boxplot() + ylab("Total mutations") + theme(text = element_text(size=20)),
  ggplot(data,aes(y=total_mutations, x=Cc)) + geom_boxplot() + ylab("Total mutations") + theme(text = element_text(size=20)),
  ggplot(data,aes(y=total_mutations, x=Cr)) + geom_boxplot() + ylab("Total mutations") + theme(text = element_text(size=20)),
  nrow=2
))
dev.off()


#
# pdf("Tumor-proportion_in_dbsnp.pdf",width=14,height=7)
# print(cowplot::plot_grid(
#   ggplot(data,aes(y=proportion_in_dbsnp, x=Tn)) + geom_boxplot() + ylab("Total mutations") + theme(text = element_text(size=20)),
#   ggplot(data,aes(y=proportion_in_dbsnp, x=Tr)) + geom_boxplot() + ylab("Total mutations") + theme(text = element_text(size=20)),
#   nrow=1
# ))
# dev.off()
#
# pdf("Normal-proportion_in_dbsnp.pdf")
# print(cowplot::plot_grid(
#   ggplot(data,aes(y=proportion_in_dbsnp, x=Cn)) + geom_boxplot() + ylab("Total mutations") + theme(text = element_text(size=20)),
#   ggplot(data,aes(y=proportion_in_dbsnp, x=Cc)) + geom_boxplot() + ylab("Total mutations") + theme(text = element_text(size=20)),
#   ggplot(data,aes(y=proportion_in_dbsnp, x=Cr)) + geom_boxplot() + ylab("Total mutations") + theme(text = element_text(size=20)),
#   nrow=2
# ))
# dev.off()


pdf("Tumor-proportion_in_dbsnp_2.pdf",width=14,height=7)
print(cowplot::plot_grid(
  ggplot(data,aes(y=proportion_in_dbsnp_2, x=Tn)) + geom_boxplot() + ylab("Proportion false positives (dbSNP)") + theme(text = element_text(size=20)),
  ggplot(data,aes(y=proportion_in_dbsnp_2, x=Tr)) + geom_boxplot() + ylab("Proportion false positives (dbSNP)") + theme(text = element_text(size=20)),
  nrow=1
))
dev.off()

pdf("Normal-proportion_in_dbsnp_2.pdf",width=12,height=12)
print(cowplot::plot_grid(
  ggplot(data,aes(y=proportion_in_dbsnp_2, x=Cn)) + geom_boxplot() + ylab("Proportion false positives (dbSNP)") + theme(text = element_text(size=20)),
  ggplot(data,aes(y=proportion_in_dbsnp_2, x=Cc)) + geom_boxplot() + ylab("Proportion false positives (dbSNP)") + theme(text = element_text(size=20)),
  ggplot(data,aes(y=proportion_in_dbsnp_2, x=Cr)) + geom_boxplot() + ylab("Proportion false positives (dbSNP)") + theme(text = element_text(size=20)),
  nrow=2
))
dev.off()
#
#
# pdf("Tn-proportion_high_quality_singleton.pdf")
# plot(proportion_high_quality_singleton~Tn,data=data)
# dev.off()
# pdf("Tr-proportion_high_quality_singleton.pdf")
# plot(proportion_high_quality_singleton~Tr,data=data)
# dev.off()
# pdf("Cn-proportion_high_quality_singleton.pdf")
# plot(proportion_high_quality_singleton~Cn,data=data)
# dev.off()
# pdf("Cc-proportion_high_quality_singleton.pdf")
# plot(proportion_high_quality_singleton~Cc,data=data)
# dev.off()
# pdf("Cr-proportion_high_quality_singleton.pdf")
# plot(proportion_high_quality_singleton~Cr,data=data)
# dev.off()
