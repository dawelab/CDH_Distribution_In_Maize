install.packages("snpR")
library(snpR)
library(ggplot2)

vcf <- read_vcf("/Volumes/Transcend/SGEOnly_Ab10BChrom.Ab10Hap.AllAb10Pos.vcf.gz") 

vcf <- calc_tajimas_d(vcf, sigma = 50, step = 25)
vcf <- calc_pi(vcf)
D <- get.snpR.stats(vcf, stats="single.window")
P <- get.snpR.stats(vcf, stats="single")

pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/8_SNPs/Nucleotide_Diversity.pdf", height = 2, width=5)
p <- ggplot(P) +
  geom_point(aes(x=position, y = pi)) +
  scale_x_continuous(breaks = c(140000000, 150000000, 160000000, 170000000, 180000000, 190000000), labels = c(140, 150, 160, 170, 180, 190)) +
  labs(x="Ab10-I (MB)", y= "Nucleotide Diversity") +
  annotate("rect", xmin=c(141115174,152050000, 158250000), xmax=c(142472000,156880000, 167721000), ymin=c(-0.2, -0.2, -0.2) , ymax=c(-0.1, -0.1, -0.1), alpha=1, color=NA, fill="darkgoldenrod") +
  annotate("rect", xmin=c(142472000,150656000, 157485200), xmax=c(146699300,153145000, 159356550), ymin=c(-0.2, -0.2, -0.2) , ymax=c(-0.1, -0.1, -0.1), alpha=1, color=NA, fill="deepskyblue") +
  annotate("rect", xmin=c(174433450), xmax=c(182846100), ymin=c(-0.2, -0.2, -0.2) , ymax=c(-0.1, -0.1, -0.1), alpha=1, color=NA, fill="darkorange3") +
geom_hline(yintercept = 0.02, color = "red")
  
p
dev.off()


pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/8_SNPs/Tajimas_D.pdf", height = 2, width=5)
p <- ggplot(D) +
  geom_segment(aes(x=start, xend = end, y = D, yend = D)) +
  scale_x_continuous(breaks = c(140000000, 150000000, 160000000, 170000000, 180000000, 190000000), labels = c(140, 150, 160, 170, 180, 190)) +
  labs(x="Ab10-I (MB)", y= "Tajima's D value") +
  annotate("rect", xmin=c(141115174,152050000, 158250000), xmax=c(142472000,156880000, 167721000), ymin=c(-3,-3,-3) , ymax=c(-3.5, -3.5, -3.5), alpha=1, color=NA, fill="darkgoldenrod") +
  annotate("rect", xmin=c(142472000,150656000, 157485200), xmax=c(146699300,153145000, 159356550), ymin=c(-3,-3,-3) , ymax=c(-3.5, -3.5, -3.5), alpha=1, color=NA, fill="deepskyblue") +
  annotate("rect", xmin=c(174433450), xmax=c(182846100), ymin=c(-3,-3,-3) , ymax=c(-3.5, -3.5, -3.5), alpha=1, color=NA, fill="darkorange3") +
  geom_hline(yintercept = 2, color = "red") +
  geom_hline(yintercept = -2, color = "red")
p
dev.off()
