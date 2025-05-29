# Estimate genetic diversity in each sample/population
# R load packages
library(adegenet)
library(ade4)
library(hierfstat)
library(RColorBrewer)

#Mywork_malsha
# Set the path to your Genepop file
file_path <- "cod_greenland-iceland-canada.gen"
# Read the first 20 lines of the file
lines <- readLines(file_path, n = 20)
# Print them to inspect
cat(lines, sep = "\n")



# import dataset.gen into R
x <- import2genind("cod_greenland-iceland-canada.gen", ncode = 3)

# 1. Compare heterozygosity as a measure of genetic diversity between populations
x.pop = seppop(x)
summary.by.pop = lapply(x.pop, summary)

# Expected heterozygosity
Hexp.ls = rep(NA, length(summary.by.pop)) 
for (i in 1:length(summary.by.pop)){ 
  Hexp.ls[i] = mean(summary.by.pop[[i]]$Hexp) 
} 
barplot(Hexp.ls, names.arg = levels(pop(x)), las = 2, main = "Expected heterozygosity", ylab = "He")


# Observed heterozygosity
Hobs.ls = rep(NA, length(summary.by.pop))
for (i in 1:length(summary.by.pop)){
  Hobs.ls[i] = mean(summary.by.pop[[i]]$Hobs)
}
barplot(Hobs.ls, names.arg = levels(pop(x)), las = 2, main =
          "Observed heterozygosity", ylab = "Ho")

#to fill the gap_Malsha
Hobs.ls = rep(NA, length(summary.by.pop))
for (i in 1:length(summary.by.pop)){
  # Use na.rm = TRUE to ignore NA values
  Hobs.ls[i] = mean(summary.by.pop[[i]]$Hobs, na.rm = TRUE)
}
barplot(Hobs.ls, names.arg = levels(pop(x)), las = 2, main =
          "Observed heterozygosity", ylab = "Ho", col = "#FD8A9F")


# 2. Determination of candidates under positive selection
# (BayeScan http://cmpg.unibe.ch/software/BayeScan/versions.html)
# Plot Bayescan result in R
#plot bayescan

# Read the Bayescan results
data <- read.table("bayes_fst_cod.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)

# Optional: set a threshold for outliers (commonly q-value < 0.05)
threshold <- 0.05

library(dplyr)
library(tidyverse)
library(ggplot2)

# Flag significant SNPs
colnames(data)

#Mywork_Malsha
str(data)
threshold <- 0.05  # or whatever your chosen threshold is

data <- data %>%
  mutate(Significant = ifelse(qval < threshold, "Outlier", "Neutral"))

ggplot(data, aes(x = -log10(qval), y = fst, color = Significant)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Neutral" = "gray", "Outlier" = "red")) +
  theme_minimal() +
  labs(
    title = "BayeScan FST Plot",
    x = "-log10(q-value)",
    y = "FST",
    color = "Category"
  ) +
  xlim(0, 5) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# save outlier snps
outlier_snps <-
  data %>% 
  filter(Significant == "Outlier")
#Maslsha_to find n code
readLines("cod_greenland-iceland-canada.gen", n = 30)

# import genotype dataset into R
#iam not sure ncode is 3 or other check again
dataset_all <- import2genind("cod_greenland-iceland-canada.gen", ncode=3)

# subset dataset for outlier snps
outlier_dataset <- x[loc=outlier_snps$snp_pos]
locNames(outlier_dataset)

pairwise.neifst <- function(dat,diploid=TRUE){
  if (is.genind(dat)) dat<-genind2hierfstat(dat)
  dat<-dat[order(dat[,1]),]
  pops<-unique(dat[,1])
  npop<-length(pops)
  fstmat <- matrix(nrow=npop,ncol=npop,dimnames=list(pops,pops))
  if (is.factor(dat[,1])) {
    dat[,1]<-as.numeric(dat[,1])
    pops<-as.numeric(pops)
  }
  for(a in 2:npop){
    for(b in 1:(a-1)){
      subdat <- dat[dat[,1] == pops[a] | dat[,1]==pops[b],]
      fstmat[a,b]<-fstmat[b,a]<- basic.stats(subdat,diploid=diploid)$overall[8]
    }
    
  }
  fstmat
}

Fst.outlier <- pairwise.neifst(outlier_dataset) 

# make a heatmap to visualize
heatmap(as.matrix(Fst.outlier), symm=T, col=brewer.pal(9,"YlGnBu"))

#malsha color heatmap
my_palette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)  # 100 smooth colors
heatmap(as.matrix(Fst.outlier), symm = TRUE, col = my_palette)


# 3. Demographic inferences from SNPs under positive selection vs. neutral loci
# PCA Analysis with FST distance matrix
library(stats)

Fst.outlier_mod <- 
  Fst.outlier %>% 
  as.data.frame() %>%
  mutate(across(everything(), ~replace_na(., 1))) %>%
  as.matrix()

outlier_pca <- princomp(as.matrix(Fst.outlier_mod))

#see variance
variance <-outlier_pca$sdev^2 / sum(outlier_pca$sdev^2)
head(variance)

#plot PCA
components <- as.data.frame(outlier_pca$scores[,1:12])

plot(components$Comp.1, components$Comp.2, xlab = "PC 1 - 95%",
     ylab= "PC 2 - 0.2")

text(components$Comp.1, components$Comp.2, row.names(components),
     cex=0.6, pos=4)

abline(h =0, lty="dashed", untf = FALSE)
abline(v =0, lty="dashed", untf = FALSE)


# repeat PCA analysis with neutral dataset
#Malsha add some to correct data
# Get all loci names
all_loci <- locNames(x)

# Get only neutral loci (those not in outlier list)
neutral_loci <- setdiff(all_loci, outlier_snps$snp_pos)

# Subset the genind object for neutral loci
neutral_dataset <- x[loc = neutral_loci]

#repeat the last steps with the neutral SNPs
neutral_dataset <- x[loc=-outlier_snps$snp_pos]
Fst.neutral <- pairwise.neifst(neutral_dataset) 

Fst.neutral_mod <- 
  Fst.neutral %>% 
  as.data.frame() %>%
  mutate(across(everything(), ~replace_na(., 1))) %>%
  as.matrix()

neutral_pca <- princomp(as.matrix(Fst.neutral_mod))

#see variance
variance <-neutral_pca$sdev^2 / sum(neutral_pca$sdev^2)
head(variance)

#plot PCA
components <- as.data.frame(neutral_pca$scores[,1:12])

plot(components$Comp.1, components$Comp.2, xlab = "PC 1 - 95%",
     ylab= "PC 2 - 0.2")

text(components$Comp.1, components$Comp.2, row.names(components),
     cex=0.6, pos=4)

abline(h =0, lty="dashed", untf = FALSE)
abline(v =0, lty="dashed", untf = FALSE)
