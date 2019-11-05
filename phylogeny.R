# Installing packages and loading libraries required

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("msa")

# Expression 'pkg::name' returns the value of the exported variable 'name' in package 'pkg' if the package has a name space.

library(Biostrings)

install.packages("seqinr")
library(seqinr)

install.packages("ape")
library(ape)

library(msa)

# Converting sequences into alignment. Sequences are stored in fasta format in "16MT7_Final.fasta" in this case. 

rRNA <- readDNAStringSet("16SMT7_Final.fasta")
rRNA_Aligned <- msa(rRNA)

rRNA_Aligned2 <- msaConvert(rRNA_Aligned, type="seqinr::alignment")

# Neighbour joining method

rRNA_dist <- dist.alignment(rRNA_Aligned2, "identity")

rRNA_Tree <- nj(rRNA_dist)

plot(rRNA_Tree)

# Maximum Parsimony 

install.packages("phangorn")
library(phangorn)

rRNA_Aligned3 <- msaConvert(rRNA_Aligned, type="phangorn::phyDat")

ParsTree <- pratchet(rRNA_Aligned3)
plot(ParsTree)

# Maximum likelihood

fit <- pml(rRNA_Tree, rRNA_Aligned3)

# There are other evolutionary models as well other than K80. Try them for comparison.

fitJC <- optim.pml(fit, model = "K80", rearrangement = "stochastic")
plot(fitJC)

# To add the bootstrapped values:

bootstrapped <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore = T, control = pml.control(trace=0))

plotBS(midpoint(fitJC$tree), bootstrapped, p = 10, type="p")

# To add a scale bar using ape: add.scale.bar(x, y, length = NULL, ask = FALSE, lwd = 1, lcol = "black", ...)

add.scale.bar(length = NULL, ask = F, cex = 0.7, font = 2, col = "black")

layout(1)






