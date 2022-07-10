## Plot the MCC tree from BEAST
## Kaichi Huang 2022 May

library(ggtree)
library(treeio)
library(tidyverse)

mu <- 6.1e-3 # mutation rate of the studied group, here we use 6.1e-3 substitution/site/Myr
my.labels <- c('H. annuus type 0','H. annuus type 2','H. argophyllus', 'H. petiolaris petiolaris','H. petiolaris fallax', 'H. niveus canescens','Perennial outgroup')
my.colors <- c("#FFCE6E","#D39117","#338732","#4BB7B7","#447499","#645084","#000000")

# Read in sample (group) information
sample_info <- read_tsv("input/sample_info.txt", col_names=F) %>% rename(sample=X1, group=X2)
sample_info$group <- factor(sample_info$group, levels=c("Ann-0", "Ann-2", "Arg", "PetPet", "PetFal", "PetCan", "outgroup"), labels=my.labels)

# Read in the MCC tree file
tree <- read.beast("output/out.mcc.tre")

pdf("output/BEAST_plot.pdf", width=7, height=5)
p <- ggtree(tree)
print(
  revts(p) %<+% sample_info +
    geom_text2(aes(subset = !isTip, label=label)) + theme_tree2() +
    geom_tippoint(aes(color=group), size=3) +
    geom_range(range='height_0.95_HPD', color='purple', alpha=.6, size=2) +
    geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
    scale_color_manual(values=my.colors, name="") +
    scale_x_continuous(limits=c(-2*2.5*mu, 0), breaks=seq(-2*2*mu, 0, 2*mu), minor_breaks=seq(-2*2.5*mu, 0, 2*mu*0.5), labels=-seq(-2,0,1)) +
    xlab("Mya") +
    theme(legend.position="right") +
    theme(panel.grid.major = element_line(color="grey", size=.2), panel.grid.minor = element_line(color="grey", size=.2), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
)
dev.off()
