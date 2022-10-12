AAI
================

### Description: R code for wrangling and plotting data outputted from CompareM

### Load dependencies

``` r
library(ggplot2)
library(ggasym)
```

### Load data

``` r
aai <- read.table("aai_summary.tsv", sep="\t", header=F)
level_order <- read.delim("order.txt", sep="\t", header=F)
```

### Wrangle data

``` r
colnames(aai) <- c("genome1","no_genes1","genome2","no_genes2","no_orthologs","AAI","AAI_SD","OF")
aai$POCP <- (2*aai$no_orthologs)/(aai$no_genes1+aai$no_genes2)*100

aai <- asymmetrise(aai, genome1, genome2)

level_order <- c(level_order$V1)
```

### Make intial plot

``` r
plot_aai <- ggplot(aai, aes(x = factor(genome1, level=rev(level_order)), y = factor(genome2, level=rev(level_order)))) +
  geom_asymmat(aes(fill_tl = AAI, fill_br = POCP, fill_diag="white")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) + labs(x="",y="", title="") + theme_bw() +
  scale_fill_tl_gradient2(midpoint = 78, low = "#1d91c0", mid="white", high = "#006d2c", limits=c(45,100)) +
  scale_fill_br_gradient2(midpoint = 45, low = "white", mid= "#fed976", high = "#800026", limits=c(0.01,100)) + 
  theme(legend.position = "right", legend.text=element_text(size=10), legend.title=element_text(size=12),
        axis.title.y = element_text(size = 12), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=12), plot.title = element_text(hjust = 0.5, size=12)) 
```

### Add labels to assymetric plot

``` r
plot_aai2 <- ggplot(aai, aes(x = factor(genome1, level=rev(level_order)), y = factor(genome2, level=rev(level_order)))) +
  geom_asymmat(aes(fill_tl=AAI, fill_br=POCP, fill_diag="white")) 

aai_list <- ggplot_build(plot_aai2)
aai_df1 <- aai_list$data[[1]]
aai_df2 <- aai_list$data[[2]]

aai_df1 <- aai_df1[c("x","y","fill_tl")]
colnames(aai_df1) <- c("x","y","value")

aai_df2 <- aai_df2[c("x","y","fill_br")]
colnames(aai_df2) <- c("x","y","value")

aai_df <- rbind(aai_df1,aai_df2)

plot_aai3 <- plot_aai + geom_text(data=aai_df, aes(x=x, y=y, label = round(value,0)), size=2.25)

ggsave(file="aai_pocp.pdf", height = 9, width = 10, useDingbats=FALSE)
```
