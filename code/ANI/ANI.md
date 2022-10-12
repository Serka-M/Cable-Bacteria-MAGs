ANI
================

### Description: R code for wrangling and plotting data outputted from pyani

### Load dependencies

``` r
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
```

### Load data

``` r
ani_ident <- read.delim("ANIb_percentage_identity.tab", sep="\t", header=T)

ani_cov <- read.delim("ANIb_alignment_coverage.tab", sep="\t", header=T)

level_order <- read.delim("order.txt", sep="\t", header=F)
```

### Wrangle data

``` r
names(ani_ident)[2:(nrow(ani_ident)+1)] <- as.vector(ani_ident[,1])
ani_ident <- ani_ident %>% pivot_longer(!X, names_to = "Y", values_to = "ident")

names(ani_cov)[2:(nrow(ani_cov)+1)] <- as.vector(ani_cov[,1])
ani_cov <- ani_cov %>% pivot_longer(!X, names_to = "Y", values_to = "cov")

ani <- cbind(ani_ident, ani_cov$cov)
colnames(ani)[4] <- "cov"

level_order <- c(level_order$V1)
```

### Plot data

``` r
plot_ident <- ggplot(data = ani, mapping = aes(x = factor(X, level=level_order), y = factor(Y, level=rev(level_order)), fill = ident*100))  +
  geom_tile(width=1) + geom_text(aes(label = round(ident*100, 0)), size=2.25) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + labs(x="",y="", title="ANIb, %") + theme_bw() +
  scale_fill_gradient2(midpoint = 80, low = "#1d91c0", mid="white", high = "#006d2c", limits=c(60,100)) +
   theme(legend.position = "none", legend.text=element_text(size=10), legend.title=element_text(size=12),
        axis.title.y = element_text(size = 12), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=12), plot.title = element_text(hjust = 0.5, size=12))

plot_cov <- ggplot(data = ani, mapping = aes(x = factor(X, level=level_order), y = factor(Y, level=rev(level_order)), fill = cov*100))  +
  geom_tile(width=1) + geom_text(aes(label = round(cov*100, 0)), size=2.25) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + labs(x="",y="", title="Alignment coverage, %") + theme_bw() +
  scale_fill_gradient2(midpoint = 55, low = "white", mid= "#fed976", high = "#800026", limits=c(0.01,100)) +
   theme(legend.position = "none", legend.text=element_text(size=10), legend.title=element_text(size=12),
        axis.title.y = element_text(size = 12), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=12), plot.title = element_text(hjust = 0.5, size=12))

plot_anib <- ggarrange(plot_ident,
                       plot_cov, 
                       nrow=2, ncol=1, align = c("hv"), legend = "none", common.legend = FALSE)

plot_anib
```

![](ANI_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
#ggsave(file="anib.pdf", width = 9, height = 17, useDingbats=FALSE)
```
