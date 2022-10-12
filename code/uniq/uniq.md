Unique gene gradient
================

### Description: R code for wrangling and plotting data outputted from Roary as the fraction of unique genes per MAG

### Load dependencies

``` r
library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
```

### Load data

#### Data is the presence/absence dataframes, where the number corresponds to the clustering identity threshold

``` r
grad_50 <- read.delim("50.csv", sep=",", header=F)
grad_60 <- read.delim("60.csv", sep=",", header=F)
grad_70 <- read.delim("70.csv", sep=",", header=F)
grad_80 <- read.delim("80.csv", sep=",", header=F)
grad_90 <- read.delim("90.csv", sep=",", header=F)
grad_95 <- read.delim("95.csv", sep=",", header=F)
grad_99 <- read.delim("99.csv", sep=",", header=F)
```

### Initial wrangling of data

``` r
grad_wrangle <- function(grad) {
names(grad) <- lapply(grad[1, ], as.character)
grad <- grad[-1,]
grad[,2:14] <- NULL
return(grad) }

grad_50 <- grad_wrangle(grad_50)
grad_60 <- grad_wrangle(grad_60)
grad_70 <- grad_wrangle(grad_70)
grad_80 <- grad_wrangle(grad_80)
grad_90 <- grad_wrangle(grad_90)
grad_95 <- grad_wrangle(grad_95)
grad_99 <- grad_wrangle(grad_99)
```

### Aggregate counts

``` r
grad_aggr <- function(grad,ident) {
grad_ <- grad %>% pivot_longer(!Gene, names_to = "bin", values_to = "id", values_drop_na = TRUE) %>% mutate_all(na_if,"")
grad_  <- na.omit(grad_)
grad_count <- aggregate(grad_$id, by=list(grad_$Gene), FUN=length)
colnames(grad_count) <- c("Gene","count")
grad_ <- merge(grad_,grad_count, by="Gene")
grad_$identity <- ident
return(grad_) }

grad <- rbind(grad_aggr(grad_50,50),grad_aggr(grad_60,60),grad_aggr(grad_70,70),
              grad_aggr(grad_80,80),grad_aggr(grad_90,90),grad_aggr(grad_95,95),grad_aggr(grad_99,99))
```

### Generate matrix: single MAGs

``` r
grad_uniq <- grad[(grad$count == 1 & (grad$bin == "ENR-cMAG" | grad$bin == "BRK-cMAG" |
                                        grad$bin == "MAR-scMAG" | grad$bin == "MAR-hqMAG" | grad$bin == "MAR-mqMAG" )), ]
grad_uniq$count_all <- str_count(grad_uniq$id, '\\w+')
grad_uniq_sum <- aggregate(grad_uniq$count_all, by=list(grad_uniq$bin,grad_uniq$identity), FUN=sum)
colnames(grad_uniq_sum) <- c("bin","identity","sum")
```

### Generate matrix: custom MAG groups

``` r
grad_substr <- function(grad, group, name) {
  grad <- grad[(grad$bin %in% group),]
  grad$Gene2 <- paste(grad$Gene,grad$identity,sep="_")
  grad$count_all <- str_count(grad$id, '\\w+')
  grad$count_all_ <- ifelse(grad$count_all >= 1,1,0)
  grad2 <- grad %>% select(bin, Gene2, count_all_) %>% pivot_wider(names_from = bin, values_from = count_all_, values_fill = 0)
  grad2 <- grad2[,c("Gene2",group)]
  
  grad2 <- grad2[(grad2[,2] == 1), ]
  grad2$total <- rowSums(grad2[,group])
  grad2 <- grad2 %>% select(Gene2, total)

  grad <- merge(grad,grad2,by="Gene2")
  grad <- grad[(grad$count == grad$total & grad$bin == group[1]), ]
  grad <- aggregate(grad$count_all, by=list(grad$identity), FUN=sum)
  colnames(grad) <- c("identity","sum")
  grad$bin <- name
  grad <- grad[, c("bin","identity","sum")]

return(grad) }

# Fist one in the list is the central one
grad_uniq_brk1 <- grad_substr(grad,c("BRK-cMAG","KV"),"BRK-cMAG (KV-excluded)") 
grad_uniq_brk2 <- grad_substr(grad,c("BRK-cMAG","ENR-cMAG","KV","SY1","UQ","GS","F5"),"BRK-cMAG (Electronema-excluded)") 
grad_uniq_sum <- rbind(grad_uniq_sum,grad_uniq_brk1,grad_uniq_brk2)
```

### Get total CDS counts and make fractions

``` r
grad_uniq_50 <- grad[(grad$identity == 50), ]
grad_uniq_50$count_all <- str_count(grad_uniq_50$id, '\\w+')
grad_cds <- aggregate(grad_uniq_50$count_all, by=list(grad_uniq_50$bin), FUN=sum)
colnames(grad_cds) <- c("bin","CDS")

# Add custom groups
grad_cds <- rbind(grad_cds,c("BRK-cMAG (KV-excluded)",grad_cds[(grad_cds$bin == "BRK-cMAG"), ]$CDS))
grad_cds <- rbind(grad_cds,c("BRK-cMAG (Electronema-excluded)",grad_cds[(grad_cds$bin == "BRK-cMAG"), ]$CDS)) 

grad_uniq_sum2 <- merge(grad_uniq_sum,grad_cds,by="bin")
grad_uniq_sum2$CDS <- as.integer(grad_uniq_sum2$CDS)
grad_uniq_sum2$fraction <- round(grad_uniq_sum2$sum / grad_uniq_sum2$CDS * 100,2)
```

### Plot matrix

``` r
level_order1 <- c("ENR-cMAG","MAR-hqMAG","MAR-mqMAG","MAR-scMAG","BRK-cMAG","BRK-cMAG (KV-excluded)","BRK-cMAG (Electronema-excluded)")
level_order2 <- c("50","60","70","80","90","95","99")

plot_uniq <- ggplot(data = grad_uniq_sum2, mapping = aes(x=factor(identity, level=level_order2),
                                                        y=factor(bin, level=rev(level_order1)), fill = fraction))  +
  geom_tile(width=1, color="black") + geom_text(aes(label = fraction), size=3.5) + 
  labs(x="Gene identity threshold (%)", y="Recovered MAGs", title="Unique gene fraction (%)") +
  scale_fill_gradient2(midpoint = 50, high = "#1c9099", mid= "#a6bddb", low = "#fff7fb", limits=c(0,100)) +
  theme_bw() + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  theme(legend.position = "none", legend.text=element_text(size=10), legend.title=element_text(size=12),
        axis.title.y = element_text(size = 12), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=12), plot.title = element_text(hjust = 0.5, size=13))

ggsave(plot_uniq, file="uniq_gene_frac.pdf",height = 6 , width = 8, useDingbats=FALSE)
```
