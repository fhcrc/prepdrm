#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(reshape2)

fisher <- read.csv('output/prep_merged_fisher.csv', as.is=TRUE)
fisher <- fisher[fisher$is_control & !fisher$is_wt_control & !is.na(fisher$p_adj), ]
logistic <- read.csv('output/prep_merged_logistic.csv', as.is=TRUE)
logistic <- logistic[logistic$is_control & !logistic$is_wt_control & !is.na(logistic$lrt_p_adj), ]

df <- rbind(with(fisher,
                 data.frame(location=location, sample_name=sample_name,
                            p_value=p_value,
                            p_adj=p_adj,
                            mut=mut,
                            nonmutant=nonmutant,
                            method='fisher')),
            with(logistic, data.frame(location=location, sample_name=sample_name,
                      p_value=lrt_p_value,
                      p_adj=lrt_p_adj,
                      mut=mut,
                      nonmutant=nonmutant,
                      method='logistic')))

tb <- xtabs(~ method + I(p_adj < 0.05) + location, df)

melted <- melt(df, measure.vars=c('p_value', 'p_adj', 'mut', 'nonmutant'))
casted <- na.omit(dcast(melted, ... ~ variable+method))
casted <- subset(casted, location != 'A98G')

# Discrepancies
cutoff <- 0.05
casted <- transform(casted,
                    enriched_fisher=p_adj_fisher < cutoff,
                    enriched_logistic=p_adj_logistic < cutoff)
discr <- casted[casted$enriched_fisher != casted$enriched_logistic, ]
