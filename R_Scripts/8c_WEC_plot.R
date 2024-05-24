# Plot Treatment effect on diversity by age group

rm(list=ls())

# libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(lme4)
library(lmerTest)
library(vegan)
library(compositions)
library(wec)
library(tm)
library(interactions)
library(broom.mixed)
library(jtools)

#install.packages("huxtable")

# data read-in

ps <- readRDS("Data/Phylo_objects/phylo-experimental-decontamPrevalence.rds")
metadata <- microbiome::meta(ps)

# set age as factor
metadata$Age.category <- factor(metadata$Age.category, levels = c("D8", "D15", "Adult"))

## WEC

#Create weighted contrast encoding, specifying "Adult" as omitted category in first variable and D8 omitted in second variable
metadata$Age.category.wec1 <- metadata$Age.category
contrasts(metadata$Age.category.wec1) <- contr.wec(metadata$Age.category, "Adult")
contrasts(metadata$Age.category.wec1)

metadata$Age.category.wec2 <- metadata$Age.category
contrasts(metadata$Age.category.wec2) <- contr.wec(metadata$Age.category, "D8")
contrasts(metadata$Age.category.wec2)


# Create new wec variable for interaction
metadata$Treatment.wec <- as.factor(metadata$Treatment)
contrasts(metadata$Treatment.wec) <- contr.wec(metadata$Treatment.wec, omitted = "Control")
contrasts(metadata$Treatment.wec)

metadata$Age.Treat.int1 <-  wec.interact(metadata$Age.category.wec1, metadata$Treatment.wec)
metadata$Age.Treat.int2 <-  wec.interact(metadata$Age.category.wec2, metadata$Treatment.wec)

## all ages chao WEC model
## age*treatment interaction

log.chao.lm.interaction.wec1 <- lmer(log(Chao1) ~ Treatment.wec + Age.category.wec1 + Age.Treat.int1 + (1|Site/Nest/Bird.ID) + (1|Plate),
                                     data = metadata, control = lmerControl(optimizer = "bobyqa"))
summary(log.chao.lm.interaction.wec1)

## save to make convenient later
w1 <- log.chao.lm.interaction.wec1 

log.chao.lm.interaction.wec2 <- lmer(log(Chao1) ~ Treatment.wec + Age.category.wec2 + Age.Treat.int2 + (1|Site/Nest/Bird.ID) + (1|Plate), data = metadata,
                                     control = lmerControl(optimizer = "bobyqa"))
summary(log.chao.lm.interaction.wec2)

## save to make convenient later
w2 <- log.chao.lm.interaction.wec2

#######
## plotting wec results

## collect all the estimates from your model in a handy table, thanks Ben Bolker
x <-broom.mixed::tidy(w1)
x2 <- broom.mixed::tidy(w2)

## unify your two sets of wec model estimates
x <- unique(rbind(x, x2))
x$estimate <- round(x$estimate, 5)
x$std.error <- round(x$std.error, 5)
x$statistic <- round(x$statistic, 5)
x$df <- round(x$df, 5)
x$p.value <- round(x$p.value, 5)
x <- unique(x)

## rename your terms so they are more meaningful and also the same across models
x2 <- unique(as.data.frame(cbind(c("(Intercept)", "Treatment.wecL. kimchicus", "Age.category.wec1D8", 
                                   "Age.category.wec1D15", "Age.Treat.int1x1D8:x2L. kimchicus", 
                                   "Age.Treat.int1x1D15:x2L. kimchicus", "sd__(Intercept)", "sd__(Intercept)", 
                                   "sd__(Intercept)", "sd__(Intercept)", "sd__Observation", "Age.category.wec2D15", 
                                   "Age.category.wec2Adult", "Age.Treat.int2x1D15:x2L. kimchicus", 
                                   "Age.Treat.int2x1Adult:x2L. kimchicus"),
                                 c("(Intercept)", "Treatment", "D8", 
                                   "D15", "TreatmentD8", 
                                   "TreatmentD15", "sd__(Intercept)", "sd__(Intercept)", 
                                   "sd__(Intercept)", "sd__(Intercept)", "sd__Observation", "D15", 
                                   "Adult", "TreatmentD15", 
                                   "TreatmentAdult"))))
names(x2) <- c("term", "name")

## quick tidy so we only have entries that we need
x <- left_join(x, x2)
x <- unique(select(x, -c("term")))
x <- x[x$effect == "fixed",]

# drop duplicated terms
x <- x[c(1,2,3,4,5,6,8,10),]
## add up your estimates so that you get the actual value that your model is assigning to each term
## helpful to make a new column (eff_sz) so you can check your maths
x$eff_sz <- x$estimate
x$eff_sz[c(3:8)] <- x$eff_sz[c(3:8)] + x$estimate[1] ## add the grand mean to everything
## add the effect of treatment and the effect of age to the treatment:age interaction 
x$eff_sz[5] <-x$eff_sz[5] + x$estimate[2] + x$estimate[3] 
x$eff_sz[6] <-x$eff_sz[6] + x$estimate[2] + x$estimate[4] 
x$eff_sz[8] <-x$eff_sz[8] + x$estimate[2] + x$estimate[7] 


# re-order x axis
x$name <- factor(x$name, levels = c("(Intercept)", "Treatment", "D8", "D15", "Adult", "TreatmentD8", "TreatmentD15", "TreatmentAdult"))

## plot these estimates
w_plot <- ggplot(x[x$name %in% c("D8", "D15", "Adult", "TreatmentD8", "TreatmentD15", "TreatmentAdult"),], aes(name, eff_sz)) + 
  geom_hline(yintercept= x$eff_sz[x$name == "(Intercept)"], alpha = 0.8) + ## adding the grand mean to the plot for visualization
  geom_line() +
  ## plot only the terms of interest, i.e. not the effect size of Treatment alone
  geom_pointrange(data = x[x$name %in% c("D8", "D15", "Adult", "TreatmentD8", "TreatmentD15", "TreatmentAdult"),], 
                  aes(name, eff_sz, ymin=eff_sz - std.error, ymax=eff_sz + std.error))

w_plot

metadata$plotVar <- ifelse(metadata$Treatment=="L. kimchicus", paste0("Treatment", as.character(metadata$Age.category.wec1)), as.character(metadata$Age.category.wec1))

w_plot.Raw <- w_plot +
  geom_point(data = metadata, aes(plotVar, log(Chao1), colour = Age.category), alpha = 0.5) +
  scale_colour_manual(values = c("D8" = "#49b7fc",
                                 "D15" = "#ff7b00",
                                 "Adult" = "#17d898"))

w_plot.Raw

# i want to add partial residuals to this model instead of the raw data
# extract and plot partial residuals
## this is making all the bird IDs CB02A for some reason
partial.resids.w1 <- partialize(w1, vars = c("Treatment.wec", "Age.category.wec1", "Age.Treat.int1"), data = w1@frame)

# add plotvar
partial.resids.w1$plotVar <- ifelse(partial.resids.w1$Treatment.wec=="L. kimchicus", paste0("Treatment", 
                                                                                            as.character(metadata$Age.category.wec1)), as.character(metadata$Age.category.wec1))
## final product
## alt plot with 95% CI's instead of SE's
## could use geom_errorbar instead of geom_pointrange to get different line type
w_plot.Final <- ggplot(x[x$name %in% c("D8", "D15", "Adult", "TreatmentD8", "TreatmentD15", "TreatmentAdult"),], aes(name, eff_sz)) + 
  geom_hline(yintercept= x$eff_sz[x$name == "(Intercept)"], alpha = 0.6) + ## adding the grand mean to the plot for visualization
  geom_line() +
  # add raw datapoints
  geom_point(data = partial.resids.w1, aes(plotVar,`log(Chao1)`, colour = Age.category.wec1), alpha = 0.5) +
  scale_colour_manual(values = c("D8" = "#49b7fc",  "D15" = "#ff7b00",  "Adult" = "#17d898")) +
  ## plot only the terms of interest, i.e. not the effect size of Treatment alone
  geom_pointrange(data = x[x$name %in% c("D8", "D15", "Adult", "TreatmentD8", "TreatmentD15", "TreatmentAdult"),], 
                  aes(name, eff_sz, ymin=eff_sz - (std.error*1.96), ymax=eff_sz + (std.error*1.96))) +
  ylab("log(Chao1)") +
  xlab("") +
  theme_classic() +
  labs(color = "Age category") +
  # add custom x axis labels
  scale_x_discrete(breaks = c("D15", "TreatmentD15"), labels = c("Control", "Treatment"))

w_plot.Final

ggsave("Figure2_WEC.pdf", plot = w_plot.Final, device = "pdf", path = "Figures/", width = 11, height = 11, units = "cm", dpi = 300)
