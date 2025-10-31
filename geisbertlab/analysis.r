#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggsurvfit))
options(stringsAsFactors=FALSE)
theme_set(ggpubr::theme_pubr() + theme(legend.position="right"))

# helper variables
cols.outcome <- c("Survived"="black",
                  "Succumbed"="#e41a1c")
cols.dose <- c("27.5"="#000000",
               "10"="#737373",
               "1"="#bdbdbd",
               "0.1"="#f0f0f0")
#sample.dates <- c(0, 3, 4, 7, 10, 14, 21, 28, 35)
breaks.dates <- c(0, 7, 14, 21, 28, 35)

# load the full Excel table 
sheets <- readxl::excel_sheets("data.xlsx")
df <- lapply(sheets, function(i) { readxl::read_excel("data.xlsx", sheet=i)})
names(df) <- sheets
rm(sheets)

## the basics: survival, scores, and clinical signs ----------------------------
# format outcome and dose
df$animals <- df$animals %>%
              rename(Dose=`Dose (PFU)`,
                     Death=`Death DPI`) %>%
              mutate(Outcome=factor(Outcome, levels=names(cols.outcome)),
                     Dose=factor(Dose, levels=names(cols.dose)))

# format for KM curve
km <- df$animals %>%
      mutate(Outcome=as.integer(Outcome)-1)

# survival: KM curve and log-rank test
pval <- km %>%
        survival::survdiff(Surv(Death, Outcome) ~ Dose, data=.)
pval <- format.pval(pval$pvalue, digits=1)
pval <- paste0("p=", pval)

# add an "overall" condition then plot it
km %>%
  mutate(Dose="Overall") %>%
  rbind(km) %>%
  mutate(Dose=factor(Dose, levels=c(levels(df$animals$Dose), "Overall"))) %>%
  survfit2(Surv(Death, Outcome) ~ Dose, data=.) %>%
  ggsurvfit(theme=ggpubr::theme_pubr(), 
            linewidth=1) +
  scale_color_manual(values=c(cols.dose, Overall="#e41a1c")) +
  # add p-value
  annotate("text", x=35, y=1.05, hjust=1, label=pval) +
  scale_y_continuous(limits=c(0, 1.1), 
                     breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks=breaks.dates) +
  labs(x="Days postinfection",
       y="Survival probability", 
       col="PFU",
       title="Survival") +
  theme(legend.position="right")
ggsave("analysis/survival.png",
       units="in", width=4, height=2.5)
rm(pval, km)

# scores
# extract the challenge time and calculate time since challenge for all obs
baseline <- df$scores %>%
            filter(DPI==0) %>%
            rename(Challenge=Datetime) %>%
            select(Tattoo, Challenge)
df$scores %>%
  # add in challenge dose
  left_join(df$animals, by="Tattoo") %>%
  # calculate time (days) since exact challenge time
  left_join(baseline, by="Tattoo") %>%
  mutate(DPI=as.integer(Datetime-Challenge),
         DPI=DPI/(60*60*24)) %>%
  arrange(Dose) %>%
  ggplot(aes(DPI, Score)) +
  geom_hline(yintercept=9, col="lightgrey", linetype=2) +
  geom_line(aes(col=Dose, group=Tattoo)) +
  geom_point(aes(col=Dose), size=0.5) +
  scale_color_manual(values=cols.dose) +
  scale_x_continuous(breaks=breaks.dates) +
  labs(x="Days postinfection",
       y="Clinical score",
       col="PFU",
       title="Clinical scores") 
ggsave("analysis/score.png",
       units="in", width=4, height=2.5)
rm(baseline)

# decreased appetite
df$observations %>%
  left_join(df$animals, by="Tattoo") %>%
  group_by(Tattoo) %>%
  mutate(MaxBis=max(`Biscuits eaten`),
         Reduction=`Biscuits eaten`/MaxBis,
         Dose=paste(Dose, "PFU")) %>%
  ungroup() %>%
  select(Tattoo, Dose, DPI, Reduction, Death) %>%
  # order NHPs by day of death
  arrange(Death) %>%
  mutate(Tattoo=factor(Tattoo, levels=unique(Tattoo)),
         Dose=factor(Dose, levels=paste(levels(df$animals$Dose), "PFU"))) %>%
  # plot it
  ggplot(aes(DPI, Tattoo)) +
  geom_tile(aes(fill=Reduction)) +
  scale_fill_gradient(low="black", high="white",
                      labels=c(0, 50, 100),
                      breaks=c(0, 0.5, 1)) +
  scale_x_continuous(breaks=breaks.dates) +
  facet_wrap(~Dose, ncol=1, scales="free_y") +
  labs(x="Days postinfection",
       y=NULL,
       fill="Appetite (%)",
       title="Appetite") +
  theme(legend.position="bottom",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave("analysis/appetite.png",
       units="in", width=2.5, height=5)

## viremia ---------------------------------------------------------------------
# add challenge dose
viremia <- df$`virology-blood` %>%
           left_join(df$animals, by="Tattoo") %>%
           select(Tattoo, DPI, Outcome, `GEq/mL`, `PFU/mL`, Dose) %>%
           arrange(Dose)
# qRT-PCR
viremia %>%
  ggplot(aes(DPI, `GEq/mL`+1e2)) +
  geom_line(aes(col=Dose, group=Tattoo)) +
  geom_point(aes(col=Dose), size=0.5) +
  scale_color_manual(values=cols.dose) +
  scale_y_continuous(limits=c(20, 1e13), 
                     breaks=c(1e2, 1e4, 1e6, 1e8, 1e10, 1e12),
                     labels=c("<LoD", "1e04", "1e+06", 
                              "1e+08", "1e+10", "1e+12"),
                     trans="log10") +
  scale_x_continuous(breaks=breaks.dates) +
  labs(x="Days postinfection",
       y="LASV GEq/mL",
       col="PFU",
       title="LASV genomes") +
  theme(legend.position="none")
ggsave("analysis/viremia-qpcr-longitudinal.png",
       units="in", width=3.25, height=2)

# plaque assay
viremia %>%
  ggplot(aes(DPI, `PFU/mL`+24)) +
  geom_line(aes(group=Tattoo, col=Dose)) +
  geom_point(aes(col=Dose), size=0.521) +
  scale_color_manual(values=cols.dose) +
  scale_y_continuous(limits=c(20, 1e13), 
                     breaks=c(24, 1e4, 1e6, 1e8, 1e10, 1e12),
                     labels=c("<LoD", "1e04", "1e+06", 
                              "1e+08", "1e+10", "1e+12"),
                     trans="log10") +
  scale_x_continuous(breaks=breaks.dates) +
  labs(x="Days postinfection",
       y="LASV PFU/mL", 
       col="PFU",
       title="LASV titer") +
  theme(legend.position=c(0.8, 0.6))
ggsave("analysis/viremia-plaqueassay-longitudinal.png",
       units="in", width=3.25, height=2)

# viremia onset hypotheses: 
# (1) time of viremia onset differs by challenge dose or outcome
# (2) viral load at onset differs by challenge dose or outcome
onset <- viremia %>%
         reshape2::melt(id.vars=c("Tattoo", "DPI", "Dose", "Outcome"),
                        variable.name="Assay",
                        value.name="Viral.load") %>%
         mutate(Assay=factor(Assay, 
                             levels=c("GEq/mL", "PFU/mL"), 
                             labels=c("qRT-PCR", "Plaque assay"))) %>%
         filter(Viral.load > 0) %>%
         group_by(Tattoo, Outcome, Assay) %>%
         top_n(n=-1, wt=DPI) %>%
         ungroup() %>%
         rename(Onset=DPI) %>%
         select(Assay, Outcome, Dose, Onset, Viral.load)
# can't test statistically; n is too small for some groups
onset %>%
  group_by(Assay, Outcome) %>%
  summarise(NHPs=n())
onset %>%
  group_by(Assay, Dose) %>%
  summarise(NHPs=n())
rm(onset)

# peak viral load hypotheses: including zeroes from survivors with no viremia
peak <- viremia %>%
        reshape2::melt(id.vars=c("Tattoo", "DPI", "Dose", "Outcome"),
                       variable.name="Assay",
                       value.name="Viral.load") %>%
        mutate(Assay=factor(Assay, 
                            levels=c("GEq/mL", "PFU/mL"), 
                            labels=c("qRT-PCR", "Plaque assay")),
               Outcome=factor(Outcome, 
                              levels=c("Succumbed", "Survived"),
                              labels=c("Succ.", "Surv."))) %>%
        group_by(Tattoo, Outcome, Assay) %>%
        slice_max(order_by=Viral.load, n=1, with_ties=FALSE) %>%
        ungroup() %>%
        select(Assay, Outcome, Dose, Viral.load)
# (1) differs by challenge dose -- test first by KW, then pairwise
peak %>%
  group_by(Assay) %>%
  rstatix::kruskal_test(Viral.load ~ Dose) %>%
  rstatix::p_format(digits=1) 
peak %>%
  filter(Assay=="qRT-PCR") %>%
  ggplot(aes(Dose, Viral.load+1e2)) +
  geom_boxplot(fill="grey", outliers=FALSE, alpha=0.5) +
  geom_jitter(height=0, width=0.2, size=0.5) +
  ggpubr::stat_pwc(label="p.adj.signif", step.increase=0.2, hide.ns=TRUE) + 
  scale_y_continuous(limits=c(20, 1e18), 
                     breaks=c(1e2, 1e6, 1e10, 1e14),
                     labels=c("<LoD", "1e+06", "1e+10", "1e+14"),
                     trans="log10") +
  labs(x="Challenge (PFU)",
       y="LASV GEq/mL",
       title="Max genomes") +
  theme(legend.position="none")
ggsave("analysis/viremia-qpcr-peak-dose.png",
       units="in", width=2.16, height=2)
peak %>%
  filter(Assay=="Plaque assay") %>%
  ggplot(aes(Dose, Viral.load+24)) +
  geom_boxplot(fill="grey", outliers=FALSE, alpha=0.5) +
  geom_jitter(height=0, width=0.2, size=0.5) +
  ggpubr::stat_pwc(label="p.adj.signif", step.increase=0.5, hide.ns=TRUE) + 
  scale_y_continuous(limits=c(20, 1e18), 
                     breaks=c(24, 1e6, 1e10, 1e14),
                     labels=c("<LoD", "1e+06", "1e+10", "1e+14"),
                     trans="log10") +
  labs(x="Challenge (PFU)",
       y="LASV PFU/mL",
       title="Max titer") +
  theme(legend.position="none")
ggsave("analysis/viremia-plaqueassay-peak-dose.png",
       units="in", width=2.16, height=2)
# (2) differs by outcome -- test first by KW, then pairwise
peak %>%
  group_by(Assay) %>%
  rstatix::kruskal_test(Viral.load ~ Outcome) %>%
  rstatix::p_format(digits=1) 
peak %>%
  filter(Assay=="qRT-PCR") %>%
  ggplot(aes(Outcome, Viral.load+1e2)) +
  geom_boxplot(fill="grey", outliers=FALSE, alpha=0.5) +
  geom_jitter(height=0, width=0.2, size=0.5) +
  ggpubr::stat_pwc(label="p.signif") + 
  scale_y_continuous(limits=c(20, 1e18), 
                     breaks=c(1e2, 1e6, 1e10, 1e14),
                     labels=c("<LoD", "1e+06", "1e+10", "1e+14"),
                     trans="log10") +
  labs(x="Outcome",
       title="") +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())
ggsave("analysis/viremia-qpcr-peak-outcome.png",
       units="in", width=1.18, height=2)
peak %>%
  filter(Assay=="Plaque assay") %>%
  ggplot(aes(Outcome, Viral.load+24)) +
  geom_boxplot(fill="grey", outliers=FALSE, alpha=0.5) +
  geom_jitter(height=0, width=0.2, size=0.5) +
  ggpubr::stat_pwc(label="p.signif") + 
  scale_y_continuous(limits=c(20, 1e18), 
                     breaks=c(1e2, 1e6, 1e10, 1e14),
                     labels=c("<LoD", "1e+06", "1e+10", "1e+14"),
                     trans="log10") +
  labs(x="Outcome",
       title="") +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())
ggsave("analysis/viremia-plaqueassay-peak-outcome.png",
       units="in", width=1.18, height=2)

# clean up
rm(viremia, peak)

## tissue virology -------------------------------------------------------------
# are all statistically significantly different by outcome? Yes.
df$`virology-tissues` %>%
  filter(Tissue != "Aqueous humor") %>%
  left_join(df$animals, by="Tattoo") %>%
  group_by(Tissue) %>%
  rstatix::wilcox_test(`GEq/g` ~ Outcome) %>%
  rstatix::p_format()
  
# qRT-PCR only; not all tissues plaque-assayed
df$`virology-tissues` %>%
  filter(Tissue != "Aqueous humor") %>%
  left_join(df$animals, by="Tattoo") %>%
  ggplot(aes(Tissue, `GEq/g`+1e2)) +
  geom_boxplot(aes(fill=Outcome, col=Outcome), alpha=0.5, na.rm=TRUE) +
  scale_fill_manual(values=cols.outcome) +
  scale_color_manual(values=cols.outcome) +
  scale_y_continuous(limits=c(20, NA), 
                     breaks=c(1e2, 1e4, 1e6, 1e8, 1e10, 1e12),
                     labels=c("<LoD", "1e+04", "1e+06", 
                              "1e+08", "1e+10", "1e+12"),
                     trans="log10") +
  labs(x=NULL,
       y="LASV GEq/g") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="bottom")
ggsave("analysis/viremia-qpcr-tissues.png",
       units="in", width=6.5, height=4)

## clinical chemistry ----------------------------------------------------------
chem <- df$chemistry %>%
        # no change/values in uric acid, so remove that
        select(-`Uric acid (UA)`) %>%
        left_join(select(df$animals, Tattoo, Outcome, Dose), 
                  by="Tattoo") %>%
        reshape2::melt(id.vars=c("Tattoo", "DPI", "Outcome", "Dose"),
                       variable.name="Analyte",
                       value.name="Concentration")

# clinical chemistry spaghetti
spaghetti <- levels(chem$Analyte) %>%
            lapply(function(i) {
              # plot it
              chem %>%
                filter(Analyte==i,
                       !is.na(Concentration)) %>%
                ggplot(aes(DPI, Concentration)) +
                geom_line(aes(group=Tattoo, 
                              col=Dose)) +
                geom_point(aes(color=Dose), 
                           size=0.5,
                           na.rm=TRUE) +
                scale_color_manual(values=cols.dose) +
                scale_x_continuous(breaks=breaks.dates) +
                labs(x=NULL,
                     y=df$units[df$units$Analyte==i, "Units"], 
                     title=str_extract(i, "(?<=\\()[A-Z]+(?=\\))")) +
                theme(legend.position="none",
                      plot.title=element_text(size=12))
            })
names(spaghetti) <- str_extract(levels(chem$Analyte), "(?<=\\()[A-Z]+(?=\\))")
# pull out the legend
legs <- spaghetti$GLU +
        theme(legend.position="bottom") +
        labs(col="PFU")
legs <- cowplot::get_plot_component(legs, "guide-box-bottom")
xlabel <- spaghetti$GLU +
          labs(x="Days postinfection")
xlabel <- cowplot::get_plot_component(xlabel, "xlab-b")
# plot all chemistry together
x <- cowplot::plot_grid(plotlist=spaghetti, labels="AUTO", ncol=3)
cowplot::plot_grid(x, xlabel, legs, ncol=1, rel_heights=c(20, 1, 1))
ggsave("analysis/clinicalchemistry-spaghetti-all.png",
       units="in", width=6.5, height=7)

# save ALT, TP, and CRP separately for the main figure
# add legend to CRP
spaghetti$CRP <- spaghetti$CRP +
                 labs(col=NULL) +
                 theme(legend.position=c(0.85, 0.65))
# add x-label to all
cowplot::plot_grid(spaghetti$CRP + labs(x="Days postinfection"), 
                   spaghetti$ALT + labs(x="Days postinfection"), 
                   spaghetti$TP + labs(x="Days postinfection"), 
                   nrow=1)
ggsave("analysis/clinicalchemistry-spaghetti-selected.png",
       units="in", width=6.5, height=2)

# clean up
rm(spaghetti, legs, xlabel, x)

# are there differences in outcome at max value, or min value? 
# pull out baseline, then calculate differences and pull the "most different"
bl <- chem %>%
      filter(DPI==0) %>%
      mutate(Baseline=Concentration) %>%
      select(Tattoo, Analyte, Baseline)
comparison <- chem %>%
              left_join(bl, by=c("Tattoo", "Analyte")) %>%
              mutate(Delta=Concentration - Baseline) %>%
              group_by(Tattoo, Outcome, Analyte)
y <- comparison %>%
     top_n(n=-1, wt=Delta) %>%
     ungroup() %>%
     mutate(Comparison="Min")
comparison <- comparison %>%
              top_n(n=1, wt=Delta) %>%
              ungroup() %>%
              mutate(Comparison="Max") %>%
              rbind(y) %>%
              select(Outcome, Comparison, Analyte, Concentration)
# create a list of plots
comparison <- levels(chem$Analyte) %>%
              lapply(function(i) {
                 comparison %>%
                  filter(Analyte==i) %>%
                  ggplot(aes(Comparison, Concentration)) +
                  geom_boxplot(aes(fill=Outcome, col=Outcome), 
                               alpha=0.5, outlier.size=0.5) +
                  scale_fill_manual(values=cols.outcome) +
                  scale_color_manual(values=cols.outcome) +
                  ggpubr::stat_compare_means(aes(group=interaction(Comparison, 
                                                                   Outcome)), 
                                             method="wilcox.test",
                                             hide.ns=TRUE, 
                                             label="p.signif") +
                  scale_y_continuous(expand=expansion(mult=c(0.1, 0.2))) +
                  labs(x=NULL,
                       y=df$units[df$units$Analyte==i, "Units"],
                       title=str_extract(i, "(?<=\\()[A-Z]+(?=\\))")) +
                  theme(legend.position="none",
                        plot.title=element_text(size=12))
              })
names(comparison) <- str_extract(levels(chem$Analyte), "(?<=\\()[A-Z]+(?=\\))")
# pull out the legend
legs <- comparison$GLU + theme(legend.position="bottom")
legs <- cowplot::get_plot_component(legs, "guide-box-bottom")
# plot all chemistry together
x <- cowplot::plot_grid(plotlist=comparison, labels="AUTO", ncol=4)
cowplot::plot_grid(x, legs, ncol=1, rel_heights=c(15, 1))
ggsave("analysis/clinicalchemistry-peak-all.png",
       units="in", width=6.5, height=6)

# save ALT, TP, and CRP separately for the main figure
# format legend to right side
legs <- comparison$GLU + theme(legend.position="right")
legs <- cowplot::get_plot_component(legs, "guide-box-right")
cowplot::plot_grid(comparison$CRP, comparison$ALT, comparison$TP, legs, nrow=1)
ggsave("analysis/clinicalchemistry-peak-selected.png",
       units="in", width=6.5, height=2)

# clean up
rm(chem, x, bl, legs, comparison, y)
  
## hematology ------------------------------------------------------------------
hema <- df$hematology %>%
        select(-starts_with("Mean"), -starts_with("Hem")) %>%
        left_join(select(df$animals, Tattoo, Outcome, Dose),
                  by="Tattoo") %>%
        reshape2::melt(id.vars=c("Tattoo", "DPI", "Outcome", "Dose"),
                       variable.name="Analyte",
                       value.name="Concentration")

# hematology spaghetti
spaghetti <- levels(hema$Analyte) %>%
            lapply(function(i) {
              hema %>%
                filter(Analyte==i) %>%
                ggplot(aes(DPI, Concentration)) +
                geom_line(aes(group=Tattoo, col=Dose)) +
                geom_point(aes(col=Dose), size=0.5) +
                scale_color_manual(values=cols.dose) +
                scale_x_continuous(breaks=breaks.dates) +
                labs(x=NULL,
                     y=df$units[df$units$Analyte==i, "Units"],
                     title=i) +
                theme(legend.position="none",
                      plot.title=element_text(size=12))
            })
names(spaghetti) <- levels(hema$Analyte)
# pull out the legend
spaghetti$legs <- spaghetti$Leukocytes +
                  theme(legend.position="right") +
                  labs(col="PFU") 
spaghetti$legs <- cowplot::get_plot_component(spaghetti$legs, "guide-box-right")
xlabel <- spaghetti$Leukocytes +
          labs(x="Days postinfection")
xlabel <- cowplot::get_plot_component(xlabel, "xlab-b")
# plot all chemistry together
x <- cowplot::plot_grid(plotlist=spaghetti, 
                        labels=LETTERS[1:length(spaghetti)-1], ncol=3)
cowplot::plot_grid(x, xlabel, ncol=1, rel_heights=c(15, 1))
ggsave("analysis/hematology-spaghetti-all.png",
       units="in", width=6.5, height=6)

# save lymphocytes, neutrophils, and platelets separately for the main figure
cowplot::plot_grid(spaghetti$Lymphocytes + labs(x="Days postinfection"),
                   spaghetti$Neutrophils + labs(x="Days postinfection"),
                   spaghetti$Platelets + labs(x="Days postinfection"),
                   nrow=1)
ggsave("analysis/hematology-spaghetti-selected.png",
       units="in", width=6.5, height=2)

# clean up
rm(spaghetti, xlabel, x)

# are there differences in outcome at max or min value?
# pull out baseline, then calculate differences and pull the "most different"
bl <- hema %>%
      filter(DPI==0) %>%
      mutate(Baseline=Concentration) %>%
      select(Tattoo, Analyte, Baseline)
comparison <- hema %>%
              left_join(bl, by=c("Tattoo", "Analyte")) %>%
              mutate(Diff=Concentration - Baseline) %>%
              group_by(Tattoo, Outcome, Analyte)
y <- comparison %>%
     top_n(n=-1, wt=Diff) %>%
     ungroup() %>%
     mutate(Comparison="Min")
comparison <- comparison %>%
              top_n(n=1, wt=Diff) %>%
              ungroup() %>%
              mutate(Comparison="Max") %>%
              # add in minimum values
              rbind(y) %>%
              select(Outcome, Comparison, Analyte, Concentration) 
comparison <- levels(hema$Analyte) %>%
              lapply(function(i) {
                comparison %>%
                  filter(Analyte==i) %>%
                  ggplot(aes(Comparison, Concentration)) +
                  geom_boxplot(aes(fill=Outcome, col=Outcome), 
                               alpha=0.5, outlier.size=0.5) +
                  scale_fill_manual(values=cols.outcome) +
                  scale_color_manual(values=cols.outcome) +
                  ggpubr::stat_compare_means(aes(group=interaction(Comparison, 
                                                                   Outcome)), 
                                             method="wilcox.test",
                                             hide.ns=TRUE, 
                                             label="p.signif") +
                  scale_y_continuous(expand=expansion(mult=c(0.1, 0.2))) +
                  labs(x=NULL,
                       y=df$units[df$units$Analyte==i, "Units"],
                       title=i) +
                  theme(plot.title=element_text(size=12),
                        legend.position="none")
              })
names(comparison) <- levels(hema$Analyte)
# pull out the legend
legs <- comparison$Leukocytes + theme(legend.position="bottom")
legs <- cowplot::get_plot_component(legs, "guide-box-bottom")
# plot all hematology together
x <- cowplot::plot_grid(plotlist=comparison, 
                        labels="AUTO", 
                        ncol=4)
cowplot::plot_grid(x, legs, ncol=1, rel_heights=c(10, 1))
ggsave("analysis/hematology-peak-all.png",
       units="in", width=6.5, height=4)

# save lymphocytes, neutrophils, and platelets separately for the main figure
legs <- comparison$Lymphocytes + theme(legend.position="right")
legs <- cowplot::get_plot_component(legs, "guide-box-right")
cowplot::plot_grid(comparison$Lymphocytes,
                   comparison$Neutrophils,
                   comparison$Platelets,
                   legs, nrow=1)
ggsave("analysis/hematology-peak-selected.png",
       units="in", width=6.5, height=2)

# clean up
rm(hema, bl, xlabel, spaghetti, comparison, y, legs)

## PRNTs -----------------------------------------------------------------------
# get the average plaques per sample
prnt <- df$PRNTs %>%
        reshape2::melt(measure.vars=c("Plaques A", "Plaques B"),
                       value.name="Plaques") %>%
        group_by(Tattoo, DPI, Batch, `Dilution factor`) %>%
        summarise(Plaques=mean(Plaques),
                  .groups="drop")
# pull out VC and calculate the baseline for each batch
vc <- prnt %>%
      filter(Tattoo=="Control") %>%
      rename(VC=Plaques) %>%
      select(Batch, VC)
# calculate % reduction
prnt <- prnt %>%
        filter(Tattoo != "Control") %>%
        # add in VC counts and calculate % reduction
        left_join(vc, by="Batch") %>%
        mutate(Reduction=VC-Plaques,
               PercentReduction=100*Reduction/VC,
               # update DPI to factor
               DPI=factor(DPI, levels=c(0, 14, "Terminal", 35))) %>%
        replace_na(list(DPI="Terminal"))
# bound the % reduction to [0, 100]
prnt$PercentReduction[prnt$PercentReduction < 0] <- 0

# plot each NHP by dose
prnt %>%
  left_join(df$animals, by="Tattoo") %>%
  mutate(Dose=factor(Dose, labels=paste(levels(Dose), "PFU"))) %>%
  arrange(desc(DPI)) %>%
  ggplot(aes(`Dilution factor`, PercentReduction)) +
  geom_hline(yintercept=50, linetype=3) +
  geom_line(aes(group=interaction(Tattoo, DPI), col=DPI)) +
  geom_point(aes(fill=DPI), col="black", size=2, pch=21) +
  scale_fill_manual(values=c("white", "darkgrey", 
                             "#e41a1c", "black")) +
  scale_color_manual(values=c("white", "darkgrey", 
                              "#e41a1c", "black")) +
  scale_x_continuous(breaks=unique(prnt$`Dilution factor`),
                     trans="log2") +
  ylim(0, 100) +
  facet_wrap(~Dose, nrow=1) +
  labs(y="Neutralization (%)") +
  guides(fill=guide_legend(override.aes=list(pch=21)),
         shape=guide_legend(override.aes=list(fill="black"))) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10))
ggsave("analysis/prnts-curves.png",
       units="in", width=6.5, height=2)

# plot neutralization for each at the lowest dilution and compare outcomes
prnt %>%
  filter(`Dilution factor`==10,
         DPI %in% c("Terminal", "35")) %>%
  group_by(Tattoo) %>%
  slice_max(order_by=PercentReduction, n=1, with_ties=FALSE) %>%
  ungroup() %>%
  left_join(df$animals, by="Tattoo") %>%
  ggplot(aes(Outcome, PercentReduction)) +
  geom_boxplot(aes(fill=Outcome), alpha=0.5, outliers=FALSE) +
  geom_jitter(aes(fill=Outcome), pch=21, height=0, width=0.2) +
  scale_fill_manual(values=cols.outcome) +
  ylim(0, 100) +
  ggpubr::stat_pwc(method="wilcox.test", label="p.signif") +
  labs(y="Neutralization (%)",
       x="Outcome",
       title="LASV neutralization") +
  theme(legend.position="none")
ggsave("analysis/prnts-max.png",
       units="in", width=2.5, height=2.5)

# clean up
rm(prnt, vc)

## ELISAs ----------------------------------------------------------------------
# over time -- very messy
df$ELISAs %>%
  reshape2::melt(measure.vars=c("IgM", "IgG"),
                 value.name="Titer") %>%
  left_join(df$animals, by="Tattoo") %>%
  ggplot(aes(DPI, Titer+100)) +
  geom_line(aes(group=Tattoo, col=Outcome)) +
  geom_point(aes(fill=Outcome), pch=21) +
  scale_fill_manual(values=cols.outcome) +
  scale_color_manual(values=cols.outcome) +
  scale_y_continuous(limits=c(NA, 30000),
                     breaks=c(100, 200, 800, 3200, 12800),
                     labels=c("<LOD", 200, 800, 3200, 12800),
                     trans="log2") +
  facet_grid(variable ~ Antigen)

# is there any difference between outcomes at endpoints?
elisas <- df$ELISAs %>%
          group_by(Tattoo, Antigen) %>%
          top_n(n=1, wt=DPI) %>%
          reshape2::melt(measure.vars=c("IgM", "IgG"),
                         variable.name="Immunoglobulin",
                         value.name="Titer") %>%
          left_join(df$animals, by="Tattoo")
# get p-value significance
psig <- elisas %>%
        mutate(Titer=log2(Titer+1)) %>%
        group_by(Immunoglobulin, Antigen) %>%
        rstatix::wilcox_test(Titer ~ Outcome) %>%
        rstatix::adjust_pvalue() %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(x="Antigen")
# add pseudocount (100) to samples <LOD
elisas$Titer[elisas$Titer==0] <- 100
# plot it and add statistical comparisons
elisas %>%
  ggplot(aes(Antigen, Titer)) +
  geom_boxplot(aes(group=interaction(Antigen, Outcome), fill=Outcome), 
               alpha=0.5, outliers=FALSE) +
  geom_point(aes(fill=Outcome, group=interaction(Antigen, Outcome)), 
             pch=21, position=position_jitterdodge(jitter.width=0.2)) +
  scale_fill_manual(values=cols.outcome) +
  scale_y_continuous(breaks=c(100, 200, 800, 3200, 12800),
                     labels=c("<LOD", 200, 800, 3200, 12800),
                     trans="log2",
                     expand = expansion(mult = c(0.05, 0.1))) +
  ggpubr::stat_pvalue_manual(psig, label="p.adj.signif") +
  facet_wrap(~Immunoglobulin) +
  labs(y="Endpoint titer",
       x="LASV antigen",
       title="ELISAs")
ggsave("analysis/elisas-endpoint.png",
       units="in", width=4, height=2.5)

# clean up
rm(psig, elisas)

## fin -------------------------------------------------------------------------
sessionInfo()
