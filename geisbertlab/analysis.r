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
sample.dates <- c(0, 3, 4, 7, 10, 14, 21, 28, 35)

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

# survival: KM curve and log-rank test
pval <- df$animals %>%
        mutate(Outcome=as.integer(Outcome)-1) %>%
        survival::survdiff(Surv(Death, Outcome) ~ Dose, data=.)
pval <- format.pval(pval$pvalue, digits=2)
pval <- paste0("p=", pval)

# plot it
df$animals %>%
  mutate(Outcome=as.integer(Outcome)-1) %>%
  survfit2(Surv(Death, Outcome) ~ Dose, data=.) %>%
  ggsurvfit::ggsurvfit(theme=ggpubr::theme_pubr(), 
                       size=1) +
  scale_color_manual(values=cols.dose) +
  # add p-value
  annotate("text", x=35, y=1.05, hjust=1, label=pval) +
  scale_y_continuous(limits=c(0, 1.1), 
                     breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks=c(0, 7, 14, 21, 28, 35)) +
  labs(x="Days postinfection",
       y="Survival probability", 
       col="Challenge (PFU)") +
  theme(legend.position=c(0.8, 0.4))
ggsave("analysis/survival.png",
       units="in", width=4, height=3)
rm(pval)

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
  geom_point(aes(fill=Dose), pch=21) +
  scale_color_manual(values=cols.dose) +
  scale_fill_manual(values=cols.dose) +
  scale_y_continuous(limits=c(0, 20)) +
  scale_x_continuous(breaks=c(0, 7, 14, 21, 28, 35)) +
  labs(x="Days postinfection",
       y="Clinical score",
       col="Challenge (PFU)",
       fill="Challenge (PFU)") +
  theme(legend.position=c(0.8, 0.8))
ggsave("analysis/score.png",
       units="in", width=6, height=4)
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
  mutate(Tattoo=factor(Tattoo, levels=unique(Tattoo))) %>%
  # plot it
  ggplot(aes(DPI, Tattoo)) +
  geom_tile(aes(fill=Reduction)) +
  scale_fill_gradient(low="black", high="white") +
  scale_x_continuous(breaks=c(0, 7, 14, 21, 28, 35)) +
  facet_wrap(~Dose, ncol=1, scales="free_y") +
  labs(x="Days postinfection",
       y=element_blank(),
       fill="Appetite")
ggsave("analysis/appetite.png",
       units="in", width=3.75, height=6)

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
  geom_point(aes(fill=Dose), pch=21) +
  scale_color_manual(values=cols.dose) +
  scale_fill_manual(values=cols.dose) +
  scale_y_continuous(limits=c(20, 1e13), 
                     breaks=c(1e2, 1e4, 1e6, 1e8, 1e10, 1e12),
                     labels=c("<LoD", "1e04", "1e+06", 
                              "1e+08", "1e+10", "1e+12"),
                     trans="log10") +
  scale_x_continuous(breaks=sample.dates) +
  labs(x="Days postinfection",
       y="LASV GEq/mL",
       col="Challenge (PFU)",
       fill="Challenge (PFU)")
ggsave("analysis/viremia-qpcr-longitudinal.png",
       units="in", width=6, height=4)

# plaque assay
viremia %>%
  ggplot(aes(DPI, `PFU/mL`+24)) +
  geom_line(aes(group=Tattoo, col=Dose)) +
  geom_point(aes(fill=Dose), pch=21) +
  scale_color_manual(values=cols.dose) +
  scale_fill_manual(values=cols.dose) +
  scale_y_continuous(limits=c(20, 1e13), 
                     breaks=c(24, 1e4, 1e6, 1e8, 1e10, 1e12),
                     labels=c("<LoD", "1e04", "1e+06", 
                              "1e+08", "1e+10", "1e+12"),
                     trans="log10") +
  scale_x_continuous(breaks=sample.dates) +
  labs(x="Days postinfection",
       y="LASV PFU/mL", 
       col="Challenge (PFU)",
       fill="Challenge (PFU)")
ggsave("analysis/viremia-plaqueassay-longitudinal.png",
       units="in", width=6, height=4)

# hypothesis: viremia onset differs by challenge dose -- no difference
# qRT-PCR
viremia %>%
  filter(`GEq/mL` > 0) %>%
  group_by(Tattoo) %>%
  top_n(n=1, wt=DPI) %>%
  ungroup() %>%
  rename(Onset=DPI) %>%
  rstatix::kruskal_test(Onset ~ Dose)
# plaque assay
viremia %>%
  filter(`PFU/mL` > 0) %>%
  group_by(Tattoo) %>%
  top_n(n=1, wt=DPI) %>%
  ungroup() %>%
  rename(Onset=DPI) %>%
  rstatix::kruskal_test(Onset ~ Dose)

# what about by outcome? still no. 
# qRT-PCR
viremia %>%
  filter(`GEq/mL` > 0) %>%
  group_by(Tattoo) %>%
  top_n(n=1, wt=DPI) %>%
  ungroup() %>%
  rename(Onset=DPI) %>%
  rstatix::wilcox_test(Onset ~ Outcome)
# plaque assay
viremia %>%
  filter(`PFU/mL` > 0) %>%
  group_by(Tattoo) %>%
  top_n(n=1, wt=DPI) %>%
  ungroup() %>%
  rename(Onset=DPI) %>%
  rstatix::kruskal_test(Onset ~ Outcome)

# hypothesis: viral load at onset differs by challenge dose -- no difference
# qRT-PCR
viremia %>%
  filter(`GEq/mL` > 0) %>%
  group_by(Tattoo) %>%
  top_n(n=1, wt=DPI) %>%
  ungroup() %>%
  rstatix::kruskal_test(`GEq/mL` ~ Dose)
# plaque assay
viremia %>%
  filter(`PFU/mL` > 0) %>%
  group_by(Tattoo) %>%
  top_n(n=1, wt=DPI) %>%
  ungroup() %>%
  rstatix::kruskal_test(`PFU/mL` ~ Dose)

# what about by outcome? Yes, but n < 3 survived with viremia, so eh...
# qRT-PCR
viremia %>%
  filter(`GEq/mL` > 0) %>%
  group_by(Tattoo) %>%
  top_n(n=1, wt=DPI) %>%
  ungroup() %>%
  rstatix::wilcox_test(`GEq/mL` ~ Outcome)
# plaque assay
viremia %>%
  filter(`PFU/mL` > 0) %>%
  group_by(Tattoo) %>%
  top_n(n=1, wt=DPI) %>%
  ungroup() %>%
  rstatix::wilcox_test(`PFU/mL` ~ Outcome)

# hypothesis: peak viral load differs by challenge dose -- absolutely
# qRT-PCR
p1 <- viremia %>%
      group_by(Tattoo) %>%
      slice_max(order_by=`GEq/mL`, n=1, with_ties=FALSE) %>%
      ungroup() %>%
      ggplot(aes(Dose, `GEq/mL`+1e2)) +
      geom_boxplot(aes(fill=Dose), outliers=FALSE, alpha=0.5) +
      geom_jitter(height=0, width=0.2) +
      scale_fill_manual(values=cols.dose) +
      ggpubr::stat_compare_means(method="kruskal.test", 
                                 label.y=17,
                                 label.x.npc=0.1) +
      ggpubr::stat_pwc(label="p.adj.signif", step.increase=0.15, hide.ns=TRUE) + 
      scale_y_continuous(limits=c(20, 1e18), 
                         breaks=c(1e2, 1e4, 1e6, 1e8, 1e10, 1e12),
                         labels=c("<LoD", "1e+04", "1e+06", 
                                  "1e+08", "1e+10", "1e+12"),
                         trans="log10") +
      facet_wrap(~ "qRT-PCR") +
      labs(x="Challenge (PFU)",
           y="Peak LASV GEq/mL") +
      theme(legend.position="none")
# plaque assay
p2 <- viremia %>%
      group_by(Tattoo) %>%
      slice_max(order_by=`PFU/mL`, n=1, with_ties=FALSE) %>%
      ungroup() %>%
      ggplot(aes(Dose, `PFU/mL`+1e2)) +
      geom_boxplot(aes(fill=Dose), outliers=FALSE, alpha=0.5) +
      geom_jitter(height=0, width=0.2) +
      scale_fill_manual(values=cols.dose) +
      ggpubr::stat_compare_means(method="kruskal.test", 
                                 label.y=17,
                                 label.x.npc=0.1) +
      ggpubr::stat_pwc(label="p.adj.signif", step.increase=0.3, hide.ns=TRUE) + 
      scale_y_continuous(limits=c(20, 1e18), 
                         breaks=c(24, 1e4, 1e6, 1e8, 1e10, 1e12),
                         labels=c("<LoD", "1e+04", "1e+06", 
                                  "1e+08", "1e+10", "1e+12"),
                         trans="log10") +
      facet_wrap(~ "Plaque assay") +
      labs(x="Challenge (PFU)",
           y="Peak LASV PFU/mL") +
      theme(legend.position="none")
cowplot::plot_grid(p1, p2, nrow=1)
ggsave("analysis/viremia-peak-dose.png",
       units="in", width=6, height=3)

# what about by outcome?
# qRT-PCR
p1 <- viremia %>%
      group_by(Tattoo) %>%
      slice_max(order_by=`GEq/mL`, n=1, with_ties=FALSE) %>%
      ungroup() %>%
      ggplot(aes(Outcome, `GEq/mL`+1e2)) +
      geom_boxplot(aes(fill=Outcome), outliers=FALSE, alpha=0.5) +
      geom_jitter(height=0, width=0.2) +
      scale_fill_manual(values=cols.outcome) +
      ggpubr::stat_compare_means(comparisons=list(c("Survived", "Succumbed")), 
                                 label="p.signif") +
      scale_y_continuous(limits=c(20, 1e13), 
                         breaks=c(1e2, 1e4, 1e6, 1e8, 1e10, 1e12),
                         labels=c("<LoD", "1e04", "1e+06", 
                                  "1e+08", "1e+10", "1e+12"),
                         trans="log10") +
      facet_wrap(~ "qRT-PCR") +
      labs(x="Outcome",
           y="Peak LASV GEq/mL") +
      theme(legend.position="none")
# plaque assay
p2 <- viremia %>%
      group_by(Tattoo) %>%
      slice_max(order_by=`PFU/mL`, n=1, with_ties=FALSE) %>%
      ungroup() %>%
      ggplot(aes(Outcome, `PFU/mL`+24)) +
      geom_boxplot(aes(fill=Outcome), outliers=FALSE, alpha=0.5) +
      geom_jitter(height=0, width=0.2) +
      scale_fill_manual(values=cols.outcome) +
      ggpubr::stat_compare_means(comparisons=list(c("Survived", "Succumbed")), 
                                 label="p.signif") +
      scale_y_continuous(limits=c(20, 1e13), 
                         breaks=c(24, 1e4, 1e6, 1e8, 1e10, 1e12),
                         labels=c("<LoD", "1e04", "1e+06", 
                                  "1e+08", "1e+10", "1e+12"),
                         trans="log10") +
      facet_wrap(~ "Plaque assay") +
      labs(x="Outcome",
           y="Peak LASV PFU/mL") +
      theme(legend.position="none")
cowplot::plot_grid(p1, p2, nrow=1)
ggsave("analysis/viremia-peak-outcome.png",
       units="in", width=6, height=3)

# clean up
rm(viremia, p1, p2)

## tissue virology -------------------------------------------------------------
# qRT-PCR only -- not all tissues plaque-assayed
df$`virology-tissues` %>%
  filter(Tissue != "Aqueous humor") %>%
  left_join(df$animals, by="Tattoo") %>%
  ggplot(aes(Tissue, `GEq/g`+1e2)) +
  geom_boxplot(aes(fill=Outcome, col=Outcome), alpha=0.5, na.rm=TRUE) +
  scale_fill_manual(values=cols.outcome) +
  scale_color_manual(values=cols.outcome) +
  ggpubr::stat_compare_means(aes(group=Outcome), na.rm=TRUE,
                             label="p.signif", method="wilcox.test") +
  scale_y_continuous(limits=c(20, NA), 
                     breaks=c(1e2, 1e4, 1e6, 1e8, 1e10, 1e12),
                     labels=c("<LoD", "1e+04", "1e+06", 
                              "1e+08", "1e+10", "1e+12"),
                     trans="log10") +
  labs(x=element_blank(),
       y="LASV GEq/g") +
  guides(fill=guide_legend(override.aes=list(pch=21)),
         shape=guide_legend(override.aes=list(fill="black"))) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="bottom")
ggsave("analysis/viremia-qpcr-tissues.png",
       units="in", width=10, height=4)

## clinical chemistry ----------------------------------------------------------
chem <- df$chemistry %>%
        # no change in uric acid, so remove that
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
                           na.rm=TRUE) +
                scale_color_manual(values=cols.dose) +
                scale_x_continuous(breaks=sample.dates) +
                labs(x=element_blank(),
                     y=df$units[df$units$Analyte==i, "Units"], 
                     title=i) +
                theme(legend.position="none",
                      plot.title=element_text(size=12))
            })
names(spaghetti) <- levels(chem$Analyte)
# pull out the legend
legs <- spaghetti$`Glucose (GLU)` +
        theme(legend.position="bottom") +
        labs(col="Challenge (PFU)",
             shape="NHP",
             linetype="NHP") +
        guides(fill=guide_legend(override.aes=list(pch=21)),
               shape=guide_legend(override.aes=list(fill="black")))
legs <- cowplot::get_plot_component(legs, "guide-box-bottom")
xlabel <- spaghetti$`Glucose (GLU)` +
          labs(x="Days postinfection")
xlabel <- cowplot::get_plot_component(xlabel, "xlab-b")
# plot all chemistry together
x <- cowplot::plot_grid(plotlist=spaghetti, labels="AUTO", ncol=3)
cowplot::plot_grid(x, xlabel, legs, ncol=1, rel_heights=c(17, 0.75, 1))
ggsave("analysis/clinicalchemistry-spaghetti.png",
       units="in", width=12, height=9)

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
              select(Outcome, Comparison, Analyte, Concentration) %>%
              # add in baseline and onset
              rbind(x) 
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
                  labs(x=element_blank(),
                       y=df$units[df$units$Analyte==i, "Units"],
                       title=str_extract(i, "(?<=\\()[A-Z]+(?=\\))")) +
                  theme(legend.position="none",
                        plot.title=element_text(size=12))
              })
names(comparison) <- levels(chem$Analyte)
# pull out the legend
legs <- comparison$`Glucose (GLU)` + theme(legend.position="right")
legs <- cowplot::get_plot_component(legs, "guide-box-right")
# plot all chemistry together
x <- cowplot::plot_grid(plotlist=comparison, labels="AUTO", ncol=4)
cowplot::plot_grid(x, legs, ncol=2, rel_widths=c(5, 1))
ggsave("analysis/clinicalchemistry-peak-outcome.png",
       units="in", width=8, height=6)

# clean up
rm(chem, x, bl, legs, xlabel, spaghetti, comparison, y)
  
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
                geom_point(aes(col=Dose)) +
                scale_color_manual(values=cols.dose) +
                scale_x_continuous(breaks=sample.dates) +
                labs(x=element_blank(),
                     y=df$units[df$units$Analyte==i, "Units"],
                     title=i) +
                theme(legend.position="none",
                      plot.title=element_text(size=12))
            })
names(spaghetti) <- levels(hema$Analyte)
# pull out the legend
spaghetti$legs <- spaghetti$`Leukocytes` +
                  theme(legend.position="right") +
                  labs(col="Challenge (PFU)") 
spaghetti$legs <- cowplot::get_plot_component(spaghetti$legs, "guide-box-right")
xlabel <- spaghetti$Leukocytes +
          labs(x="Days postinfection")
xlabel <- cowplot::get_plot_component(xlabel, "xlab-b")
# plot all chemistry together
x <- cowplot::plot_grid(plotlist=spaghetti, 
                        labels=LETTERS[1:length(spaghetti)-1], ncol=3)
cowplot::plot_grid(x, xlabel, ncol=1, rel_heights=c(15, 1))
ggsave("analysis/hematology-spaghetti.png",
       units="in", width=12, height=9)

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
                  labs(x=element_blank(),
                       y=df$units[df$units$Analyte==i, "Units"],
                       title=i) +
                  theme(plot.title=element_text(size=12),
                        legend.position="none")
              })
names(comparison) <- levels(hema$Analyte)
# pull out the legend
comparison$legs <- comparison$Leukocytes + theme(legend.position="right")
comparison$legs <- cowplot::get_plot_component(comparison$legs, 
                                               "guide-box-right")
# plot all hematology together
cowplot::plot_grid(plotlist=comparison, 
                   labels=LETTERS[1:length(comparison)-1], 
                   ncol=3)
ggsave("analysis/hematology-peak-outcome.png",
       units="in", width=6, height=6)

# clean up
rm(hema, bl, xlabel, spaghetti, comparison, y)

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

# plot each NHP
prnt %>%
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
  facet_wrap(~Tattoo) +
  labs(y="Neutralization (%)") +
  guides(fill=guide_legend(override.aes=list(pch=21)),
         shape=guide_legend(override.aes=list(fill="black"))) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10))
ggsave("analysis/prnts-curves.png",
       units="in", width=7.5, height=6)

# clean up
rm(prnt, vc)

## fin -------------------------------------------------------------------------
sessionInfo()
