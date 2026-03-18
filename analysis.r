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
cols.dose <- c("27.5"="black",
               "10"="grey35",
               "1"="grey53",
               "0.1"="grey80")
cols.group <- c(Survived="darkgrey",
                Succumbed="#e41a1c",
                "1-C"="black")
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
                     Dose=factor(Dose, levels=names(cols.dose)),
                     NHP=factor(NHP, levels=NHP),
                     NHP=relevel(NHP, "1-C"))

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
            linewidth=1, lineend="round") +
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

# percent succumbed per dose
df$animals %>%
  reshape2::dcast(Dose ~ Outcome, fun.aggregate=length, value.var="NHP") %>%
  mutate(PercentLethality=Succumbed/(Survived+Succumbed),
         PercentLethality=round(100*PercentLethality, digits=1))

# mean, median, and range of presented doses
df$animals %>%
  group_by(Dose) %>%
  summarise(Mean=mean(`Presented dose (DP)`))

# percent succumbed overall
df$animals %>%
  group_by(Outcome) %>%
  summarise(NHPs=n())

# mean TTD per dose
df$animals %>%
  filter(Outcome=="Succumbed") %>%
  group_by(Dose) %>%
  summarise(Mean=mean(Death))

# mean and 95% CI TTD in 27.5 and 10 PFU groups (n=1 for 1 PFU = meaningless)
df$animals %>%
  filter(Dose==27.5) %>%
  select(Death) %>%
  unlist() %>%
  t.test()
df$animals %>%
  filter(Dose==10) %>%
  select(Death) %>% 
  unlist() %>%
  t.test()

# any difference in TTD between 10 and 27.5 PFU groups?
df$animals %>%
  filter(Dose %in% c(10, 27.5)) %>%
  rstatix::wilcox_test(Death ~ Dose)

# scores
# extract the challenge time and calculate time since challenge for all obs
baseline <- df$scores %>%
            filter(DPI==0) %>%
            rename(Challenge=Datetime) %>%
            select(NHP, Challenge)
scores <- df$scores %>%
          # add in challenge dose
          left_join(df$animals, by="NHP") %>%
          # calculate time (days) since exact challenge time
          left_join(baseline, by="NHP") %>%
          mutate(DPI=as.integer(Datetime-Challenge),
                 DPI=DPI/(60*60*24)) %>%
          arrange(Dose)
scores %>%
  ggplot(aes(DPI, Score)) +
  geom_hline(yintercept=9, col="lightgrey", linetype=2) +
  geom_line(aes(col=Dose, group=NHP), linewidth=0.5, lineend="round") +
  scale_color_manual(values=cols.dose) +
  scale_x_continuous(breaks=breaks.dates) +
  labs(x="Days postinfection",
       y="Clinical score",
       col="PFU",
       title="Clinical scores") 
ggsave("analysis/score.png",
       units="in", width=4, height=2.5)

# who scored at least once?
scores %>%
  filter(Score > 0) %>%
  select(NHP, Dose, Outcome) %>%
  distinct()

# score starts per group (10 and 27.5 PFU)
scorestart <- scores %>%
              filter(Score > 0,
                     Dose %in% c(10, 27.5)) %>%
              group_by(NHP) %>%
              top_n(n=-1, wt=DPI) %>%
              ungroup()
# mean, median, and range
scorestart %>%
  group_by(Dose) %>%
  summarise(Mean=mean(DPI),
            Median=median(DPI),
            Min=min(DPI),
            Max=max(DPI))
# get 95% CI for these
scorestart %>%
  filter(Dose==27.5) %>%
  select(DPI) %>%
  unlist() %>%
  t.test()
scorestart %>%
  filter(Dose==10) %>%
  select(DPI) %>%
  unlist() %>%
  t.test()
# any difference?
scorestart %>%
  rstatix::wilcox_test(DPI ~ Dose)
# clean up
rm(scores, scorestart)

# decreased appetite
appetite <- df$observations %>%
            select(NHP, DPI, `Biscuits eaten`) %>%
            left_join(df$animals, by="NHP") %>%
             # remove D0 since we don't count biscuits that day
            filter(DPI > 0)
# for baseline, take D1-D4 average
baseline <- appetite %>%
            filter(DPI %in% 1:4) %>%
            group_by(NHP) %>%
            summarise(Baseline=mean(`Biscuits eaten`),
                      .groups="drop")
# calculate appetite relative to baseline
appetite <- appetite %>%
            left_join(baseline, by="NHP") %>%
            mutate(Appetite=`Biscuits eaten`/Baseline) 
# max appetite capped at 1 for visualization purposes
appetite$Appetite[appetite$Appetite > 1] <- 1
# plot it
appetite %>%
  # order NHPs by day of death
  arrange(Death) %>%
  mutate(NHP=factor(NHP, levels=unique(NHP)),
         Dose=factor(Dose, labels=paste(levels(appetite$Dose), "PFU"))) %>%
  # plot it
  ggplot(aes(DPI, NHP)) +
  geom_tile(aes(alpha=Appetite, fill=Outcome), col="white") +
  scale_fill_manual(values=cols.outcome) +
  scale_alpha(range=c(1, 0)) +
  scale_x_continuous(breaks=breaks.dates) +
  facet_wrap(~Dose, ncol=1, scales="free_y") +
  labs(x="Days postinfection",
       y=NULL,
       fill=NULL,
       title="Appetite") +
  guides(alpha="none") +
  theme(legend.position="bottom",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave("analysis/appetite.png",
       units="in", width=2.5, height=5)

# define the first day of decrease without subsequent recovery in succumbing 
# NHPs (1 PFU succumbed = 10 DPI)
# doing this manually since I'm not sure the easiest way to do it with code
appetite <- data.frame(NHP=c("27.5-F", "27.5-E", "27.5-D", 
                                "27.5-C", "27.5-B", "27.5-A",
                                "10-F", "10-E", "10-D", 
                                "10-C", "10-B", "10-A"),
                       DPI=c(7, 8, 8, 
                             8, 8, 7,
                             9, 8, 8, 
                             7, 8, 7)) %>%
            left_join(select(appetite, NHP, Dose, Outcome, Death),
                      by="NHP") %>%
            distinct() %>%
            mutate(Duration=Death-DPI+1) #need the inclusive number of days
# mean/median/range/95% CI per group relative to D0
appetite %>%
  group_by(Dose) %>%
  summarise(Median=median(DPI),
            Mean=mean(DPI),
            Min=min(DPI),
            Max=max(DPI))
appetite %>%
  filter(Dose==27.5) %>%
  select(DPI) %>%
  unlist() %>%
  t.test()
appetite %>%
  filter(Dose==10) %>%
  select(DPI) %>%
  unlist() %>%
  t.test()
appetite %>%
  rstatix::wilcox_test(DPI ~ Dose)
# mean/median/range/95% CI per group relative to terminal
appetite %>%
  group_by(Dose) %>%
  summarise(Median=median(Duration),
            Mean=mean(Duration),
            Min=min(Duration),
            Max=max(Duration))
appetite %>%
  filter(Dose==27.5) %>%
  select(Duration) %>%
  unlist() %>%
  t.test()
appetite %>%
  filter(Dose==10) %>%
  select(Duration) %>%
  unlist() %>%
  t.test()
appetite %>%
  rstatix::wilcox_test(Duration ~ Dose)

# clean up
rm(appetite, baseline)

## viremia ---------------------------------------------------------------------
# add challenge dose
viremia <- df$virology.blood %>%
           left_join(df$animals, by="NHP") %>%
           select(NHP, DPI, Outcome, `GEq/mL`, `PFU/mL`, Dose) %>%
           arrange(Dose)

# when was virus first detectable?
viremia %>%
  reshape2::melt(id.vars=c("NHP", "DPI", "Dose", "Outcome")) %>%
  filter(value > 0) %>%
  group_by(NHP, variable) %>%
  top_n(n=-1, wt=DPI) %>%
  reshape2::dcast(NHP + Dose + Outcome ~ variable, value.var="DPI") %>%
  arrange(Dose, Outcome)

# qRT-PCR
viremia %>%
  ggplot(aes(DPI, `GEq/mL`+1e2)) +
  geom_line(aes(col=Dose, group=NHP), linewidth=0.5, lineend="round") +
  scale_color_manual(values=cols.dose) +
  scale_y_continuous(limits=c(20, 1e13), 
                     breaks=c(1e2, 1e4, 1e8, 1e12),
                     labels=c("<LoD", "1e04", "1e+08", "1e+12"),
                     trans="log10") +
  scale_x_continuous(breaks=breaks.dates) +
  labs(x="Days postinfection",
       y="N gene eq/mL",
       col="PFU",
       title="LASV RNA") +
  theme(legend.position="none")
ggsave("analysis/viremia-qpcr-longitudinal.png",
       units="in", width=3.25, height=2)

# plaque assay
viremia %>%
  ggplot(aes(DPI, `PFU/mL`+24)) +
  geom_line(aes(group=NHP, col=Dose), linewidth=0.5, lineend="round") +
  scale_color_manual(values=cols.dose) +
  scale_y_continuous(limits=c(20, 1e13), 
                     breaks=c(24, 1e4, 1e8, 1e12),
                     labels=c("<LoD", "1e04",  "1e+08", "1e+12"),
                     trans="log10") +
  scale_x_continuous(breaks=breaks.dates) +
  labs(x="Days postinfection",
       y="LASV PFU/mL", 
       col="PFU",
       title="LASV titer") +
  theme(legend.position=c(0.8, 0.6))
ggsave("analysis/viremia-plaqueassay-longitudinal.png",
       units="in", width=3.25, height=2)

# is there a difference in viral load at onset? 
# Only considering 27.5 and 10 PFU groups because of numbers
viremia %>%
  reshape2::melt(id.vars=c("NHP", "Outcome", "Dose", "DPI"),
                 variable.name="Assay") %>%
  filter(value > 0,
         Dose %in% c(27.5, 10)) %>%
  group_by(NHP, Assay) %>%
  top_n(n=-1, wt=DPI) %>%
  ungroup() %>%
  group_by(Assay) %>%
  rstatix::wilcox_test(value ~ Dose)

# is there a difference in peak viral load?
# using all groups here since 0 is a valid value
viro.peaks <- viremia %>%
              select(-DPI) %>%
              reshape2::melt(id.vars=c("NHP", "Outcome", "Dose"),
                             variable.name="Assay") %>%
              group_by(NHP, Assay) %>%
              slice_max(value, n=1, with_ties=FALSE) %>%
              ungroup() 
# kw test first -- yes to both assays
viro.peaks %>%
  group_by(Assay) %>%
  rstatix::kruskal_test(value ~ Dose)
# dunn post-hoc
pval <- viro.peaks %>%
        group_by(Assay) %>%
        rstatix::dunn_test(value ~ Dose) %>%
        rstatix::add_xy_position(y.trans=log10, scales="free_y")
# plot these
viro.peaks %>%
  filter(Assay=="GEq/mL") %>%
  ggplot(aes(Dose, value+1e2)) +
  geom_boxplot(fill="grey", outliers=FALSE, alpha=0.5) +
  geom_jitter(aes(col=NHP, size=NHP), height=0, width=0.2) +
  ggpubr::stat_pvalue_manual(data=filter(pval, Assay=="GEq/mL"), hide.ns=TRUE,
                             label="p.adj.signif", step.increase=0.12,
                             bracket.nudge.y=1) +
  scale_color_manual(values=c("1-A"="#377eb8", "1-C"="#4daf4a"), 
                     na.value="black") +
  scale_size_manual(values=c("1-A"=2, "1-C"=2), na.value=1) +
  scale_y_continuous(limits=c(NA, 1e16), 
                     breaks=c(1e2, 1e4, 1e8, 1e12),
                     labels=c("<LoD", "1e+04", "1e+08", "1e+12"),
                     trans="log10") +
  labs(x="Challenge dose (PFU)",
       y="N gene eq/mL",
       col="AGM",
       title="Max LASV RNA") +
  theme(legend.position="none")
ggsave("analysis/viremia-qpcr-peak-dose.png",
       units="in", width=3, height=2)
viro.peaks %>%
  filter(Assay=="PFU/mL") %>%
  ggplot(aes(Dose, value+24)) +
  geom_boxplot(fill="grey", outliers=FALSE, alpha=0.5) +
  geom_jitter(aes(col=NHP, size=NHP), height=0, width=0.2) +
  ggpubr::stat_pvalue_manual(data=filter(pval, Assay=="PFU/mL"), hide.ns=TRUE,
                             label="p.adj.signif", step.increase=0.2, 
                             bracket.nudge.y=1) +
  scale_color_manual(values=c("1-A"="#377eb8", "1-C"="#4daf4a"), 
                     na.value="black") +
  scale_size_manual(values=c("1-A"=2, "1-C"=2), na.value=1) +
  scale_y_continuous(limits=c(NA, 1e16), 
                     breaks=c(24, 1e4, 1e8, 1e12),
                     labels=c("<LoD", "1e+04", "1e+08", "1e+12"),
                     trans="log10") +
  labs(x="Challenge dose (PFU)",
       y="LASV PFU/mL",
       col="AGM",
       size="AGM",
       title="Max titer") +
  guides(color=guide_legend(override.aes=list(size=3)))
ggsave("analysis/viremia-plaqueassay-peak-dose.png",
       units="in", width=3.5, height=2)

# now test peak again but by outcome
pval <- viro.peaks %>%
        group_by(Assay) %>%
        rstatix::wilcox_test(value ~ Outcome) %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(y.trans=log10, scales="free_y")
viro.peaks %>%
  filter(Assay=="GEq/mL") %>%
  ggplot(aes(Outcome, value+1e2)) +
  geom_boxplot(fill="grey", outliers=FALSE, alpha=0.5) +
  geom_jitter(aes(col=NHP, size=NHP), height=0, width=0.2) +
  ggpubr::stat_pvalue_manual(data=filter(pval, Assay=="GEq/mL"), 
                             bracket.nudge.y=1) +
  scale_color_manual(values=c("1-A"="#377eb8", "1-C"="#4daf4a"),
                     na.value="black") +
  scale_size_manual(values=c("1-A"=2, "1-C"=2), na.value=1) +
  scale_y_continuous(limits=c(NA, 1e16), 
                     breaks=c(1e2, 1e4, 1e8, 1e12),
                     labels=c("<LoD", "1e+04", "1e+08", "1e+12"),
                     trans="log10") +
  labs(x="Outcome",
       y="N gene eq/mL",
       title="Max LASV RNA") +
  theme(legend.position="none")
ggsave("analysis/viremia-qpcr-peak-outcome.png",
       units="in", width=3, height=2.5)
viro.peaks %>%
  filter(Assay=="PFU/mL") %>%
  ggplot(aes(Outcome, value+24)) +
  geom_boxplot(fill="grey", outliers=FALSE, alpha=0.5) +
  geom_jitter(aes(col=NHP, size=NHP), height=0, width=0.2) +
  ggpubr::stat_pvalue_manual(data=filter(pval, Assay=="PFU/mL"), 
                             bracket.nudge.y=1) + 
  scale_color_manual(values=c("1-A"="#377eb8", "1-C"="#4daf4a"),
                     na.value="black") +
  scale_size_manual(values=c("1-A"=2, "1-C"=2),
                    na.value=1) +
  scale_y_continuous(limits=c(NA, 1e16), 
                     breaks=c(1e2, 1e4, 1e8, 1e12),
                     labels=c("<LoD", "1e+04", "1e+08", "1e+12"),
                     trans="log10") +
  labs(x="Outcome",
       y="LASV PFU/mL",
       col="AGM",
       size="AGM",
       title="Max titer")  +
  guides(color=guide_legend(override.aes=list(size=3)))
ggsave("analysis/viremia-plaqueassay-peak-outcome.png",
       units="in", width=3.5, height=2.5)

# check: what is min/max each outcome?
viro.peaks %>%
  group_by(Assay, Outcome) %>%
  summarise(Max=max(value),
            Min=min(value))

# clean up
rm(viremia, viro.peaks, pval)

# tissue virology
df$virology.tissues %>%
  filter(Tissue != "Aqueous humor") %>%
  left_join(df$animals, by="NHP") %>%
  ggplot(aes(Tissue, `GEq/g`+1e4)) +
  geom_boxplot(aes(fill=Outcome, col=Outcome), alpha=0.5, na.rm=TRUE,
               outlier.size=0.5) +
  ggpubr::stat_compare_means(aes(group=interaction(Outcome, Tissue)), 
                             method="wilcox.test", label="p.signif",
                             angle=90, vjust=1, hjust=0) +
  scale_fill_manual(values=cols.outcome) +
  scale_color_manual(values=cols.outcome) +
  scale_y_continuous(breaks=c(1e4, 1e8, 1e12),
                     labels=c("<LoD", "1e+08", "1e+12"),
                     trans="log10", 
                     expand=expansion(mult=c(0.1, 0.2))) +
  labs(x=NULL,
       y="N gene eq/g",
       title="LASV RNA in tissues") +
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position="bottom")
ggsave("analysis/viremia-qpcr-tissues.png",
       units="in", width=6.5, height=4)

## clinical pathology ----------------------------------------------------------
# format clinical chemistry
chem <- df$chemistry %>%
        left_join(select(df$animals, NHP, Outcome, Dose, Death), 
                  by="NHP") %>%
        reshape2::melt(id.vars=c("NHP", "DPI", "Outcome", "Dose", "Death"),
                       variable.name="Analyte",
                       value.name="Concentration") %>%
        mutate(NHP=factor(NHP, levels=rev(levels(df$animals$NHP))),
               Type="chem")
# format hematology
hema <- df$hematology %>%
        select(-starts_with("Mean"), -starts_with("Hem")) %>%
        left_join(select(df$animals, NHP, Outcome, Dose, Death),
                  by="NHP") %>%
        reshape2::melt(id.vars=c("NHP", "DPI", "Outcome", "Dose", "Death"),
                       variable.name="Analyte",
                       value.name="Concentration") %>%
        mutate(Type="hema")
# combine and format
clinpath <- rbind(chem, hema)
rm(chem, hema)
# add group for plotting
clinpath$Group <- factor(clinpath$Outcome, levels=names(cols.group))
clinpath$Group[clinpath$NHP=="1-C"] <- "1-C"
# date a "datetime" variable for terminals
clinpath$Datetime <- clinpath$DPI
clinpath$Datetime[(clinpath$DPI != 35) & (clinpath$DPI==clinpath$Death)] <- "Terminal"
clinpath$Datetime <- factor(clinpath$Datetime, levels=unique(clinpath$Datetime))

# when do we see changes per analyte in succumbing animals?
# not all animals have D14, so using KW instead of Friedman
analytes <- clinpath %>%
            filter(Outcome=="Succumbed") %>%
            mutate(Datetime=droplevels(Datetime)) %>%
            group_by(Analyte) %>%
            rstatix::kruskal_test(Concentration ~ Datetime) %>%
            rstatix::p_format(new.col=TRUE, digits=1) %>%
            filter(p < 0.05) 
# run Dunn post-hoc and filter for comparisons relative to D0/baseline
clinpath %>%
  filter(Outcome=="Succumbed",
         Analyte %in% analytes$Analyte) %>%
  mutate(Datetime=droplevels(Datetime)) %>%
  group_by(Analyte) %>%
  rstatix::dunn_test(Concentration ~ Datetime) %>%
  rstatix::p_format(new.col=TRUE, digits=1) %>%
  filter(group1==0,
         p.adj < 0.05)

# what about aviremic surviving animals?
# we could use friedman, but need to compare to succumbed AGMs -- 
# using KW to be consistent
analytes <- clinpath %>%
            filter(Outcome=="Survived",
                   NHP != "D8993") %>%
            mutate(Datetime=droplevels(Datetime)) %>%
            group_by(Analyte) %>%
            rstatix::kruskal_test(Concentration ~ Datetime) %>%
            rstatix::p_format(new.col=TRUE, digits=1) %>%
            filter(p < 0.05) 
# run Dunn post-hoc and filter for comparisons relative to D0/baseline
clinpath %>%
  filter(Outcome=="Survived",
         NHP != "D8993",
         Analyte %in% analytes$Analyte) %>%
  mutate(Datetime=droplevels(Datetime)) %>%
  group_by(Analyte) %>%
  rstatix::dunn_test(Concentration ~ Datetime) %>%
  rstatix::p_format(new.col=TRUE, digits=1) %>%
  filter(group1==0,
         p.adj < 0.05)
rm(analytes)

# plot all clinical chemistry spaghetti
analytes <- clinpath %>%
            filter(Type=="chem") %>%
            select(Analyte) %>%
            unlist() %>%
            droplevels() %>%
            levels()
chem <- analytes %>%
        lapply(function(i) {
        # plot it
        clinpath %>%
          filter(Analyte==i,
                 !is.na(Concentration)) %>%
          ggplot(aes(DPI, Concentration)) +
          geom_line(aes(group=NHP, 
                        col=Group),
                    linewidth=0.5, lineend="round") +
          scale_color_manual(values=cols.group) +
          scale_x_continuous(breaks=breaks.dates) +
          labs(x=NULL,
               y=df$units[df$units$Analyte==i, "Units"], 
               title=str_extract(i, "(?<=\\()[A-Z]+(?=\\))")) +
          theme(legend.position="none")
        })
names(chem) <- str_extract(analytes, "(?<=\\()[A-Z]+(?=\\))")
# pull out the legend
legs <- chem$GLU +
        theme(legend.position="bottom") +
        labs(col=NULL)
legs <- cowplot::get_plot_component(legs, "guide-box-bottom")
xlabel <- chem$GLU +
          labs(x="Days postinfection")
xlabel <- cowplot::get_plot_component(xlabel, "xlab-b")
# plot all chemistry together
x <- cowplot::plot_grid(plotlist=chem, labels="AUTO", ncol=3)
cowplot::plot_grid(x, xlabel, legs, ncol=1, rel_heights=c(20, 1, 1))
ggsave("analysis/clinpath-spaghetti-chemistry.png",
       units="in", width=6.5, height=7)
rm(legs, xlabel, x, analytes)

# hematology spaghetti
analytes <- clinpath %>%
            filter(Type=="hema") %>%
            select(Analyte) %>%
            unlist() %>%
            droplevels() %>%
            levels()
hema <- analytes %>%
        lapply(function(i) {
          clinpath %>%
            filter(Analyte==i) %>%
            ggplot(aes(DPI, Concentration)) +
            geom_line(aes(group=NHP, col=Group), 
                      linewidth=0.5, lineend="round") +
            scale_color_manual(values=cols.group) +
            scale_x_continuous(breaks=breaks.dates) +
            labs(x=NULL,
                 y=df$units[df$units$Analyte==i, "Units"],
                 title=i) +
            theme(legend.position="none")
        })
names(hema) <- analytes
# pull out the legend
hema$legs <- hema$Leukocytes +
             theme(legend.position="right") +
             labs(col=NULL) 
hema$legs <- cowplot::get_plot_component(hema$legs, "guide-box-right")
xlabel <- hema$Leukocytes +
          labs(x="Days postinfection")
xlabel <- cowplot::get_plot_component(xlabel, "xlab-b")
# plot all chemistry together
x <- cowplot::plot_grid(plotlist=hema, 
                        labels=LETTERS[1:length(hema)-1], ncol=3)
cowplot::plot_grid(x, xlabel, ncol=1, rel_heights=c(15, 1))
ggsave("analysis/clinpath-spaghetti-hematology.png",
       units="in", width=6.5, height=6)
rm(analytes, xlabel, x)

# plot selected spaghetti
cowplot::plot_grid(chem$ALB + labs(x="Days postinfection"), 
                   hema$Platelets + labs(x="Days postinfection"), 
                   hema$Neutrophils + labs(x="Days postinfection"),
                   chem$CRP + labs(x="Days postinfection"), 
                   chem$ALT + labs(x="Days postinfection"), 
                   hema$legs,
                   nrow=2, labels=LETTERS[1:5])
ggsave("analysis/clinpath-spaghetti-selected.png",
       units="in", width=6.5, height=4)
rm(hema, chem)

# calculate and compare baseline, max, and min
comparison <- clinpath %>%
              filter(DPI==0) %>%
              rename(Baseline=Concentration) %>%
              select(NHP, Analyte, Baseline)
comparison <- clinpath %>%
              left_join(comparison, by=c("NHP", "Analyte")) %>%
              mutate(Delta=Concentration-Baseline) %>%
              group_by(NHP, Outcome, Analyte, Group)
# get min and max
conc.min <- comparison %>%
            top_n(n=-1, wt=Delta) %>%
            mutate(Comparison="Min") %>%
            select(NHP, Outcome, Analyte, Concentration, 
                   Group, Type, Comparison)
conc.max <- comparison %>%
            top_n(n=1, wt=Delta) %>%
            mutate(Comparison="Max") %>%
            select(colnames(conc.min))
# get baseline
baseline <- clinpath %>%
            filter(DPI==0) %>%
            mutate(Comparison="D0") %>%
            select(colnames(conc.min))
# compile all
comparison <- baseline %>%
              rbind(conc.min) %>%
              rbind(conc.max)
rm(conc.min, conc.max, baseline)

# plot all chemistry
analytes <- comparison %>%
            filter(Type=="chem") %>%
            select(Analyte) %>%
            unlist() %>%
            droplevels() %>%
            levels()
# create a list of plots
chem <- analytes %>%
        lapply(function(i) {
          comparison %>%
            filter(Analyte==i) %>%
            ggplot(aes(Comparison, Concentration)) +
            geom_boxplot(aes(group=interaction(Outcome, Comparison)),
                         fill="lightgrey", col="black",
                         alpha=0.5, outliers=FALSE) +
            geom_point(aes(group=interaction(Outcome, Comparison),
                           col=Group, size=Group), 
                       position=position_jitterdodge(jitter.height=0)) +
            scale_color_manual(values=cols.group) +
            scale_size_manual(values=c("1-C"=1), na.value=0.5) +
            ggpubr::stat_compare_means(aes(group=interaction(Comparison, 
                                                             Outcome)), 
                                       method="wilcox.test",
                                       hide.ns=TRUE, 
                                       label="p.signif") +
            scale_y_continuous(expand=expansion(mult=c(0.1, 0.2))) +
            labs(x=NULL,
                 y=df$units[df$units$Analyte==i, "Units"],
                 title=str_extract(i, "(?<=\\()[A-Z]+(?=\\))")) +
            theme(legend.position="none")
        })
names(chem) <- str_extract(analytes, "(?<=\\()[A-Z]+(?=\\))")
# pull out the legend
legs <- chem$GLU + 
        labs(col=NULL) + 
        guides(col=guide_legend(override.aes=list(size=3)),
               size="none") +
        theme(legend.position="bottom")
legs <- cowplot::get_plot_component(legs, "guide-box-bottom")
# plot all chemistry together
x <- cowplot::plot_grid(plotlist=chem, labels="AUTO", ncol=3)
cowplot::plot_grid(x, legs, ncol=1, rel_heights=c(15, 1))
ggsave("analysis/clinpath-peak-chemistry.png",
       units="in", width=6.5, height=7)
rm(legs, x, analytes)

# plot all hematology
analytes <- comparison %>%
            filter(Type=="hema") %>%
            select(Analyte) %>%
            unlist() %>%
            droplevels() %>%
            levels()
# create a list of plots
hema <- analytes %>%
        lapply(function(i) {
          comparison %>%
            filter(Analyte==i) %>%
            ggplot(aes(Comparison, Concentration)) +
            geom_boxplot(aes(group=interaction(Outcome, Comparison)),
                         fill="lightgrey", col="black",
                         alpha=0.5, outliers=FALSE) +
            geom_point(aes(group=interaction(Outcome, Comparison),
                           col=Group, size=Group), 
                       position=position_jitterdodge(jitter.height=0)) +
            scale_color_manual(values=cols.group) +
            scale_size_manual(values=c("1-C"=1), na.value=0.5) +
            ggpubr::stat_compare_means(aes(group=interaction(Comparison, 
                                                             Outcome)), 
                                       method="wilcox.test",
                                       hide.ns=TRUE, 
                                       label="p.signif") +
            scale_y_continuous(expand=expansion(mult=c(0.1, 0.2))) +
            labs(x=NULL,
                 y=df$units[df$units$Analyte==i, "Units"],
                 title=i) +
            theme(legend.position="none")
        })
names(hema) <- analytes
# pull out the legend
hema$legs <- hema$Leukocytes +
             labs(col=NULL) + 
             guides(col=guide_legend(override.aes=list(size=3)),
                    size="none") +
             theme(legend.position="right") 
hema$legs <- cowplot::get_plot_component(hema$legs, "guide-box-right")
# plot all chemistry together
cowplot::plot_grid(plotlist=hema, labels=LETTERS[1:length(hema)-1], ncol=3)
ggsave("analysis/clinpath-peak-hematology.png",
       units="in", width=6.5, height=6)

# plot selected peaks
cowplot::plot_grid(chem$ALB, hema$Platelets, hema$Neutrophils, 
                   chem$CRP, chem$ALT, hema$legs,
                   nrow=2, labels=LETTERS[6:10])
ggsave("analysis/clinpath-peak-selected.png",
       units="in", width=6.5, height=4)

# clean up
rm(clinpath, chem, hema, comparison)

## ELISAs ----------------------------------------------------------------------
# over time -- very messy
df$ELISAs %>%
  reshape2::melt(measure.vars=c("IgM", "IgG"),
                 value.name="Titer") %>%
  left_join(df$animals, by="NHP") %>%
  ggplot(aes(DPI, Titer+100)) +
  geom_line(aes(group=NHP, col=Outcome)) +
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
          group_by(NHP, Antigen) %>%
          top_n(n=1, wt=DPI) %>%
          reshape2::melt(measure.vars=c("IgM", "IgG"),
                         variable.name="Immunoglobulin",
                         value.name="Titer") %>%
          left_join(df$animals, by="NHP")
elisas$Group <- factor(elisas$Outcome, levels=names(cols.group))
elisas$Group[elisas$NHP=="1-C"] <- "1-C"

# get p-value significance
psig <- elisas %>%
        mutate(Titer=log2(Titer+1)) %>%
        group_by(Immunoglobulin, Antigen) %>%
        rstatix::wilcox_test(Titer ~ Outcome) %>%
        rstatix::adjust_pvalue() %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(x="Antigen")
# add pseudocount to samples <LOD
elisas$Titer[elisas$Titer==0] <- 50
# plot it and add statistical comparisons
elisas %>%
  ggplot(aes(Antigen, Titer)) +
  geom_boxplot(aes(group=interaction(Antigen, Outcome)), 
               fill="lightgrey", col="black",
               alpha=0.5, outliers=FALSE) +
  geom_point(aes(group=interaction(Antigen, Outcome), size=Group, col=Group), 
             position=position_jitterdodge()) +
  scale_color_manual(values=cols.group) +
  scale_size_manual(values=c("1-C"=1), na.value=0.5) +
  scale_y_continuous(breaks=c(50, 200, 800, 3200, 12800),
                     labels=c("<LOD", 200, 800, 3200, 12800),
                     trans="log2",
                     expand = expansion(mult = c(0.05, 0.1))) +
  ggpubr::stat_pvalue_manual(psig, label="p.adj.signif") +
  facet_wrap(~Immunoglobulin) +
  labs(y="Endpoint titer",
       x="LASV antigen",
       col=NULL,
       title="ELISAs") +
  theme(legend.position="none")
ggsave("analysis/elisas-endpoint.png",
       units="in", width=2.5, height=2.5)

# clean up
rm(psig, elisas)

## PRNTs -----------------------------------------------------------------------
# get the average plaques per sample
prnt <- df$PRNTs %>%
        reshape2::melt(measure.vars=c("Plaques A", "Plaques B"),
                       value.name="Plaques") %>%
        group_by(NHP, DPI, Batch, `Dilution factor`) %>%
        summarise(Plaques=mean(Plaques),
                  .groups="drop")
# pull out VC and calculate the baseline for each batch
vc <- prnt %>%
      filter(NHP=="Control") %>%
      rename(VC=Plaques) %>%
      select(Batch, VC)
# calculate % reduction
prnt <- prnt %>%
        filter(NHP != "Control") %>%
        # add in VC counts and calculate % reduction
        left_join(vc, by="Batch") %>%
        mutate(Reduction=VC-Plaques,
               PercentReduction=100*Reduction/VC,
               # update DPI to factor
               DPI=factor(DPI, levels=c(0, 14, "Terminal", 35))) %>%
        replace_na(list(DPI="Terminal")) %>%
        # add NHP info
        left_join(df$animals, by="NHP")
# bound the % reduction to [0, 100]
prnt$PercentReduction[prnt$PercentReduction < 0] <- 0

# add group
prnt$Group <- factor(prnt$Outcome, levels=names(cols.group))
prnt$Group[prnt$NHP=="1-C"] <- "1-C"

# plot each NHP by dose
prnt %>%
  mutate(Dose=factor(Dose, labels=paste(levels(Dose), "PFU"))) %>%
  arrange(desc(DPI)) %>%
  ggplot(aes(`Dilution factor`, PercentReduction)) +
  geom_hline(yintercept=50, linetype=3) +
  geom_line(aes(group=interaction(NHP, DPI), col=DPI), linewidth=0.5) +
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
       units="in", width=6.5, height=3)

# plot neutralization for each at the lowest dilution and compare outcomes
prnt <- prnt %>%
        filter(`Dilution factor`==10,
               DPI %in% c("Terminal", "35")) %>%
        group_by(NHP) %>%
        slice_max(order_by=PercentReduction, n=1, with_ties=FALSE) %>%
        ungroup()
pval <- prnt %>%
        rstatix::wilcox_test(PercentReduction ~ Outcome) %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(x="Outcome")
pval
prnt %>%
  ggplot(aes(Outcome, PercentReduction)) +
  geom_hline(yintercept=50, linetype=3) +
  geom_boxplot(alpha=0.5, outliers=FALSE, col="black", fill="lightgrey") +
  geom_jitter(aes(col=Group, size=Group), height=0, width=0.2) +
  scale_color_manual(values=cols.group) +
  scale_size_manual(values=c("1-C"=1), na.value=0.5) +
  ylim(0, 100) +
  ggpubr::stat_pvalue_manual(data=pval, label="p.signif", bracket.nudge.y=5) +
  labs(y="Neutralization (%)",
       x="Outcome",
       col=NULL,
       title="LASV neutralization") +
  guides(fill="none",
         size="none",
         col=guide_legend(override.aes=list(size=3)))
ggsave("analysis/prnts-max.png",
       units="in", width=4, height=2.5)

# clean up
rm(prnt, vc, pval)

## fin -------------------------------------------------------------------------
sessionInfo()
