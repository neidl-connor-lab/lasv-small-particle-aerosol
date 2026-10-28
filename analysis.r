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

# figures
fig1 <- list()
fig3 <- list()
fig6 <- list()
sup1 <- list()
sup3 <- list()
sup4 <- list()
sup5 <- list()
sup6 <- list()
sup8 <- list()

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
fig1$A <- km %>%
          mutate(Dose="Overall") %>%
          rbind(km) %>%
          mutate(Dose=factor(Dose, 
                             levels=c(levels(df$animals$Dose), "Overall"))) %>%
          survfit2(Surv(Death, Outcome) ~ Dose, data=.) %>%
          ggsurvfit(theme=ggpubr::theme_pubr(), 
                    linewidth=1, lineend="round") +
          scale_color_manual(values=c(cols.dose, Overall="#e41a1c")) +
          # add p-value
          annotate("text", x=35, y=1.1, hjust=1, label=pval) +
          scale_y_continuous(limits=c(0, 1.2), 
                             breaks=c(0, 0.25, 0.5, 0.75, 1)) +
          scale_x_continuous(breaks=breaks.dates) +
          labs(x="Days postinfection",
               y="Probability", 
               col=NULL) +
          theme(legend.position="right")
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


# set up actual (presented) and target doses. no need to take log?
ld50 <- df$animals %>%
        mutate(Target=as.numeric(as.character(Dose)),
               Presented=`Presented dose (DP)`,
               Outcome=as.integer(Outcome)-1)

# using target dose
drc::drm(select(ld50, Outcome, Target),
         fct=drc::LL.2(), 
         type="binomial") %>%
  drc::ED(c(10, 50, 90), interval="delta", display=FALSE) %>%
  round(digits=1)

# using presented dose
drc::drm(select(ld50, Outcome, Presented),
         fct=drc::LL.2(), 
         type="binomial") %>%
  drc::ED(c(10, 50, 90), interval="delta", display=FALSE) %>%
  round(digits=1)
rm(ld50)

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
                 DPI=DPI/(60*60*24),
                 Score=if_else(Score > 9, 9, Score)) %>%
          arrange(Dose)
fig1$B <- scores %>%
          ggplot(aes(DPI, Score)) +
          geom_hline(yintercept=9, col="lightgrey", linetype=2) +
          geom_line(aes(col=Dose, group=NHP), linewidth=0.5, lineend="round") +
          scale_color_manual(NULL, values=cols.dose) +
          scale_y_continuous("Clinical score", breaks=c(0, 3, 6, 9),
                             labels=c("0", "3", "6", "≥9")) +
          scale_x_continuous("Days postinfection", breaks=breaks.dates) +
          guides(col=guide_legend(override.aes=list(linewidth=1))) +
          theme(legend.position=c(0.9, 0.7))


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
fig1$C <- appetite %>%
          # order NHPs by day of death
          arrange(Death) %>%
          mutate(NHP=factor(NHP, levels=unique(NHP)),
                 Dose=factor(Dose, 
                             labels=paste(levels(appetite$Dose), "PFU")),
                 DPI=DPI-1) %>%
          # plot it
          ggplot(aes(DPI, NHP)) +
          geom_tile(aes(alpha=Appetite, fill=Outcome), col="black") +
          scale_fill_manual(NULL, values=cols.outcome) +
          scale_alpha(NULL, range=c(1, 0)) +
          scale_x_continuous("Days postinfection", breaks=breaks.dates,
                             expand=expansion(c(0, 0))) +
          scale_y_discrete("AGM", expand=expansion(c(0, 0))) +
          facet_wrap(~Dose, ncol=1, scales="free_y") +
          guides(alpha="none") +
          theme(legend.position="bottom",
                axis.text.y=element_blank(),
                axis.ticks.y=element_line(color="transparent"))

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
fig1$D <- viremia %>%
          mutate(`GEq/mL`=if_else(`GEq/mL`==0, 99, `GEq/mL`)) %>%
          ggplot(aes(DPI, `GEq/mL`)) +
          geom_line(aes(col=Dose, group=NHP), linewidth=0.5, lineend="round") +
          scale_color_manual(NULL, values=cols.dose) +
          scale_y_continuous("N gene eq/mL",
                             limits=c(99, NA), 
                             breaks=c(99, 1e4, 1e8, 1e12),
                             labels=c("<LoD", "1e04", "1e+08", "1e+12"),
                             trans="log10") +
          scale_x_continuous("Days postinfection", breaks=breaks.dates) +
          guides(color=guide_legend(override.aes=list(linewidth=1))) +
          theme(legend.position=c(0.7, 0.6),
                legend.background=element_blank())

# plaque assay
fig1$E <- viremia %>%
          mutate(`PFU/mL`=if_else(`PFU/mL` < 25, 24, `PFU/mL`)) %>%
          ggplot(aes(DPI, `PFU/mL`)) +
          geom_line(aes(group=NHP, col=Dose), linewidth=0.5, lineend="round") +
          scale_color_manual(NULL, values=cols.dose) +
          scale_y_continuous("LASV PFU/mL",
                             limits=c(24, 1e7), 
                             breaks=c(24, 1e3, 1e5, 1e7),
                             labels=c("<LoD", "1e03",  "1e+05", "1e+07"),
                             trans="log10") +
          scale_x_continuous("Days postinfection", breaks=breaks.dates) +
          guides(color=guide_legend(override.aes=list(linewidth=1))) +
          theme(legend.position=c(0.9, 0.6))

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
# plot with dunn post-hoc
fig1$F <- viro.peaks %>%
          filter(Assay=="GEq/mL") %>%
          mutate(value=if_else(value==0, 99, value)) %>%
          ggplot(aes(Dose, value)) +
          geom_boxplot(fill="grey", outliers=FALSE, alpha=0.5) +
          geom_jitter(aes(col=NHP, size=NHP), height=0, width=0.2) +
          ggpubr::stat_pwc(method="dunn.test", 
                           hide.ns=TRUE, 
                           label="p.adj.signif") +
          scale_color_manual(NULL,
                             values=c("1-A"="#377eb8", "1-C"="#4daf4a"), 
                             na.value="black") +
          scale_size_manual(NULL, values=c("1-A"=2, "1-C"=2), na.value=1) +
          scale_y_continuous("N gene eq/mL",
                             limits=c(99, NA), 
                             breaks=c(99, 1e4, 1e8, 1e12),
                             labels=c("<LoD", "1e04", "1e+08", "1e+12"),
                             expand=expansion(c(0.05, 0.1)), 
                             trans="log10") +
          labs(x="Challenge dose (PFU)") +
          guides(color=guide_legend(override.aes=list(size=3)))
fig1$G <- viro.peaks %>%
          filter(Assay=="PFU/mL") %>%
          mutate(value=if_else(value < 25, 24, value)) %>%
          ggplot(aes(Dose, value)) +
          geom_boxplot(fill="grey", outliers=FALSE, alpha=0.5) +
          geom_jitter(aes(col=NHP, size=NHP), height=0, width=0.2) +
          ggpubr::stat_pwc(method="dunn.test", 
                           hide.ns=TRUE, 
                           label="p.adj.signif") +
          scale_color_manual(NULL,
                             values=c("1-A"="#377eb8", "1-C"="#4daf4a"), 
                             na.value="black") +
          scale_size_manual(NULL, values=c("1-A"=2, "1-C"=2), na.value=1) +
          scale_y_continuous("LASV PFU/mL",
                             limits=c(24, NA), 
                             breaks=c(24, 1e3, 1e5, 1e7),
                             labels=c("<LoD", "1e03",  "1e+05", "1e+07"),
                             expand=expansion(c(0.05, 0.1)),
                             trans="log10") +
          labs(x="Challenge dose (PFU)") +
          guides(color=guide_legend(override.aes=list(size=3)))

# now test peak again but by outcome
sup1$A <- viro.peaks %>%
          filter(Assay=="GEq/mL") %>%
          mutate(value=if_else(value==0, 99, value)) %>%
          ggplot(aes(Outcome, value)) +
          geom_boxplot(fill="grey", outliers=FALSE, alpha=0.5) +
          geom_jitter(aes(col=NHP, size=NHP), height=0, width=0.2) +
          ggpubr::stat_pwc(label="p.signif") +
          scale_color_manual(NULL, values=c("1-A"="#377eb8", "1-C"="#4daf4a"),
                             na.value="black") +
          scale_size_manual(NULL, values=c("1-A"=2, "1-C"=2), na.value=1) +
          scale_y_continuous("N gene eq/mL",
                             limits=c(99, NA), 
                             breaks=c(99, 1e4, 1e8, 1e12),
                             labels=c("<LoD", "1e04", "1e+08", "1e+12"),
                             expand=expansion(c(0.05, 0.1)), 
                             trans="log10") +
          labs(x=NULL) +
          guides(color=guide_legend(override.aes=list(size=3))) +
          theme(axis.text.x=element_text(angle=30, hjust=1))
sup1$B <- viro.peaks %>%
          filter(Assay=="PFU/mL") %>%
          mutate(value=if_else(value==0, 24, value)) %>%
          ggplot(aes(Outcome, value)) +
          geom_boxplot(fill="grey", outliers=FALSE, alpha=0.5) +
          geom_jitter(aes(col=NHP, size=NHP), height=0, width=0.2) +
          ggpubr::stat_pwc(label="p.signif") +
          scale_color_manual(NULL, values=c("1-A"="#377eb8", "1-C"="#4daf4a"),
                             na.value="black") +
          scale_size_manual(NULL, values=c("1-A"=2, "1-C"=2),
                            na.value=1) +
          scale_y_continuous("LASV PFU/mL",
                             limits=c(24, NA), 
                             breaks=c(24, 1e3, 1e5, 1e7),
                             labels=c("<LoD", "1e03",  "1e+05", "1e+07"),
                             expand=expansion(c(0.05, 0.1)),
                             trans="log10") +
          labs(x=NULL)  +
          guides(color=guide_legend(override.aes=list(size=3)))  +
          theme(axis.text.x=element_text(angle=30, hjust=1))

# check: what is min/max each outcome?
viro.peaks %>%
  group_by(Assay, Outcome) %>%
  summarise(Max=max(value),
            Min=min(value))

# clean up
rm(viremia, viro.peaks)

# tissue virology
sup1$C <- df$virology.tissues %>%
          filter(Tissue != "Aqueous humor") %>%
          left_join(df$animals, by="NHP") %>%
          mutate(`GEq/g`=if_else(`GEq/g`==0, 1e4-1, `GEq/g`)) %>%
          ggplot(aes(Tissue, `GEq/g`)) +
          geom_boxplot(aes(fill=Outcome, col=Outcome), alpha=0.5, na.rm=TRUE,
                       outlier.size=0.5) +
          ggpubr::stat_compare_means(aes(group=interaction(Outcome, Tissue)), 
                                     method="wilcox.test", label="p.signif",
                                     angle=90, vjust=1, hjust=0) +
          scale_fill_manual(NULL, values=cols.outcome) +
          scale_color_manual(NULL, values=cols.outcome) +
          scale_y_continuous("N gene eq/g",
                             breaks=c(1e4-1, 1e8, 1e12),
                             labels=c("<LoD", "1e+08", "1e+12"),
                             trans="log10", 
                             expand=expansion(mult=c(0.1, 0.2))) +
          labs(x=NULL) +
          theme(axis.text.x=element_text(angle=60, hjust=1),
                legend.position="bottom")

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
# add a "datetime" variable and merge terminal < D14 with D14
clinpath$Datetime <- clinpath$DPI
clinpath$Datetime[(clinpath$DPI != 35) & 
                  (clinpath$DPI==clinpath$Death) &
                  (clinpath$Death < 14)] <- 14
clinpath$Datetime <- factor(clinpath$Datetime, 
                            levels=c(0, 3, 4, 5, 10, 14, 21, 28, 35))

# when do we see changes per analyte in succumbing animals?
# using paired (Freidman pipeline). Exclude 1-C as well so survived=aviremic
clinpath %>%
  filter(!is.na(Datetime),
         NHP != "1-C") %>%
  group_by(Analyte, Outcome) %>%
  rstatix::friedman_test(Concentration ~ Datetime | NHP) %>%
  filter(p < 0.05) %>%
  select(Analyte, Outcome) %>% 
  # add back data for outcome/analytes that pass friedman
  left_join(filter(clinpath, !is.na(Datetime)), by=colnames(.)) %>%
  group_by(Analyte, Outcome) %>%
  rstatix::friedman_nemenyi_test(Concentration ~ Datetime | NHP) %>%
  # only keep comparisons that are significantly different from baseline
  filter(group1==0,
         p.adj < 0.05) %>%
  rstatix::p_format(accuracy=1e-4, digits=3) %>% 
  select(Outcome, Analyte, group2, p.adj) %>%
  filter(Outcome=="Survived")

# plot all clinical chemistry spaghetti
analytes <- clinpath %>%
            filter(Type=="chem") %>%
            select(Analyte) %>%
            unlist() %>%
            droplevels() %>%
            levels()
sup4 <- analytes %>%
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
          labs(x="Days postinfection",
               y=df$units[df$units$Analyte==i, "Units"], 
               title=i) +
          theme(legend.position="none")
        })
names(sup4) <- analytes

# hematology spaghetti
analytes <- clinpath %>%
            filter(Type=="hema") %>%
            select(Analyte) %>%
            unlist() %>%
            droplevels() %>%
            levels()
sup3 <- analytes %>%
        lapply(function(i) {
          clinpath %>%
            filter(Analyte==i) %>%
            ggplot(aes(DPI, Concentration)) +
            geom_line(aes(group=NHP, col=Group), 
                      linewidth=0.5, lineend="round") +
            scale_color_manual(values=cols.group) +
            scale_x_continuous(breaks=breaks.dates) +
            labs(x="Days postinfection",
                 y=df$units[df$units$Analyte==i, "Units"],
                 title=i) +
            theme(legend.position="none")
        })
names(sup3) <- analytes

# subset selected spaghetti
fig3$A <- sup4$ALB
fig3$B <- sup3$Platelets
fig3$C <- sup3$Neutrophils
fig3$D <- sup4$CRP
fig3$E <- sup4$ALT

# get baseline, min, and max per animal per analyte
baseline <- clinpath %>%
            filter(DPI==0) %>%
            mutate(Comparison="D0")
conc.min <- clinpath %>%
            group_by(Analyte, NHP) %>%
            slice_min(order_by=Concentration, n=1, with_ties=FALSE) %>%
            ungroup() %>%
            mutate(Comparison="Min") 
conc.max <- clinpath %>%
            group_by(Analyte, NHP) %>%
            slice_max(order_by=Concentration, n=1, with_ties=FALSE) %>%
            ungroup() %>%
            mutate(Comparison="Max") 
clinpath <- baseline %>%
            rbind(conc.min) %>%
            rbind(conc.max)
rm(conc.min, conc.max, baseline)


# compare min/max to baseline (FN pipeline)
p.fn <- clinpath %>%
        group_by(Outcome, Analyte) %>%
        rstatix::friedman_test(Concentration ~ Comparison | NHP) %>%
        filter(p < 0.05) %>%
        select(Outcome, Analyte) %>%
        left_join(clinpath, by=colnames(.)) %>%
        group_by(Outcome, Analyte) %>%
        rstatix::friedman_nemenyi_test(Concentration ~ Comparison | NHP) %>%
        filter(p.adj < 0.05, 
               group1=="D0") %>%
        # add x positions
        mutate(xmin=1, 
               xmax=if_else(group2=="Max", 2, 3),
               xnudge=factor(Outcome, labels=c(-0.2, 0.2)),
               xnudge=as.numeric(as.character(xnudge)),
               xmax=xmax+xnudge,
               xmin=xmin+xnudge) %>%
        arrange(xmax)
# compare between outcomes (MWU test)
p.mw <- clinpath %>%
        group_by(Analyte, Comparison) %>%
        rstatix::wilcox_test(Concentration ~ Outcome) %>%
        filter(p < 0.05) %>%
        rstatix::add_significance() %>%
        rstatix::add_y_position(step.increase=0) %>%
        mutate(xmin=as.numeric(factor(Comparison)),
               xmax=xmin)

# plot all analytes
analytes <- unique(clinpath$Analyte) %>%
            lapply(function(i) {
              # subset data
              cmat <- filter(clinpath, Analyte==i)
              
              # init increment and ypos
              increment <- 0.15*diff(range(cmat$Concentration))
              ypos <- max(cmat$Concentration)+increment
              
              # init plot
              p <- cmat %>%
                    ggplot(aes(Comparison, Concentration)) +
                    geom_boxplot(aes(group=interaction(Outcome, Comparison)),
                                 fill="lightgrey", col="black",
                                 alpha=0.5, outliers=FALSE) +
                    geom_point(aes(group=interaction(Outcome, Comparison),
                                   col=Group, size=Group), 
                               position=position_jitterdodge(jitter.height=0)) +
                    scale_color_manual(values=cols.group) +
                    scale_size_manual(values=c("1-C"=1), na.value=0.5) +
                    scale_y_continuous(df$units[df$units$Analyte==i, "Units"],
                                       expand=expansion(mult=c(0.1, 0.2))) +
                    labs(x=NULL, title=i) +
                    theme(legend.position="none")
              
              # add outcome comparisons
              pval <- filter(p.mw, Analyte==i)
              if(nrow(pval) > 0) {
                pval$y.position <- pval$y.position+increment
                p <- p + ggpubr::stat_pvalue_manual(data=pval, 
                                                    remove.bracket=TRUE)
                # increment ypos, if necessary
                ypos <- max(c(ypos, pval$y.position))
              }
              
              # add baseline comparisons
              pval <- filter(p.fn, Analyte==i)
              if(nrow(pval) > 0) {
                pval$y.position <- increment*(1:nrow(pval))
                pval$y.position <- pval$y.position + ypos
                p <- p + ggpubr::stat_pvalue_manual(data=pval, color="Outcome")
              }
              
              # fin
              return(p)
            })
names(analytes) <- unique(clinpath$Analyte)
rm(p.fn, p.mw)

# plot all chemistry
sup6 <- clinpath %>%
        filter(Type=="chem") %>%
        select(Analyte) %>%
        unlist() %>%
        droplevels() %>%
        levels()
sup6 <- analytes[sup6]
# plot all hematology
sup5 <- clinpath %>%
        filter(Type=="hema") %>%
        select(Analyte) %>%
        unlist() %>%
        droplevels() %>%
        levels()
sup5 <- analytes[sup5]

# subset selected peaks
fig3$F <- sup6$ALB
fig3$G <- sup5$Platelets
fig3$H <- sup5$Neutrophils
fig3$I <- sup6$CRP
fig3$J <- sup6$ALT

# clean up
rm(clinpath, analytes)

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
# add pseudocount to samples <LOD
elisas$Titer[elisas$Titer==0] <- 50

# plot it and add statistical comparisons
fig6$A <- elisas %>%
          ggplot(aes(Antigen, Titer)) +
          geom_boxplot(aes(group=interaction(Antigen, Outcome)), 
                       fill="lightgrey", col="black",
                       alpha=0.5, outliers=FALSE) +
          geom_point(aes(group=interaction(Antigen, Outcome), 
                         size=Group, col=Group), 
                     position=position_jitterdodge()) +
          scale_color_manual(NULL, values=cols.group) +
          scale_size_manual(NULL, values=c("1-C"=1), na.value=0.5) +
          ggpubr::stat_pwc(aes(group=Outcome), label="p.signif") +
          scale_y_continuous("Endpoint titer",
                             breaks=c(50, 200, 1000, 10000),
                             labels=c("<LOD", 200, 1000, "10,000"),
                             trans="log10",
                             expand = expansion(mult = c(0.05, 0.1))) +
          facet_wrap(~Immunoglobulin) +
          labs(x="LASV antigen", title="ELISAs") +
          guides(fill="none", size="none", 
                 col=guide_legend(override.aes=list(size=3)))
  
# clean up
rm(elisas)

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
sup8 <- prnt %>%
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

# plot neutralization for each at the lowest dilution and compare outcomes
fig6$B <- prnt %>%
          filter(`Dilution factor`==10,
                 DPI %in% c("Terminal", "35")) %>%
          group_by(NHP) %>%
          slice_max(order_by=PercentReduction, n=1, with_ties=FALSE) %>%
          ungroup() %>%
          ggplot(aes(Outcome, PercentReduction)) +
          geom_hline(yintercept=50, linetype=3) +
          geom_boxplot(alpha=0.5, outliers=FALSE, 
                       col="black", fill="lightgrey") +
          geom_jitter(aes(col=Group, size=Group), height=0, width=0.2) +
          ggpubr::stat_pwc(aes(group=Outcome), label="p.signif") +
          scale_color_manual(NULL, values=cols.group) +
          scale_size_manual(NULL, values=c("1-C"=1), na.value=0.5) +
          scale_y_continuous("LASV neutralization", limits=c(0, 100)) +
          labs(x="Outcome",
               title="LASV neutralization") +
          theme(legend.position="none")

# clean up
rm(prnt, vc)

## assemble figures ------------------------------------------------------------
# figure 1
x <- cowplot::plot_grid(fig1$A, fig1$B, ncol=1, labels=c("A", "B"))
x <- cowplot::plot_grid(x, fig1$C, ncol=2, labels=c(NA, "C"))
y <- cowplot::plot_grid(plotlist=fig1[4:7], labels=LETTERS[4:7])
cowplot::plot_grid(x, y, ncol=1)
ggsave("analysis/figure1.pdf", units="in", width=6.5, height=9)
rm(x, y, fig1)

# figure 3
# get line and point legends
x <- fig3$A + 
     labs(col=NULL) + 
     guides(color=guide_legend(override.aes=list(linewidth=1))) +
     theme(legend.position="right")
x <- cowplot::get_legend(x)
x <- cowplot::ggdraw(x)
x <- cowplot::plot_grid(plotlist=append(fig3[1:5], x), labels=LETTERS[1:5])
y <- fig3$F + 
     labs(col=NULL, size=NULL) + 
     guides(color=guide_legend(override.aes=list(size=3)),
            size="none") +
     theme(legend.position="right")
y <- cowplot::get_legend(y)
y <- cowplot::ggdraw(y)
y <- cowplot::plot_grid(plotlist=append(fig3[6:10], y), labels=LETTERS[6:10])
cowplot::plot_grid(x, y, ncol=1)
ggsave("analysis/figure3.pdf", units="in", width=6.5, height=8)
rm(x, y, fig3)

# figure 6
cowplot::plot_grid(plotlist=fig6, nrow=1, rel_widths=c(4, 2.5), labels="AUTO")
ggsave("analysis/figure6.pdf", units="in", width=6.5, height=2.5)
rm(fig6)

# supplemental 1
x <- cowplot::plot_grid(plotlist=sup1[1:2], nrow=1, labels="AUTO")
cowplot::plot_grid(x, sup1$C, labels=c(NA, "C"), ncol=1, rel_heights=c(1, 2))
ggsave("analysis/supplemental1.png", units="in", width=6.5, height=7)
rm(x, sup1)

# supplemental 3
sup3$legs <- sup3$Leukocytes + 
             labs(col=NULL, size=NULL) + 
             guides(color=guide_legend(override.aes=list(linewidth=1)),
                    size="none") +
             theme(legend.position="right")
sup3$legs <- cowplot::get_legend(sup3$legs)
sup3$legs <- cowplot::ggdraw(sup3$legs)
cowplot::plot_grid(plotlist=sup3, ncol=3, labels=LETTERS[1:length(sup3)-1])
ggsave("analysis/supplemental3.png", units="in", width=6.5, height=6)
rm(sup3)

# supplemental 4
legs <- sup4$GLU + 
        labs(col=NULL, size=NULL) + 
        guides(color=guide_legend(override.aes=list(linewidth=1)),
               size="none") +
        theme(legend.position="bottom")
legs <- cowplot::get_legend(legs)
legs <- cowplot::ggdraw(legs)
sup4 <- cowplot::plot_grid(plotlist=sup4, ncol=3, labels="AUTO")
cowplot::plot_grid(sup4, legs, ncol=1, rel_heights=c(20, 1))
ggsave("analysis/supplemental4.png", units="in", width=6.5, height=7)
rm(sup4, legs)

# supplemental 5
sup5$legs <- sup5$Leukocytes + 
             labs(col=NULL, size=NULL) + 
             guides(color=guide_legend(override.aes=list(size=3)),
                    size="none") +
             theme(legend.position="right")
sup5$legs <- cowplot::get_legend(sup5$legs)
sup5$legs <- cowplot::ggdraw(sup5$legs)
cowplot::plot_grid(plotlist=sup5, ncol=3, labels=LETTERS[1:length(sup5)-1])
ggsave("analysis/supplemental5.png", units="in", width=6.5, height=6)
rm(sup5)

# supplemental 6
legs <- sup6$GLU + 
        labs(col=NULL, size=NULL) + 
        guides(color=guide_legend(override.aes=list(size=3)),
               size="none") +
        theme(legend.position="bottom")
legs <- cowplot::get_legend(legs)
legs <- cowplot::ggdraw(legs)
sup6 <- cowplot::plot_grid(plotlist=sup6, ncol=3, labels="AUTO")
cowplot::plot_grid(sup6, legs, ncol=1, rel_heights=c(20, 1))
ggsave("analysis/supplemental6.png", units="in", width=6.5, height=7)
rm(sup6, legs)

# supplemental 8
ggsave("analysis/supplemental8.png", sup8,
       units="in", width=6.5, height=3)
rm(sup8)

## fin -------------------------------------------------------------------------
sessionInfo()
