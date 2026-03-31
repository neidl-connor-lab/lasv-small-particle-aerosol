#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
theme_set(ggpubr::theme_pubr() + theme(legend.position="right"))

# helper variables
cols.stage <- c(Baseline="white", 
                Early="#e9c46a", 
                Middle="#2a9d8f", 
                Late="#225ea8", 
                Recovered="grey50")
cols.reg <- c(Up="#e41a1c", Down="#377eb8")
shapes.outcome <- c(Survived=22, Succumbed=21)
alpha.group <- c(HD=0.5, "LD-V"=1, "LD-NV"=0.25)
breaks.dates <- c(0, 7, 14, 21, 28, 35)
fig4 <- list()

# helper functions
format.counts <- function(counts.mat, n=3, threshold=20) {
  # get the set of genes with >= n samples with >threshold counts
  keepgenes <- rowSums(counts.mat > threshold)
  keepgenes <- (keepgenes >= n)
  # return the subsetted matrix
  counts.mat[keepgenes, ]
}
run.pca <- function(counts.mat, meta.mat) {
  # format counts to log2 CPM, then run PCA
  pca <- counts.mat %>%
         DGEList() %>%
         calcNormFactors() %>%
         cpm(log=TRUE) %>%
         t() %>%
         prcomp()
  # get PCs
  pcs <- summary(pca)$importance["Proportion of Variance", 1:2]
  pcs <- round(100*pcs)
  pcs <- paste0(names(pcs), " (", pcs, "%)")
  # get contributions
  contribs <- pca$rotation %>%
              as.data.frame() %>%
              rownames_to_column("Gene") %>%
              select(Gene, PC1, PC2)
  # format PCA output for plotting
  pca <- pca$x %>%
         as.data.frame() %>%
         rownames_to_column("ID") %>%
         select(ID, PC1, PC2) %>%
         left_join(meta.mat, by="ID")
  # return outputs
  list(pca=pca, contribs=contribs, pcs=pcs)
}
run.de <- function(counts.mat, model.mat, contrasts.obj) {
  counts.mat %>%
    DGEList() %>%
    calcNormFactors() %>%
    voom(model.mat) %>%
    lmFit() %>%
    contrasts.fit(contrasts.obj) %>%
    eBayes()
}
format.results <- function(coef, ebayes.output, psig=0.05, fsig=1) {
  # get the full toptable and format
  degs <- topTable(ebayes.output, coef, number=Inf) %>%
          rownames_to_column("Gene") %>%
          rename(lfc=logFC, padj=adj.P.Val) %>%
          mutate(Stage=str_extract(coef, "^[A-z]+"),
                 psig=(padj < psig),
                 fsig=(abs(lfc) > fsig),
                 Significant=(psig & fsig)) %>%
          select(Gene, lfc, padj, Stage, Significant)
  # format regulation column
  degs$Regulation <- "None"
  degs$Regulation[degs$Significant] <- "Up"
  degs$Regulation[degs$Significant & (degs$lfc < 0)] <- "Down"
  degs$Regulation <- factor(degs$Regulation, levels=c("Up", "None", "Down"))
  # return the formatted matrix
  return(degs)
}

## load data -------------------------------------------------------------------
# sample metadata
meta <- readxl::read_excel("data.xlsx", "nanostring.samplesheet") %>%
        filter(QC) %>%
        mutate(Stage=factor(Stage, levels=names(cols.stage)),
               Outcome=factor(Outcome, levels=names(shapes.outcome)),
               Group=factor(Group, levels=names(alpha.group))) %>%
        as.data.frame() 
rownames(meta) <- meta$ID

# thresholded counts
cmat <- readxl::read_excel("data.xlsx", "nanostring.thresholded") %>%
        # remove positive and negative control probes
        filter(!str_detect(Gene, "^POS_|^NEG_")) %>%
        # format as a counts matrix
        column_to_rownames("Gene") %>%
        as.matrix()

# align rows and columns 
x <- intersect(rownames(meta), colnames(cmat))
meta <- meta[x, ]
cmat <- cmat[, x]
rm(x)

## sample timeline -------------------------------------------------------------
fig4$A <- meta %>%
          select(DPI, Stage) %>%
          distinct() %>%
          ggplot(aes(DPI, 0)) +
          annotate("segment", x=0, xend=35, y=0, linewidth=1) +
          geom_point(aes(fill=Stage), pch=21, size=3) +
          scale_fill_manual(NULL, values=cols.stage) +
          geom_text(aes(label=DPI), y=1.75, size=3) +
          ylim(-1, 2) +
          theme_void() +
          theme(legend.position="top")

## PCAs ------------------------------------------------------------------------
# pca of all samples
pca <- cmat %>%
       format.counts() %>%
       run.pca(meta)

# make labels for 1-C and 1-A, then plot
labelpts <- pca$pca %>%
            filter(Group=="LD-V",
                   !(Stage %in% c("Baseline", "Early", "Recovered"))) %>%
            mutate(Label=paste0(NHP, ", ", DPI, " DPI")) 
fig4$B <- pca$pca %>%
          arrange(Group=="LD-V") %>% 
          ggplot(aes(PC1, PC2, fill=Stage, shape=Outcome, alpha=Group)) +
          geom_hline(yintercept=0, linetype=3, col="lightgrey") +
          geom_vline(xintercept=0, linetype=3, col="lightgrey") +
          # all points
          geom_point() +
          scale_fill_manual(values=cols.stage, guide="none") +
          scale_shape_manual(values=shapes.outcome) +
          scale_alpha_manual(values=alpha.group, guide="none") +
          # add labels
          ggrepel::geom_text_repel(data=labelpts, aes(label=Label), 
                                   size=3, fontface="bold",
                                   min.segment.length=2) +
          # format axes and legends
          labs(x=pca$pcs[1],
               y=pca$pcs[2],
               shape=NULL) +
          guides(shape=guide_legend(override.aes=list(fill="lightgrey", 
                                                      size=3))) +
          theme(legend.position="none")
# save pca legend
fig4$leg <- fig4$B + theme(legend.position="bottom")
fig4$leg <- cowplot::get_legend(fig4$leg)
fig4$leg <- cowplot::ggdraw(fig4$leg)
rm(labelpts)

# contributors plot: plot top 30 using vector length
fig4$C <- pca$contribs %>%
          mutate(Magnitude=sqrt(PC1^2 + PC2^2)) %>%
          slice_max(n=30, order_by=Magnitude, with_ties=FALSE) %>%
          ggplot() +
          geom_vline(xintercept=0, linetype=3, col="lightgrey") +
          geom_hline(yintercept=0, linetype=3, col="lightgrey") +
          geom_segment(aes(x=0, y=0, xend=PC1, yend=PC2), col="darkgrey",
                       alpha=0.9, arrow=arrow(length=unit(0.05, "in"))) +
          ggrepel::geom_text_repel(aes(PC1, PC2, label=Gene), 
                                   size=2, fontface="bold") +
          labs(x=pca$pcs[1],
               y=pca$pcs[2]) +
          theme(legend.position="none")
# clean up
rm(pca)

## viremic/symptomatic longitudinal --------------------------------------------
# subset data
m <- meta[meta$Group != "LD-NV", ]

# define our model and contrasts
model <- model.matrix(~0 + Stage, data=m)
colnames(model) <- levels(m$Stage)
contrs <- makeContrasts(Early - Baseline,
                        Middle - Baseline,
                        Late - Baseline,
                        levels=colnames(model))

# make DE model
res.viro <- cmat[, rownames(m)] %>%
            format.counts() %>%
            run.de(model, contrs)

# extract results table, format, and save
res.viro <- colnames(contrs) %>%
            lapply(format.results, ebayes.out=res.viro) %>%
            do.call(rbind, .) %>%
            mutate(Stage=factor(Stage, levels=names(cols.stage)))
write.csv(res.viro, "analysis/nanostring-viremic.csv", row.names=FALSE)

# save for IPA 
res.viro %>%
  select(Gene, Stage, lfc, padj) %>%
  reshape2::melt(id.vars=c("Gene", "Stage")) %>%
  mutate(Group=paste(Stage, variable, sep=".")) %>%
  reshape2::dcast(Gene ~ Group, value.var="value") %>%
  write.csv("analysis/ipa-input-viremic.csv", row.names=FALSE)

# plot volcano
topgenes <- res.viro %>%
            filter(Regulation=="Up") %>%
            group_by(Stage) %>%
            top_n(n=8, wt=lfc) %>%
            ungroup()
fig4$D <- res.viro %>%
          arrange(Significant) %>%
          ggplot(aes(lfc, -log10(padj))) +
          geom_point(aes(size=Regulation, col=Regulation), alpha=0.6) +
          scale_color_manual(values=cols.reg, na.value="lightgrey") +
          scale_size_manual(values=c(Up=1, Down=1), na.value=0.5) +
          ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), 
                                   size=2, fontface="bold", force=50,
                                   nudge_y=1, min.segment.length=2) +
          scale_x_continuous(limits=c(-3, 8),
                             breaks=c(-3, 0, 3, 6)) +
          scale_y_continuous(limits=c(NA, 30)) +
          facet_wrap(~Stage) +
          labs(x="Fold change (log2)",
               y="padj (-log10)") +
          theme(legend.position="none")
# clean up
rm(topgenes, m, contrs, model, res.viro)

## aviremic longitudinal -------------------------------------------------------
# subset data
m <- meta[meta$Group=="LD-NV", ]

# define our model and contrasts
model <- model.matrix(~0 + Stage, data=m)
colnames(model) <- levels(m$Stage)
contrs <- makeContrasts(Early - Baseline,
                        Middle - Baseline,
                        Late - Baseline,
                        levels=colnames(model))

# make DE model
res.avir <- cmat[, rownames(m)] %>%
            format.counts() %>%
            run.de(model, contrs)

# extract results table, format, and save
res.avir <- colnames(contrs) %>%
            lapply(format.results, ebayes.out=res.avir) %>%
            do.call(rbind, .) %>%
            mutate(Stage=factor(Stage, levels=names(cols.stage)))
write.csv(res.avir, "analysis/nanostring-aviremic.csv", row.names=FALSE)

# save for IPA 
res.avir %>%
  select(Gene, Stage, lfc, padj) %>%
  reshape2::melt(id.vars=c("Gene", "Stage")) %>%
  mutate(Group=paste(Stage, variable, sep=".")) %>%
  reshape2::dcast(Gene ~ Group, value.var="value") %>%
  write.csv("analysis/ipa-input-aviremic.csv", row.names=FALSE)

# plot volcano
topgenes <- res.avir %>%
            filter(Regulation=="Up") %>%
            group_by(Stage) %>%
            top_n(n=5, wt=lfc) %>%
            ungroup()
fig4$E <- res.avir %>%
          arrange(Significant) %>%
          ggplot(aes(lfc, -log10(padj))) +
          geom_point(aes(size=Regulation, col=Regulation), alpha=0.6) +
          scale_color_manual(values=cols.reg, na.value="lightgrey") +
          scale_size_manual(values=c(Up=1, Down=1), na.value=0.5) +
          ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), 
                                   size=2, fontface="bold", force=50,
                                   nudge_y=1, min.segment.length=2) +
          scale_x_continuous(limits=c(-3, 8),
                             breaks=c(-3, 0, 3, 6)) +
          scale_y_continuous(limits=c(NA, 30)) +
          facet_wrap(~Stage) +
          labs(x="Fold change (log2)",
               y="padj (-log10)") +
          theme(legend.position="none")
# clean up
rm(topgenes, m, contrs, model)

## IPA, viremic ----------------------------------------------------------------
# load shortlist
ipa.pathlist <- read.csv("analysis/ipa-shortlist.csv")

# load IPA data and subset
ipa.pathlist <- read.delim("analysis/ipa-output.tsv") %>%
                right_join(ipa.pathlist, by="Canonical.Pathways") %>%
                mutate(Nickname=factor(Nickname, 
                                       levels=ipa.pathlist$Nickname)) %>%
                arrange(Nickname) %>%
                select(-Canonical.Pathways) %>%
                replace_na(list(Late=0, Early=0)) %>%
                column_to_rownames("Nickname") %>%
                t()
# define min and max of each column, then format dataframe 
ipa.pathlist <- data.frame(Max=rep(3, dim(ipa.pathlist)[2]),
                           Min=0) %>%
                t() %>%
                rbind(ipa.pathlist)
ipa.pathlist <- ipa.pathlist[c("Max", "Min", "Early", "Late"), ]
ipa.pathlist <- as.data.frame(ipa.pathlist)

# format pathway names so they will render newlines
varnames <- str_replace(colnames(ipa.pathlist), "\\\\n", "\n")

# plot it
fig4$F <- ggplotify::as.ggplot(
            ~fmsb::radarchart(ipa.pathlist,
                             seg=3, 
                             pcol=cols.stage[c("Early", "Late")],
                             plty=1, plwd=1,
                             cglty=1, cglcol="lightgrey", 
                             vlabels=varnames, vlcex=0.75))

# clean up
rm(ipa.pathlist, varnames)

## heatmap showing phenotypes --------------------------------------------------
# gene subset
genelist <- c("CXCL10", "IFIT2", "MX1", "IFIT3", "IFI44", "OAS2", "OAS3", 
              "OASL", "GBP1", "SERPING1","GZMA")

# subset results from LD-NV
ldnv <- res.avir %>%
        filter(Gene %in% genelist) %>%
        mutate(Group="LD-NV") %>%
        select(Gene, lfc, Group, Stage)

# run DE for just LD-V (no real significance testing, just getting lfc)
m <- filter(meta, Group=="LD-V")
model <- model.matrix(~0 + Stage, data=m)
colnames(model) <- levels(m$Stage)
contrs <- makeContrasts(Early - Baseline,
                        Middle - Baseline,
                        Late - Baseline,
                        levels=colnames(model))
ldv <- cmat[, rownames(m)] %>%
       format.counts() %>%
       run.de(model, contrs)
ldv <- colnames(contrs) %>%
       lapply(format.results, ebayes.out=ldv) %>%
       do.call(rbind, .) %>%
       filter(Gene %in% genelist) %>%
       mutate(Group="LD-V") %>%
       select(Gene, lfc, Group, Stage)
rm(m, model, contrs)

# run DE for HD 
m <- meta %>%
     filter(Group=="HD") %>%
     mutate(Stage=droplevels(Stage))
model <- model.matrix(~0 + Stage, data=m)
colnames(model) <- levels(m$Stage)
contrs <- makeContrasts(Early - Baseline,
                        Middle - Baseline,
                        Late - Baseline,
                        levels=colnames(model))
hd <- cmat[, rownames(m)] %>%
      format.counts() %>%
      run.de(model, contrs)
hd <- colnames(contrs) %>%
      lapply(format.results, ebayes.out=hd) %>%
      do.call(rbind, .) %>%
      filter(Gene %in% genelist) %>%
      mutate(Group="HD") %>%
      select(Gene, lfc, Group, Stage)
rm(m, model, contrs)

# make a helper variable to control splitting up the columns
csplit <- c(rep("Early", 3), rep("Middle", 3), rep("Late", 3))
csplit <- factor(csplit, levels=names(cols.stage))

# make a helper color palette variable
colfun <- circlize::colorRamp2(c(-7, 0, 7), 
                               c(cols.reg["Down"], "white", cols.reg["Up"]))

# combine tables and plot
fig4$G <- hd %>%
          rbind(ldv) %>%
          rbind(ldnv) %>%
          mutate(Group=factor(Group, levels=c("HD", "LD-V", "LD-NV")),
                 Stage=factor(Stage, levels=names(cols.stage)),
                 Label=paste(Group, tolower(Stage), sep=", ")) %>%
          arrange(Stage, Group) %>%
          mutate(Label=factor(Label, levels=unique(Label)))%>%
          reshape2::dcast(Gene ~ Label, value.var="lfc", fill=0) %>%
          column_to_rownames("Gene") %>%
          as.matrix() %>% 
          ComplexHeatmap::Heatmap(name="lfc", 
                                  rect_gp=grid::gpar(col="black", lwd=1), 
                                  col=colfun,
                                  cluster_columns=FALSE, 
                                  column_split=csplit,
                                  column_labels=rep(c("HD", "LD-V", "LD-NV"), 3),
                                  column_names_gp=grid::gpar(fontsize=8),
                                  row_names_gp=grid::gpar(fontsize=8),
                                  column_title_gp=grid::gpar(fontsize=10)) %>%
          ComplexHeatmap::draw() %>%
          grid::grid.grabExpr() %>%
          ggplotify::as.ggplot()

# clean up
rm(hd, ldv, ldnv, csplit, colfun, genelist, res.avir)

# make lfc plots for selected genes --------------------------------------------
genelist <- c("CXCL10", "GZMA")

# extract counts and merge with metadata
cmat <- cmat[genelist, ] %>%
        as.data.frame() %>%
        rownames_to_column("Gene") %>%
        reshape2::melt(id.vars="Gene",
                       variable.name="ID", 
                       value.name="Counts") %>%
        left_join(meta, by="ID") %>%
        # only want HD and LD-V
        filter(Group != "LD-NV")

# extract and add back baseline to calculate fold change (log2)
cmat <- cmat %>%
        filter(DPI==0) %>%
        rename(Baseline=Counts) %>%
        select(Gene, NHP, Baseline) %>%
        # add baseline to counts
        right_join(cmat, by=c("NHP", "Gene")) %>%
        arrange(Group!="HD") %>%
        # calculate lfc
        mutate(lfc=log2(Counts/Baseline),
               # also set up a new group factor for colors
               Group=factor(NHP, levels=c("HD", "1-A", "1-C")),
               NHP=factor(NHP, levels=unique(NHP))) %>%
        replace_na(list(Group="HD"))

# plot CXCL10
fig4$H <- cmat %>%
          filter(Gene=="CXCL10") %>%
          arrange(NHP, Group!="HD") %>% 
          ggplot(aes(DPI, lfc)) +
          geom_line(aes(group=NHP, col=Group, linewidth=Group), 
                    lineend="round") +
          scale_linewidth_manual(values=c("1-A"=1, "1-C"=1), 
                                 na.value=0.5, guide="none") +
          scale_color_manual(NULL, 
                             values=c(HD="lightgrey", 
                                      "1-A"="#377eb8", 
                                      "1-C"="#4daf4a")) +
          scale_x_continuous(breaks=breaks.dates) +
          labs(x="Days postinfection",
               y="Fold change (log2)",
               title="CXCL10/IP-10") +
          guides(col=guide_legend(override.aes=list(linewidth=1))) +
          theme(legend.position=c(0.9, 0.7),
                legend.background=element_blank())

# plot GZMA
fig4$I <- cmat %>%
          filter(Gene=="GZMA") %>%
          arrange(NHP, Group!="HD") %>% 
          ggplot(aes(DPI, lfc)) +
          geom_line(aes(group=NHP, col=Group, linewidth=Group), 
                    lineend="round") +
          scale_linewidth_manual(values=c("1-A"=1, "1-C"=1), 
                                 na.value=0.5, guide="none") +
          scale_color_manual(NULL, 
                             values=c(HD="lightgrey", 
                                      "1-A"="#377eb8", 
                                      "1-C"="#4daf4a")) +
          scale_x_continuous(breaks=breaks.dates) +
          labs(x="Days postinfection",
               y="Fold change (log2)",
               title="GZMA") +
          guides(col=guide_legend(override.aes=list(linewidth=2))) +
          theme(legend.position="none")

# clean up
rm(genelist)

## build figure ----------------------------------------------------------------
# put B-D in column
x <- cowplot::plot_grid(plotlist=fig4[c("leg", "B", "C", "D", "E")],
                        ncol=1, rel_heights=c(0.3, 2, 2, 2, 2),
                        labels=c(NA, "B", "C", "D", "E"))

# put F-I in a column
y <- cowplot::plot_grid(plotlist=fig4[c("F", "G", "H", "I")],
                        ncol=1, rel_heights=c(2, 2.3, 2, 2), 
                        labels=c("F", "G", "H", "I"))

# combine columns then add timeline on top
x <- cowplot::plot_grid(x, y, ncol=2)
cowplot::plot_grid(fig4$A, x, ncol=1, 
                   labels=c("A", NA), 
                   rel_heights=c(0.65, 8.3))
ggsave("analysis/figure4.pdf", units="in", width=6.5, height=9)

# clean up
rm(x, y, fig4)

## fin -------------------------------------------------------------------------
sessionInfo()

