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
fig5 <- list()
sup7 <- list()

# helper functions
run.pca <- function(logconcentration.mat, meta.mat) {
  # format counts to log2 CPM, then run PCA
  pca <- logconcentration.mat %>%
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
run.de <- function(logconcentration.mat, model.mat, contrasts.obj) {
  logconcentration.mat %>%
    lmFit(model.mat) %>%
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
meta <- readxl::read_excel("data.xlsx", "legendplex.samplesheet") %>%
        mutate(Stage=factor(Stage, levels=names(cols.stage)),
               Outcome=factor(Outcome, levels=names(shapes.outcome)),
               Group=factor(Group, levels=names(alpha.group)),
               Batch=factor(Group, labels=c("A", "B", "B"))) %>%
        as.data.frame() 
rownames(meta) <- meta$ID

# concentration data; remove poor-performing analytes
discard.pile <- c("IFN-y_cc", "IL-1β_cc", "IL-6_cc", "MCP-1_cc", "TNF-α_cc", 
                  "IL-8_cc", "IL-6_th", "IL-8_th", "IP-10_inf")
cmat <- readxl::read_excel("data.xlsx", "legendplex.concentration") %>%
        # merge replicates 
        reshape2::melt(id.vars=c("ID", "Replicate"),
                       variable.name="Analyte",
                       value.name="Concentration") %>%
        group_by(ID, Analyte) %>%
        summarise(Concentration=mean(Concentration),
                  .groups="drop") %>%
        # remove selected analytes, then remove panel names since no more dups
        filter(!(Analyte %in% discard.pile)) %>%
        mutate(Analyte=str_remove(Analyte, "_.+")) %>%
        # format with analytes in rows and samples in columns
        reshape2::dcast(Analyte ~ ID, value.var="Concentration") %>%
        column_to_rownames("Analyte") %>%
        as.matrix()
rm(discard.pile)

# align rows and columns 
x <- intersect(rownames(meta), colnames(cmat))
meta <- meta[x, ]
cmat <- cmat[, x]
rm(x)

## PCAs ------------------------------------------------------------------------
# pca on all samples, remove batch effect
pca <- model.matrix(~0 + Stage, data=meta)
pca <- cmat %>%
       log2() %>%
       removeBatchEffect(batch=meta$Batch, design=pca) %>%
       run.pca(meta)

# make labels for 1-C and 1-A, then plot
labelpts <- pca$pca %>%
            filter(Group=="LD-V",
                   !(Stage %in% c("Baseline", "Early", "Recovered"))) %>%
            mutate(Label=paste0(NHP, ", ", DPI, " DPI")) 
sup7$A <- pca$pca %>%
          arrange(Group=="LD-V") %>% 
          ggplot(aes(PC1, PC2, fill=Stage, shape=Outcome, alpha=Group)) +
          geom_hline(yintercept=0, linetype=3, col="lightgrey") +
          geom_vline(xintercept=0, linetype=3, col="lightgrey") +
          # all points
          geom_point(size=3) +
          scale_fill_manual(values=cols.stage) +
          scale_shape_manual(values=shapes.outcome) +
          scale_alpha_manual(values=alpha.group, guide="none") +
          # add labels
          ggrepel::geom_text_repel(data=labelpts, aes(label=Label), 
                                   size=3, fontface="bold", nudge_y=1,
                                   min.segment.length=2) +
          # format axes and legends
          labs(x=pca$pcs[1],
               y=pca$pcs[2],
               shape=NULL,
               fill=NULL) +
          guides(shape=guide_legend(override.aes=list(fill="lightgrey", size=5)),
                 fill=guide_legend(override.aes=list(size=5, pch=21)))
# now subset to HD and LD-V
m <- filter(meta, Group!="LD-NV")
pca <- model.matrix(~0 + Stage, data=m)
pca <- cmat[, m$ID] %>%
       log2() %>%
       removeBatchEffect(batch=m$Batch, design=pca) %>%
       run.pca(m)

# make labels for 1-C and 1-A, then plot
labelpts <- pca$pca %>%
            filter(Group=="LD-V",
                   !(Stage %in% c("Baseline", "Early", "Recovered"))) %>%
            mutate(Label=paste0(NHP, ", ", DPI, " DPI")) 
fig5$A <- pca$pca %>%
          arrange(Group=="LD-V") %>% 
          ggplot(aes(PC1, PC2, fill=Stage, shape=Outcome, alpha=Group)) +
          geom_hline(yintercept=0, linetype=3, col="lightgrey") +
          geom_vline(xintercept=0, linetype=3, col="lightgrey") +
          # all points
          geom_point(size=3) +
          scale_fill_manual(values=cols.stage) +
          scale_shape_manual(values=shapes.outcome) +
          scale_alpha_manual(values=alpha.group, guide="none") +
          # add labels
          ggrepel::geom_text_repel(data=labelpts, aes(label=Label), 
                                   size=3, fontface="bold",
                                   min.segment.length=2, nudge_y=1) +
          # format axes and legends
          labs(x=pca$pcs[1],
               y=pca$pcs[2],
               shape=NULL,
               fill=NULL) +
          guides(shape=guide_legend(override.aes=list(fill="lightgrey", size=5)),
                 fill=guide_legend(override.aes=list(size=5, pch=21)))
# now plot contributors
fig5$B <- pca$contribs %>%
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
rm(pca, m, labelpts)

## viremic/symptomatic longitudinal --------------------------------------------
# subset data
m <- meta[meta$Group != "LD-NV", ]

# define our model and contrasts
model <- model.matrix(~0 + Stage + Batch, data=m)
colnames(model) <- c(levels(m$Stage), "Batch")
contrs <- makeContrasts(Early - Baseline,
                        Middle - Baseline,
                        Late - Baseline,
                        levels=colnames(model))

# make DE model
res.viro <- cmat[, rownames(m)] %>%
            log2() %>%
            run.de(model, contrs)

# extract results table, format, and save
res.viro <- colnames(contrs) %>%
            lapply(format.results, ebayes.out=res.viro) %>%
            do.call(rbind, .) %>%
            mutate(Stage=factor(Stage, levels=names(cols.stage)))
write.csv(res.viro, "analysis/legendplex-viremic.csv", row.names=FALSE)

# plot volcano
topgenes <- res.viro %>%
            filter(Regulation=="Up") %>%
            group_by(Stage) %>%
            top_n(n=10, wt=lfc) %>%
            ungroup()
fig5$C <- res.viro %>%
          arrange(Significant) %>%
          ggplot(aes(lfc, -log10(padj))) +
          geom_point(aes(size=Regulation, col=Regulation), alpha=0.6) +
          scale_color_manual(values=cols.reg, na.value="lightgrey") +
          scale_size_manual(values=c(Up=1, Down=1), na.value=0.5) +
          ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), 
                                   size=2, fontface="bold", force=50,
                                   nudge_y=1, min.segment.length=2) +
          scale_x_continuous(limits=c(-4, 8),
                             breaks=c(-4, 0, 4, 8)) +
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

# define our model and contrasts (no batch for this one)
model <- model.matrix(~0 + Stage, data=m)
colnames(model) <- levels(m$Stage)
contrs <- makeContrasts(Early - Baseline,
                        Middle - Baseline,
                        Late - Baseline,
                        levels=colnames(model))

# make DE model
res.avir <- cmat[, rownames(m)] %>%
            log2() %>%
            run.de(model, contrs)

# extract results table, format, and save
res.avir <- colnames(contrs) %>%
            lapply(format.results, ebayes.out=res.avir) %>%
            do.call(rbind, .) %>%
            mutate(Stage=factor(Stage, levels=names(cols.stage)))
write.csv(res.avir, "analysis/legendplex-aviremic.csv", row.names=FALSE)

# plot volcano
topgenes <- res.avir %>%
            filter(Significant) %>%
            group_by(Stage, Regulation) %>%
            top_n(n=10, wt=abs(lfc)) %>%
            ungroup()
sup7$B <- res.avir %>%
          arrange(Significant) %>%
          ggplot(aes(lfc, -log10(padj))) +
          geom_point(aes(size=Regulation, col=Regulation), alpha=0.6) +
          scale_color_manual(values=cols.reg, na.value="lightgrey") +
          scale_size_manual(values=c(Up=1, Down=1), na.value=0.5) +
          ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), 
                                   size=2, fontface="bold", force=50,
                                   nudge_y=1, min.segment.length=2) +
          scale_x_continuous(limits=c(-4, 8),
                             breaks=c(-4, 0, 4, 8)) +
          scale_y_continuous(limits=c(NA, 30)) +
          facet_wrap(~Stage) +
          labs(x="Fold change (log2)",
               y="padj (-log10)") +
          theme(legend.position="none")

# clean up
rm(topgenes, m, contrs, model, res.avir)

## plot selected lfcs ----------------------------------------------------------
selected <- c("IP-10", "MCP-1", "MIG", "Eotaxin")

# extract counts and merge with metadata
cmat <- cmat[selected, ] %>%
        as.data.frame() %>%
        rownames_to_column("Gene") %>%
        reshape2::melt(id.vars="Gene",
                       variable.name="ID", 
                       value.name="Counts") %>%
        left_join(meta, by="ID") 

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
               Group=factor(Group, levels=c("HD", "LD-NV", "1-A", "1-C")),
               NHP=factor(NHP, levels=unique(NHP)))
cmat$Group[cmat$NHP=="1-A"] <- "1-A"
cmat$Group[cmat$NHP=="1-C"] <- "1-C"

plotlist <- selected %>%
            lapply(function(i) {
              cmat %>%
                filter(Gene==i) %>%
                arrange(NHP, Group!="HD") %>% 
                ggplot(aes(DPI, lfc)) +
                geom_line(aes(group=NHP, col=Group, linewidth=Group), 
                          lineend="round") +
                scale_linewidth_manual(values=c("1-A"=1, "1-C"=1), 
                                       na.value=0.5, guide="none") +
                scale_color_manual(NULL, 
                                   values=c(HD="lightgrey", 
                                            "LD-NV"="grey30",
                                            "1-A"="#377eb8", 
                                            "1-C"="#4daf4a")) +
                scale_x_continuous(breaks=breaks.dates) +
                labs(x="Days postinfection",
                     y="Fold change (log2)",
                     title=i) +
                theme(legend.position="none")
            })
fig5$D <- plotlist[[1]] + 
          guides(col=guide_legend(override.aes=list(linewidth=1))) +
          labs(title="CXCL10/IP-10") +
          theme(legend.position=c(0.8, 0.9), 
                legend.background=element_blank())
fig5$E <- plotlist[[2]] + labs(title="CCL2/MCP-1")
fig5$F <- plotlist[[3]] + labs(title="CXCL9/MIG")
fig5$G <- plotlist[[4]] + labs(title="CCL11/Eotaxin") 

# clean up
rm(plotlist, selected)

## assemble figures ------------------------------------------------------------
# figure 5
x <- cowplot::plot_grid(plotlist=fig5[-1], ncol=2, labels=LETTERS[2:7])
cowplot::plot_grid(fig5$A, x, ncol=1, rel_heights=c(3, 6), labels=c("A", NA))
ggsave("analysis/figure5.pdf", units="in", width=6.5, height=9)
rm(x, fig5)

# supplemental 7
cowplot::plot_grid(plotlist=sup7, labels="AUTO", ncol=1)
ggsave("analysis/supplemental7.png", units="in", width=6.5, height=6)
rm(sup7)

## fin -------------------------------------------------------------------------
sessionInfo()
