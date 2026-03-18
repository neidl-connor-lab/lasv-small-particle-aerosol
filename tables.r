#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
options(stringsAsFactors=FALSE)

## inputs ----------------------------------------------------------------------
# load the data 
sheets <- readxl::excel_sheets("data.xlsx")
df <- lapply(sheets, function(i) { readxl::read_excel("data.xlsx", sheet=i)})
names(df) <- sheets
rm(sheets)

## clinical illness ------------------------------------------------------------
# temperature
baseline <- df$WT %>%
            filter(DPI==0) %>%
            mutate(Baseline=Temperature) %>%
            select(NHP, Baseline)
clinical <- df$WT %>%
            left_join(baseline, by="NHP") %>%
            mutate(Delta=Temperature-Baseline,
                   Fever=(Delta > 2.5) | (Delta > 1.5 & Temperature > 103.5),
                   Hypothermia=(Delta <= -3.5)) %>%
            filter(Fever | Hypothermia) %>%
            select(NHP, DPI, Fever, Hypothermia) %>%
            # melt to long
            reshape2::melt(id.vars=c("NHP", "DPI"),
                           variable.name="Change",
                           value.name="Present") %>%
            filter(Present) %>%
            # compile by group
            group_by(NHP, Change) %>%
            summarise(DPI=paste(DPI, collapse=", "),
                      .groups="drop") %>%
            mutate(Change=paste0(Change, " (d", DPI, ")")) %>%
            # compile by NHP
            group_by(NHP) %>%
            summarise(Temperature=paste(Change, collapse="; "),
                      .groups="drop") %>%
            mutate(Temperature=tolower(Temperature))

# appetite & anorexia
clinical <- df$observations %>%
            group_by(NHP) %>%
            mutate(MaxB=max(`Biscuits eaten`),
                   FoodConsumed=`Biscuits eaten`/MaxB) %>%
            ungroup() %>%
            # decreased appetite <= 60% food consumed
            # anorexia = 0% food consumed 
            mutate(`decreased appetite`=(FoodConsumed <= 0.6),
                   anorexia=(FoodConsumed==0)) %>%
            select(NHP, DPI, `decreased appetite`, anorexia) %>%
            # melt to long
            reshape2::melt(id.vars=c("NHP", "DPI"),
                           variable.name="Appetite",
                           value.name="Present") %>%
            filter(Present) %>%
            mutate(Appetite=factor(Appetite, 
                                   levels=c("decreased appetite", 
                                            "anorexia"))) %>%
            # remove decreased appetite when anorexia is present
            group_by(NHP, DPI) %>%
            top_n(n=1, wt=Appetite) %>%
            ungroup() %>%
            # compile by appetite label 
            group_by(NHP, Appetite) %>%
            summarise(DPI=paste(DPI, collapse=", "),
                      .groups="drop") %>%
            mutate(Appetite=paste0(Appetite, " (d", DPI, ")")) %>%
            # compile by NHP
            group_by(NHP) %>%
            summarise(Appetite=paste(Appetite, collapse="; "),
                      .groups="drop") %>%
            # add to clinical
            full_join(clinical, by="NHP")

# observations
clinical <- df$observations %>%
            select(-`Biscuits eaten`) %>%
            # melt to long format
            reshape2::melt(id.vars=c("NHP", "DPI"),
                           variable.name="Observation",
                           value.name="Present") %>%
            filter(!is.na(Present)) %>%
            select(-Present) %>% 
            # compile by observation type
            group_by(NHP, Observation) %>%
            summarise(DPI=paste(DPI, collapse=", "),
                      .groups="drop") %>%
            mutate(Observation=paste0(tolower(Observation), 
                                      " (d", DPI, ")")) %>%
            # compile by NHP
            group_by(NHP) %>%
            summarise(Observation=paste(Observation, collapse="; "),
                      .groups="drop") %>%
            # add to clinical illness
            full_join(clinical, by="NHP")

# survival
endpoint <- df$animals %>%
            select(NHP, `Death DPI`, Outcome) %>%
            # format outcome for human readability
            mutate(Outcome=factor(Outcome, 
                                  levels=c("Succumbed", 
                                           "Survived", 
                                           "Sacrificed"),
                                  labels=c("Succumbed on d", 
                                           "Survived to study endpoint (d",
                                           "Scheduled euthanasia at d")),
                   Outcome=paste0(Outcome, `Death DPI`)) 
# add the last parenthesis for "survived"
x <- str_detect(endpoint$Outcome, "Survived")
endpoint$Outcome[x] <- paste0(endpoint$Outcome[x], ")")
# add to clinical illness
clinical <- endpoint %>%
            select(NHP, Outcome) %>%
            mutate(Outcome=paste0(Outcome, ".")) %>%
            full_join(clinical, by="NHP")

# assemble the clinical illness column
clinical <- clinical %>%
            # compile by NHP. Keep outcome separate for now.
            reshape2::melt(id.vars=c("NHP", "Outcome"),
                           variable.name="Metric",
                           value.name="Clinical")
illness <- clinical %>%
           # remove metrics with no values
           na.omit() %>%
           # compile by NHP and outcome
           group_by(NHP, Outcome) %>%
           summarise(Clinical=paste(Clinical, collapse="; "),
                     .groups="drop") %>% 
           mutate(Clinical=str_to_sentence(Clinical))
# add formatted clinical illness column back to the data frame. 
clinical <- clinical %>%
            group_by(NHP, Outcome) %>%
            summarize(.groups="drop") %>%
            left_join(illness, by=c("NHP", "Outcome")) %>%
            # fill in missing values (NHPs with no clinical illness)
            replace_na(list(Clinical="None")) %>%
            # compile by NHP
            group_by(NHP) %>%
              summarise(`Clinical illness`=paste0(Clinical, ". ", Outcome),
                        .groups="drop")
# clean up
rm(baseline, x, endpoint, illness)

## clinical pathology ----------------------------------------------------------
# hematology
baseline <- df$hematology %>%
            filter(DPI==0) %>%
            select(NHP, Leukocytes, Lymphocytes, Monocytes, Erythrocytes,
                   Platelets, Neutrophils, Eosinophils, Basophils) %>%
            reshape2::melt(id.vars="NHP",
                           variable.name="Celltype",
                           value.name="Baseline")
clinpath <- df$hematology %>%
            reshape2::melt(id.vars=c("NHP", "DPI"),
                           variable.name="Celltype",
                           value.name="Concentration") %>%
            filter(Celltype %in% unique(baseline$Celltype)) %>%
            # add in baseline and calculate fraction
            # "-penia" = 0.65x baseline or lower
            # "-tosis" = 2x baseline or higher
            left_join(baseline, by=c("NHP", "Celltype")) %>%
            mutate(Fraction=(Concentration/Baseline)) %>%
            mutate(Penia=(Fraction <= 0.65),
                   Tosis=(Fraction >= 2)) %>%
            filter(Penia | Tosis) %>%
            select(NHP, DPI, Celltype, Penia, Tosis) %>%
            # reshape again and modulate to get names
            reshape2::melt(id.vars=c("NHP", "DPI", "Celltype"),
                           variable.name="Change",
                           value.name="Present") %>%
            filter(Present) %>%
            # smash the names together and re-level
            mutate(Change=factor(paste0(Celltype, Change),
                                 levels=c("LeukocytesPenia", 
                                          "LymphocytesPenia",
                                          "MonocytesPenia",
                                          "NeutrophilsPenia",
                                          "EosinophilsPenia",
                                          "BasophilsPenia",
                                          "ErythrocytesPenia",
                                          "PlateletsPenia",
                                          "LeukocytesTosis",
                                          "LymphocytesTosis",
                                          "MonocytesTosis",
                                          "ErythrocytesTosis",
                                          "PlateletsTosis",
                                          "NeutrophilsTosis",
                                          "EosinophilsTosis",
                                          "BasophilsTosis"),
                                 labels=c("leukopenia", 
                                          "lymphopenia",
                                          "monocytopenia",
                                          "neutropenia",
                                          "eosinopenia",
                                          "basopenia",
                                          "erythrocytopenia",
                                          "thrombocytopenia",
                                          "leukocytosis",
                                          "lymphocytosis",
                                          "monocytosis",
                                          "erythrocytosis",
                                          "thrombocytosis",
                                          "neutrophilia",
                                          "eosinophilia",
                                          "basophilia"))) %>%
            # compile by NHP and type of change
            group_by(NHP, Change) %>% 
            summarise(DPI=paste(DPI, collapse=", "), 
                      .groups="drop") %>%
            mutate(Change=paste0(Change, " (d", DPI, ")")) %>%
            group_by(NHP) %>%
            summarise(Hematology=paste(Change, collapse="; "),
                      .groups="drop")

# chemistry
baseline <- df$chemistry %>%
            left_join(df$hematology,
                      by=c("NHP", "DPI")) %>%
            filter(DPI==0) %>%
            select(NHP, `Glucose (GLU)`, Hematocrit, Hemoglobin, Erythrocytes,
                   `Albumin (ALB)`, `Total protein (TP)`, `Amylase (AMY)`, 
                   `Calcium (CA)`, `Alanine aminotransferase (ALT)`, 
                   `Aspartate transferase (AST)`, `Alkaline phosphatase (ALP)`, 
                   `Creatinine (CRE)`, `C-reactive protein (CRP)`) %>%
            reshape2::melt(id.vars="NHP",
                           variable.name="Analyte",
                           value.name="Baseline")
# make a temporary DF to keep for fold-change/arrows calculations
chem <- df$chemistry %>%
        left_join(df$hematology,
                  by=c("NHP", "DPI")) %>%
        select(NHP, DPI, unique(baseline$Analyte)) %>%
        # melt for baseline calculations
        reshape2::melt(id.vars=c("NHP", "DPI"),
                       variable.name="Analyte",
                       value.name="Concentration") %>%
        left_join(baseline, by=c("NHP", "Analyte")) %>%
        mutate(Change=(Concentration/Baseline)) %>%
        select(NHP, DPI, Analyte, Change) %>%
        # expand for terminology calculations
        reshape2::dcast(NHP + DPI ~ Analyte, value.var="Change")
clinpath <- chem %>%
            mutate(Hyperglycemia=(`Glucose (GLU)` >= 2),
                   Hypoglycemia=(`Glucose (GLU)` <= 0.75),
                   Hypoalbuminemia=(`Albumin (ALB)` <= 0.75),
                   Hypoproteinemia=(`Total protein (TP)` <= 0.75),
                   Hyperamylasemia=(`Amylase (AMY)` >= 2),
                   Hypoamylasemia=(`Amylase (AMY)` <= 0.75),
                   Hypocalcemia=(`Calcium (CA)` <= 0.75),
                   Erythrocytes=(Erythrocytes <= 0.7),
                   Hematocrit=(Hematocrit <= 0.7),
                   Hemoglobin=(Hemoglobin <= 0.7),
                   Anemia=(Erythrocytes & Hematocrit & Hemoglobin)) %>%
            select(NHP, DPI, Anemia, 
                   starts_with("Hyper"), 
                   starts_with("Hypo")) %>%
            # melt for filtering and formatting
            reshape2::melt(id.vars=c("NHP", "DPI"),
                           variable.name="Analyte",
                           value.name="Present") %>%
            filter(Present) %>%
            # compile by NHP and analyte
            group_by(NHP, Analyte) %>%
            summarise(DPI=paste(DPI, collapse=", "),
                      .groups="drop") %>%
            mutate(Analyte=paste0(tolower(Analyte), " (d", DPI, ")")) %>%
            # compile by NHP
            group_by(NHP) %>%
            summarise(Chemistry=paste(Analyte, collapse="; "),
                      .groups="drop") %>%
            # add to clinical pathology
            full_join(clinpath, by="NHP")

# arrows. Use temporary chemistry data as a starting point
# ↑ is [2-5), ↑↑ is [5-10), ↑↑↑ is [10-20), ↑↑↑↑ is >= 20, and ↓ is <= 0.5
chem <- chem %>%
        select(NHP, DPI, `Alanine aminotransferase (ALT)`, 
               `Aspartate transferase (AST)`, `Albumin (ALB)`, 
               `Creatinine (CRE)`, `C-reactive protein (CRP)`, 
               Hematocrit, Hemoglobin) %>%
        # transform to long
        reshape2::melt(id.vars=c("NHP", "DPI"),
                       variable.name="Analyte",
                       value.name="Change")
chem$Arrows[chem$Change >= 2] <- "↑"
chem$Arrows[chem$Change >= 5] <- "↑↑"
chem$Arrows[chem$Change >= 10] <- "↑↑↑"
chem$Arrows[chem$Change >= 20] <- "↑↑↑↑"
chem$Arrows[chem$Change <= 0.5] <- "↓"
# format and add to clinical pathology
clinpath <- chem %>%
            na.omit() %>%
            mutate(Analyte=str_extract(Analyte, "(?<=\\()[A-Z]+(?=\\))")) %>%
            # first compile by number of arrows
            group_by(NHP, Analyte, Arrows) %>%
            summarise(DPI=paste(DPI, collapse=", "),
                      .groups="drop") %>%
            # combine with arrows
            mutate(Arrows=paste0(Arrows, " (d", DPI, ")")) %>%
            # next compile by analyte, since there can be multiple 
            # arrows per analyte
            group_by(NHP, Analyte) %>%
            summarise(Arrows=paste(Arrows, collapse=", "),
                      .groups="drop") %>%
            mutate(Analyte=paste(Analyte, Arrows)) %>%
            # last, compile by NHP
            group_by(NHP) %>%
            summarise(Arrows=paste(Analyte, collapse="; "),
                      .groups="drop") %>%
            # add to the rest of the data
            full_join(clinpath, by="NHP")

# assemble the clinical pathology column
clinpath <- clinpath %>%
            # melt to long format and remove NA
            reshape2::melt(id.vars="NHP",
                           variable.name="Metric",
                           value.name="clinpath") %>%
            na.omit() %>%
            mutate(Metric=factor(Metric, 
                                 levels=c("Hematology", 
                                          "Chemistry", 
                                          "Arrows"))) %>%
            arrange(Metric) %>% 
            # compile by NHP
            group_by(NHP) %>%
            summarise(`Clinical pathology`=paste(clinpath, collapse="; "),
                      .groups="drop") %>%
            mutate(`Clinical pathology`=DescTools::StrCap(`Clinical pathology`))

# clean up
rm(baseline, chem)

## assemble the full table and save --------------------------------------------
df$animals %>%
  mutate(Sex=factor(Sex, levels=c("Female", "Male"), labels=c("F", "M"))) %>%
  select(NHP, Sex, `Dose (PFU)`, `Presented dose (DP)`) %>%
  # add clinical illness
  left_join(clinical, by="NHP") %>%
  # add clinical pathology
  left_join(clinpath, by="NHP") %>%
  # arrange by dose
  arrange(desc(`Dose (PFU)`)) %>%
  openxlsx::write.xlsx("analysis/clinicaltables.xlsx")

## fin -------------------------------------------------------------------------
sessionInfo()
