# 16S rRNA analysis for Vitor's column experiment

# --- ENVIRONMENT SETUP ---
# Ensure a writable user-specific library exists and is added to .libPaths()
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "" || !dir.exists(user_lib)) {
  user_lib <- file.path(Sys.getenv("USERPROFILE"), "Documents", "R", "win-library", "4.5")
}
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(user_lib, .libPaths()))

# Set CRAN mirror to avoid interactive prompt
options(repos = c(CRAN = "https://cloud.r-project.org"))

# --- PACKAGE INSTALLATION ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Fix rbibutils if needed (for FSA)
if (requireNamespace("rbibutils", quietly = TRUE)) {
  if (packageVersion("rbibutils") < "2.4.1") {
    BiocManager::install("rbibutils", update = TRUE, ask = FALSE, force = TRUE)
  }
} else {
  BiocManager::install("rbibutils", update = FALSE, ask = FALSE)
}

# Install microViz specially
if (!requireNamespace("microViz", quietly = TRUE)) {
  install.packages("microViz", repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos")))
}

# Helper to install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }, error = function(e) {
      message(paste("Failed to install", pkg, ":", e$message))
    })
  }
}

pkgs <- c("DECIPHER", "phytools", "lubridate", "microeco", "phangorn", "DESeq2", 
          "FSA", "rcompanion", "ggsci", "microbiome", "phyloseq", "picante", 
          "vegan", "ape", "nlme", "ggplot2", "reshape2", "tidyverse", 
          "RColorBrewer", "patchwork", "file2meco", "GUniFrac", "hillR", 
          "ggridges", "tibble", "scales", "decontam", "ggpubr", "knitr", "dplyr")

for (pkg in pkgs) {
  install_if_missing(pkg)
}

# --- LOAD LIBRARIES ---
library(vegan)
library(picante)
library(ape)
library(nlme)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(phytools)
library(phyloseq)
library(RColorBrewer)
library(microbiome)
library(microeco)
library(microViz)
library(DECIPHER)
library(phangorn)
library(decontam)
library(patchwork)
library(file2meco)
library(DESeq2)
library(GUniFrac)
library(hillR)
library(FSA)
library(rcompanion)
library(dplyr)
library(ggridges)
library(tibble)
library(scales)
library(knitr)

# --- THEME DEFINITION (Makie-inspired) ---
# colors = [:crimson, :steelblue, :forestgreen, :darkorange, :darkgrey]
makie_colors <- c(
  "Column 1" = "#DC143C", # Crimson hex
  "Column 2" = "#4682B4", # Steelblue hex
  "Column 3" = "#228B22", # Forestgreen hex
  "Column 4" = "#FF8C00", # Darkorange hex
  "Control"  = "#A9A9A9"  # Darkgrey hex
)

theme_makie <- function() {
  theme_classic(base_size = 11) +
    theme(
      text = element_text(family = "Helvetica"),
      plot.title = element_text(size = 10, face = "plain", hjust = 0),
      plot.title.position = "plot",
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      axis.ticks.length = unit(3, "pt"),
      axis.line = element_line(linewidth = 1.2),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 9, face = "bold"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 9),
      legend.key.size = unit(10, "pt"),
      legend.margin = margin(2, 2, 2, 2)
    )
}



####LOAD DATA####
# asv_table (16S rDNA marker gene count data)
asv_table <- read.delim("data/raw_lab_data/biomolecular/DADA2_counts_as_matrix.tsv",
                        check.names = FALSE, # so it does not change "-" to "." in JMF names 
                        row.names=1, # make rownames out of the first column 
                        header=TRUE, # column names are defined
                        sep="\t") 

# tax_table (taxonomy of each OTU)
# will need to be further cleaned later on
tax_table <- read.delim("data/raw_lab_data/biomolecular/DADA2_ASVs.rRNA_SSU.SILVA_reference.DADA2_classified.tsv", 
                        header=TRUE, 
                        sep=",")

# metadata (data where the first column has the same identifier as the column names of asv_table - sample names)
metadata_16S <- read.table("data/raw_lab_data/biomolecular/vitor_kassel_metadata.csv", 
                           row.names=1, # make rownames out of the first column 
                           header=TRUE, 
                           sep=",",
                           #  your decimal is likely written as "," instead of the usual "." 
                           # if your laptop language is german
                           # if that is not the case then change the dec = "."
                           dec = ".") 

#####CLEAN DATA####
#clean asv_table
# remove unnecessary suffix
names(asv_table) <- sub("A", "", names(asv_table))

# clean taxonomy
tax_table <- tax_table %>%
  tidyr::separate_wider_delim(
    `Sequence_ID.Domain.Phylum.Class.Order.Family.Genus.Species`,
    delim = "\t",
    names = c(
      "Sequence_ID", "Domain", "Phylum", "Class", "Order",
      "Family", "Genus", "Species"
    ),
    too_few = "align_start",
    too_many = "drop"
  ) %>%
  tibble::column_to_rownames("Sequence_ID")

# with the separate_wider_column the first NA in line end up just empty
## so make that an NA as well:
tax_table[tax_table == ''] <- NA

# some ASVs are classified as "uncultured" at Genus level
# this does not help us 
### use this: tax_table[tax_table == 'uncultured'] <- NA
### use this: tax_table[tax_table == 'Incertae Sedis'] <- NA

## check first columns of the dataset
asv_table[1:10,7:11] %>%  knitr::kable()
tax_table[1:5,1:5] %>%  knitr::kable()
metadata_16S[2:2,1:2] %>%  knitr::kable()

# inspect all of our columns
names(metadata_16S)

# make sure numeric variables are numeric
### not possible since there's no chem data: metadata_16S[,5:14] <- as.numeric(unlist(metadata_16S[,5:14]))

# make sure categorical variables are categorical
metadata_16S %>% mutate_at(c("site",
                             "sample_type") , as.factor ) -> metadata_16S
# check
str(metadata_16S[1:2, 1:1]) %>%  knitr::kable()

#Final cleaning before proceeding with analyses - make sure that row names and column names are present and identical in all three datasets
identical(rownames(tax_table),rownames(asv_table))
identical(colnames(asv_table), rownames(metadata_16S))

####DECONTAMINATION####
# first inspect the number of taxa before decontamination
nrow(asv_table) # 2853 taxa present

# make vector for blanks
# i.e. an object containing information whether sample is a blank
# if it is blank it will say TRUE, otherwise it will be FALSE
vector_for_decontam <- print(metadata_16S$well_id == "blank") 

# identify blanks
# create TRUE/FALSE vector identifying blanks
vector_for_decontam <- metadata_16S$sample_type == "blank"

# run decontam without batch column
contam_df <- isContaminant(
  t(asv_table),
  neg = vector_for_decontam,
  method = "prevalence",
  threshold = 0.1
)

# to see how many are contaminants
table(contam_df$contaminant) 

# identified 1  as contaminant, find the names
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ]) %>% as.data.frame()
contam_asvs[,1]

# identify which taxa are contaminants
# connect names to the taxonomy table to identify them 
cont_identified <- tax_table %>% 
  rownames_to_column(".") %>% 
  # left_join/right_join/full_join are popular ways of gluing datasets together
  # for details try the help manual
  left_join (contam_asvs, ., by=".") 

# to see which are the most frequent contaminants arrange the data frame from previous command
# in this case not needed since there is only one contaminant (it is Deinococci)
# but in your project you could have many
table(cont_identified$Class) %>% 
  as.data.frame() %>%  
  arrange(desc(Freq)) %>% 
  knitr::kable()

# remove contaminants asvs from the asv_table
# the ! symbol means "not the same/different than"
asv_table_no_cont <- asv_table[!row.names(asv_table) %in% contam_asvs$., ]


# removes blank samples from the asv_table
asv_table_no_cont <- asv_table_no_cont[ , !vector_for_decontam]

# removes contaminants from taxonomy table
tax_table_no_cont <- tax_table[!row.names(tax_table) %in% contam_asvs$., ]

# final check - have we removed the contaminant taxa
nrow(asv_table_no_cont)



####OBTAIN PHYLOSEQ####
# Filter metadata to remove blanks (to match the filtered asv_table_no_cont)
metadata_16S_no_blanks <- metadata_16S[!vector_for_decontam, ]

# Now create phyloseq with matching filtered data
ps_16S_sediments <- phyloseq(
  otu_table(as.matrix(asv_table_no_cont), taxa_are_rows=T),
  tax_table(as.matrix(tax_table_no_cont[,-7])), 
  sample_data(metadata_16S_no_blanks)  # Use filtered metadata
) %>% 
  tax_fix()

# Check - should now have 198 samples and 8134 taxa
ps_16S_sediments

####CHECK PHYLOSEQ OBJECT####
# For prevalence checking, use UNFILTERED data (before removing contaminants/blanks)
ps_for_prev <- phyloseq(
  otu_table(as.matrix(asv_table), taxa_are_rows=T),
  tax_table(as.matrix(tax_table[,-7])), 
  sample_data(metadata_16S)  # Original metadata with blanks
)

# Transform into presence-absence
ps_pa <- transform_sample_counts(ps_for_prev, 
                                 function(abund) 1*(abund>0))

# Separate blanks and samples
ps_pa_neg <- prune_samples(sample_data(ps_pa)$sample_type == "blank", ps_pa)
ps_pa_pos <- prune_samples(sample_data(ps_pa)$sample_type != "blank", ps_pa)  # FIXED

# Create prevalence dataframe
df_pa <- data.frame(
  pa_pos = taxa_sums(ps_pa_pos), 
  pa_neg = taxa_sums(ps_pa_neg),
  contaminant = contam_df$contaminant
)

# Plot
ggplot(data = df_pa,
       aes(x = pa_neg, y = pa_pos, 
           color = contaminant, size = contaminant)) + 
  geom_point(alpha = 0.4) +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)") +
  theme_bw() + 
  theme(legend.position = c(0.8, 0.8))

# Check - should now have 198 samples and 8134 taxa
ps_16S_sediments
ps_for_prev 


####FILTERING DATASET####
ps_16S_sediments <-  subset_taxa(ps_16S_sediments, Domain %in% c("Bacteria", "Archaea"),
                                 !is.na(Phylum) & 
                                   !Phylum %in% c("", "uncharacterized", 
                                                  "Bacteria Kingdom", 
                                                  "Archaea Kingdom") & 
                                   !Order %in% c("Chloroplast", "Incertae Sedis") &
                                   !Family %in% c("Mitochondria"))

#how many taxa we have now 
ps_16S_sediments


####Dealing with uneven sequencing counts####
# look at numbers
as.data.frame(sample_sums(ps_16S_sediments) ) %>% summary()

# plotting
as.data.frame(sample_sums(ps_16S_sediments)) %>% 
  mutate(site = as.data.frame(sample_data(ps_16S_sediments))$site) %>% 
  ggplot(aes(x = sample_sums(ps_16S_sediments), y = site, fill = site)) + 
  
  # density ridges (removed bins parameter)
  ggridges::geom_density_ridges(alpha = 0.7, 
                                position = 'identity', 
                                quantile_lines = TRUE) + 
  
  # add title with mean and sd
  ggtitle(paste0("Total number of DNA reads", 
                 " (",
                 as.data.frame(sample_sums(ps_16S_sediments))[,1] %>% mean() %>% round(., 2),
                 " ± ",
                 as.data.frame(sample_sums(ps_16S_sediments))[,1] %>% sd() %>% round(., 2),
                 ")")) +
  
  labs(y = "Frequency of sequence count", 
       x = "Sequence count in sample", 
       fill = "Site") +
  
  # Either remove scale_fill_manual() to use default colors:
  # (delete the line below)
  
  # OR provide specific colors for your 8 sites:
  scale_fill_brewer(palette = "Set2") +  # Use this instead
  
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        ggside.panel.scale.x = 0.3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0))


####RAREFRACTION####
# 1st ps object
# even though we are not rarefying we will remove the samples that have an extremely low read count
ps_16S_sediments

# when removing samples with less than 20 count 
ps_16S_sediments <- prune_samples(sample_sums(ps_16S_sediments) >= 20, ps_16S_sediments) 

# 2nd ps object
ps_16S_sediments_rar <- rarefy_even_depth(ps_16S_sediments, sample.size = 1500,
                                          rngseed = 123, replace = TRUE, 
                                          trimOTUs = TRUE, verbose = TRUE)

# 3rd ps object
ps_16S_sediments_rel <- transform_sample_counts(ps_16S_sediments, function(x){x / sum(x)})

ps_16S_sediments
ps_16S_sediments_rar


####Downstream analysis####
# look at all the samples of our dataset
ps_16S_sediments %>%
  comp_barplot( tax_level = "Phylum", # choose the taxonomic level you want
                n_taxa = 15,          # choose number of taxa you want to show
                facet_by = "site", # group by well id
                other_name = "Other phyla", # set custom name for the "other" category
                merge_other = FALSE, # split the "Other" category to display alpha diversity
                bar_width = 0.9, # reduce the bar width to 70% of one row
                bar_outline_colour = "grey5" ) + 
  theme_minimal() + 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        panel.grid = element_blank())  +
  scale_y_continuous(expand = c(0,0),
                     labels = function(x) x * 100) +
  labs(y = "Relative abundance (%)") + 
  guides(fill = guide_legend(reverse = TRUE)) -> plot_1.1

# once you give it a name (plot_1.1) you can call on this plot whenever you want
plot_1.1



ps_16S_sediments_rar %>%
  comp_barplot( tax_level = "Phylum", # choose the taxonomic level you want
                n_taxa = 15,          # choose number of taxa you want to show
                facet_by = "site", # group by well id
                other_name = "Other phyla", # set custom name for the "other" category
                merge_other = FALSE, # split the "Other" category to display alpha diversity
                bar_width = 0.9, # reduce the bar width to 70% of one row
                bar_outline_colour = "grey5" ) + 
  theme_minimal() + 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        panel.grid = element_blank())  +
  scale_y_continuous(expand = c(0,0),
                     labels = function(x) x * 100) +
  labs(y = "Relative abundance (%)") + 
  guides(fill = guide_legend(reverse = TRUE)) -> plot_1.2


plot_1.2


#########################################
### A1. FILTER FOR Fuhrberg columns ONLY
#########################################
# Keep only Fuhrberg (drops sonja + anything else)
ps_fuhrberg <- prune_samples(sample_data(ps_16S_sediments)$site == "fuhrberg",
                             ps_16S_sediments)
sample_sums(ps_fuhrberg)
# Total ASVs present (non-zero after subsetting)
ntaxa(ps_fuhrberg)

read_counts <- sample_sums(ps_fuhrberg)

mean_reads <- mean(read_counts)
sd_reads   <- sd(read_counts)
min_reads  <- min(read_counts)
max_reads  <- max(read_counts)

mean_reads
sd_reads
min_reads
max_reads

## --- Fuhrberg ASV abundance table (counts) ---
asv_fuhrberg_mat <- as(otu_table(ps_fuhrberg), "matrix")
if (!taxa_are_rows(ps_fuhrberg)) asv_fuhrberg_mat <- t(asv_fuhrberg_mat)

# dimensions
dim(asv_fuhrberg_mat)          # ASVs x samples
n_asv_fuhrberg <- nrow(asv_fuhrberg_mat)

# total reads per sample
reads_fuhrberg <- colSums(asv_fuhrberg_mat)
summary(reads_fuhrberg)

# ASVs present across Fuhrberg (non-zero across all samples pooled)
n_asv_present <- sum(rowSums(asv_fuhrberg_mat) > 0)

n_asv_fuhrberg
n_asv_present

#OTU for family level
ps_fuhrberg_family <- tax_glom(ps_fuhrberg, taxrank = "Family", NArm = TRUE)

otu_family_mat <- as(otu_table(ps_fuhrberg_family), "matrix")
if (!taxa_are_rows(ps_fuhrberg_family)) otu_family_mat <- t(otu_family_mat)

dim(otu_family_mat)   # families x samples
sum(rowSums(otu_family_mat) > 0)


# (optional) drop taxa that become all-zero after subsetting
ps_fuhrberg <- prune_taxa(taxa_sums(ps_fuhrberg) > 0, ps_fuhrberg)

ps_fuhrberg

# 1) filter first
ps_fuhrberg <- prune_samples(sample_data(ps_16S_sediments)$site == "fuhrberg",
                             ps_16S_sediments)
ps_fuhrberg <- prune_taxa(taxa_sums(ps_fuhrberg) > 0, ps_fuhrberg)

# 2) rarefy (choose sample.size based on Fuhrberg min depth)
min(sample_sums(ps_fuhrberg))
ps_fuhrberg_rar <- rarefy_even_depth(ps_fuhrberg, sample.size = 3000,
                                     rngseed = 123, replace = TRUE,
                                     trimOTUs = TRUE, verbose = TRUE)

# create relative-abundance object
ps_fuhrberg_rel <- transform_sample_counts(
  ps_fuhrberg,
  function(x) x / sum(x)
)
# quick check
range(otu_table(ps_fuhrberg_rel))


# 3) relative abundance
plot_fuhrberg <- ps_fuhrberg %>%
  comp_barplot(
    tax_level = "Phylum",
    n_taxa = 15,
    facet_by = "site",
    other_name = "Other phyla",
    merge_other = FALSE,
    bar_width = 0.9,
    bar_outline_colour = "grey5"
  ) +
  theme_minimal() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0, 0),
                     labels = function(x) x * 100) +
  labs(y = "Relative abundance (%)") +
  guides(fill = guide_legend(reverse = TRUE))

plot_fuhrberg

#Quick sanity check (so you don’t accidentally use rel-abund data for count-required steps)
# should be TRUE for counts
all(otu_table(ps_fuhrberg) == round(otu_table(ps_fuhrberg)))

# should be FALSE for rel abund
all(otu_table(ps_fuhrberg_rel) == round(otu_table(ps_fuhrberg_rel)))


## =========================
## 0) Add labels to samples (Fuhrberg only)
## =========================
sd <- data.frame(sample_data(ps_fuhrberg))
sd$sample_id <- rownames(sd)

sd <- sd %>%
  mutate(sample_label = case_when(
    sample_id == "JMF-2511-04-0034" ~ "Control – Inlet (0-4cm)",
    sample_id == "JMF-2511-04-0039" ~ "Control – Outlet (4-8cm)",

    sample_id == "JMF-2511-04-0031" ~ "Column 1 – Inlet (0-4cm)",
    sample_id == "JMF-2511-04-0032" ~ "Column 2 – Inlet (0-4cm)",
    sample_id == "JMF-2511-04-0033" ~ "Column 3 – Inlet (0-4cm)",

    sample_id == "JMF-2511-04-0036" ~ "Column 1 – Outlet (4-8cm)",
    sample_id == "JMF-2511-04-0037" ~ "Column 2 – Outlet (4-8cm)",
    sample_id == "JMF-2511-04-0038" ~ "Column 3 – Outlet (4-8cm)",  # <-- change this to the correct ID for Column 3 outlet

    TRUE ~ as.character(sd$sample_type)
  ))
sample_data(ps_fuhrberg)$sample_label <- sd$sample_label

sample_data(ps_fuhrberg)$sample_label <- factor(
  sample_data(ps_fuhrberg)$sample_label,
  levels = c(
    "Control – Inlet (0-4cm)", "Control – Outlet (4-8cm)",
    "Column 1 – Inlet (0-4cm)", "Column 1 – Outlet (4-8cm)",
    "Column 2 – Inlet (0-4cm)", "Column 2 – Outlet (4-8cm)",
    "Column 3 – Inlet (0-4cm)", "Column 3 – Outlet (4-8cm)"
  )
)

# Optional helper: Inlet vs Outlet grouping
sample_data(ps_fuhrberg)$position <- ifelse(
  grepl("Inlet", as.character(sample_data(ps_fuhrberg)$sample_label)),
  "Inlet", "Outlet"
)
sample_data(ps_fuhrberg)$position <- factor(sample_data(ps_fuhrberg)$position, levels = c("Inlet", "Outlet"))

## =========================
## 1) Phylum composition (counts preserved; percent axis)
## =========================
# Uses comp_barplot(); if your comp_barplot has no `x=` argument, see fallback below.
pA <- ps_fuhrberg %>%
  comp_barplot(
    tax_level = "Phylum",
    n_taxa = 15,
    x = "sample_label",
    other_name = "Other phyla",
    merge_other = FALSE,
    bar_width = 0.9,
    bar_outline_colour = "grey20"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Relative abundance (%)", title = "Community composition (Phylum)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  ) +
  guides(fill = guide_legend(reverse = TRUE))

pA


## =========================
## Family composition (same order as heatmap, no Mitochondria)
## =========================

# enforce the exact order used in the heatmap
sample_data(ps_fuhrberg)$sample_label <- factor(
  sample_data(ps_fuhrberg)$sample_label,
  levels = c(
    "Control – Inlet (0-4cm)", "Column 1 – Inlet (0-4cm)", "Column 2 – Inlet (0-4cm)", "Column 3 – Inlet (0-4cm)",
    "Control – Outlet (4-8cm)", "Column 1 – Outlet (4-8cm)", "Column 2 – Outlet (4-8cm)", "Column 3 – Outlet (4-8cm)"
  )
)

# remove unwanted taxa first
ps_fuhrberg_no_mito <- subset_taxa(
  ps_fuhrberg,
  Domain == "Bacteria" &
    !Order %in% c("Chloroplast") &
    !Family %in% c(
      "Bacteria Domain", "Sericytochromatia Class", "Kapabacteriales Order",
      "mle1-27 Order", "053A03-B-DI-P58 Order", "Lineage IV Order"
    )
)

# make sure sample metadata has the needed variables
sd <- data.frame(sample_data(ps_fuhrberg_no_mito))
sd$sample_id <- rownames(sd)

# derivation of column and position from sample_label
sd$column <- dplyr::case_when(
  grepl("^Control",  sd$sample_label) ~ "Control",
  grepl("^Column 1", sd$sample_label) ~ "Column 1",
  grepl("^Column 2", sd$sample_label) ~ "Column 2",
  grepl("^Column 3", sd$sample_label) ~ "Column 3",
  TRUE ~ NA_character_
)

sd$position <- ifelse(grepl("Inlet", sd$sample_label), "Inlet (0-4cm)", "Outlet (4-8cm)")

sd$column <- factor(sd$column, levels = c("Control", "Column 1", "Column 2", "Column 3"))
sd$position <- factor(sd$position, levels = c("Inlet (0-4cm)", "Outlet (4-8cm)"))

sample_data(ps_fuhrberg_no_mito)$column <- sd$column
sample_data(ps_fuhrberg_no_mito)$position <- sd$position
sample_data(ps_fuhrberg_no_mito)$sample_label <- factor(
  sd$sample_label,
  levels = c(
    "Control – Inlet (0-4cm)", "Column 1 – Inlet (0-4cm)", "Column 2 – Inlet (0-4cm)", "Column 3 – Inlet (0-4cm)",
    "Control – Outlet (4-8cm)", "Column 1 – Outlet (4-8cm)", "Column 2 – Outlet (4-8cm)", "Column 3 – Outlet (4-8cm)"
  )
)

# reorder samples in the phyloseq object based on sample_label
sample_order <- rownames(data.frame(sample_data(ps_fuhrberg_no_mito)))[
  order(sample_data(ps_fuhrberg_no_mito)$sample_label)
]

ps_fuhrberg_ordered <- prune_samples(sample_order, ps_fuhrberg_no_mito)

# quick check of valid sample variables
sample_variables(ps_fuhrberg_ordered)

# family barplot with inlet/outlet separation
pB <- ps_fuhrberg_ordered %>%
  comp_barplot(
    tax_level = "Family",
    n_taxa = 30,
    x = "column",
    facet_by = "position",
    other_name = "Other family",
    merge_other = FALSE,
    bar_width = 0.9,
    bar_outline_colour = "grey20"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = NULL,
    y = "Relative abundance (%)",
    title = "a. Family Composition"
  ) +
  theme_makie() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.text = element_text(size = 7, face = "italic"),
    legend.key.size = unit(8, "pt")
  ) +
  guides(fill = guide_legend(reverse = TRUE, ncol = 1))

pB


## =========================
## Mean and SD of Family relative abundance
## =========================

# Agglomerate to Family level
ps_fuhrberg_family_rel <- tax_glom(ps_fuhrberg_rel, taxrank = "Family", NArm = TRUE)

# Convert to long format
fam_rel_df <- psmelt(ps_fuhrberg_family_rel) %>%
  filter(!is.na(Family))

# Calculate mean and SD across all samples
fam_mean_sd <- fam_rel_df %>%
  group_by(Family) %>%
  summarise(
    mean_relative_abundance = mean(Abundance, na.rm = TRUE),
    sd_relative_abundance   = sd(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mean_percent = mean_relative_abundance * 100,
    sd_percent   = sd_relative_abundance * 100
  ) %>%
  arrange(desc(mean_percent))

# View top families
head(fam_mean_sd, 30)

############## Heatmap #############

stopifnot("sample_label" %in% sample_variables(ps_fuhrberg))

# Inlet vs Outlet grouping (optional but matches your request)
sample_data(ps_fuhrberg)$io_group <- ifelse(
  grepl("Inlet", as.character(sample_data(ps_fuhrberg)$sample_label)),
  "Inlet (0-4cm)", "Outlet (4-8cm)"
)
sample_data(ps_fuhrberg)$io_group <- factor(sample_data(ps_fuhrberg)$io_group, levels = c("Inlet (0-4cm)", "Outlet (4-8cm)"))

# 1) Agglomerate to Family
ps_fam <- tax_glom(ps_fuhrberg, taxrank = "Family", NArm = TRUE)

# 2) Relative abundance
ps_fam_rel <- transform_sample_counts(ps_fam, function(x) x / sum(x))

# 3) Long format + attach sample_label + io_group
df <- psmelt(ps_fam_rel) %>%
  mutate(
    sample_label = sample_data(ps_fuhrberg)$sample_label[match(Sample, sample_names(ps_fuhrberg))],
    io_group     = sample_data(ps_fuhrberg)$io_group[match(Sample, sample_names(ps_fuhrberg))]
  ) %>%
  filter(Family != "Mitochondria")

# 4) Top N families + collapse others
topN <- 40
top_fams <- df %>%
  group_by(Family) %>%
  summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = topN) %>%
  pull(Family)

df2 <- df %>%
  mutate(Family2 = ifelse(Family %in% top_fams, as.character(Family), "Other family")) %>%
  group_by(sample_label, io_group, Family2) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate(
    sample_label = factor(sample_label, levels = levels(sample_data(ps_fuhrberg)$sample_label)),
    pct = Abundance * 100,
    pct_lab = ifelse(pct < 0, "0", sprintf("%.2f", pct)),  # hide tiny values; adjust threshold if you want
  )

# 5) Order families by mean abundance (nice reading order)
fam_order <- df2 %>%
  group_by(Family2) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  pull(Family2)

df2$Family2 <- factor(df2$Family2, levels = rev(fam_order))

# 6) Heatmap + numbers
p_heat_vals <- ggplot(df2, aes(x = sample_label, y = Family2, fill = pct)) +
  geom_tile(color = "grey85", linewidth = 0.2) +
  geom_text(aes(label = pct_lab), size = 3, color = "black") +
  facet_grid(. ~ io_group, scales = "free_x", space = "free_x") +
  scale_fill_gradient(
    low = "white",
    high = "purple",
    name = "Relative Abundance (%)"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Family"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.15),
    plot.title.position = "plot",   # ← THIS is key
    legend.position = "right"
  )

p_heat_vals

# -----------------------------
# 7) Add sample_label, position, io_group, column  (CRITICAL for microeco)
# -----------------------------
sd <- data.frame(sample_data(ps_fuhrberg))
sd$sample_id <- rownames(sd)

sd <- sd %>%
  mutate(
    sample_label = case_when(
          sample_id == "JMF-2511-04-0034" ~ "Control – Inlet (0-4cm)",
          sample_id == "JMF-2511-04-0039" ~ "Control – Outlet (4-8cm)",

          sample_id == "JMF-2511-04-0031" ~ "Column 1 – Inlet (0-4cm)",
          sample_id == "JMF-2511-04-0032" ~ "Column 2 – Inlet (0-4cm)",
          sample_id == "JMF-2511-04-0033" ~ "Column 3 – Inlet (0-4cm)",

          sample_id == "JMF-2511-04-0036" ~ "Column 1 – Outlet (4-8cm)",
          sample_id == "JMF-2511-04-0037" ~ "Column 2 – Outlet (4-8cm)",
          sample_id == "JMF-2511-04-0038" ~ "Column 3 – Outlet (4-8cm)",

          TRUE ~ as.character(sample_type)
        )
,
    position = factor(ifelse(grepl("Inlet", sample_label), "Inlet (0-4cm)", "Outlet (4-8cm)"), levels = c("Inlet (0-4cm)","Outlet (4-8cm)")),
    io_group = position,
    column = case_when(
      grepl("^Control",  sample_label) ~ "Control",
      grepl("^Column 1", sample_label) ~ "Column 1",
      grepl("^Column 2", sample_label) ~ "Column 2",
      grepl("^Column 3", sample_label) ~ "Column 3",
      TRUE ~ NA_character_
    )
  )

sd$sample_label <- factor(
  sd$sample_label,
  levels = c(
    "Control – Inlet (0-4cm)", "Control – Outlet (4-8cm)",
    "Column 1 – Inlet (0-4cm)", "Column 1 – Outlet (4-8cm)",
    "Column 2 – Inlet (0-4cm)", "Column 2 – Outlet (4-8cm)",
    "Column 3 – Inlet (0-4cm)", "Column 3 – Outlet (4-8cm)"
  )
)
sd$position <- factor(sd$position, levels = c("Inlet (0-4cm)","Outlet (4-8cm)"))
sd$io_group <- factor(sd$io_group, levels = c("Inlet (0-4cm)","Outlet (4-8cm)"))
sd$column   <- factor(sd$column, levels = c("Control","Column 1","Column 2","Column 3"))

# write back into phyloseq
sample_data(ps_fuhrberg)$sample_label <- sd$sample_label
sample_data(ps_fuhrberg)$position     <- sd$position
sample_data(ps_fuhrberg)$io_group     <- sd$io_group
sample_data(ps_fuhrberg)$column       <- sd$column

# sanity
stopifnot("sample_label" %in% sample_variables(ps_fuhrberg))
stopifnot("position" %in% sample_variables(ps_fuhrberg))
stopifnot("column" %in% sample_variables(ps_fuhrberg))

# -----------------------------
# 8) Rarefy + relative abundance objects (FROM LABELED ps_fuhrberg)
# -----------------------------
set.seed(123)
ps_fuhrberg_rar <- rarefy_even_depth(
  ps_fuhrberg, sample.size = 3000,
  rngseed = 123, replace = TRUE,
  trimOTUs = TRUE, verbose = TRUE
)

ps_fuhrberg_rel <- transform_sample_counts(ps_fuhrberg, function(x) x / sum(x))

# -----------------------------
# 9) Convert to microeco (now column exists!)
# -----------------------------
library(microeco)
library(file2meco)

meco_rar <- phyloseq2meco(ps_fuhrberg_rar)
meco_rel <- phyloseq2meco(ps_fuhrberg_rel)

meco_rar$tidy_dataset()
meco_rel$tidy_dataset()

print(names(meco_rar$sample_table))  # should include column, position, sample_label

# ============================================================
# 10) Alpha diversity (Inlet vs Outlet; points colored by column)
# ============================================================
library(knitr)
library(ggplot2)
library(patchwork)

t_alpha <- trans_alpha$new(dataset = meco_rar, group = "position")
t_alpha$cal_diff(method = "wilcox")
t_alpha$res_diff %>% knitr::kable()

# palettes
col_palette <- c(
  "Control"  = "red",
  "Column 1" = "blue",
  "Column 2" = "gold",
  "Column 3" = "darkgreen"
)

# box/outline same for both positions (and hide legend)
pos_outline <- c("Inlet (0-4cm)" = "grey30", "Outlet (4-8cm)" = "grey30")

# ---- Build base plots FIRST (so $data exists) ----
p_alpha_chao_base <- t_alpha$plot_alpha(
  measure = "Chao1",
  add = NULL,
  size = 1
)

p_alpha_shan_base <- t_alpha$plot_alpha(
  measure = "Shannon",
  add = NULL,
  size = 1
)

# ---- Now add the jitter points using the base plot data ----
p_alpha_chao <- p_alpha_chao_base +
  geom_jitter(
    data = p_alpha_chao_base$data,
    aes(x = position, y = Value, fill = column),
    width = 0.12,
    shape = 21,
    color = "black",
    size = 2.8,
    stroke = 0.25,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = makie_colors, name = NULL) +
  scale_color_manual(values = pos_outline, guide = "none") +
  labs(title = "Chao1 Index") +
  theme_makie()

p_alpha_shan <- p_alpha_shan_base +
  geom_jitter(
    data = p_alpha_shan_base$data,
    aes(x = position, y = Value, fill = column),
    width = 0.12,
    shape = 21,
    color = "black",
    size = 2.8,
    stroke = 0.25,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = makie_colors, name = NULL) +
  scale_color_manual(values = pos_outline, guide = "none") +
  labs(title = "Shannon Index") +
  theme_makie()

p_alpha <- (p_alpha_shan + p_alpha_chao) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p_alpha

# ============================================================
# 11) Beta diversity: PERMANOVA + PERMDISP + PCoA (paired segments)
# ============================================================
library(vegan)

otu_rel <- as(otu_table(ps_fuhrberg_rel), "matrix")
if (taxa_are_rows(ps_fuhrberg_rel)) otu_rel <- t(otu_rel)

meta <- data.frame(sample_data(ps_fuhrberg_rel))
meta$Sample <- rownames(meta)

bray <- vegdist(otu_rel, method = "bray")

metaMDS(bray)

adon_col <- adonis2(bray ~ column, data = meta, permutations = 9999)
print(adon_col)

bd <- betadisper(bray, meta$column)
print(anova(bd))
print(permutest(bd, permutations = 999))

pcoa <- cmdscale(bray, eig = TRUE, k = 2)

ord_df <- data.frame(
  Sample = rownames(pcoa$points),
  PCoA1  = pcoa$points[, 1],
  PCoA2  = pcoa$points[, 2]
) %>%
  left_join(meta, by = "Sample")

var_expl <- pcoa$eig / sum(pcoa$eig)
var_expl <- var_expl[1:2]

ord_df2 <- ord_df %>%
  mutate(position = factor(position, levels = c("Inlet (0-4cm)", "Outlet (4-8cm)"))) %>%
  arrange(column, position)

p_beta_pcoa <- ggplot(ord_df2, aes(PCoA1, PCoA2)) +
  geom_path(aes(group = column, color = column),
            linewidth = 0.8, alpha = 0.6, show.legend = FALSE) +
  geom_point(aes(fill = column, shape = position),
             size = 3.2, color = "black", stroke = 0.3) +
  scale_fill_manual(values = makie_colors, name = NULL) +
  scale_color_manual(values = makie_colors, guide = "none") +
  scale_shape_manual(values = c("Inlet (0-4cm)" = 21, "Outlet (4-8cm)" = 24), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, color = NA))) +
  labs(
    title = "Beta diversity (Bray–Curtis PCoA)",
    x = paste0("PCoA1 (", round(100 * var_expl[1], 1), "%)"),
    y = paste0("PCoA2 (", round(100 * var_expl[2], 1), "%)")
  ) +
  theme_makie() +
  theme(legend.position = "bottom")

p_beta_pcoa


############################################################
## DESeq2 analysis for Fuhrberg columns
## Family level
############################################################

library(phyloseq)
library(DESeq2)
library(dplyr)
library(tibble)
library(knitr)

############################################################
## 1) Start from Fuhrberg phyloseq object
##    Assumes ps_fuhrberg already exists
############################################################

# Agglomerate to Family level
ps_family <- tax_glom(ps_fuhrberg, taxrank = "Family", NArm = TRUE)

# Optional: remove very low-count families
ps_family <- prune_taxa(taxa_sums(ps_family) > 10, ps_family)

############################################################
## 2) Build metadata variables from sample_label
############################################################

meta <- data.frame(sample_data(ps_family))
meta$Sample <- rownames(meta)

# If column is not already present, derive it
if (!"column" %in% colnames(meta)) {
  meta$column <- case_when(
    grepl("^Control", meta$sample_label)  ~ "Control",
    grepl("^Column 1", meta$sample_label) ~ "Column 1",
    grepl("^Column 2", meta$sample_label) ~ "Column 2",
    grepl("^Column 3", meta$sample_label) ~ "Column 3",
    TRUE ~ NA_character_
  )
}

# If treatment is not already present, derive it
if (!"treatment" %in% colnames(meta)) {
  meta$treatment <- ifelse(meta$column == "Control", "Control", "Nitrate-fed")
}

# Make factors with explicit order
meta$column <- factor(
  meta$column,
  levels = c("Control", "Column 1", "Column 2", "Column 3")
)

meta$treatment <- factor(
  meta$treatment,
  levels = c("Control", "Nitrate-fed")
)

meta$position <- factor(
  meta$position,
  levels = c("Inlet (0-4cm)", "Outlet (4-8cm)")
)

# Write back into phyloseq
sample_data(ps_family)$column <- meta$column
sample_data(ps_family)$treatment <- meta$treatment
sample_data(ps_family)$position <- meta$position

# Check
sample_variables(ps_family)
table(data.frame(sample_data(ps_family))$treatment, useNA = "ifany")
table(data.frame(sample_data(ps_family))$position, useNA = "ifany")
table(data.frame(sample_data(ps_family))$column, useNA = "ifany")

############################################################
## 3) DESeq2: Control vs Nitrate-fed
############################################################

dds_treatment <- phyloseq_to_deseq2(ps_family, ~ treatment)

dds_treatment <- DESeq(dds_treatment, fitType = "parametric")

res_treatment <- results(
  dds_treatment,
  contrast = c("treatment", "Nitrate-fed", "Control")
)

res_treatment <- as.data.frame(res_treatment)
res_treatment$taxon_id <- rownames(res_treatment)

############################################################
## 4) DESeq2: Inlet vs Outlet
############################################################

dds_position <- phyloseq_to_deseq2(ps_family, ~ position)

dds_position <- DESeq(dds_position, fitType = "parametric")

res_position <- results(
  dds_position,
  contrast = c("position", "Outlet (4-8cm)", "Inlet (0-4cm)")
)

res_position <- as.data.frame(res_position)
res_position$taxon_id <- rownames(res_position)

############################################################
## 5) Add taxonomy to results
############################################################

tax_df <- as.data.frame(tax_table(ps_family))
tax_df$taxon_id <- rownames(tax_df)

res_treatment_annot <- left_join(res_treatment, tax_df, by = "taxon_id")
res_position_annot  <- left_join(res_position,  tax_df, by = "taxon_id")

############################################################
## 6) Significant results
############################################################

sig_treatment <- res_treatment_annot %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj)

sig_position <- res_position_annot %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj)

############################################################
## 7) View top results
############################################################

# Top taxa enriched in nitrate-fed (log2FoldChange)
sig_treatment %>%
  arrange(desc(log2FoldChange)) %>%
  select(Family, log2FoldChange, pvalue, padj) %>%
  head(40) %>%
  knitr::kable()


# Top taxa enriched in outlet (log2FoldChange)
sig_position %>%
  arrange(desc(log2FoldChange)) %>%
  select(Family, log2FoldChange, pvalue, padj) %>%
  head(40) %>%
  knitr::kable()

############################################################
## 8) export results
############################################################
write.csv(sig_treatment,
          "data/processed_results/DESeq2_Fuhrberg_treatment_family_significant.csv",
          row.names = FALSE)

write.csv(sig_position,
          "data/processed_results/DESeq2_Fuhrberg_position_family_significant.csv",
          row.names = FALSE)

############################################################
## 9) Deseq plot
############################################################


deseq_plot_top <- sig_treatment %>%
  filter(!is.na(Family), !Family %in% c("Bacteria Domain", "Sericytochromatia Class", "Kapabacteriales Order", "mle1-27 Order",
                                        "053A03-B-DI-P58 Order", "053A03-B-DI-P58 Order", "Lineage IV Order", "MA-28-I98C", "Candidatus Lloydbacteria Order",
                                        "Bacteriovoracaceae", "Nocardiaceae", "Bacteroidetes BD2-2",
                                        "Lentimicrobiaceae")) %>%
  slice_max(order_by = abs(log2FoldChange), n = 30) %>%
  mutate(
    direction = ifelse(log2FoldChange > 0, "Nitrate-fed", "Control")
  ) %>%
  arrange(log2FoldChange)

deseq_plot_top$Family <- factor(deseq_plot_top$Family, levels = deseq_plot_top$Family)

# Shorten long family names for the y-axis
shorten_label <- function(x, n = 10) {
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "."), x)
}

p_deseq_bar_top <- ggplot(deseq_plot_top, aes(x = log2FoldChange, y = Family, fill = direction)) +
  geom_col(color = "black", width = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  scale_fill_manual(
    values = c("Control" = "#4d9078", "Nitrate-fed" = "#4d9078"),
  ) +
  scale_x_continuous(limits = c(-10, 10)) +
  scale_y_discrete(labels = shorten_label) +
  labs(
    x = "Log2 Fold Change",
    y = NULL,
    title = "b. Differential Abundance (Nitrate-fed vs Control)"
  ) +
  theme_makie() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "italic")
  )

p_deseq_bar_top


p_combined <- pB / p_deseq_bar_top +
  plot_layout(heights = c(1.8, 1))

p_combined

#### SAVE PLOTS ####

# p_alpha (Double column width, shorter height)
ggsave("figs/p_alpha.pdf", plot = p_alpha, width = 174, height = 100, units = "mm", device = cairo_pdf)
ggsave("figs/p_alpha.png", plot = p_alpha, width = 174, height = 100, units = "mm", dpi = 600)

# p_beta_pcoa (1.5 column width)
ggsave("figs/p_beta_pcoa.pdf", plot = p_beta_pcoa, width = 129, height = 100, units = "mm", device = cairo_pdf)
ggsave("figs/p_beta_pcoa.png", plot = p_beta_pcoa, width = 129, height = 100, units = "mm", dpi = 600)

# p_combined (Double column width, adjusted height)
ggsave("figs/p_combined.pdf", plot = p_combined, width = 174, height = 180, units = "mm", device = cairo_pdf)
ggsave("figs/p_combined.png", plot = p_combined, width = 174, height = 180, units = "mm", dpi = 600)