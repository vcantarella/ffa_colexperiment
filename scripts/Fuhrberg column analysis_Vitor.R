#16S rRNA analysis for Vitor's column experiment

#Load packages
# Set CRAN mirror to avoid interactive prompt
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

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

# Helper to install
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }, error = function(e) {
      message(paste("Failed to install", pkg, ":", e$message))
    })
  }
}

# List of packages
pkgs <- c("DECIPHER", "phytools", "lubridate", "microeco", 
          "phangorn", "DESeq2", "FSA", 
          "rcompanion", "ggsci", "microbiome", "phyloseq",
          "picante", "vegan", "ape", "nlme", "ggplot2", 
          "reshape2", "tidyverse", "RColorBrewer", "patchwork", 
          "file2meco", "GUniFrac", "hillR", "ggridges", "tibble", "scales",
          "decontam", "ggpubr")

for (pkg in pkgs) {
  install_if_missing(pkg)
}


library(vegan)        # ecological functions
library(picante)      # ecological functions
library(ape)          # phylogenetic functions
library(nlme)
library(ggplot2)      # plotting package
library(reshape2)     # data manipulation
library(tidyverse)    # data manipulation
library(phytools)     # phylogenetic functions
library(phyloseq)     # microbial analyses
library(RColorBrewer) # plotting package
library(microbiome)   # microbial analyses
library(microeco)     # microbial analyses
library(microViz)     # microbial analyses
library(DECIPHER)     # phylogenetic functions
library(phangorn)     # phylogenetic functions
library(decontam)     # decontamination
library(patchwork)    # plotting package
library(file2meco)    # switching between micro packages
library(DESeq2)       # normalization
library(picante)      # phylogenetic diversity
library(GUniFrac)     # phylogenetic diversity
library(hillR)        # hill numbers
library(FSA)          # KW_dunn letters
library(rcompanion)    # KW_dunn letters
library(dplyr)
library(ggridges)
library(tibble)
library(scales)

# setwd("~/Library/CloudStorage/OneDrive-UniversityofthePhilippines/Vitor Column Sequences")


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

# once you give it a name (plot_1.1) you can call on this plot whenever you want
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
ps_fuhrberg_rar <- rarefy_even_depth(ps_fuhrberg, sample.size = 1500,
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
    sample_id == "JMF-2511-04-0034" ~ "Control – Inlet",
    sample_id == "JMF-2511-04-0039" ~ "Control – Outlet",
    
    sample_id == "JMF-2511-04-0031" ~ "Column 1 – Inlet",
    sample_id == "JMF-2511-04-0032" ~ "Column 2 – Inlet",
    sample_id == "JMF-2511-04-0033" ~ "Column 3 – Inlet",
    
    sample_id == "JMF-2511-04-0036" ~ "Column 1 – Outlet",
    sample_id == "JMF-2511-04-0037" ~ "Column 2 – Outlet",
    sample_id == "JMF-2511-04-0038" ~ "Column 3 – Outlet",  # <-- change this to the correct ID for Column 3 outlet
    
    TRUE ~ as.character(sd$sample_type)
  ))

sample_data(ps_fuhrberg)$sample_label <- sd$sample_label

sample_data(ps_fuhrberg)$sample_label <- factor(
  sample_data(ps_fuhrberg)$sample_label,
  levels = c(
    "Control – Inlet", "Control – Outlet",
    "Column 1 – Inlet", "Column 1 – Outlet",
    "Column 2 – Inlet", "Column 2 – Outlet",
    "Column 3 – Inlet", "Column 3 – Outlet"
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
## 1) Family composition (counts preserved; percent axis)
## =========================
# Uses comp_barplot(); if your comp_barplot has no `x=` argument, see fallback below.

pB <- ps_fuhrberg %>%
  comp_barplot(
    tax_level = "Family",
    n_taxa = 30,
    x = "sample_label",
    other_name = "Other family",
    merge_other = FALSE,
    bar_width = 0.9,
    bar_outline_colour = "grey20"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Relative abundance (%)", title = "Community composition (Family)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  ) +
  guides(fill = guide_legend(reverse = TRUE))

pB



library(phyloseq)
library(dplyr)
library(tibble)
library(ggplot2)

stopifnot("sample_label" %in% sample_variables(ps_fuhrberg))

# -----------------------------
# Move Metadata Enrichment Here (before Heatmap)
# -----------------------------
sd <- data.frame(sample_data(ps_fuhrberg))
sd$sample_id <- rownames(sd)

sd <- sd %>%
  mutate(
    sample_label = case_when(
      sample_id == "JMF-2511-04-0034" ~ "Column 4 – 0-4 cm",
      sample_id == "JMF-2511-04-0039" ~ "Column 4 – 4-8 cm",
      
      sample_id == "JMF-2511-04-0031" ~ "Column 1 – 0-4 cm",
      sample_id == "JMF-2511-04-0032" ~ "Column 2 – 0-4 cm",
      sample_id == "JMF-2511-04-0033" ~ "Column 3 – 0-4 cm",
      
      sample_id == "JMF-2511-04-0036" ~ "Column 1 – 4-8 cm",
      sample_id == "JMF-2511-04-0037" ~ "Column 2 – 4-8 cm",
      sample_id == "JMF-2511-04-0038" ~ "Column 3 – 4-8 cm",
      
      TRUE ~ as.character(sample_type)
    ),
    position = ifelse(grepl("0-4 cm", sample_label), "0-4 cm", "4-8 cm"),
    io_group = position,
    column = case_when(
      grepl("^Column 1", sample_label) ~ "Column 1",
      grepl("^Column 2", sample_label) ~ "Column 2",
      grepl("^Column 3", sample_label) ~ "Column 3",
      grepl("^Column 4", sample_label) ~ "Column 4",
      TRUE ~ NA_character_
    )
  )

sd$sample_label <- factor(
  sd$sample_label,
  levels = c(
    "Column 1 – 0-4 cm", "Column 1 – 4-8 cm",
    "Column 2 – 0-4 cm", "Column 2 – 4-8 cm",
    "Column 3 – 0-4 cm", "Column 3 – 4-8 cm",
    "Column 4 – 0-4 cm", "Column 4 – 4-8 cm"
  )
)
sd$position <- factor(sd$position, levels = c("0-4 cm","4-8 cm"))
sd$io_group <- factor(sd$io_group, levels = c("0-4 cm","4-8 cm"))
sd$column   <- factor(sd$column, levels = c("Column 1","Column 2","Column 3","Column 4"))

# write back into phyloseq
sample_data(ps_fuhrberg)$sample_label <- sd$sample_label
sample_data(ps_fuhrberg)$position     <- sd$position
sample_data(ps_fuhrberg)$io_group     <- sd$io_group
sample_data(ps_fuhrberg)$column       <- sd$column


# 1) Agglomerate to Family
ps_fam <- tax_glom(ps_fuhrberg, taxrank = "Family", NArm = TRUE)

# 2) Relative abundance
ps_fam_rel <- transform_sample_counts(ps_fam, function(x) x / sum(x))

# 3) Long format + attach sample_label + io_group + column
df <- psmelt(ps_fam_rel) %>%
  mutate(
    sample_label = sample_data(ps_fuhrberg)$sample_label[match(Sample, sample_names(ps_fuhrberg))],
    io_group     = sample_data(ps_fuhrberg)$io_group[match(Sample, sample_names(ps_fuhrberg))],
    column       = sample_data(ps_fuhrberg)$column[match(Sample, sample_names(ps_fuhrberg))]
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
  group_by(sample_label, column, Family2) %>%
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
  facet_grid(. ~ column, scales = "free_x", space = "free_x") +
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
# ALREADY DONE ABOVE - logic moved before heatmap

# -----------------------------
# 8) Rarefy + relative abundance objects (FROM LABELED ps_fuhrberg)
# -----------------------------
set.seed(123)
ps_fuhrberg_rar <- rarefy_even_depth(
  ps_fuhrberg, sample.size = 1500,
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
  "Column 4" = "red",
  "Column 1" = "blue",
  "Column 2" = "gold",
  "Column 3" = "darkgreen"
)

# box/outline same for both positions (and hide legend)
pos_outline <- c("0-4 cm" = "grey30", "4-8 cm" = "grey30")

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
  scale_color_manual(values = pos_outline, guide = "none") +
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
  scale_fill_manual(values = col_palette, name = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    axis.line = element_line(color = "black"),
    plot.title = element_text(face = "bold")
  )

p_alpha_shan <- p_alpha_shan_base +
  scale_color_manual(values = pos_outline, guide = "none") +
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
  scale_fill_manual(values = col_palette, name = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    axis.line = element_line(color = "black"),
    plot.title = element_text(face = "bold")
  )

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
  mutate(position = factor(position, levels = c("0-4 cm", "4-8 cm"))) %>%
  arrange(column, position)

p_beta_pcoa <- ggplot(ord_df2, aes(PCoA1, PCoA2)) +
  geom_path(aes(group = column, color = column),
            linewidth = 0.8, alpha = 0.9, show.legend = FALSE) +
  geom_point(aes(fill = column, shape = position),
             size = 3.2, color = "black", stroke = 0.3) +
  scale_fill_manual(values = col_palette, name = NULL) +
  scale_color_manual(values = col_palette, guide = "none") +
  scale_shape_manual(values = c("0-4 cm" = 21, "4-8 cm" = 24), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, color = NA))) +
  labs(
    title = "Beta diversity (Bray–Curtis PCoA)",
    x = paste0("PCoA1 (", round(100 * var_expl[1], 1), "%)"),
    y = paste0("PCoA2 (", round(100 * var_expl[2], 1), "%)")
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    axis.line = element_line(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.title.position = "plot"
  )

p_beta_pcoa

#### SAVE PLOTS ####

# plot_1.1
ggsave("figs/plot_1.1.png", plot = plot_1.1, width = 10, height = 8)
ggsave("figs/plot_1.1.pdf", plot = plot_1.1, width = 10, height = 8)

# plot_1.2
ggsave("figs/plot_1.2.png", plot = plot_1.2, width = 10, height = 8)
ggsave("figs/plot_1.2.pdf", plot = plot_1.2, width = 10, height = 8)

# plot_fuhrberg
ggsave("figs/plot_fuhrberg.png", plot = plot_fuhrberg, width = 10, height = 8)
ggsave("figs/plot_fuhrberg.pdf", plot = plot_fuhrberg, width = 10, height = 8)

# pA
ggsave("figs/pA.png", plot = pA, width = 10, height = 8)
ggsave("figs/pA.pdf", plot = pA, width = 10, height = 8)

# pB
ggsave("figs/pB.png", plot = pB, width = 12, height = 8)
ggsave("figs/pB.pdf", plot = pB, width = 12, height = 8)

# p_heat_vals
# Saved with high quality (300 dpi) and specified dimensions (max 174mm width)
# 174 mm = 6.85 inches
# 234 mm = 9.21 inches
ggsave("figs/p_heat_vals.png", plot = p_heat_vals, width = 6.85, height = 9, dpi = 300)
ggsave("figs/p_heat_vals.pdf", plot = p_heat_vals, width = 6.85, height = 9)

# p_alpha
ggsave("figs/p_alpha.png", plot = p_alpha, width = 10, height = 6)
ggsave("figs/p_alpha.pdf", plot = p_alpha, width = 10, height = 6)

# p_beta_pcoa
ggsave("figs/p_beta_pcoa.png", plot = p_beta_pcoa, width = 8, height = 6)
ggsave("figs/p_beta_pcoa.pdf", plot = p_beta_pcoa, width = 8, height = 6)

