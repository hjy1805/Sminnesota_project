# Load necessary libraries
library(tidyverse)
library(ape)
library(treedataverse)
library(ggnewscale)
library(rhierbaps)
library(phytools)

## Construct the neighbour-joining tree

# Read DNA sequences in FASTA format
dna <- read.dna("./sminnesota_final_all_clincal.snp_sites.aln", format = "fasta")

# Compute distance matrix using nucleotide model (N) and pairwise deletion
distdna <- dist.dna(dna, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)

# Construct a Neighbor-Joining (NJ) tree from the distance matrix
njtree <- nj(distdna)

# Root the NJ tree at its midpoint
njtree <- midpoint.root(njtree)

# Plot the NJ tree
plot(njtree)

# Save the NJ tree to a file
write.tree(njtree, file = "clincal_njtree.tree")

# Run RhierBAPS analysis
# Load the SNP matrix from the FASTA file
snp.matrix <- load_fasta(dna)

# Perform hierarchical BAPS analysis
hb.results <- hierBAPS(snp.matrix, max.depth = 2, n.pops = 10, quiet = TRUE)

# Display the top rows of the partition dataframe from RhierBAPS results
head(hb.results$partition.df)

# Extract Isolate information from RhierBAPS results
hb.results$partition.df$Isolate

# Read metadata and continent data
metadata <- read_tsv("../Sminnesota_all_metadata_28Feb.tsv")
continent <- read_csv("/Users/huanj0f/Documents/E.coli_selected/continent.csv")

# Merge metadata with continent data
metadata <- left_join(metadata, continent, by = c("Country" = 'country'))

# Merge RhierBAPS results with metadata and rename BAPS column
metadata_hirebaps <- left_join(hb.results$partition.df, metadata, by = c("Isolate" = 'Run')) %>%
  rename(BAPS = `level 2`) %>%
  separate(collection_date, into = c("collection_year", "collection_month", "collection_day"), sep = "/", remove = FALSE)

# Merge NJ tree with metadata
njtree <- left_join(njtree, metadata, by = c("label" = "Isolate"))

# Save the merged metadata with BAPS results to a file
write_tsv(metadata_hirebaps, "./Sminnesota_all_metadata_baps.tsv")

# Create the base tree plot with tips colored by BAPS clusters
p1 <- ggtree(njtree) +
  geom_tippoint(aes(color = factor(BAPS)), size = 1.5, alpha = 1) +
  theme_tree2(legend.position = 'right') +
  scale_colour_discrete("BAPS")
p1

# Prepare metadata for continent heatmap
tree_continent <- metadata %>% select(continent)
tree_continent <- as.data.frame(tree_continent)
rownames(tree_continent) <- metadata$Isolate

# Define colors for continents and add 'grey' to the color palette
colors <- c(colors, "grey")

# Add the continent heatmap to the tree plot
p2 <- gheatmap(p1, tree_continent, offset = 0.05, width = 0.1, font.size = 3, colnames = FALSE, hjust = 0, color = NA) +
  scale_x_ggtree() +
  scale_y_continuous(expand = c(0, 0.00003)) +
  scale_fill_manual(breaks = c(unique(tree_continent$continent)), values = colors, name = "Continent")
p2

# Prepare metadata for country heatmap
tree_country <- metadata_hirebaps %>% select(Country)
tree_country <- as.data.frame(tree_country)
rownames(tree_country) <- metadata_hirebaps$Isolate

# Display the number of unique countries
length(unique(tree_country$Country))

# Add a new scale for the country heatmap
p3 <- p2 + new_scale_fill()

# Add the country heatmap to the tree plot
p4 <- gheatmap(p3, tree_country, offset = 0.07, width = 0.1, font.size = 3, colnames = FALSE, hjust = 0, color = NA) +
  scale_x_ggtree() +
  scale_y_continuous(expand = c(0, 0.00003)) +
  scale_fill_manual(values = c("#DD6E42", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#00FFFF", "#FF00FF", "#800000", "#008000", "#000080", "#808000", "#800080", "#008080", "#B1740F", "#FFD07B", "#FFA500", "#FFFF80", "#80FF00", "#80FFFF", "#FF80C0", "#FF0080", "#8000FF"), name = "Continent")
p4

# Prepare metadata for focused country heatmap
tree_country_f <- metadata_hirebaps %>% select(Country_focused)
tree_country_f <- as.data.frame(tree_country_f)
rownames(tree_country_f) <- metadata_hirebaps$Isolate

# Add a new scale for the focused country heatmap
p5 <- p4 + new_scale_fill()

# Add the focused country heatmap to the tree plot
p6 <- gheatmap(p5, tree_country_f, offset = 0.09, width = 0.1, font.size = 3, colnames = FALSE, hjust = 0, color = NA) +
  scale_x_ggtree() +
  scale_y_continuous(expand = c(0, 0.00003)) +
  scale_fill_manual(breaks = c("Brazil", "Saudi Arabia", "Others"), values = c('#432772', '#8FCF63', '#C0C0C0'), name = "Country focused")
p6

# Prepare metadata for sample origin heatmap
tree_origin <- metadata_hirebaps %>% select(Origin)
tree_origin <- as.data.frame(tree_origin)
rownames(tree_origin) <- metadata_hirebaps$Isolate

# Add a new scale for the sample origin heatmap
p7 <- p6 + new_scale_fill()

# Add the sample origin heatmap to the tree plot
p8 <- gheatmap(p7, tree_origin, offset = 0.11, width = 0.1, font.size = 3, colnames = FALSE, hjust = 0, color = NA) +
  scale_x_ggtree() +
  scale_y_continuous(expand = c(0, 0.00003)) +
  scale_fill_manual(breaks = c("External", "Internal"), values = c('blue', '#90323D'), name = "Sample Origin")
p8

# Prepare metadata for collection year heatmap
tree_year <- metadata_hirebaps %>% select(collection_year)
tree_year <- as.data.frame(tree_year)
rownames(tree_year) <- metadata_hirebaps$Isolate

# Display the number of unique collection years
length(unique(tree_year$collection_year))

# Add a new scale for the collection year heatmap
p9 <- p8 + new_scale_fill()

# Add the collection year heatmap to the tree plot
p10 <- gheatmap(p9, tree_year, offset = 0.13, width = 0.1, font.size = 3, colnames = FALSE, hjust = 0, color = NA) +
  scale_x_ggtree() +
  scale_y_continuous(expand = c(0, 0.00003)) +
  scale_fill_viridis_d(option = "C", name = "Collection year") +
  coord_cartesian(clip = "off")
p10

# Prepare metadata for ONT sequencing heatmap
tree_ont <- metadata_hirebaps %>% select(ONT_seq)
tree_ont <- as.data.frame(tree_ont)
rownames(tree_ont) <- metadata_hirebaps$Isolate

# Add a new scale for the ONT sequencing heatmap
p11 <- p10 + new_scale_fill()

# Add the ONT sequencing heatmap to the tree plot
p12 <- gheatmap(p11, tree_ont, offset = 0.15, width = 0.1, font.size = 3, colnames = FALSE, hjust = 0, color = NA) +
  scale_x_ggtree() +
  scale_y_continuous(expand = c(0, 0.00003)) +
  scale_fill_manual(breaks = c("Yes", "No"), values = c('#FF6663', '#BFD7EA'), name = "ONT sequencing")
p12


## Bubble plot for the AMR & virulence factor detected from AMRFinderPlus

# Read in metadata with BAPS assignments
metadata_hirebaps <- read_tsv("./Sminnesota_all_metadata_baps.tsv")

# Filter out specific BAPS groups from the metadata
meta_baps_filtered <- metadata_hirebaps %>% filter(!BAPS %in% c(10, 11, 14)) 

# Read in AMRFinder results and select relevant columns
amrfinder <- read_tsv("./amrfinder_result.tsv") %>% select(c("Sample_name", "Gene_symbol", "Element_type"))

# Separate AMR and virulence genes
amr <- amrfinder %>% filter(Element_type == "AMR") %>% select(-c("Element_type")) 
vf <- amrfinder %>% filter(Element_type == "VIRULENCE") %>% select(-c("Element_type"))

# Count the number of samples in each BAPS group
BAPS_count <- meta_baps_filtered %>% group_by(BAPS) %>% summarise(sample_count = n()) 

# Calculate frequency of AMR genes in each BAPS group
amr_freq <- left_join(amr, meta_baps_filtered, by = c("Sample_name" = "Isolate")) %>% 
  filter(!is.na(BAPS)) %>% 
  group_by(BAPS, Gene_symbol) %>% 
  summarise(count = n())
amr_freq_2 <- left_join(amr_freq, BAPS_count, by = c("BAPS" = "BAPS")) %>% 
  mutate(frequency = count / sample_count)

# Create bubble plot for AMR gene frequencies
p13 <- ggplot(amr_freq_2, aes(x = Gene_symbol, y = factor(BAPS), size = frequency))
p13 + geom_point(alpha = 0.7, color = "#107E7D", fill = "#107E7D", shape = 21) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1, face = "bold"), 
        plot.margin = margin(t = 10, r = 25, b = 10, l = 30)) +
  scale_size(range = c(0.1, 8), name = "Frequency") +
  labs(x = "AMR gene", y = "BAPS") + 
  theme(text = element_text(size = 18, face = "bold"), 
        plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold"), 
        axis.text.y = element_text(size = 15, color = "black", face = "bold"), 
        legend.title = element_text(size = 15, face = "bold"), 
        legend.text = element_text(size = 18, face = "bold"))

# Calculate frequency of virulence genes in each BAPS group
vf_freq <- left_join(vf, meta_baps_filtered, by = c("Sample_name" = "Isolate")) %>% 
  filter(!is.na(BAPS)) %>% 
  group_by(BAPS, Gene_symbol) %>% 
  summarise(count = n())
vf_freq_2 <- left_join(vf_freq, BAPS_count, by = c("BAPS" = "BAPS")) %>% 
  mutate(frequency = count / sample_count)

# Create bubble plot for virulence gene frequencies
p14 <- ggplot(vf_freq_2, aes(x = Gene_symbol, y = factor(BAPS), size = frequency))
p14 + geom_point(alpha = 0.7, color = "#E3B505", fill = "#E3B505", shape = 21) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1, face = "bold"), 
        plot.margin = margin(t = 10, r = 25, b = 10, l = 30)) +
  scale_size(range = c(0.1, 8), name = "Frequency") +
  labs(x = "Virulence gene", y = "BAPS") + 
  theme(text = element_text(size = 18, face = "bold"), 
        plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold"), 
        axis.text.y = element_text(size = 15, color = "black", face = "bold"), 
        legend.title = element_text(size = 15, face = "bold"), 
        legend.text = element_text(size = 18, face = "bold"))


## Skygrowth plot for the effective population size
tree <- read.nexus("./BAPS1/BAPS1_city_mascot.tree")
fit <- skygrowth.map( tree
                      , tau0 = 0.05    # Smoothing parameter. If prior is not specified, this will also set the scale of the prior
)
p15<-plot(fit)+theme_bw()+
  theme(axis.text.x = element_text(size=30, hjust = 1),
        axis.text.y = element_text(size=30, hjust = 1),
        axis.title.x = element_text(color="black", size=35, face="bold"),
        axis.title.y = element_text(color="black", size=35, face="bold"),
        strip.text.x = element_text(color="black", size=35, face = "bold")
  )
p15

## BEAST analysis result

# Compute mutation rate and its confidence interval
beast_result_saudi <- beast_result_saudi %>% 
  mutate(mutation_rate = clockrate * beast_alignment) %>% 
  mutate(mt_lower = cr_lower * beast_alignment) %>% 
  mutate(mt_upper = cr_upper * beast_alignment)

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Plot 1: Root age by BAPS
p1 <- ggplot(beast_result_saudi, aes(x = BAPS, y = age)) +
  geom_bar(stat = "identity", fill = "#FF6666") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15, vjust = 1, hjust = 1, face = "bold"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 35)) + 
  coord_cartesian(ylim = c(0, 15)) + 
  scale_y_continuous(breaks = seq(0, 15, 5)) +
  labs(x = "BAPS", y = "Age") + 
  theme(text = element_text(size = 35, face = "bold", color = "black"),
        plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold", color = "black"), 
        axis.text.x = element_text(size = 35, color = "black", face = "bold"),
        axis.text.y = element_text(size = 35, color = "black", face = "bold"), 
        legend.title = element_text(size = 35, face = "bold", color = "black"), 
        legend.text = element_text(size = 35, face = "bold", color = "black")) + 
  geom_errorbar(aes(x = BAPS, ymin = age_lower, ymax = age_upper), width = 0.4, colour = "black", alpha = 0.9, size = 1)

print(p1)

# Plot 2: Migration rate by BAPS and route
p2 <- ggplot(beast_migration, aes(x = BAPS, y = rate, fill = route)) +
  geom_bar(position = "dodge", stat = "identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1, face = "bold"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 35)) + 
  coord_cartesian(ylim = c(0, 5)) + 
  scale_y_continuous(breaks = seq(0, 5, 0.5)) +
  labs(x = "BAPS", y = "Migration rate", fill = "Route") + 
  theme(text = element_text(size = 18, face = "bold", color = "black"), 
        plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold", color = "black"), 
        axis.text.y = element_text(size = 18, color = "black", face = "bold"), 
        legend.title = element_text(size = 18, face = "bold", color = "black"), 
        legend.text = element_text(size = 18, face = "bold", color = "black")) + 
  geom_errorbar(aes(x = BAPS, ymin = lower, ymax = upper), width = 0.4, colour = "black", alpha = 0.9, size = 1, position = position_dodge(.9)) +
  scale_fill_manual(values = c("#3C5488FF","#4DBBD5FF","#E64B35FF","#00A087FF","#9D1E3E","#F9E20D", "#7D60D5"))

print(p2)

# Plot 3: Effective population size by BAPS and city
p3 <- ggplot(beast_ne, aes(x = BAPS, y = ne, fill = city)) +
  geom_bar(position = "dodge", stat = "identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1, face = "bold"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 35)) + 
  coord_cartesian(ylim = c(0, 20)) + 
  scale_y_continuous(breaks = seq(0, 20, 5)) +
  labs(x = "BAPS", y = "Effective population size", fill = "City") + 
  theme(text = element_text(size = 18, face = "bold", color = "black"), 
        plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold", color = "black"), 
        axis.text.y = element_text(size = 18, color = "black", face = "bold"), 
        legend.title = element_text(size = 18, face = "bold", color = "black"), 
        legend.text = element_text(size = 18, face = "bold", color = "black")) + 
  geom_errorbar(aes(x = BAPS, ymin = lower, ymax = upper), width = 0.4, colour = "black", alpha = 0.9, size = 1, position = position_dodge(.9)) +
  scale_fill_manual(values = c("#48A9A6","#E4DFDA","#D4B483"))

print(p3)

# Plot 4: Mutation rate by BAPS
p4 <- ggplot(beast_result_saudi, aes(x = BAPS, y = mutation_rate_genome)) +
  geom_bar(stat = "identity", fill = "#73AB84") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 25, vjust = 1, hjust = 1, face = "bold"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 35)) + 
  coord_cartesian(ylim = c(0, 0.000004)) + 
  scale_y_continuous(breaks = seq(0, 0.000004, 0.000001)) +
  labs(x = "BAPS", y = "Mutation rate") + 
  theme(text = element_text(size = 25, face = "bold", color = "black"), 
        plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold", color = "black"), 
        axis.text.y = element_text(size = 25, color = "black", face = "bold"), 
        legend.title = element_text(size = 25, face = "bold", color = "black"), 
        legend.text = element_text(size = 25, face = "bold", color = "black")) + 
  geom_errorbar(aes(x = BAPS, ymin = mt_lower_genome, ymax = mt_upper_genome), width = 0.4, colour = "black", alpha = 0.9, size = 1, position = position_dodge(.9))

print(p4)











