# Load necessary library
library(tidyverse)


## Processing Internal sequencing metadata
# Read and process forward size data
forward_size <- read_tsv("./forward_size.tsv") %>%
  separate(file, into = c("text1", "tag", "text3"), sep = "/", remove = TRUE) %>%
  select(-c("text1", "text3")) %>%
  rename("forward_size" = "size") %>%
  relocate(tag)

# Read and process reverse size data
reverse_size <- read_tsv("./reverse_size.tsv") %>%
  separate(file, into = c("text1", "tag", "text3"), sep = "/", remove = TRUE) %>%
  select(-c("text1", "text3")) %>%
  rename("reverse_size" = "size")

# Read other necessary datasets
map_coverage <- read_tsv("./coverage.tsv")
assembly_qc <- read_tsv("./qc_table.tsv")

# Read and process MLST result data
mlst <- read_tsv("./mlst_result.txt") %>%
  separate(file, into = c("text1", "tag", "text3", "test4"), sep = "/", remove = TRUE) %>%
  select(-c("text1", "text3", "test4"))

# Combine all the processed data into a single dataset
id_label_seq_meta <- left_join(id_label_meta, forward_size, by = c("sample_tag" = "tag")) %>%
  left_join(reverse_size, by = c("sample_tag" = "tag")) %>%
  left_join(map_coverage, by = c("sample_tag" = "tag")) %>%
  left_join(assembly_qc, by = c("sample_tag" = "tag")) %>%
  left_join(mlst, by = c("sample_tag" = "tag"))

# Write the combined dataset to a file
write_tsv(id_label_seq_meta, "Sminnesota_seq_meta_new.tsv")

# Read and process CARD genes data
card_genes <- read_tsv("./CARD_result_clincal.tsv") %>%
  select(c("Sample", "gene"))

# Pivot the CARD genes data to wider format
card_genes_pivot <- card_genes %>%
  pivot_wider(names_from = gene, values_from = gene, values_fn = length, values_fill = 0)

# Calculate the sum of AMR genes
colnms <- colnames(card_genes_pivot)
card_genes_pivot$AMR_sum <- rowSums(card_genes_pivot[, colnms[2:25]])
card_genes_pivot <- card_genes_pivot %>% relocate(AMR_sum)

# Combine CARD genes data with the previously combined dataset
seq_card_gene_meta <- left_join(id_label_seq_meta, card_genes_pivot, by = c("sample_tag" = "Sample"))

# Read and process VFDB genes data
vfdb_genes <- read_tsv("./VFDB_result_clincal.tsv") %>%
  select(c("Sample", "gene"))

# Pivot the VFDB genes data to wider format
vfdb_genes_pivot <- vfdb_genes %>%
  pivot_wider(names_from = gene, values_from = gene, values_fn = length, values_fill = 0)

# Read additional VFDB data for matching gene names
vfdb_name <- read_tsv("/Users/huanj0f/Downloads/srst_database/VFDB.txt")
vfdb_numname <- read_tsv("/Users/huanj0f/Downloads/srst_database/VFDB_am.txt")
vfdb_match <- data.frame(Name = vfdb_name$name, Number_name = vfdb_numname$num_name)

# Rename the columns in VFDB genes data
for (x in colnames(vfdb_genes_pivot[, -1])) {
  new_name <- vfdb_match$Name[vfdb_match$Number_name == x]
  colnames(vfdb_genes_pivot)[which(names(vfdb_genes_pivot) == x)] <- new_name
}

# Calculate the sum of VFDB genes
colnms <- colnames(vfdb_genes_pivot)
vfdb_genes_pivot$VF_sum <- rowSums(vfdb_genes_pivot[, colnms[2:178]])
vfdb_genes_pivot <- vfdb_genes_pivot %>% relocate(VF_sum)

# Combine VFDB genes data with the previously combined dataset
seq_card_vfdb_gene_meta <- left_join(seq_card_gene_meta, vfdb_genes_pivot, by = c("sample_tag" = "Sample"))

# Write the combined dataset to a file
write_csv(seq_card_vfdb_gene_meta, "Sminnesota_seq_meta_card_vfdb.csv")

# Read and process plasmid genes data
plasmid_genes <- read_tsv("./plasmidfinder_result.txt") %>%
  select(c("Sample", "gene"))

# Pivot the plasmid genes data to wider format
plasmid_genes_pivot <- plasmid_genes %>%
  pivot_wider(names_from = gene, values_from = gene, values_fn = length, values_fill = 0)

# Calculate the sum of plasmid genes
colnms <- colnames(plasmid_genes_pivot)
plasmid_genes_pivot$plasmid_sum <- rowSums(plasmid_genes_pivot[, colnms[2:31]])
plasmid_genes_pivot <- plasmid_genes_pivot %>% relocate(plasmid_sum)

# Combine plasmid genes data with the previously combined dataset
seq_card_vfdb_plasmid_gene_meta <- left_join(seq_card_vfdb_gene_meta, plasmid_genes_pivot, by = c("sample_tag" = "Sample"))

# Write the final combined dataset to a file
write_csv(seq_card_vfdb_plasmid_gene_meta, "Sminnesota_seq_meta_card_vfdb_plasmid.csv")




## Process the external sequencing metadata

# Read and process Entrobase data
entrobase <- read_tsv("./sminnesota_entrobase.tsv")
entrobase <- entrobase %>%
  separate(`Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Bases;Average Length;Status)`,
           into = c("Accession", "Sequencing_Platform", "Sequencing_Library", "Insert_Size", "Experiment", "Bases", "Average_Length", "Status"), 
           sep = ";", remove = TRUE) %>%
  select(-c("Sequencing_Platform", "Sequencing_Library", "Insert_Size"))

# Clean Entrobase data by filtering out unwanted rows and columns
entrobase_clean <- entrobase %>%
  filter(!Experiment == "NA") %>%
  filter(!Accession == "GCF_000486855") %>%
  filter(!Accession == "SRR1177262") %>%
  drop_na('Collection Year') %>%
  drop_na('Country')

# Write cleaned Entrobase data to a file
write_tsv(entrobase_clean, "./sminnesota_entrobase_clean.tsv")

# Generate a list of external projects from the cleaned Entrobase data
external_project_list <- entrobase_clean %>%
  group_by(`Bio Project ID`) %>%
  summarise() %>%
  rename(Bioproject_ID = `Bio Project ID`)

# Write the external project list to a file
write_tsv(external_project_list, "./external_project_list.tsv")

# Read and process NCBI data
pathogen_sra <- read_tsv("./PDG000000002.2847.metadata.tsv")

# Clean NCBI data by filtering out unwanted rows and columns
pathogen_sra_clean <- pathogen_sra %>%
  filter(computed_types == "serotype=Minnesota,antigen_formula=21:b:e,n,x") %>%
  filter(!Run %in% c("NULL")) %>%
  filter(!collection_date %in% c("NULL", "missing")) %>%
  filter(!geo_loc_name %in% c("NULL"))

# Combine cleaned NCBI data with cleaned Entrobase data
combined_ena_list <- unique(c(pathogen_sra_clean$Run, entrobase_clean$Accession))

# Write cleaned NCBI data to a file
write_tsv(pathogen_sra_clean, "./sminnesota_pathogen_clean.tsv")

# Read the cleaned Entrobase and NCBI data for further processing
entrobase_clean_2 <- read_tsv("./sminnesota_entrobase_clean.tsv")
pathogen_clean_2 <- read_tsv("./sminnesota_pathogen_clean.tsv")

# Further clean and process the Entrobase data
entrobase_clean_2 <- entrobase_clean_2 %>%
  mutate(collection_day = ifelse(is.na(collection_month) & is.na(collection_day), "30", collection_day),
         collection_month = ifelse(is.na(collection_month) & collection_day == "30", "6", collection_month)) %>%
  mutate(collection_day = ifelse(is.na(collection_day), "15", collection_day)) %>%
  unite(collection_date, collection_year, collection_month, collection_day, sep = "/", remove = TRUE)

# Further clean and process the NCBI data
pathogen_clean_2 <- pathogen_clean_2 %>%
  separate(geo_loc_name, c("Country", "rest"), sep = ":") %>%
  select(-c(rest))

# Combine the cleaned Entrobase and NCBI data into a single dataset
all_external <- rbind(entrobase_clean_2, pathogen_clean_2) %>%
  distinct(Run, .keep_all = TRUE)

# Write the combined dataset to a file
write_tsv(all_external, "./Sminnesota_external_metadata.tsv")















