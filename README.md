# Large-scale Genomic Survey of Non-typhoidal *Salmonella enterica* serovar Minnesota Strains in Chicken Products Reveals the Emergence of Multidrug Resistant Clones

## Abstract

### Background:
*Salmonella enterica* serovar Minnesota (*S. Minnesota*) is an emerging serovar of non-typhoidal Salmonella, known to persist in the food chain and distribution systems, potentially leading to outbreaks of Salmonella infections in human settings. Understanding the population dynamics and dynamics of genomic determinants is key to designing preventive measures and containing the spread of the pathogen in the poultry production chain.

### Methods:
In this study, we conducted a large-scale population-level study on *S. Minnesota* by fully characterizing population diversity and dynamics of a systematically consistent collection from the poultry production chain and one clinical strain from a patient with Salmonellosis in the Kingdom of Saudi Arabia. We sequenced 240 S. Minnesota strains from the western, eastern, and central regions of the country. We analyzed short-read and long-read sequencing data to decipher the population diversity and dynamics of resistance, as well as to contextualize the collections. Long-read sequencing was employed to retrieve and analyze the population diversity of plasmids carrying antimicrobial resistance and virulence factors.

### Results:
Our results indicate the rise of four clones (BAPS groups) in Saudi Arabia, three of which were mixed with global strains and originated from Brazil. The clones emerged over the past five to ten years and exhibited circulation between countries. The transmission analysis shows evidence of the spread of strains across cities, between countries, and mixing of strains from different suppliers. The clones with Saudi strains harbored a higher resistance and virulence level, owing to the acquisition of multiple plasmids, most importantly the IncC plasmid. The IncC plasmid was a mosaic plasmid, which carried antimicrobial resistance islands with *blaCMY-2*, ESBL *blaCTX-M*, aminoglycoside, and tetracycline resistance genes, as well as hyperpathogenicity island virulence island with yersiniabactin genes, providing evidence for the convergence of resistance and virulence. The plasmidome analysis revealed a high level of dynamics in the IncC plasmid structures with various configurations of resistance genes within the stains in each clone. The strains from humans showed high identity (10 SNPs) with the strains from food products, showing potential clinical importance of the strains.

### Conclusion:
Taken together, our results demonstrate a dynamic population and the emergence of multidrug-resistant clones for *Salmonella Minnesota*. They also highlight the variety of plasmids carrying antimicrobial resistance genes and the genomic contexts for resistance genes, as well as the genomic rearrangements in the antimicrobial resistance (AMR) regions in IncC plasmids. These changes occur in response to antimicrobial therapy in the poultry sector.


## Files Descriptions

| File & folder Name                  | Description                                                                                  |
| ----------------------------------- | ---------------------------------------------------------------------------------------------|
| `Annotated_plasmid`                 | Annotated plasmid sequence in gbk format                                                     |
| `plasmid_amr_vf_amrfinderplus.tsv`  | AMR genes and virulence factors detected in plasmid by AMRFinderPlus                         |
| `plasmid_amr_vf_blast.tsv`          | AMR genes and virulence factors detected in plasmid by BLASTN against CARD & VFDB database   |
| `External_samples_metadata.tsv`     | Metadata of external samples fetched from public database in this study                      |
| `Internal_samples_metadata.csv`     | Metadata of internal samples that were in-house sequenced in this study                      |
| `assembly_unicycler.sh`             | short-read genome assembly bash script                                                       |
| `assembly_unicycler_hybird.sh`      | Long-read hybrid genome assembly bash script                                                 |
| `assembly_qc_quast.sh`              | Genome assembly quality accession bash script                                                |
| `call_snp_sites.sh`                 | Bash script for calling SNP sites from multiple sequence alignment                           |
| `mapping_snippy.sh`                 | Bash script for mapping short reads against reference genome                                 |
| `mlst_mlst.sh`                      | In silico multi-locus sequence typing (MLST) bash script                                     |
| `run_gubbins.sh`                    | Bash script for filtering out polymorphic sites                                              |
| `run_beast2.sh`                     | Bash script for running Bayesian Evolutionary Analysis Sampling Trees 2 (BEAST2)             |
| `scan_amrfinder.sh`                 | Bash script for AMR gene & virulence factor detection by AMRFinderPlus                       |
| `typing_srst.sh`                    | Bash script for gene detection against the relevant database by srst2                        |
| `figure_plot.R`                     | R script for the relevant figure plot                                                        |
| `metadata_processing.R`             | R script for the relevant metadata processing                                                |
| `pvgenomic.ipynb`                   | Python script for the plasmid alignment and visualisation by pyGenomeViz                     |


| `K_158_3.gbff`     | Annotated plasmid sequence of K_158_3             |
| `K_160_2.gbff`     | Annotated plasmid sequence of K_160_2             |
| `K_195_2.gbff`     | Annotated plasmid sequence of K_195_2             |
| `K_229_3.gbff`     | Annotated plasmid sequence of K_229_3             |



















