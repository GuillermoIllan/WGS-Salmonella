# WGS-Salmonella
This report presents the results of the phylogenetic analysis conducted on presumptive *Salmonella* spp. isolates obtained from environmental samples. The objective of this analysis is to confirm the taxonomic identity of the isolates, characterize their genetic diversity, and establish their evolutionary relationships. These insights are essential for epidemiological surveillance and for understanding potential sources of contamination.

# Material and Methods

Raw sequencing reads were trimmed and filtered using **Fastp v0.24.3**, followed by quality assessment with **FastQC v0.12**. High-quality reads were then assembled using the pipeline implemented in **Shovill v1.1.0** (<https://github.com/tseemann/shovill>), a method previously applied in *Salmonella* spp. genomic studies.

To identify potential contamination, **Kraken2** was used to detect contigs originating from non-target organisms. Draft genome assemblies were evaluated for quality using **QUAST v5.3.0** and **CheckM2**.

**Multilocus sequence typing (MLST)** was performed using the Achtman 7-gene scheme through the **mlst** tool(<https://github.com/tseemann/mlst>). **In silico** prediction of *Salmonella* serovar and core genome MLST (cgMLST) type was conducted with **SISTR v1.1.1**, a method considered one of the most reliable for *Salmonella* subtyping.

To identify clusters for high-quality SNP (hqSNP) analysis, a reference-free SNP analysis was initially performed using **kSNP4**. Isolate pairs showing fewer than 100 SNP differences were subsequently analyzed using the **CFSAN SNP Pipeline**, with the reference assembly selected as previously described. The reference assembly for the CFSAN SNP pipeline was selected as previously described .

Maximum likelihood phylogenetic trees based on hqSNP data were constructed using **RAxML v8.2.13**, employing the GTRGAMMA nucleotide substitution model and 1,000 bootstrap replicates. An SNP matrix was generated using **snp-sites** and **snp-dists** (<https://github.com/tseemann/snp-dists>) to visualize SNP variation among isolates.
