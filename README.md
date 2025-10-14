# GEMap_ms
This repository contains code associated with the manuscript describing GEMap, which is a method capable of quantifying the differential dependence of stratified subgroups (such as KRAS-mutant cancers), therapeutic target prioritization, tissue-specific targeting, and drug candidate discovery. GEMap mainly includes four modules: Essentiality Analysis, Essentiality2Target, Essentiality2Drug and DIY Tools.

# Data
GEMap consolidates CRISPR screening data from Broad Achilles 23Q2, Sanger, Avana, GeCKO and supplements these with manually curated literature-derived datasets from PubMed. Drug data were sourced from LINCS 2020 and DGIdb datasets，including nearly 20,000 drug entries, enabling users to easily identify candidate drugs for genes of interest.
Processing was applied to files larger than 25 MB (CRISPR_DATA.Rdata, Mut_Data.rds, drug_data_list_20.RData, Avana_DATA.Rdata, Sanger_DATA.Rdata). These data only show a portion of the information to ensure the transparency and reproducibility.

# Essentiality Analysis
The main scripts are concentrated in 01-08. The Essentiality Analysis module facilitates exploration of multimodal gene data across diverse cell lines and tissue types, including Dependency, Expression, CNV and Mutation.

# Essentiality2Target
The main scripts are concentrated in 09-13. Essentiality2Target prioritizes context-specific therapeutic targets by quantifying differential genetic dependencies in molecularly stratified cancers, such KRAS mutation-cancers. The Pancancer section of the Essentiality2Target module enables users to predict and prioritize therapeutic targets based on the mutation or copy number variation (CNV) of a gene of interest. The Tissue section was able to identify which cancer type with KRAS mutation exhibits greater sensitivity to this therapeutic target.

# Essentiality2Drug
The main scripts are concentrated in 14-15. The Essentiality2Drug module facilitates drug discovery by allowing users to identify potential therapeutic compounds through two key drug-target relationships:  Drug Perturbational Profiles (DPP) and Drug Target Interactions (DTI). Discover tailored treatments and drugs candidate targeting particular genes using large-scale drug perturbation data and drug physical targets.

# DIY Tools
The main scripts are concentrated in 16-20. DIY Tools can be categorized into Correlation and Visualize Your Way. The Correlation includes: Co-Dependency, Dependency/Expression, Dependency/CNV, Dependency/Mutation. It has been reported genes with correlated knockout fitness profiles tend to operate in the same biological process or pathway. Additionally, the Visualize You Way module can assess both dependency and expression patterns across multiple genes.
