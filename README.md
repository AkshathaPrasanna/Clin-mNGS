# Clin-mNGS : Automated pipeline for pathogen detection from clinical metagenomic data


The Clin-mNGS pipeline is an integrated, open-source, scalable, reproducible, and user-friendly framework scripted using
the Snakemake workflow management software. The implementation avoids the hassle of manual installation and configuration of the multiple
command-line tools and dependencies and can be deployed on most Linux workstations and clusters. The versions of the implemented tools are made user modifiable. The approach directly screens pathogens from clinical raw reads and generates consolidated reports for each sample.

The pipeline is currently automated to perform in quality check, filtering, host subtraction, assembly of reads into contigs, assembly metrics, relative abundances of bacterial species, antimicrobial resistance genes, plasmid finding, and virulence factors identification. 


# Quick Start

Commands to start analysing your data:

    conda create -c conda-forge -c bioconda -n snakemake snakemake      
    conda activate snakemake
    snakemake --use-conda --cores 16
    
    
# Citation

DOI : 10.2174/1574893615999200608130029
    
    



