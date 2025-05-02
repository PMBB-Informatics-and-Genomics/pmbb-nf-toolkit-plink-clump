
Documentation for PLINK Clump
=============================

# Module Overview


Plink clumping is a procedure used to determine which variants belong to a GWAS signal via LD. It starts by identifying the most significant SNP and designates any variants in LD with that lead SNP as part of its “clump.” We have added additional functionality to return “loci” made up of clumps that are physically overlapping.

[Paper Link for Reference](https://academic.oup.com/gigascience/article/4/1/s13742-015-0047-8/2707533)

[Tool Documentation Link](https://www.cog-genomics.org/plink/2.0/postproc#clump)

[Example Module Config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/plink_clump.config)

[Example nextflow.config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/nextflow.config)
## Cloning Github Repository:


* Command: `git clone https://github.com/PMBB-Informatics-and-Genomics/geno_pheno_workbench.git`

* Navigate to relevant workflow directory run commands (our pipelines assume all of the nextflow files/scripts are in the current working directory)
## Software Requirements:


* [Nextflow version 23.04.1.5866](https://www.nextflow.io/docs/latest/cli.html)

* [Singularity version 3.8.3](https://sylabs.io/docs/) OR [Docker version 4.30.0](https://docs.docker.com/)
## Commands for Running the Workflow


* Singularity Command: `singularity build plink_clump.sif docker://guarelin/plink_clump:latest`

* Docker Command: `docker pull guarelin/plink_clump:latest`

* Command to Pull from Google Container Registry: `docker pull guarelin/plink_clump:latest`

* Run Command: `nextflow run /path/to/toolkit/module/plink_clump.nf`

* Common `nextflow run` flags:

    * `-resume` flag picks up the workflow where it left off, otherwise, the workflow will rerun from the beginning

    * `-stub` performs a sort of dry run of the whole workflow, checks channels without executing any code

    * `-profile` selects the compute profiles we set up in nextflow.config (see nextflow.config file below)

    * `-profile` selects the compute profiles we set up in nextflow.config (see nextflow.config file below)

    * `-profile standard` uses the docker image to executes the processes

    * `-profile cluster` uses the singularity container and submits processes to a queue- optimal for HPC or LPC computing systems

    * `-profile all_of_us` uses the docker image to execute pipelines on the All of Us Researcher Workbench

* for more information visit the [Nextflow documentation](https://www.nextflow.io/docs/latest/cli.html)
# Configuration Parameters and Input File Descriptions

## Workflow


* `my_python` (Type: File Path)

    * Path to the python executable to be used for python scripts - often it comes from the docker/singularity container (/opt/conda/bin/python)
## Pre-Processing


* `input_colnames` (Type: Map (Dictionary))

    * This map indicates the column names in the GWAS summary stats. All summary stats files should have the same columns. The keys required for this map are chr_col, pos_col, id_col, p_col, A1_col, and A2_col.

    * Corresponding Input File: GWAS Summary Stats (specified with input table)

        * GWAS summary statistics - required columns are chromosome, position, alleles, and p-value. Column names are specified in the parameters.

        * Type: Summary Statistics

        * Format: tsv.gz

        * Input File Header:





        ```
        #CHROM  BP      A1      A2      BETA    OR      SE      P       A1_FREQ N       ID_CHR:POS      ID_CHR:POS:REF:ALT
        1       1       C       A       -0.1119171733535359     0.8941183112547663      0.032313165627030196    0.0019817724363402917   0.12782190647715028     10105   chr1:1  chr1:1:A:C
        1       2       C       A       -0.0864716424473868     0.9171615568730449      0.12186030955155278     0.6202982240603924      0.42131284017764215     10105   chr1:2  chr1:2:A:C
        1       3       C       A       -0.20796068564107204    0.8122389687385936      0.15240824661729327     0.3145166720046993      0.2654984327348345      10105   chr1:3  chr1:3:A:C
        1       4       C       A       -0.0010440676231661449  0.9989564772257986      0.0014464012295794873   0.6148871476723993      0.4516872588226362      10105   chr1:4  chr1:4:A:C
        
        ```

* `input_descriptor_table_filename (Plink Clump)` (Type: File Path)

    * Path to a .csv file with rows describing each set of GWAS summary statistics to be clumped. The columns (in order) should be: analysis,file,ref_prefix,ref_suffix,ref_plink_flag. 

    * Corresponding Input File: GWAS Summary Stats (specified with input table)

        * GWAS summary statistics - required columns are chromosome, position, alleles, and p-value. Column names are specified in the parameters.

        * Type: Summary Statistics

        * Format: tsv.gz

        * Input File Header:





        ```
        #CHROM  BP      A1      A2      BETA    OR      SE      P       A1_FREQ N       ID_CHR:POS      ID_CHR:POS:REF:ALT
        1       1       C       A       -0.1119171733535359     0.8941183112547663      0.032313165627030196    0.0019817724363402917   0.12782190647715028     10105   chr1:1  chr1:1:A:C
        1       2       C       A       -0.0864716424473868     0.9171615568730449      0.12186030955155278     0.6202982240603924      0.42131284017764215     10105   chr1:2  chr1:2:A:C
        1       3       C       A       -0.20796068564107204    0.8122389687385936      0.15240824661729327     0.3145166720046993      0.2654984327348345      10105   chr1:3  chr1:3:A:C
        1       4       C       A       -0.0010440676231661449  0.9989564772257986      0.0014464012295794873   0.6148871476723993      0.4516872588226362      10105   chr1:4  chr1:4:A:C
        
        ```

    * Corresponding Input File: PLINK Clump Input Descriptor Table

        * The ‘analysis’ column should be a nickname for that set of summary stats (often a population and a phenotype). The ‘file’ columns should contain the path to that set of summary stats. The columns ‘ref_prefix’ and ‘ref_suffix’ should correspond to the chromosome-separated plink reference files you want to use. Each plink prefix will be computed as: ref_prefix{chromosome}ref_suffix. Finally, the ‘ref_plink_flag’ columns should have either --bfile for plink 1.9 inputs and --pfile for plink 2.0 inputs.

        * Type: Data Table

        * Format: csv

        * Input File Header:





        ```
        analysis,file,ref_prefix,ref_suffix,ref_plink_flag
        POP1_Q1,${launchDir}/Input/POP1.quant_PHENO1.sumstats.sim.txt.gz,${launchDir}/Input/genotype_100markers_2chr.chr,,--bfile
        POP1_B1,${launchDir}/Input/POP1.bin_PHENO1.sumstats.sim.txt.gz,${launchDir}/Input/genotype_100markers_2chr.chr,,--pfile
        POP1_Q2,${launchDir}/Input/POP1.quant_PHENO2.sumstats.sim.txt.gz,${launchDir}/Input/genotype_100markers_2chr.chr,,--pfile
        POP1_B2,${launchDir}/Input/POP1.bin_PHENO2.sumstats.sim.txt.gz,${launchDir}/Input/genotype_100markers_2chr.chr,,--bfile
        ```

    * Corresponding Input File: Chr-Separated Plink Files (specified with input table)

        * Chromosome-separated plink file sets

        * Type: Plink Set

        * Format: plink binary
## PLINK


* `clump_options` (Type: Map (Dictionary))

    * The four required keys for this map are clump_p1, clump_p2, clump_r2, and clump_kb.  The PLINK docs have more details, but briefly: clump_p1 lead SNP p-value threshold, clump_p2 is the threshold for other SNPs in the clump, clump_r2 is the min R^2 for being clumped with a lead SNP, and clump_kb is the max kb distance for being clumped with a lead SNP

* `my_plink2` (Type: File Path)

    * Path to the PLINK2 executable to be used for PLINK2 score - often it comes from the docker or singularity container (plink2) (is on the path in the container
## Post-Processing


* `biofilter_close_dist` (Type: Float)

    * The distance in bp for something to be considered “close” vs “far” with respect to nearest gene annotation. Value is often 5E4

* `biofilter_script` (Type: File Path)

    * The path to the biofilter script to use. If using the singularity container, should be ‘/app/biofilter.py’

* `biofilter_loki` (Type: File Path)

    * The path to a loki.db file to be used for nearest gene annotation

* `biofilter_build` (Type: String)

    * The build to pass to biofilter - can be 19 or 38

* `annotate` (Type: Bool (Java: true or false))

    * Whether or not to annotate results with the RSIDs and nearest genes for plotting and summary files.
# Output Files from PLINK_Clump


* Loci for each Analysis

    * `Clump/Results/{analysis}.loci.csv`

    * Merged “loci” for one set of GWAS summary statistics. A “locus” is one or more clumps which were combined because they physically overlap.

    * Type: Data Table

    * Format: csv

    * Parallel By: Analysis

* Clumps for each Analysis

    * `Clump/Results/{analysis}.clumps.csv`

    * Merged clumps for one set of GWAS summary statistics

    * Type: Data Table

    * Format: csv

    * Parallel By: Analysis

* Multi-locus Manhattan plot

    * `Plots/{analysis}.loci.png`

    * For each analysis, a plot is generated that zooms in on each locus and presents it in a vertical Manhattan plot

    * Type: Manhattan Plot

    * Format: png

    * Parallel By: Analysis

* Multi-clump Manhattan plot

    * `Plots/{analysis}.clumps.png`

    * For each analysis, a plot is generated that zooms in on each clump and presents it in a vertical Manhattan plot

    * Type: Manhattan Plot

    * Format: png

    * Parallel By: Analysis

* Locus Summary

    * `Summary/all_loci_with_annot.csv`

    * Merged “loci” across all summary stats. A “locus” is one or more clumps which were combined because they physically overlap.

    * Type: Summary Table

    * Format: csv

    * Output File Header:





    ```
    Lead_SNP,#CHROM,POS,P,TOTAL,NONSIG,S0.05,S0.01,S0.001,S0.0001,SP2,MIN_POS,MAX_POS,ANALYSIS,Lead_SNP_Nearest_Gene,Lead_SNP_RSID
    chr15:57726260:C:T,15,57726260,2.02e-08,5,0,4,1,0,0,"chr15:57729010:G:A,chr15:57734220:C:G,chr15:57734776:A:G,chr15:57734983:T:C,chr15:57741156:G:A",57726260.0,57741156.0,AFR_Endo,POLR2M/LOC105370834,
    chr1:22139327:T:C,1,22139327,2.88e-18,103,0,2,6,4,91,"chr1:22020377:G:C,chr1:22020903:A:G,chr1:22021916:G:A,chr1:22022094:A:C,chr1:22022893:C:T,chr1:22023122:A:G,chr1:22024341:C:T,chr1:22025547:T:G,ch
    chr2:11591720:C:T,2,11591720,1.7e-10,55,4,2,6,10,33,"chr2:11562535:A:C,chr2:11564710:G:T,chr2:11571949:G:T,chr2:11578465:T:C,chr2:11578732:T:C,chr2:11581409:T:G,chr2:11581886:C:G,chr2:11581956:T:C,chr
    chr2:49092263:T:C,2,49092263,9.37e-10,143,0,3,1,15,124,"chr2:49018269:C:T,chr2:49022873:C:T,chr2:49032563:G:T,chr2:49033096:G:A,chr2:49035532:A:C,chr2:49038693:C:A,chr2:49039715:T:C,chr2:49041339:C:G,
    
    ```

* Clump Summary

    * `Summary/all_clumps_with_annot.csv`

    * Merged clumps across all summary stat / analyses. These can also be annotated with the RSID of the lead SNP and nearest gene.

    * Type: Summary Table

    * Format: csv

    * Output File Header:





    ```
    Lead_SNP,#CHROM,POS,P,TOTAL,NONSIG,S0.05,S0.01,S0.001,S0.0001,SP2,MIN_POS,MAX_POS,ANALYSIS,Lead_SNP_Nearest_Gene,Lead_SNP_RSID
    chr15:57726260:C:T,15,57726260,2.02e-08,5,0,4,1,0,0,"chr15:57729010:G:A,chr15:57734220:C:G,chr15:57734776:A:G,chr15:57734983:T:C,chr15:57741156:G:A",57726260.0,57741156.0,AFR_Endo,POLR2M/LOC105370834,
    chr1:22139327:T:C,1,22139327,2.88e-18,103,0,2,6,4,91,"chr1:22020377:G:C,chr1:22020903:A:G,chr1:22021916:G:A,chr1:22022094:A:C,chr1:22022893:C:T,chr1:22023122:A:G,chr1:22024341:C:T,chr1:22025547:T:G,ch
    chr2:11591720:C:T,2,11591720,1.7e-10,55,4,2,6,10,33,"chr2:11562535:A:C,chr2:11564710:G:T,chr2:11571949:G:T,chr2:11578465:T:C,chr2:11578732:T:C,chr2:11581409:T:G,chr2:11581886:C:G,chr2:11581956:T:C,chr
    chr2:49092263:T:C,2,49092263,9.37e-10,143,0,3,1,15,124,"chr2:49018269:C:T,chr2:49022873:C:T,chr2:49032563:G:T,chr2:49033096:G:A,chr2:49035532:A:C,chr2:49038693:C:A,chr2:49039715:T:C,chr2:49041339:C:G,
    
    ```
# Current Dockerfile for the Container/Image


```docker
FROM continuumio/miniconda3

WORKDIR /app

# biofilter version argument
ARG BIOFILTER_VERSION=2.4.3

RUN apt-get update \
    # install packages needed for PLINK, NEAT plots and biofilter installation
    && apt-get install -y --no-install-recommends libz-dev g++ gcc git wget tar unzip make \
    && rm -rf /var/lib/apt/lists/* \
    # install PLINK
    && wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_amd_avx2_20250129.zip \
    && unzip plink2_linux_amd_avx2_20250129.zip \
    && rm -rf plink2_linux_amd_avx2_20250129.zip \
    # move plink2 executable to $PATH
    && mv plink2 /usr/bin \
    # install python packages needed for pipeline
    && conda install -y -n base -c conda-forge python=3.10 scikit-learn dominate wget libtiff conda-build scipy pandas seaborn matplotlib numpy apsw sqlite \
    && conda clean --all --yes \
    # install biofilter
    && wget https://github.com/RitchieLab/biofilter/releases/download/Biofilter-${BIOFILTER_VERSION}/biofilter-${BIOFILTER_VERSION}.tar.gz -O biofilter.tar.gz \
    && tar -zxvf biofilter.tar.gz --strip-components=1 -C /app/ \
    && /opt/conda/bin/python setup.py install \
    && chmod a+rx /app/biofilter.py \
    # install NEAT plots
    && git clone https://github.com/PMBB-Informatics-and-Genomics/NEAT-Plots.git \
    && mv NEAT-Plots/manhattan-plot/ /app/ \
    && conda develop /app/manhattan-plot/ \
    # remove NEAT-plots directory and biofilter tarball
    && rm -R NEAT-Plots biofilter.tar.gz

USER root
```