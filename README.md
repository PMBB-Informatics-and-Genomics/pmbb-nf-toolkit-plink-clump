
Documentation for PLINK Clump
=============================

# Module Overview


Plink clumping is a procedure used to determine which variants belong to a GWAS signal via LD. It starts by identifying the most significant SNP and designates any variants in LD with that lead SNP as part of its “clump.” We have added additional functionality to return “loci” made up of clumps that are physically overlapping.
- [Tool Paper Link for Reference](https://academic.oup.com/gigascience/article/4/1/s13742-015-0047-8/2707533)
- [Tool Documentation Link for Reference](https://www.cog-genomics.org/plink/2.0/postproc#clump)
- [Example Config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/plink_clump.config)
- [Example nextflow.config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/nextflow.config)

## Software Requirements


* [Nextflow version 24.04.3](https://www.nextflow.io/docs/latest/cli.html)

* [Singularity 3.8.3](https://sylabs.io/docs/) OR [Docker 4.30.0](https://docs.docker.com/)
## Commands for Running the Workflow


* Singularity Command: `singularity build plink_clump.sif docker://guarelin/plink_clump:latest`

* Docker Command: `docker pull guarelin/plink_clump:latest`

* Pull from Google Container Registry: `docker pull guarelin/plink_clump:latest`

* Run Command: `nextflow run /path/to/toolkit/module/plink_clump.nf`

* Common `nextflow run` flags:

    * `-resume` flag picks up workflow where it left off

    * `-stub` performs a dry run, checks channels without executing code

    * `-profile` selects the compute profiles in nextflow.config

    * `-profile standard` uses the Docker image to execute processes

    * `-profile cluster` uses the Singularity container and submits processes to a queue

    * `-profile all_of_us` uses the Docker image on All of Us Workbench

* More info: [Nextflow documentation](https://www.nextflow.io/docs/latest/cli.html)
# Detailed Pipeline Steps

## Part I: Setup


1. Start your own tools directory and go there. You may do this in your project analysis directory, but it often makes sense to clone into a general `tools` location

```sh
# Make a directory to clone the pipeline into
TOOLS_DIR="/path/to/tools/directory"
mkdir $TOOLS_DIR
cd $TOOLS_DIR
```

2. Download the source code by cloning from git

```sh
git clone None
cd $TOOLS_DIR/pmbb-nf-toolkit-plink-clump
```

3. Build the singularity image
    - you may call the image whatever you like, and store it wherever you like. Just make sure you specify the name in `nextflow.conf`
    - this does NOT have to be done for every saige-based analysis, but it is good practice to re-build every so often as we update regularly.


```sh
cd $TOOLS_DIR/pmbb-nf-toolkit-plink-clump
singularity build plink_clump.sif docker://guarelin/plink_clump:latest
```
## Part II: Configure your run


1. Make a separate analysis/run/working directory.
    - The quickest way to get started, is to run the analysis in the folder the pipeline is run. However, subsequent analyses will over-write results from previous analyses.
    - ❗This step is optional, but We Highly recommend making a `tools` directory separate from your `run` directory. We recommend storing the `nextflow.conf` in here as it shouldn't change between runs.


```sh
WDIR="/path/to/analysis/run1"
mkdir -p $WDIR
cd $WDIR
```

2. Fill out the `nextflow.config` file for your system.
    - See [Nextflow configuration documentation](https://www.nextflow.io/docs/latest/config.html) for information on how to configure this file. An example can be found on our GitHub: [Nextflow Config](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/blob/main/Example_Configs/nextflow.config).
    - ❗IMPORTANTLY, you must configure a user-defined profile for your run environments (local, docker, saige, cluster, etc.). If multiple profiles are specified, run with a specific profile using `nextflow run -profile $MY_PROFILE`.
    - For singularity, The profile's attribute `process.container` should be set to `'/path/to/plink_clump.sif'` (replace `/path/to` with the location where you built the image above). See [Nextflow Executor Information](https://www.nextflow.io/docs/latest/executor.html) for more details.
    - ⚠️As this file remains mostly unchanged for your system, We recommend storing this file in the `tools/pipeline` directory and passing it to the pipeline with `-c /path/to/nextflow.config`.


3. Create a pipeline-specific `.config` file specifying your run parameters and input files. See Below for workflow-specific parameters and what they mean.
    - Everything in here can be configured in `nextflow.config`, however we find it easier to separate the system-level profiles from the individual run parameters.
    - Examples can be found in our Pipeline-Specific [Example Config Files](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs).
    - you can compartamentalize your config file as much as you like by passing
    - There are 2 ways to specify the config file during a run:

        - with the `-c` option on the command line: `nextflow run -c Plink_Clump/plink_clump.config`
        - in the `nextflow.config`: at the top of the file add: `includeConfig Plink_Clump/plink_clump.config`

## Part III: Run your analysis


❗We HIGHLY recommend doing a STUB run to test the analysis using the `-stub` flag. This is a dry run to make sure your environment, parameters, and input_files are specified and formatted correctly.❗We also HIGHLY recommend doing a TEST run with the included test data in `$TOOLS_DIR/pmbb-nf-toolkit-plink-clump/test_data`we have several pre-configured analyses runs with input data and fully-specified config files.

```sh
# run an exwas stub
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-plink-clump/plink_clump.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c Plink_Clump/plink_clump.config \
   -stub

# run an exwas for real
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-plink-clump/plink_clump.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c Plink_Clump/plink_clump.config

# resume an exwas run if it was interrupted or ran into an error
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-plink-clump/plink_clump.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c Plink_Clump/plink_clump.config \
   -resume
```
# Pipeline Parameters

## Input Files for PLINK_Clump


* GWAS Summary Stats (specified with input table)

    * GWAS summary statistics - required columns are chromosome, position, alleles, and p-value. Column names are specified in the parameters.

    * Type: Summary Statistics

    * Format: tsv.gz

    * File Header:


    ```
    #CHROM  BP      A1      A2      BETA    OR      SE      P       A1_FREQ N       ID_CHR:POS      ID_CHR:POS:REF:ALT
    1       1       C       A       -0.1119171733535359     0.8941183112547663      0.032313165627030196    0.0019817724363402917   0.12782190647715028     10105   chr1:1  chr1:1:A:C
    1       2       C       A       -0.0864716424473868     0.9171615568730449      0.12186030955155278     0.6202982240603924      0.42131284017764215     10105   chr1:2  chr1:2:A:C
    1       3       C       A       -0.20796068564107204    0.8122389687385936      0.15240824661729327     0.3145166720046993      0.2654984327348345      10105   chr1:3  chr1:3:A:C
    1       4       C       A       -0.0010440676231661449  0.9989564772257986      0.0014464012295794873   0.6148871476723993      0.4516872588226362      10105   chr1:4  chr1:4:A:C
    
    ```

* PLINK Clump Input Descriptor Table

    * The ‘analysis’ column should be a nickname for that set of summary stats (often a population and a phenotype). The ‘file’ columns should contain the path to that set of summary stats. The columns ‘ref_prefix’ and ‘ref_suffix’ should correspond to the chromosome-separated plink reference files you want to use. Each plink prefix will be computed as: ref_prefix{chromosome}ref_suffix. Finally, the ‘ref_plink_flag’ columns should have either --bfile for plink 1.9 inputs and --pfile for plink 2.0 inputs.

    * Type: Data Table

    * Format: csv

    * File Header:


    ```
    analysis,file,ref_prefix,ref_suffix,ref_plink_flag
    POP1_Q1,${launchDir}/Input/POP1.quant_PHENO1.sumstats.sim.txt.gz,${launchDir}/Input/genotype_100markers_2chr.chr,,--bfile
    POP1_B1,${launchDir}/Input/POP1.bin_PHENO1.sumstats.sim.txt.gz,${launchDir}/Input/genotype_100markers_2chr.chr,,--pfile
    POP1_Q2,${launchDir}/Input/POP1.quant_PHENO2.sumstats.sim.txt.gz,${launchDir}/Input/genotype_100markers_2chr.chr,,--pfile
    POP1_B2,${launchDir}/Input/POP1.bin_PHENO2.sumstats.sim.txt.gz,${launchDir}/Input/genotype_100markers_2chr.chr,,--bfile
    ```

* Chr-Separated Plink Files (specified with input table)

    * Chromosome-separated plink file sets

    * Type: Plink Set

    * Format: plink binary
## Output Files for PLINK_Clump


* Clumps for each Analysis

    * `Clump/Results/{analysis}.clumps.csv`

    * Merged clumps for one set of GWAS summary statistics

    * Type: Data Table

    * Format: csv

        * Parallel By: Analysis

* Locus Summary

    * `Summary/all_loci_with_annot.csv`

    * Merged “loci” across all summary stats. A “locus” is one or more clumps which were combined because they physically overlap.

    * Type: Summary Table

    * Format: csv

    * File Header:


    ```
    Lead_SNP,#CHROM,POS,P,TOTAL,NONSIG,S0.05,S0.01,S0.001,S0.0001,SP2,MIN_POS,MAX_POS,ANALYSIS,Lead_SNP_Nearest_Gene,Lead_SNP_RSID
    chr15:57726260:C:T,15,57726260,2.02e-08,5,0,4,1,0,0,"chr15:57729010:G:A,chr15:57734220:C:G,chr15:57734776:A:G,chr15:57734983:T:C,chr15:57741156:G:A",57726260.0,57741156.0,AFR_Endo,POLR2M/LOC105370834,
    chr1:22139327:T:C,1,22139327,2.88e-18,103,0,2,6,4,91,"chr1:22020377:G:C,chr1:22020903:A:G,chr1:22021916:G:A,chr1:22022094:A:C,chr1:22022893:C:T,chr1:22023122:A:G,chr1:22024341:C:T,chr1:22025547:T:G,ch
    chr2:11591720:C:T,2,11591720,1.7e-10,55,4,2,6,10,33,"chr2:11562535:A:C,chr2:11564710:G:T,chr2:11571949:G:T,chr2:11578465:T:C,chr2:11578732:T:C,chr2:11581409:T:G,chr2:11581886:C:G,chr2:11581956:T:C,chr
    chr2:49092263:T:C,2,49092263,9.37e-10,143,0,3,1,15,124,"chr2:49018269:C:T,chr2:49022873:C:T,chr2:49032563:G:T,chr2:49033096:G:A,chr2:49035532:A:C,chr2:49038693:C:A,chr2:49039715:T:C,chr2:49041339:C:G,
    
    ```

* Multi-clump Manhattan plot

    * `Plots/{analysis}.clumps.png`

    * For each analysis, a plot is generated that zooms in on each clump and presents it in a vertical Manhattan plot

    * Type: Manhattan Plot

    * Format: png

        * Parallel By: Analysis

* Clump Summary

    * `Summary/all_clumps_with_annot.csv`

    * Merged clumps across all summary stat / analyses. These can also be annotated with the RSID of the lead SNP and nearest gene.

    * Type: Summary Table

    * Format: csv

    * File Header:


    ```
    Lead_SNP,#CHROM,POS,P,TOTAL,NONSIG,S0.05,S0.01,S0.001,S0.0001,SP2,MIN_POS,MAX_POS,ANALYSIS,Lead_SNP_Nearest_Gene,Lead_SNP_RSID
    chr15:57726260:C:T,15,57726260,2.02e-08,5,0,4,1,0,0,"chr15:57729010:G:A,chr15:57734220:C:G,chr15:57734776:A:G,chr15:57734983:T:C,chr15:57741156:G:A",57726260.0,57741156.0,AFR_Endo,POLR2M/LOC105370834,
    chr1:22139327:T:C,1,22139327,2.88e-18,103,0,2,6,4,91,"chr1:22020377:G:C,chr1:22020903:A:G,chr1:22021916:G:A,chr1:22022094:A:C,chr1:22022893:C:T,chr1:22023122:A:G,chr1:22024341:C:T,chr1:22025547:T:G,ch
    chr2:11591720:C:T,2,11591720,1.7e-10,55,4,2,6,10,33,"chr2:11562535:A:C,chr2:11564710:G:T,chr2:11571949:G:T,chr2:11578465:T:C,chr2:11578732:T:C,chr2:11581409:T:G,chr2:11581886:C:G,chr2:11581956:T:C,chr
    chr2:49092263:T:C,2,49092263,9.37e-10,143,0,3,1,15,124,"chr2:49018269:C:T,chr2:49022873:C:T,chr2:49032563:G:T,chr2:49033096:G:A,chr2:49035532:A:C,chr2:49038693:C:A,chr2:49039715:T:C,chr2:49041339:C:G,
    
    ```

* Multi-locus Manhattan plot

    * `Plots/{analysis}.loci.png`

    * For each analysis, a plot is generated that zooms in on each locus and presents it in a vertical Manhattan plot

    * Type: Manhattan Plot

    * Format: png

        * Parallel By: Analysis

* Loci for each Analysis

    * `Clump/Results/{analysis}.loci.csv`

    * Merged “loci” for one set of GWAS summary statistics. A “locus” is one or more clumps which were combined because they physically overlap.

    * Type: Data Table

    * Format: csv

        * Parallel By: Analysis
## Other Parameters for PLINK_Clump

### PLINK


* `my_plink2` (Type: File Path)

    * Path to the PLINK2 executable to be used for PLINK2 score - often it comes from the docker or singularity container (plink2) (is on the path in the container

* `clump_options` (Type: Map (Dictionary))

    * The four required keys for this map are clump_p1, clump_p2, clump_r2, and clump_kb.  The PLINK docs have more details, but briefly: clump_p1 lead SNP p-value threshold, clump_p2 is the threshold for other SNPs in the clump, clump_r2 is the min R^2 for being clumped with a lead SNP, and clump_kb is the max kb distance for being clumped with a lead SNP
### Post-Processing


* `biofilter_close_dist` (Type: Float)

    * The distance in bp for something to be considered “close” vs “far” with respect to nearest gene annotation. Value is often 5E4
### Pre-Processing


* `input_colnames` (Type: Map (Dictionary))

    * This map indicates the column names in the GWAS summary stats. All summary stats files should have the same columns. The keys required for this map are chr_col, pos_col, id_col, p_col, A1_col, and A2_col.

    * Corresponding Input File: GWAS Summary Stats (specified with input table)

        * GWAS summary statistics - required columns are chromosome, position, alleles, and p-value. Column names are specified in the parameters.

        * Type: Summary Statistics

        * Format: tsv.gz

        * File Header:


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

        * File Header:


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

        * File Header:


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
### Workflow


* `my_python` (Type: File Path)

    * Path to the python executable to be used for python scripts - often it comes from the docker/singularity container (/opt/conda/bin/python)
# Configuration and Advanced Workflow Files

## Example Config File Contents (From Path)


```
params {
    // python version (default is for the container)
    my_python = '/opt/conda/bin/python'

    // Must be plink 2.0
    // my_plink2 = '/usr/bin/plink2'
    my_plink2 = '/appl/plink2-20240804/plink2'

    // chromsome list
    chromosome_list = [6, 12]
    // chromosome_list = 1..22

    // plink clump parameters (PLINK docs have more details)
    // clump_p1 lead SNP p-value threshold
    // clump_p2 is the threshold for other SNPs in the clump
    // clump_r2 is the min R^2 for being clumped with a lead SNP
    // clump_kb is the max kb distance for being clumped with a lead SNP
    clump_options = [
        clump_p1: 0.00000005,
        clump_p2: 1,
        clump_r2: 0.25,
        clump_kb: 1000
    ]

    // Each summary stats file gets an analysis nickname
    // (the key such as Meta_AFR or Meta_ALL or GWAS_EUR)
    // Columns of this .csv include:
    // analysis, file, ref_prefix, ref_suffix, ref_plink_flag
    input_descriptor_table_filename = 'sumstats_manifest.csv'

    // Input column names
    // All must match in your summary stats provided
    // chr_col, pos_col are for coordinates
    // id_col is for the variant IDs
    // p_col is for p values
    // A1 is the effect allele and A2 is the other allele
    input_colnames = [
        'chr_col': 'CHR',
        'pos_col': 'POS',
        'id_col': 'variant_id',
        'p_col': 'p-value',
        'A1_col': 'reference_allele',
        'A2_col': 'other_allele'
    ]
    
    // parameters for getting RSIDs and nearest genes
    annotate = true // whether or not to annotate with biofilter
    biofilter_build = '38' // can be 19 or 38
    biofilter_loki = '/path/to/data/loki.db'
    biofilter_script = '/home/guarelin/mambaforge/envs/py38/bin/biofilter.py' // Must be an executable python file
    biofilter_close_dist = 5E4 // How "close" to a gene a variant should be to be considered "close" vs "far"
}
```
## Current `nextflow.config` contents


```
includeConfig 'plink_clump.config'

profiles {
    non_docker_dev {
        process.executor = awsbatch-or-lsf-or-slurm-etc
        process.queue = 'epistasis_normal'
        process.memory = '15GB'
    }

    standard {
        process.executor = awsbatch-or-lsf-or-slurm-etc
    }

    cluster {
        process.executor = awsbatch-or-lsf-or-slurm-etc
        process.queue = 'epistasis_normal'
        process.memory = '15GB'
        process.container = 'plink_clump.sif'
        singularity.enabled = true
        singularity.runOptions = '-B /root/,/directory/,/names/'
    }
}
```