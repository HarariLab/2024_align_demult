# Snakemake workflow: `align-and-demultiplex`

A Snakemake workflow for Multiome ATAC + Gene Expression sequencing alignment with cellranger and demultiplex with demuxlet.

## Usage

This pipeline needs the following inputs:  
    - A list of multiplexed fastq files (`libraries.csv` file);  
    - VCF files with genotypes for each sample in the pool (`data/genotypes`)  

### Inputs

- The `libraries.csv` file:

It has the following columns: `pool_id`, `fastqs`, `sample`, and `library_type`. It is similar to the one required by cellranger-arc, with the addition of the pool_id column. Here's an example:

```
pool_id,fastqs,sample,library_type
pool1,/fs/ess/PAS2694/Data/Karch_multiome_SR004606_10X/FASTQ,NDRI_ADRC_Frontal_and_Cerebellum_Pool1,Gene Expression
pool1,/fs/ess/PAS2694/Data/Karch_multiome_SR004606_10X/FASTQ,NDRI_ADRC_Frontal_and_Cerebellum_Pool1_atac,Chromatin Accessibility
pool2,/fs/ess/PAS2694/Data/Karch_multiome_SR004606_10X/FASTQ,NDRI_ADRC_Frontal_and_Cerebellum_Pool2,Gene Expression
pool2,/fs/ess/PAS2694/Data/Karch_multiome_SR004606_10X/FASTQ,NDRI_ADRC_Frontal_and_Cerebellum_Pool2_atac,Chromatin Accessibility
```

- Genotypes: 

The `populate_data.sh` script will create soft links for the VCFs on the `data/genotypes` directory. Make sure all VCFs are following the sample naming convention. 

###  Pipeline configuration

The config file for this pipeline is at `config/config.yaml`. These are some parameters that need attention:

* **Cellranger-arc**:
    - `libraries`: Path to the libraries.csv file. I usually put it at `data/libraries.csv`.
    - `cellrangerarc_path`: Path to the cellranger executable. E.g. `/users/PAS2713/iaradsouza1/cellranger-arc-2.0.2/cellranger-arc`
    - `refdir`: Path to the reference directory for cellranger-arc. Check cellranger-arc docs to download it. 
    - `jobmode`: Path to the slurm template to be used by cellranger-arc to trigger jobs on slurm. An example can be found at `assets/slurm.template`. 

* **Skip vcf sorting?**
    - skip_vcf_sort: Default set to `False`. It has an optional step to normalize VCF files based on the bam files produced by cellranger-arc. 

* **Demuxlet**
- `genotypes`: Path to the directory containing the VCF files. E.g. `data/genotypes`.
- `conda_env`: Conda environment name for demuxlet and picard. The environment file is at `workflow/envs/demuxlet.yml`.


## Run pipeline on the OSC

1. Create conda environment
```
conda env create -f workflow/envs/demuxlet.yml
```

2. Create symbolic links for VCF files

Edit the `incoming_dir` variable on the `populate_data.sh` to the path where VCFs are located. Then run ` bash populate_data.sh`.

3. Set up the `slurm.template` for cellranger

An example is provided on `assets/slurm.template`. Set the path for your template on the `config/config.yaml` file.

*Optional: Create a snakemake environment and install slurm executor plugin. This pipeline considers snakemake version >= 8*
```
conda create -n snakemake -c bioconda snakemake
conda activate snakemake
pip install snakemake-executor-plugin-slurm
```

4. Run snakemake

On the pipeline root directory, run:

```
conda activate snakemake
snakemake --use-conda --rerun-incomplete --scheduler=greedy --dry-run  ## dry-run the pipeline
snakemake --use-conda --rerun-incomplete --scheduler=greedy --executor slurm --jobs <njobs>
```

## Profile configuration

The slurm configuration for the profile can be found at `profiles/default/config.yaml`.

# TODO

* Add cellranger-arc to a container
* Automatically download and set up cellranger-arc reference. 
* Add report with some statistics for the demultiplexing step. 
* Test other configurations for slurm profiles based on this: https://github.com/jdblischak/smk-simple-slurm


