# Parse single cell sequencing analysis

- Stephen Doyle

# https://support.parsebiosciences.com/hc/en-us/articles/17166220335636



#working dir: /nfs/users/nfs_s/sd21/lustre_link/parse_single_cell 

# setting up the analysis tools
```bash
conda create -n spipe python=3.10

conda activate spipe

```



mkdir analysis expdata genomes


```bash
split-pipe \
--mode mkref \
--genome_name hg38 \
--fasta /newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
--genes /newvolume/genomes/Homo_sapiens.GRCh38.109.gtf.gz \
--output_dir /newvolume/genomes/hg38

```



```bash
# Example with explicit sample definitions
split-pipe \
    --mode all \
    --chemistry v2 \     
    --genome_dir /newvolume/genomes/hg38_mm10/ \
    --fq1 /newvolume/expdata/202306_ex1/s1_S1_R1_001.fastq.gz \
    --fq2 /newvolume/expdata/202306_ex1/s1_S1_R2_001.fastq.gz \
    --output_dir /newvolume/analysis/202306_ex1_an1 \
    --sample xcond_1 A1-A6 \
    --sample xcond_2 C1-C6 \
    --sample xcond_3 D1-D6


```

## Running the test data to ensure the pipeline works.
```bash
## download some test data to test pipeline
wget https://www.dropbox.com/sh/r33i9enwz1rnhn4/AABDWXd1o7Swcq9-RgSXY0Kra/hg-ref/genome.fas.gz
wget  https://www.dropbox.com/sh/r33i9enwz1rnhn4/AABGoHYFlrkw69PW2RkZ1vRDa/hg-ref/genes.gtf.gz

# make the test reference
split-pipe \
--mode mkref \
--genome_name hg38 \
--fasta /nfs/users/nfs_s/sd21/lustre_link/parse_single_cell/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
--genes /nfs/users/nfs_s/sd21/lustre_link/parse_single_cell/genomes/Homo_sapiens.GRCh38.109.gtf.gz \
--output_dir /nfs/users/nfs_s/sd21/lustre_link/parse_single_cell/genomes/hg38


mkdir /nfs/users/nfs_s/sd21/lustre_link/parse_single_cell/expdata/test

wget https://www.dropbox.com/sh/r33i9enwz1rnhn4/AABTIISx9pIZcqZlm_r6YM1Ka/pbmc1M_R1.fastq.gz
wget https://www.dropbox.com/sh/r33i9enwz1rnhn4/AAAE_gZo_dD49jt5eEVUAmxKa/pbmc1M_R2.fastq.gz

split-pipe \
--mode all \
--chemistry v2 \
--fq1 /nfs/users/nfs_s/sd21/lustre_link/parse_single_cell/expdata/test/pbmc1M_R1.fastq.gz \
--fq2 /nfs/users/nfs_s/sd21/lustre_link/parse_single_cell/expdata/test/pbmc1M_R2.fastq.gz \
--output_dir /nfs/users/nfs_s/sd21/lustre_link/parse_single_cell/analysis/test-out \
--genome_dir /nfs/users/nfs_s/sd21/lustre_link/parse_single_cell/genomes/hg38 \
--sample donor1 A1-A12 \
--sample donor2 B1-B12 \
--sample donor3 C1-C12 \
--sample donor_4 D1-D12

```


