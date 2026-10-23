# Parse single cell sequencing analysis

- Stephen Doyle


## Installed Split Pipe
- https://support.parsebiosciences.com/hc/en-us/articles/23060102930580-Pipeline-Installation-Current-Version 
```bash
cd /nfs/users/nfs_s/sd21/lustre_link/software/

# downloaded and transferred ParseBiosciences-Pipeline.1.9.0.zip to directory
unzip ParseBiosciences-Pipeline.1.9.0.zip
rm ParseBiosciences-Pipeline.1.9.0.zip

conda create -n spipe_1.9.0 conda-forge::python==3.12.8
conda activate spipe_1.9.0

cd ParseBiosciences-Pipeline.1.9.0
bash ./install_dependencies_conda.sh -i -y
pip install --no-cache-dir ./

# test
split-pipe

```


## Experiment 1

| Sample number | Species                | Sex/stage    | Sample_name             | Plate location |
|---------------|------------------------|--------------|-------------------------|----------------|
| 1             | Caenorhabditis elegans | mixed adults | celegans_adult_mixed    | A1,A7          |
| 2             | Trichuris muris        | adult male   | tmuris_adult_male       | A2,A8          |
| 3             | Trichuris muris        | adult female | tmuris_adult_female     | A3,A9          |
| 4             | Haemonchus contortus   | adult male   | hcontortus_adult_male   | A4,A10         |
| 5             | Haemonchus contortus   | adult female | hcontortus_adult_female | A5,A11         |
| 6             | Haemonchus contortus   | L3           | hcontortus_L3           | A6,A12         |


## Experiment 2

| Sample number | Species                | Sex/stage    | Sample_name             | Plate location |
|---------------|------------------------|--------------|-------------------------|----------------|
| 1             | Haemonchus contortus   | adult female | hcontortus_adult_female | A1,A2          |
| 2             | Haemonchus contortus   | adult male   | hcontortus_adult_male   | A3,A4          |
| 3             | Trichuris muris        | adult female | tmuris_adult_female     | A5,A6          |
| 4             | Trichuris muris        | adult male   | tmuris_adult_male       | A7,A8,A9       |
| 5             | Caenorhabditis elegans | mixed adults | celegans_adult_mixed    | A10,A11        |
| 6             | Haemonchus contortus   | L3           | hcontortus_L3           | A12            |





# Working directory and setup
```bash
cd /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation

```

## Setting up the reference files
- Haemonchus contortus
- Trichuris muris

### Haemonchus
```bash
# reference directory
cd /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/01_references

#------ Haemonchus
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS19.genomic.fa.gz

wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS19.annotations.gff3.gz

gunzip haemonchus_contortus.PRJEB506.WBPS19.annotations.gff3.gz



# convert the gff to gtf
module load agat/1.2.0--pl5321hdfd78af_0

# need to make a fix - split-pipe needs a "gene_biotype"
sed -i 's/; biotype/; gene_biotype/' haemonchus_contortus.PRJEB506.WBPS19.annotations.gtf

agat_convert_sp_gff2gtf.pl --gff haemonchus_contortus.PRJEB506.WBPS19.annotations.gff3 --gtf_version 3 --output  haemonchus_contortus.PRJEB506.WBPS19.annotations.gtf

agat_sp_manage_attributes.pl -f haemonchus_contortus.PRJEB506.WBPS19.annotations.gtf -att description,description_source,description_source_acc,Name,transcript_id -o haemonchus_contortus.PRJEB506.WBPS19.annotations.minimal.gff

cat haemonchus_contortus.PRJEB506.WBPS19.annotations.minimal.gff | awk '{if($3=="gene" || $3=="exon" || $3=="RNA") print}' |  sed -e 's/exon://g' -e 's/transcript://g' -e 's/gene://g' > haemonchus_contortus.PRJEB506.WBPS19.annotations.minimal_clean.gff

agat_convert_sp_gff2gtf.pl --gff haemonchus_contortus.PRJEB506.WBPS19.annotations.minimal_clean.gff --gtf_version 3 --output  haemonchus_contortus.PRJEB506.WBPS19.annotations.minimal_clean.gtf




# need to make a fix - split-pipe needs a "gene_biotype"
sed -i 's/; biotype/; gene_biotype/' haemonchus_contortus.PRJEB506.WBPS19.annotations.gtf

grep -v "agat" haemonchus_contortus.PRJEB506.WBPS19.annotations.genes_exons.gtf | grep -v "rRNA" > haemonchus_contortus.PRJEB506.WBPS19.annotations.genes_exons.no-agat.gtf

agat_sp_manage_attributes.pl -f haemonchus_contortus.PRJEB506.WBPS18.annotations.genes_exons.no-agat.gtf -att description,description_source,description_source_acc,Name,transcript_id -o haemonchus_contortus.PRJEB506.WBPS18.annotations.minimal.gff

cat haemonchus_contortus.PRJEB506.WBPS18.annotations.minimal.gff | awk '{if($3=="gene" || $3=="exon" || $3=="RNA") print}' |  sed -e 's/exon://g' -e 's/transcript://g' -e 's/gene://g' > haemonchus_contortus.PRJEB506.WBPS18.annotations.minimal_clean.gff

agat_convert_sp_gff2gtf.pl --gff haemonchus_contortus.PRJEB506.WBPS18.annotations.minimal_clean.gff --gtf_version 3 --output  haemonchus_contortus.PRJEB506.WBPS18.annotations.minimal_clean.gtf

# run spipe make reference
conda activate spipe_1.9.0

bsub.py 10 spipe_mkref_hcon "
    split-pipe 
        --mode mkref 
        --genome_name hcontortus 
        --fasta haemonchus_contortus.PRJEB506.WBPS19.genomic.fa.gz 
        --genes haemonchus_contortus.PRJEB506.WBPS19.annotations.minimal_clean.gtf 
        --output_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/01_references/hcontortus"

```



### Trichuris muris
```bash
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/trichuris_muris/PRJEB126/trichuris_muris.PRJEB126.WBPS19.genomic.fa.gz

wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/trichuris_muris/PRJEB126/trichuris_muris.PRJEB126.WBPS19.annotations.gff3.gz

cat trichuris_muris.PRJEB126.WBPS19.annotations.gff3 | awk '{if($3=="gene" || $3=="exon" || $3=="mRNA") print}' |  sed -e 's/exon://g' -e 's/Transcript://g' -e 's/Gene://g' > trichuris_muris.PRJEB126.WBPS19.annotations.minimal_clean.gff3

# need to make a fix - split-pipe needs a "gene_biotype"
sed -i 's/;biotype/;gene_biotype/' trichuris_muris.PRJEB126.WBPS19.annotations.minimal_clean.gff3


# gff to gtf
agat_convert_sp_gff2gtf.pl --gff trichuris_muris.PRJEB126.WBPS19.annotations.minimal_clean.gff3 --gtf_version 3 --output  trichuris_muris.PRJEB126.WBPS19.annotations.gtf

bsub.py 10 spipe_mkref_tmur "split-pipe \
        --mode mkref \
        --genome_name tmuris \
        --fasta trichuris_muris.PRJEB126.WBPS19.genomic.fa.gz \
        --genes trichuris_muris.PRJEB126.WBPS19.annotations.gtf \
        --output_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/01_references/tmuris"




# Combined reference
bsub.py 10 spipe_mkref_species-combo "split-pipe \
--mode mkref \
--genome_name hcontortus tmuris \
--fasta haemonchus_contortus.PRJEB506.WBPS19.genomic.fa.gz  \
        trichuris_muris.PRJEB126.WBPS19.genomic.fa.gz \
--genes haemonchus_contortus.PRJEB506.WBPS19.annotations.minimal_clean.gtf \
        trichuris_muris.PRJEB126.WBPS19.annotations.gtf \
--output_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/01_references/hcontortus_tmuris_combined"
```







## Get the data

```bash
# data dir 
cd /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/02_raw-data

module load irods_extractor/v4.0.1

bsub.py 10 extract_run1 "irods_extractor —-studyid 7067 --runid 48163 --laneid 1"
bsub.py 10 extract_run2 "irods_extractor —-studyid 7067 --runid 48240 --laneid 1"



Run: 48240
Individual tags
Tag Index	Tag	Reads, count	Reads, %
1	CAGATC	153,808,413	42.1
2	ACTTGA	144,960,947	39.6
0	 	66,897,611	18.3

Run: 48163
Tag Index	Tag	Reads, count	Reads, %
1	CAGATC	317,343,483	58.7
888	TGTGCA	211,409,771	39.1
0	 	11,602,344	2.1

#- 888 - there was 40% PhIX in the run - not sure why it goes to an index?



# check read quality with fastqc and multiqc

bsub.py 10 fastqc "~sd21/bash_scripts/run_fastqc"

scp sd21@farm5-head2:~/lustre_link/parse_single_cell/expdata/skb_run1/multiqc_report.html skb_run1_multiqc_report.html
```

[Multiqc report](../04_analysis/skb_run1_multiqc_report.html) 





## Run the pipeline
- https://support.parsebiosciences.com/hc/en-us/articles/23060220494228-Running-the-Pipeline-Current-Version#h_01GT2VJQS0Y83BN1N7J4R1DMDE 
```bash
# working directory
cd /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/04_analysis 

# Generate parameter file
echo "post_min_map_frac 0.01" > parfile.txt

# Add this parameter to your spilt-pipe "--mode all" call
# --parfile parfile.txt

# load spipe
conda activate spipe_1.9.0

# hcontortus - experiment 1
bsub.py --queue hugemem --threads 30 200 split-pipe_hcontortus_48163_1_1 "\
split-pipe \
--mode all \
--chemistry v2 \
--kit WT_mini \
--nthreads 30 \
--parfile parfile.txt \
--fq1 /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/02_raw-data/results/raw_fastqs/48163_1#1_raw_1.fastq.gz \
--fq2 /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/02_raw-data/results/raw_fastqs/48163_1#1_raw_2.fastq.gz \
--output_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/04_analysis/hcontortus_48163_1_1 \
--genome_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/01_references/hcontortus \
--sample hcontortus_adult_male A4,A10 \
--sample hcontortus_adult_female A5,A11 \
--sample hcontortus_L3 A6,A12"



# tmuris - experiment 1
bsub.py --queue hugemem --threads 30 200 split-pipe_tmuris_48163_1_1 "\
split-pipe \
--mode all \
--chemistry v2 \
--kit WT_mini \
--nthreads 30 \
--parfile parfile.txt \
--fq1 /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/02_raw-data/results/raw_fastqs/48163_1#1_raw_1.fastq.gz \
--fq2 /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/02_raw-data/results/raw_fastqs/48163_1#1_raw_2.fastq.gz \
--output_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/04_analysis/tmuris_48163_1_1 \
--genome_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/01_references/tmuris \
--sample tmuris_adult_male A2,A8 \
--sample tmuris_adult_female A3,A9"



# hcontortus - experiment 2 - sublibrary 1
bsub.py --queue hugemem --threads 30 200 split-pipe_hcontortus_48240_1_1 "\
split-pipe \
--mode all \
--chemistry v2 \
--kit WT_mini \
--nthreads 30 \
--parfile parfile.txt \
--fq1 /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/02_raw-data/results/raw_fastqs/48240_1#1_raw_1.fastq.gz \
--fq2 /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/02_raw-data/results/raw_fastqs/48240_1#1_raw_2.fastq.gz \
--output_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/04_analysis/hcontortus_48240_1_1 \
--genome_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/01_references/hcontortus \
--sample hcontortus_adult_female A1,A2 \
--sample hcontortus_adult_male A3,A4 \
--sample hcontortus_L3 A12"

# hcontortus - experiment 2 - sublibrary 2
bsub.py --queue hugemem --threads 30 200 split-pipe_hcontortus_48240_1_2 "\
split-pipe \
--mode all \
--chemistry v2 \
--kit WT_mini \
--nthreads 30 \
--parfile parfile.txt \
--fq1 /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/02_raw-data/results/raw_fastqs/48240_1#2_raw_1.fastq.gz \
--fq2 /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/02_raw-data/results/raw_fastqs/48240_1#2_raw_2.fastq.gz \
--output_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/04_analysis/hcontortus_48240_1_2 \
--genome_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/01_references/hcontortus \
--sample hcontortus_adult_female A1,A2 \
--sample hcontortus_adult_male A3,A4 \
--sample hcontortus_L3 A12"



# tmuris - experiment 2 - sublibrary 1
bsub.py --queue hugemem --threads 30 200 split-pipe_tmuris_48240_1_1 "\
split-pipe \
--mode all \
--chemistry v2 \
--kit WT_mini \
--nthreads 30 \
--parfile parfile.txt \
--fq1 /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/02_raw-data/results/raw_fastqs/48240_1#1_raw_1.fastq.gz \
--fq2 /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/02_raw-data/results/raw_fastqs/48240_1#1_raw_2.fastq.gz \
--output_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/04_analysis/tmuris_48240_1_1 \
--genome_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/01_references/tmuris \
--sample tmuris_adult_female A5,A6 \
--sample tmuris_adult_male A7,A8,A9"

# tmuris - experiment 2 - sublibrary 2
bsub.py --queue hugemem --threads 30 200 split-pipe_tmuris_48240_1_2 "\
split-pipe \
--mode all \
--chemistry v2 \
--kit WT_mini \
--nthreads 30 \
--parfile parfile.txt \
--fq1 /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/02_raw-data/results/raw_fastqs/48240_1#2_raw_1.fastq.gz \
--fq2 /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/02_raw-data/results/raw_fastqs/48240_1#2_raw_2.fastq.gz \
--output_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/04_analysis/tmuris_48240_1_2 \
--genome_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/01_references/tmuris \
--sample tmuris_adult_female A5,A6 \
--sample tmuris_adult_male A7,A8,A9"
```


## Combine sublibraries from experiment 2
```bash
conda activate spipe_1.9.0

# hcontortus
bsub.py --queue hugemem --threads 30 200 hcontortus_48240_combine "\
    split-pipe \
    --mode comb \
    --sublibraries /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/04_analysis/hcontortus_48240_1_1 \
                   /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/04_analysis/hcontortus_48240_1_2 \
    --output_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/04_analysis/hcontortus_48240_1_combined"

#tmuris
bsub.py --queue hugemem --threads 30 200 tmuris_48240_combine "\
split-pipe \
    --mode comb \
    --sublibraries /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/04_analysis/tmuris_48240_1_1 \
                   /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/04_analysis/tmuris_48240_1_2 \
    --output_dir /nfs/users/nfs_s/sd21/lustre_link/single-nuclei-rnaseq/protocol_optimisation/04_analysis/tmuris_48240_1_combined"

