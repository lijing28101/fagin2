# Fagin2

A pipeline for the classification of orphans into origin classes using a syntenic filter.


## Funding

This work is funded by the National Science Foundation grant:

[NSF-IOS 1546858](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1546858)
Orphan Genes: An Untapped Genetic Reservoir of Novel Traits


## Installation

```
# Clone the repository
git clone https://github.com/lijing28101/fagin2.git

cd fagin2

# Create a new environment for Fagin2
conda create --name fagin2 --file requirements.txt

# Activate the new environment
conda activate fagin2

# Go back to the parent dir and build fagin2
cd ..
R CMD build fagin2

# Install
R CMD INSTALL fagin2_2.0.0.tar.gz

```

**Potential error for installation**
```
Error in dyn.load(file, DLLpath = DLLpath, ...) :
  unable to load shared object '/work/LAS/mash-lab/jing/bin/Anaconda3/envs/fagin2/lib/R/library/L1pack/libs/L1pack.so':
  libRlapack.so: cannot open shared object file: No such file or directory
Calls: <Anonymous> ... asNamespace -> loadNamespace -> library.dynam -> dyn.load
Execution halted
ERROR: lazy loading failed for package ‘fagin2’
```
**Method to solve the problem (not the best way)**
```
# change to your lib path
cd Anaconda3/envs/fagin2/lib/

ls -ltrh 'liblapack.so'
ls -ltrh 'libblas.so'
mv liblapack.so libRlapack.so
mv libblas.so libRblas.so
```

## Input

1. Phylogeny tree for all included species, e.g.
   ```
   ((Saccharomyces_cerevisiae,Saccharomyces_kudriavzevii)s3,Saccharomyces_arboricola)s2;
   ```

2. Genome sequence for each species

3. GFF file for each species
   - Must include "gene", "mRNA", "CDS", and "exon" in column 3
   - Must include "Parent", "ID" or "Name" (or both) in column 9

  ```
  # sample GFF file
  NC_001133.9     RefSeq  gene    1807    2169    .       -       .       ID=PAU8;Parent=-
  NC_001133.9     RefSeq  mRNA    1807    2169    .       -       .       ID=NM_001180043.1;Parent=PAU8
  NC_001133.9     RefSeq  exon    1807    2169    .       -       .       ID=id3;Parent=NM_001180043.1
  NC_001133.9     RefSeq  CDS     1807    2169    .       -       0       ID=NP_009332.1;Parent=NM_001180043.1
  NC_001133.9     RefSeq  gene    2480    2707    .       +       .       ID=YAL067W-A;Parent=-
  NC_001133.9     RefSeq  mRNA    2480    2707    .       +       .       ID=NM_001184582.1;Parent=YAL067W-A
  NC_001133.9     RefSeq  exon    2480    2707    .       +       .       ID=id4;Parent=NM_001184582.1
  NC_001133.9     RefSeq  CDS     2480    2707    .       +       0       ID=NP_878038.1;Parent=NM_001184582.1
  NC_001133.9     RefSeq  gene    7235    9016    .       -       .       ID=SEO1;Parent=-
  NC_001133.9     RefSeq  mRNA    7235    9016    .       -       .       ID=NM_001178208.1;Parent=SEO1
  NC_001133.9     RefSeq  exon    7235    9016    .       -       .       ID=id5;Parent=NM_001178208.1
  NC_001133.9     RefSeq  CDS     7235    9016    .       -       0       ID=NP_009333.1;Parent=NM_001178208.1
  ```

4. Synteny map for the focal species versus each other species
   - Suggest obtain from mummer4
   ```bash
   # Name of query and target species
   query=saccharomyces_cerevisiae
   target=saccharomyces_arboricola
   prefix=${query}_${target}
   #
   # ${query}_genome.fna and ${target}_genome.fna are the path of genome sequences
   nucmer --mum -t 36 -p ${prefix} ${query}_genome.fna ${target}_genome.fna
   show-coords -T -r -c -l ${prefix}.delta  > ${prefix}.coords
   #
   # Convert coords format to the syn format required for fagin
   awk 'BEGIN{FS="\t"; OFS="\t"} {if ($4>$3) $8="+"; else $8="-"; print  $12,$1,$2,$13,$3,$4,$7,$8}' ${prefix}.coords | awk 'BEGIN{OFS="\t"}{if($5<$6){start=$5;stop=$6}else{start=$6;stop=$5};print $1,$2,$3,$4,start,stop,$7,$8}' | sed '1,4d' > ${prefix}.syn
   ```

   - Synteny format as below (without header):

   ```diff
   ! start coordinate should be always smaller than end even if the orientation is -
   ! use tab as delimiter
   - query_chromosome   query_start   query_end    focal_chromosome   focal_start  focal_end  orientation
   ```
   ```
   NC_001133.9     1783    2167    CM001577.1      1020875 1021256 90.65   -
   NC_001133.9     1801    2167    CM001574.1      205548  205924  85.68   +
   NC_001133.9     1804    2167    CM001576.1      770149  770512  89.29   -
   NC_001133.9     1805    2167    CM001568.1      106576  106941  90.44   +
   NC_001133.9     1805    2165    CM001574.1      491     860     87.57   +
   NC_001133.9     1806    2167    CM001574.1      84073   84443   79.78   -
   ```

5. A gene list for each focal species which will compare to target species

   - must use the Name but not the ID of mRNA in GFF file if both Name and ID existing on column 9


```diff
! Important notes:
! Chromosome ID in GFF, fna, and syn files should be consistent.
! The number of chromosomes or scaffolds in GFF, fna, and syn files should be same.
! You may need to request large memory to run the program if the synteny map is large
```

## Configuration

To run and configure `fagin`, you need to set paths to your data in
a configuration object. The default configuration can be generated

```R
config()
```

Example to generate a configuration object:

```R
library(fagin2)

get_yeast_config <- function(){

  # create a default configuration
  con <- config()

  # set parameter
  con@synder@offsets = c(1L,1L) # offsets for mummer, c(0L,1L) for satsuma
  con@synder@trans = "p" # percent identity transform (0-100), or "d" for proportion transform (0-1)
  # alnrate is the proportion of target sequence match to query sequence for alignment
  # default cutoff for alnrate are 0
  con@alignment@alnrate@prot2prot=0.5 # protein vs protein
  con@alignment@alnrate@prot2allorf=0.5 # protein vs all ORF
  con@alignment@alnrate@prot2transorf=0.5 #protein vs ORF on the mRNA
  con@alignment@alnrate@dna2dna=0.5 # DNA vs DNA

  # focal species can be one or more species (max: the total number of species in tree file)
  con@input@focal_species = c("Saccharomyces_cerevisiae","Saccharomyces_arboricola")

  # set path for input files
  # name of list should be same as the species in tree
  con@input@gff <- list(
    "Saccharomyces_arboricola"   = "saccharomyces_arboricola_annotation.gff"
    , "Saccharomyces_cerevisiae"   = "saccharomyces_cerevisiae_annotation.gff"
    , "Saccharomyces_kudriavzevii" = "saccharomyces_kudriavzevii_annotation.gff"
  )
  con@input@fna <- list(
    "Saccharomyces_arboricola"   = "saccharomyces_arboricola_genome.fna"
    , "Saccharomyces_cerevisiae"   = "saccharomyces_cerevisiae_genome.fna"
    , "Saccharomyces_kudriavzevii" = "saccharomyces_kudriavzevii_genome.fna"
  )

  # syn path is a two-dimensional list
  # first level: focal species; second level: target species
  con@input@syn <- list(
    Saccharomyces_cerevisiae = list(
      "Saccharomyces_arboricola"   = "saccharomyces_cerevisiae.vs.saccharomyces_arboricola.syn"
      , "Saccharomyces_kudriavzevii" = "saccharomyces_cerevisiae.vs.saccharomyces_kudriavzevii.syn"
    ),
    Saccharomyces_arboricola = list(
      "Saccharomyces_cerevisiae"   = "saccharomyces_arboricola_saccharomyces_cerevisiae.syn"
      , "Saccharomyces_kudriavzevii" = "saccharomyces_arboricola_saccharomyces_kudriavzevii.syn"
    )
  )
  con@input@tree <- "tree"
  con@input@gene_list <- list(
    "Saccharomyces_cerevisiae" = "sc.list.txt"
    , "Saccharomyces_arboricola"   =  "sa.list.txt"
  )
  con@archive = "yeast-archive" # dir for output files

  validate_config(con)
  con
}

con <- get_yeast_config() # create revised configuration
```

## Running fagin

1. Load input genome sequence and GFF file for each species, then processing and summarize the input data (save as `[species_name]_data.rds` for reuse).

2. Do alignment and comparison for each pair of focal and target species (save as `[focal_species]-[target_species].rds` for reuse)

3. Summarize the comparison and export as xlsx file for each focal species (also save as `[focal_species]_result.rds`)


```R
# con is your configuration object
# overwrite.result: whether rerun the summarize step for each focal species
# overwrite.result=TRUE, when you add new target species in con for comparison
# overwrite.result=FALSE, only when you add new focal species, but not add new target species for previous focal
# cores: number of cores used for parallel steps, recommend for large data

m <- run_fagin(con, overwrite.result=T) # without parallel, 30 min for sample data

or

m <- run_fagin_parallel(con, cores=16) # run loading data and comparison in parallel, 13 min for sample data
```

**Note: Fagin will save the input data for each individual species, and pair of focal-target comparison as rds file for reuse. `run_fagin` and `run_fagin_parallel` will skip to load input data or pairwise comparison if the rds file already existed in the archive dir. You should delete the rds file for a specific species or pairwise comparison if you want to re-analysis that step.**

## Output

The output m is a multi-dimensional list. First level is focal species. Second level include 5 data frames:

1. `feature_focal`: main result for homolog classes
2. `aatab_focal`: Significant hits for protein vs protein
3. `orftab_focal`: Significant hits for protein vs all ORFs
4. `transtab_focal`: Significant hits for protein vs ORFs on the mRNA
5. `gentab_focal`: Significant hits for DNA vs DNA

All of these result also write into xlsx files for each focal species in the archive dir
