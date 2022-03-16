## Bioinformatics pipeline overview

Due to the large size of sequence files, much of the following was performed on the High Performance Computing Cluster at Northern Arizona University. Many of the scipts can therefore not be compiled here and instead serve as a description of the steps taken in processing the metagenomes. The prefix `pre_*`  designates scripts part of this upstream analysis.

### QC filtering


[Using the program Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to trim adapter sequences and remove low quality and short reads. Output for paired end reads will end with `*PE_1.fastq` and `*PE_2.fastq` for forward and reverse reads respectively. 

```
java -jar $PATH/trimmomatic-0.39.jar ILLUMINACLIP:$PATH/adapters/TruSeq2-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:10:20 MINLEN:50
```


### Average Genome Size (AGS)

AGS was calculated using [MicrobeCensus](https://github.com/snayfach/MicrobeCensus), which uses the abundance of single copy genes to get an estimate of how many individuals are in a metagenome. It then divides this by the total number of base-pairs to get an AGS of the whole metagenome. 

```
run_microbe_census.py -t8 Seq_PE_1.fastq,Seq_PE_2.fastq ${1/_1.fastq/_census.txt}
```



### Assembling contigs

MEGAHIT was used to assemble QC filtered reads into contigs. 

```shell
megahit -1 Seq_PE_1.fastq -2 Seq_PE_2.fastq --k-list 21,29,39,59,79 --min-contig-len 400 -o megahit_output -t 16
```

QC filtered reads were then mapped against contigs using BBMap using the following

```shell
srun bbmap.sh in=Seq_PE_1.fastq in2=Seq_PE_1.fastq ref=contigs.fa covstats=constats.txt bincov=bincov.txt
```


### Gene calling

 Open Reading Frames (ORFs) were predicted usiung[Prodigal](https://github.com/hyattpd/Prodigal). 
```shell
prodigal -i contigs.fa -o genes -a genes.faa -d genes.fna -p meta
```

Taxonomy predicted from Kaiju and bacterial reads were identified using the following commands. These commands create 3 files: 
1) A line-by-line annotation of each gene
2) A summary file for the abundance of each taxa
3) A list of scaffolds with bacterial annnotations 

```shell

kaiju -v -o $kaijuout -t $PATH/nodes.dmp -f $PATH/kaiju_db_nr_euk.fmi -z 4 -i $1

kaiju2table -t $PATH/nodes.dmp -n $PATH/names.dmp -r genus -p -o $kaijusummary $kaijuout

kaiju-addTaxonNames -t $PATH/nodes.dmp -n $PATH/names.dmp -p -i $kaijuout -o $kaijuoutnamess

grep "Bacteria" $kaijuoutnamess > $kaijubact
```


### Gene Annotations

#### Function

KEGG functional genes were identified using [kofamscan](https://taylorreiter.github.io/2019-05-11-kofamscan/):

```shell
srun $PATH/kofam_scan-1.3.0/exec_annotation -f mapper -o result.txt genes.faa
```

Functional gene annotations were not included in this analysis, but we have included them here for public use (see `Kegg_output.csv`). We ask that any use of these data references the associated paper.

#### Taxonomy 

Kaiju was used to annotate the taxonomy of ORFs identified via prodigal. These annotations were then used to isolate bacterial genes and contigs.
The following commands create 3 files: 
1) A line-by-line annotation of each gene
2) A summary file for the abundance of each taxa
3) A list of contigs with bacterial annnotations 

```shell

kaiju -v -o $kaijuout -t /scratch/pfc25/kaijudb/nodes.dmp -f /scratch/pfc25/kaijudb/nr_euk/kaiju_db_nr_euk.fmi -z 4 -i $1

kaiju2table -t /scratch/pfc25/kaijudb/nodes.dmp -n /scratch/pfc25/kaijudb/names.dmp -r genus -p -o $kaijusummary $kaijuout

kaiju-addTaxonNames -t /scratch/pfc25/kaijudb/nodes.dmp -n /scratch/pfc25/kaijudb/names.dmp -p -i $kaijuout -o $kaijuoutnamess

grep "Bacteria" $kaijuoutnamess > $kaijubact
```

#### Codon and AA frequency

From the ORFs, we gathered the codon (`codon_parser.py`) and amino acid (`AA_parser.py` ) frequency for each metagenome. `AA_parser.py` also uses formulas for each amino acid to calculate the total elemental composition of the amino acids for a metagenome.  

### Combining and merging data

The script `merge_output.py` was used to merge files, adjust for read depth, filter bacterial contigs, and generate the following files:
* GC_output.txt - The adjusted GC and C:N of amino acids of bacterial readsd
* AA_codon_output.csv - The adjusted codon and amino acid frequencies, as well as the elemental ratios of the AA
* kegg_output.csv - Depth adjusted gene calls
* census.txt - the census output

The jupyter notebook files `Output_reformat.ipynb` and `kegg_output_parser.ipynb` were used to merge multiple samples into one dataframe. 

### Second analysis

We also determined the average GC content for each domain and for each gene in a metagenome using the script `remerge.py`. The output of these files was merged using `GC_by_kegg.ipynb` and `Second_output_reformat.ipynb`. These data were not included in this manuscript but have been included here for public use. Again, we ask that any use of these data references the associated paper. 
