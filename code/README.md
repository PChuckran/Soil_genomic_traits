## Code 

Scripts which are not stand-alone functions have prefixes which indicate the order in which they should be run. Functions are located in the functions directory. Naming/numbering scheme is as follows: 

**pre_** : Upstream bioinformatics code which cannot be run as part of this github repo. Details about this code can be found in the `upstream_bioinformatics` folder and markdown

**1[a-e]_** : Code which creates a sample list and downloads corresponding nongenomic data. Much of this code takes a considerable amount of time to run, so running this code is therefore only recommended if nongenomic data needs updating - such as if the metagenomic sample list has been updated. 

**2_** : Currently consists of one file - `2_construct_full_dataframe.R`, which merges genomic and nongenomic data into one dataframe, and coalesces portions of the soil data according to geographic distance. Should be run if changes to either the nongenomic or genomic data are made. Produces a csv fille - `data/input/derived/full_df.csv`.

**3_** : Analysis files. This includes Random Forest Analysis, table and figure generation, descriptive statistics. 

**4_** : Scripts (if any) dependent on output from **3_**. 
