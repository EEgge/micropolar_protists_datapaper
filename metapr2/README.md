# MetaPR2 files

## Data sets
* 207:	MicroPolar_Hiseq A
* 208:	MicroPolar_Hiseq B
* 209:	MicroPolar_Miseq

## Dada2 processing scripts
Scripts are in directory `/scripts`
The script is run on the Roscoff ABIMS server as follows where $DIR is the directory where scripts are and $DATASET_ID is the dataset ID (207, 208, 209)

```
cd $DIR

source $CONDA3/activate cutadapt-2.8

/opt/6.x/R-3.5.1/bin/Rscript --no-save --no-restore script_dada2.R -d $DATASET_ID > script_dada2_$1_$DATE.out
```


## Instructions to regenerate the data (no need to do a priori)
* Do not change the structure of the directories because the R markdown script will not work
* Launch the R project in the dada2 directory "metaPR2 micropolar.Rproj"

## ASVs
* ASVs from different datasets have been merged using the sequence hashtag (2 identical sequences have the same hashtag).  The hashtag is used as the asv name.
* The following filtering has been used
  * At least 10 reads for any asv 
  * Bootstrap value at class level >= 90


## Directories and files

* functions_db.R: database function
* functions_metapr2.R: metapr2 function
* metaPR2 micropolar.Rmd: run to regenerate data (no need)
* /mysql : Configuration file to read from the MySQL dadabase
* /export : Files produced by R markdown file
    * **metapr2_asv_set_xxx_Eukaryota.fasta**: fasta file with all asvs
    * **metapr2_long_asv_set_xxx_Eukaryota.xlsx**: long form of the dataset - each line corresponds to a different sample and different asv
    * **metapr2_samples_asv_set_xxx.xlsx**: list of all samples with metadata
    * **metapr2_wide_asv_set_xxx_Eukaryota.xlsx**: wide form of the dataset - each line corresponds to a different asv and columns correspond to samples
    * **phyloseq_metapr2_asv_set_xxx_Eukaryota.rds**: phyloseq file for all samples and asvs

## Structure of the phyloseq object
* otu_table()   OTU Table:         [ 7036 taxa and 196 samples ]  
* sample_data() Sample Data:       [ 196 samples by 55 sample variables ]  
* tax_table()   Taxonomy Table:    [ 7036 taxa by 8 taxonomic ranks ]  

## Sample variables
 [1] "sample_id"                   "file_name"                   "sample_name"                 "sample_code"                 "metadata_code"              
 [6] "replicate"                   "DNA_RNA"                     "fraction_name"               "fraction_min"                "fraction_max"               
[11] "reads_total"                 "sample_remark"               "metadata_id"                 "project"                     "cruise"                     
[16] "station_id"                  "year"                        "date"                        "season"                      "depth_level"                
[21] "depth"                       "substrate"                   "latitude"                    "longitude"                   "country"                    
[26] "oceanic_region"              "bottom_depth"                "temperature"                 "salinity"                    "O2"                         
[31] "fluorescence"                "Chla"                        "NO2"                         "NO3"                         "PO4"                        
[36] "Si"                          "bact_ml"                     "syn_ml"                      "peuk_ml"                     "neuk_ml"                    
[41] "crypto_ml"                   "virus_small_ml"              "virus_large_ml"              "dataset_code"                "dataset_name"               
[46] "processing_pipeline_metapr2" "processing_date"             "sequencing_technology"       "sequencing_type"             "sequencing_company"         
[51] "region"                      "ecosystem"                   "substrate_type"              "data_available"              "contact_name" 
