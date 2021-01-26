# =================================
# Swtiches to adjust the processing
# =================================

  do_cutadapt     <- TRUE
  do_summary      <- TRUE
  do_plot_quality <- TRUE
  do_filtering <- TRUE # If TRUE and primers are present must also do_cutadapt = TRUE
  do_dada2    <- TRUE
  do_taxo <- TRUE
  multithread <- TRUE
  multithread_filter <- TRUE
  bigdata <- TRUE
  
# =================================
# Read the parameters 
# =================================

dataset_code <- "MicroPolar_HiSeqB"
dataset_id <- 208
dataset_path <- "/shared/projects/metapr2/208_MicroPolar_HiSeqB/" # Server

pr2_file <- "/projet/umr7144/dipo/vaulot/metabarcodes/database/pr2_4.12.0/pr2_version_4.12.0_18S_dada2.fasta.gz"

# -- File structure

  paired_reads = TRUE
  
  file_identifier = ".fastq"  # String to identify the files to be processed.
  R1_identifier = "1.fastq"
  R2_identifier = "2.fastq"
  file_name_separator = ".[12].fastq"
# This the first character of the file name to consider to extract the sample name (usually = 1)
  sample.names_first_character =  1


# --  Other parameters
# target
  gene = "18S rRNA"
  gene_region = "V4"
  organelle = "nucleus"


# primer sets (TAReuk454FWD1	CCAGCASCYGCGGTAATTCC	V4 18S Next.Rev	ACTTTCGTTCTTGATYRATGA) 

  FWD = "CCAGCASCYGCGGTAATTCC"
  REV = "ACTTTCGTTCTTGATYRATGA"

  anchor = ""  # Put in front of primer to anchor primer at start of sequence
  
# parameters for filterAndTrim
  sequencer = "Illumina"

# parameters for filterAndTrim
  truncLen = c(240,200) # This influences the number of ASVs and the percent of asv recovered (need to remove 20 and 21)
  minLen = c(240,200)
  truncQ = 2         
  maxEE = c(10, 10) 
  maxLen = 400  # This is for 454 to remove long and bad reads

# Reduce the number of asvs for problematic cases
  max_number_asvs = 0

# parameters for removeBimeraDenovo
  method_chimera = "pooled"

