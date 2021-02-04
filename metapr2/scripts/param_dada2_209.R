# =================================
# Swtiches to adjust the processing
# =================================
  remove_primers <- TRUE
  do_cutadapt     <- FALSE
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

dataset_code <- "MicroPolar_MiSeq"
dataset_id <- 209
dataset_path <- "/shared/projects/metapr2/209_MicroPolar_MiSeq/" # Server

pr2_file <- "/projet/umr7144/dipo/vaulot/metabarcodes/database/pr2_4.12.0/pr2_version_4.12.0_18S_dada2.fasta.gz"

# -- File structure

  paired_reads = TRUE
  
  file_identifier = ".fastq"  # String to identify the files to be processed.
  R1_identifier = "_R1"
  R2_identifier = "_R2"
  file_name_separator = "_L001"
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
  truncLen = c(270,240) # This influences the number of ASVs and the percent of asv recovered (need to remove 20 and 21)
  minLen = c(270,240)
  truncQ = 2         
  maxEE = c(10, 10) 
  maxLen = 400  # This is for 454 to remove long and bad reads

# Reduce the number of asvs for problematic cases
  max_number_asvs = 0

# parameters for removeBimeraDenovo
  method_chimera = "pooled"

# With qsub

# Job 2272076 (qsub_dada2.sh) Complete
#  User             = vaulot
#  Queue            = short.q@n116.sb-roscoff.fr
#  Host             = n116.sb-roscoff.fr
#  Start Time       = 12/17/2020 01:58:39
#  End Time         = 12/17/2020 03:02:32
#  User Time        = 20:38:35
#  System Time      = 00:31:32
#  Wallclock Time   = 01:03:53
#  CPU              = 21:10:07
#  Max vmem         = 183.943G
#  Exit Status      = 0

# usage    1:                 cpu=20:58:38, mem=151363.54358 GBs, io=4.70839, vmem=3.475G, maxvmem=183.943G
# Thu Dec 17 03:02:32 CET 2020 job qsub_dada2.sh done

# With SLURM (64 proc 8 GB each) - died

#SBATCH --cpus-per-task 64
#SBATCH --mem 8GB                    # mémoire vive pour l'ensemble des cœurs

# slurmstepd-n96: error: Job 233059 exceeded memory limit (238 083 160 > 8 388 608), being killed
# slurmstepd-n96: error: Exceeded job memory limit
