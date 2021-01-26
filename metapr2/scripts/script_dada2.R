# =================================
# Name: script_dada2.R
# Object: dada2 processing script for 18S metabarcodes
# Author: Daniel Vaulot

# Changes

# version 1.0.0 - 2020-11-24

# version 1.0.1 - 2020-11-25
#  * Add max_number_asvs (if = 0 then no filtering)
#  * do_filter = FALSE, removes any filtering
#  * Taxonomy is now assigned by slices of 5000 (need to turn Multithread to FALSE if problem)

# version 2.0.0 - 2020-11-30
#  * Use the optparse library for parameters
#  * Now parameters are in file param_dada2_xxx and read from the main script
#  * Testing option do not perform dada2 and taxo
#  * Use cutadapt 2.8 (needs conda before: source $CONDA3/activate cutadapt-2.8)
#  * Print all parameters in the output files
#
# version 2.0.1 - 2020-12-08
# * Add a multithread_filter parameter because FilterandTrim use mcapply function which causes problem with slurm 
# * Add "--cores=0" to cutadapt for automatic detection of multicores
# * Two new parameters
#   * do_cutadapt - runs the cutadapat parameter (this is set true the first time use_cutadapt is set true)
#   * remove_primers - use cutadapt files (if existing)
# * If testing, multithread is false


# To run

# cd /projet/umr7144/dipo/vaulot/metabarcodes/R
# /opt/6.x/R-3.5.1/bin/Rscript --no-save --no-restore script_dada2.R -d 070 > script_dada2_070.out

# For testing
# /opt/6.x/R-3.5.1/bin/Rscript --no-save --no-restore script_dada2.R -d 070 -t > script_dada2_070_test.out


# ================================= 

# Packages ---------------------------------------------


suppressPackageStartupMessages({
  library(dada2) # Must use version >= 1.12
  library(Biostrings)
  library(ShortRead)
  library(stringr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(purrr)
  
  library("optparse")
})

# Parameters ---------------------------------------------

## First read in the arguments listed at the command line
option_list = list(
  make_option(c("-d", "--dataset"), type="character", action = "store", default=006, 
              help="ID of dataset to process [default= %default]", metavar="number"),
  make_option(c("-t", "--test"), type="character", action = "store_true", default=FALSE, 
              help="Test [default= %default]", metavar="logical")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dataset_id  = opt$dataset
testing = opt$test

## Read the parameter file and print it in output file  

file_param <- str_c("param_dada2_",dataset_id,".R")

source(file_param)

system(str_c("cat ",file_param))

## When testing is on, disable dada2 and below  

if (testing) {
  do_dada2    <- FALSE
  do_taxo <- FALSE
  bigdata <- FALSE
  multithread <- FALSE
  multithread_filter <- FALSE  
}

## Take care of novel parameters 

if(!exists("multithread_filter")) multithread_filter <- FALSE # To prevent problem with SLURM
if(!exists("remove_primers") & do_cutadapt == TRUE) remove_primers <- TRUE # For older param versions that did not contain remove_primers


print(sessionInfo())
# =================================
#   Define variables
# =================================

tax_levels <- c("kingdom", "supergroup","division", "class", "order", "family", "genus", "species")

# Maximum number of quality plots to print
max_plot_quality = 6

# For assigning taxonomy by chunks
taxo_slice_size = 1000

# Primer information
primer_length_fwd <- str_length(FWD)
primer_length_rev <- str_length(REV)


# =================================
#    Create directories 
# =================================


scratch_path <- stringr::str_replace(dataset_path, "projet", "scratch")
scratch_path <- stringr::str_replace(scratch_path, "sbr/ccebarcodep1408", "umr7144/dipo/vaulot")
if(!dir.exists(scratch_path)) dir.create(scratch_path)

path_dataset <- function(file_name) str_c(dataset_path, file_name)
path_scratch <- function(file_name) str_c(scratch_path, file_name)
file_dataset <- function(file_end)  str_c(dataset_path, "dada2/", dataset_code, file_end)

dir_fastq <-    path_dataset("fastq")

# dir_fastqN <-   path_scratch("fastqN/") # for filtered files removing reads with N
dir_cutadapt <- path_scratch("cutadapt/") # files after cutadapt
dir_filtered <- path_scratch("fastq_filtered/")
dir_qual <-     path_dataset("qual_pdf/")
dir_dada2 <-    path_dataset("dada2/")
dir_blast <-    path_dataset("blast/")

# if(!dir.exists(dir_fastqN)) dir.create(dir_fastqN)
if(!dir.exists(dir_cutadapt)) dir.create(dir_cutadapt)
if(!dir.exists(dir_filtered)) dir.create(dir_filtered)
if(!dir.exists(dir_qual)) dir.create(dir_qual)
if(!dir.exists(dir_dada2)) dir.create(dir_dada2)
if(!dir.exists(dir_blast)) dir.create(dir_blast)


# =================================
#    Get the file names 
# =================================

fns <- sort(list.files(dir_fastq, full.names = TRUE))
fns <- fns[str_detect( basename(fns),file_identifier)]
if (testing) fns <- fns[1:6]

# print(fns)

# fastq
if (paired_reads){
  
  fns_R1 <- fns[str_detect( basename(fns),R1_identifier)]
  fns_R2 <- fns[str_detect( basename(fns),R2_identifier)]
  
  fns_R1.fastq <- fns_R1
  fns_R2.fastq <- fns_R2
  
  # filters with reads with N removed
  # fns_R1.filtN <- str_c(dir_fastqN, basename(fns_R1))  # Put N-filterd files in filtN/ subdirectory
  # fns_R2.filtN <- str_c(dir_fastqN, basename(fns_R2))
  
  # after cutadapt
  fns_R1.cut <- str_c(dir_cutadapt, basename(fns_R1))  # Put files in /cutadapt subdirectory
  fns_R2.cut <- str_c(dir_cutadapt, basename(fns_R2))
  
  # after FilterAndTrim
  fns_R1.filt <- str_c(dir_filtered, basename(fns_R1))
  fns_R2.filt <- str_c(dir_filtered, basename(fns_R2))
  
  sample.names <- str_split(basename(fns_R1), pattern = file_name_separator, simplify = TRUE)
  
}  else {
  if (sequencer == "Illumina") {
    fns_R1 <- fns[str_detect( basename(fns),R1_identifier)]
    fns <- fns_R1
  }
  
  fns.fastq <- fns
  #  fns.filtN <- str_c(dir_fastqN, basename(fns))
  fns.cut <- str_c(dir_cutadapt, basename(fns))
  fns.filt <- str_c(dir_filtered, basename(fns))
  # print(fns)
  sample.names <- str_split(basename(fns), pattern = file_name_separator, simplify = TRUE)
}
# print(sample.names)
sample.names <- sample.names[,1]
sample.names <- str_sub(sample.names, start=sample.names_first_character)
sample.names

# =================================
#   Get the number of reads in each file (just R1 if paired_reads)
# =================================

if (do_summary){
  
  if (paired_reads){
    fns.summary <- fns_R1 
  } else {
    fns.summary <- fns
  }
  
  summary <- data.frame()
  for(i in 1:length(fns.summary)) {
    # For the next line to work needs to install the latest Biostrings lib (see https://github.com/benjjneb/dada2/issues/766)
    geom <- fastq.geometry(fns.summary[i])
    summary_one_row <- data.frame (n_seq=geom[1], file_name=basename(fns.summary[i]))
    summary <- bind_rows(summary, summary_one_row)
    print(paste("Finished with file:",i ,fns.summary[i], summary_one_row$n_seq,"sequences", sep=" "))
  }
  
  write_tsv(summary, file_dataset( "_summary_raw_files.tsv"))
}


if (remove_primers){
  
  # =================================
  #   Remove primers using cutadapt 
  # =================================
  # following: https://benjjneb.github.io/dada2/ITS_workflow.html
  
  # ==========
  # Step # 1   - Create all orientations of the input primers
  # ==========
  
  allOrients <- function(primer) {
    
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
  }
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  FWD.orients
  REV.orients
  
  # ==========
  # Step # 2 - Remove any reads that contain N - This is now done by cutadapt
  # ==========
  # if (paired_reads){
  #   out_N <- filterAndTrim(fns_R1.fastq, fns_R1.filtN, fns_R2.fastq, fns_R2.filtN, 
  #                        maxN = 0, minQ = -10, multithread = multithread)
  # } else {
  #   out_N <- filterAndTrim(fns.fastq, fns.filtN,  
  #                        maxN = 0, minQ = -10, multithread = multithread)    
  # }
  # out_N
  
  # ==========
  # Step # 3 - Check primers in one sample
  # ==========
  primerHits <- function(primer, fn) {
    # Exist if primer is empty
    if (primer == "") return(0)
    # Counts number of reads in which the primer is found 
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  
  if (paired_reads){  
    primer_test <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fns_R1.fastq[[1]]), 
                         FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fns_R2.fastq[[1]]), 
                         REV.ForwardReads = sapply(REV.orients, primerHits, fn = fns_R1.fastq[[1]]), 
                         REV.ReverseReads = sapply(REV.orients, primerHits, fn = fns_R2.fastq[[1]]))
  } else {
    primer_test <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fns.fastq[[1]]),  
                         REV.ForwardReads = sapply(REV.orients, primerHits, fn = fns.fastq[[1]]))
  }
  print(primer_test) 
  
  # ==========
  # Step # 4 - Run cutadapt - This does not work under R studio....
  # ==========
  
  # CONDA cutadapt 2.3 must be activated BEFORE launching R and R cannot be launched using only R which is a CONDA environement so
  
  # $ source $CONDA3/activate cutadapt-2.3
  # $ /opt/6.x/R-3.5.1/bin/R 
  # 
  # > system2("cutadapt", args = "--version") 
  # 2.3 
  
  # system2("source", args = c("$CONDA3/activate", "cutadapt-2.3")) # does not work
  
  # reticulate::conda_list()
  # reticulate::use_condaenv(condaenv = "cutadapt-2.3", required = TRUE) # does not work
  
  cutadapt <- "cutadapt" 
  system2(cutadapt, args = "--version")
  
  FWD.RC <- dada2:::rc(FWD)
  REV.RC <- dada2:::rc(REV)
  
  if (paired_reads){   
        # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
        # R1.flags <- paste("-g", FWD, "-a", REV.RC) 
        R1.flags <- str_c("-g ", anchor , FWD) # just the FWD
        
        # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
        # R2.flags <- paste("-G", REV, "-A", FWD.RC) 
        R2.flags <- str_c("-G ", anchor, REV) 
        
        # Run Cutadapt
        
        if(do_cutadapt) {  
          for(i in seq_along(fns_R1)) {
          system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,              # -n 2 required to remove FWD and REV from reads
                                     "--trim-n",                               # Remove any N present (some Illumina have N at the end)    
                                     "-o", fns_R1.cut[i], "-p", fns_R2.cut[i], # output files
                                     fns_R1.fastq[i], fns_R2.fastq[i],         # input files
                                     "--cores=0",                              # automatic detection of number of cores
                                     "--discard-untrimmed",                    # remove all reads where primer not found
                                     "--minimum-length 50"))                   # removal reads that are too short (will cause an error in Plot Quality)   
          } 
        }
      
      # =======================  
      # Using now the cutadapt files as the starting files
      # =======================
      fns_R1 <- fns_R1.cut
      fns_R2 <- fns_R2.cut
      fns <- c(fns_R1, fns_R2)
      
      
    } else {
      # Assume all reads to be in the good orientation
      # Trim FWD and the reverse-complement of REV off of reads.  The reverse cannot be anchored)
      flags_FWD = ""
      flags_REV = ""
      if (FWD != "") flags_FWD <- str_c("-g ", anchor, FWD) 
      if (REV != "") flags_REV <- str_c("-a ", REV.RC) 
      
      # # Must run cutadapt in 2 steps to make sure both forward and reverse are removed and that all files
      # temp_file <- str_c(dir_cutadapt, "temp.fastq.gz")
      
      # Run Cutadapt
      
      if(do_cutadapt) {  
      for(i in seq_along(fns)) {
        system2(cutadapt, args = c(flags_FWD, 
                                   flags_REV,
                                   "-o",  fns.cut[i],           # output files
                                   fns.fastq[i],                # input files
                                   "--cores=0",                 # automatic detection of number of cores
                                   "--max-n=0",                 # Remove any N
                                   "--trim-n",                  # Remove N at end of reads
                                   "--discard-untrimmed"))      # remove all reads where primer not found   
        # system2(cutadapt, args = c(flags_REV, 
        #                    "-o", fns.cut[i],            # output files
        #                     temp_file))      # Keep the reads even if reverse primer not found  
        }
      }
      
      # =======================  
      # Using now the cutadapt files as the starting files
      # =======================
      fns <- fns.cut
    } 
    
  }  
  
  
  # =================================
  #   Get the number of reads in each file (just R1)
  # =================================
  
  if (do_summary){
    
    print("=== Number of reads in each file ===")  
    
    summary <- data.frame()
    
    if (paired_reads){ 
      fns.summary <- fns_R1
    } else {
      fns.summary <- fns
    }
    
    for(i in 1:length(fns.summary)) {
      geom <- fastq.geometry(fns.summary[i])
      summary_one_row <- data.frame (n_seq=geom[1], file_name=basename(fns.summary[i]))
      summary <- bind_rows(summary, summary_one_row)
      print(paste("Finished with file", fns.summary[i], ". ", round(i/length(fns.summary)*100, 2), "%", sep=""))
    }
    
    write_tsv(summary, file_dataset( "_summary_after_cutadapt.tsv"))
  }
  
  # =================================
  #   Plot quality after cutadapt
  # =================================
  if (do_plot_quality) {
    
    print("=== Plot quality ===")   
    
    for(i in 1:min(length(fns), max_plot_quality)) {
      print(str_c("i = ", i))
      p1 <- plotQualityProfile(fns[i])
      # if (i <= 2) {print(p1)}
      p1_file <- paste0(dir_qual, basename(fns[i]),".qual.pdf")
      ggsave( plot=p1, filename= p1_file, device = "pdf", width = 15, height = 15, scale=1, units="cm")
      
      read_length <- data.frame(length = width(ShortRead::readFastq(fns[i]))) # Read the fastQ file and get the length of all the reads...
      print(str_c("File before filtration", fns[i], "- Read length min=", min(read_length$length),"max=", max(read_length$length), "mean=", mean(read_length$length, na.rm=TRUE),  sep=" "))
      
      print(paste("Finished with file", fns[i], ". ", sep=""))
    }
    
  }
  
  
  # =================================
  #   Filtering
  # =================================
  if (do_filtering) {
    
    print("=== Filtering ===") 
    
    if (paired_reads){
      out <- filterAndTrim(fns_R1, fns_R1.filt, fns_R2, fns_R2.filt, 
                           maxN=0, rm.phix=TRUE,
                           truncLen=truncLen, minLen=minLen, maxEE=maxEE,truncQ=truncQ,
                           compress=TRUE, multithread = multithread_filter)
      
      fns.filt <- c(fns_R1.filt, fns_R2.filt)
      
    } else { 
      out <- filterAndTrim(fns, fns.filt, 
                           maxN=0, rm.phix=TRUE,
                           truncLen=truncLen[1], maxLen = maxLen, minLen=minLen[1], maxEE=maxEE[1],truncQ=truncQ,
                           compress=TRUE, multithread = multithread_filter)
    }  
    
    print(out)
    write_tsv(data.frame(out), file_dataset( "_summary_filtered_files.tsv")) 
    
  } else {
    if (paired_reads){
      fns_R1.filt <- fns_R1
      fns_R2.filt <- fns_R2
      fns.filt <- c(fns_R1.filt, fns_R2.filt)
      
    } else { 
      fns.filt <- fns
    }  
    
  }
  
  # =================================
  #   Dada2 block
  # =================================  
  
  if (do_dada2) {  
    
    print("=== Start Dada2 ===") 
    
    # =================================
    #   Error
    # =================================
    
    print("=== Error computing ===") 
    
    if (paired_reads){ 
      err_R1 <- learnErrors(fns_R1.filt, multithread = multithread)
      p <- plotErrors(err_R1, nominalQ=TRUE)
      p_file <- file_dataset("_LearnErrors_R1.pdf")
      ggsave( plot=p, filename= p_file, device = "pdf", 
              width = 15, height = 15, scale=1, units="cm")
      
      
      
      err_R2 <- learnErrors(fns_R2.filt, multithread = multithread)
      p <- plotErrors(err_R2, nominalQ=TRUE)
      p_file <- file_dataset("_LearnErrors_R2.pdf")
      ggsave( plot=p, filename= p_file, device = "pdf", 
              width = 15, height = 15, scale=1, units="cm")
    } else {
      err_R1 <- learnErrors(fns.filt, multithread = multithread)
      p <- plotErrors(err_R1, nominalQ=TRUE)
      p_file <- file_dataset("_LearnErrors_R1.pdf")
      ggsave( plot=p, filename= p_file, device = "pdf", 
              width = 15, height = 15, scale=1, units="cm")
    }
    
    # ================================
    # Big data set - files are processed one by one
    # See: https://benjjneb.github.io/dada2/bigdata_paired.html
    # ================================
    
    if (bigdata && paired_reads) {
      print("=== dereplicate big data sets ===") 
      print(str_c("logical expression: ", bigdata, paired_reads, bigdata && paired_reads, sep = "-" ))
      mergers <- vector("list", length(sample.names))
      names(mergers) <- sample.names 
      
      for(i in 1:length(fns_R1)) {
        print(cat("Processing file # :", i, "\n"))
        derep_R1 <- derepFastq(fns_R1.filt[i], verbose=T)
        dada_R1 <- dada(derep_R1, err=err_R1, multithread=multithread)
        derep_R2 <- derepFastq(fns_R2.filt[i], verbose=T)
        dada_R2 <- dada(derep_R2, err=err_R2, multithread=multithread)
        merger <- mergePairs(dada_R1, derep_R1,dada_R2, derep_R2)
        mergers[[i]] <- merger
      }
      
      rm(derep_R1)
      rm(derep_R2)
      
      # Construct sequence table and remove chimeras
      seqtab <- makeSequenceTable(mergers)
      
    } else { 
      # ================================
      # Small data set _ files processed together
      # ================================
      # =================================
      #   Dereplicate
      # =================================
      
      print("=== dereplicate small data sets ===") 
      
      if (paired_reads){ 
        derep_R1 <- derepFastq(fns_R1.filt, n = 1e+05, verbose=T)
        derep_R2 <- derepFastq(fns_R2.filt, n = 1e+05, verbose=T)
        
        names(derep_R1) <- sample.names
        names(derep_R2) <- sample.names
      } else {
        derep_R1 <- derepFastq(fns.filt, n = 1e+05, verbose=T)
        
        names(derep_R1) <- sample.names
      }
      
      
      # =================================
      #   dada
      # ================================= 
      
      print("=== dada ===")  
      
      
      if (paired_reads){   
        dada_R1 <- dada(derep_R1, err=err_R1, multithread = multithread)
        dada_R2 <- dada(derep_R2, err=err_R2, multithread = multithread)
        
        # dada_R1[[1]]
        # dada_R2[[1]]
      }  else {
        # following https://benjjneb.github.io/dada2/faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data
        if (sequencer == "454") {
          dada_R1 <- dada(derep_R1, err=err_R1, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32, multithread = multithread)
        } else{
          dada_R1 <- dada(derep_R1, err=err_R1, multithread = multithread)
        }
        seqtab <- makeSequenceTable(dada_R1)
      }
      # =================================
      #  Merge pairs
      # =================================  
      
      print("=== Merging pairs ===")  
      
      if (paired_reads){ 
        mergers <- mergePairs(dada_R1, derep_R1, dada_R2, derep_R2, verbose=TRUE)
        
        seqtab <- makeSequenceTable(mergers)
      } 
    }
    
    # ================================
    # End of small vs big data set
    # ================================
    
    t_seqtab <- t(seqtab)
    
    # Only takes the first max_number_asvs rows
    if(exists("max_number_asvs")) {
      if(max_number_asvs > 0) {
        t_seqtab <- head(t_seqtab, max_number_asvs)
      }
    }
    
    print(table(nchar(getSequences(seqtab))))
    
    print(sprintf("Mean asv length : %.2f", mean(nchar(getSequences(seqtab)))))
    
    
    # =================================
    #  Remove chimera
    # =================================  
    
    print("=== Remove Chimera ===")  
    
    seqtab.nochim <- removeBimeraDenovo(seqtab, method=method_chimera, multithread=multithread, verbose=TRUE)
    
    p <- ggplot(data.frame(seq_length=nchar(getSequences(seqtab.nochim)))) +
      geom_histogram(aes(x=seq_length)) +
      ggtitle(str_c("Number of asv: ", ncol(seqtab.nochim)))
    p_file <- file_dataset("_asv_length_hist.pdf")
    ggsave( plot=p, filename= p_file, device = "pdf", 
            width = 15, height = 15, scale=1, units="cm")  
    #
    # =================================
    #  saveRDS
    # ================================= 
    
    saveRDS(seqtab.nochim, file_dataset("_seqtab.nochim.rds")) 
    
    # =================================
    #  Compile number of reads at each step
    # ================================= 
    
    print("=== Compile number of reads at each step ===")  
    
    getN <- function(x) sum(getUniques(x))
    
    if(paired_reads){
      if (bigdata) {
        # if big data dada_R1 cannot have the stats for dada_R1
        track <- cbind(  sapply(mergers, getN), 
                         rowSums(seqtab), 
                         rowSums(seqtab.nochim))
        colnames(track) <- c("merged", "tabled", "nonchim")
      } else {
        track <- cbind(sapply(dada_R1, getN), 
                       sapply(mergers, getN), 
                       rowSums(seqtab), 
                       rowSums(seqtab.nochim))
        colnames(track) <- c("denoised", "merged", "tabled", "nonchim")
      }
    } else {
      track <- cbind(sapply(dada_R1, getN), 
                     rowSums(seqtab), 
                     rowSums(seqtab.nochim))    
      colnames(track) <- c("denoised", "tabled", "nonchim")
    }
    
    
    track <- data.frame(track) %>% 
      mutate(file_code = sample.names)
    
    write_tsv(track, file_dataset("_summary_dada2.txt"))
    
    # =================================
    #    Write fasta file without taxo
    # =================================   
    
    seqtab.nochim_trans <- as.data.frame(t(seqtab.nochim)) %>% 
      rownames_to_column(var = "sequence") %>%
      rowid_to_column(var = "asv_number") %>%
      mutate(asv_code = sprintf("asv_%03d_%05d", dataset_id, asv_number)) %>% 
      mutate(sequence = str_replace_all(sequence, "(-|\\.)",""))
    
    seq_out <- Biostrings::DNAStringSet(seqtab.nochim_trans$sequence)
    names(seq_out) <- seqtab.nochim_trans$asv_code
    Biostrings::writeXStringSet(seq_out, file_dataset("_no_taxo.fasta"), 
                                compress=FALSE, width = 20000)
  }
  
  
  
  if (do_taxo) {
    
    print("=== Assigning Taxonomy ===")
    
    # =================================
    #  Reload the files
    # =================================
    
    seqtab.nochim <- readRDS(file_dataset("_seqtab.nochim.rds"))
    
    # =================================
    #   Assign taxonomy - do by 5000 to prevent memory explosion
    # =================================  
    
    taxa_list <- list()
    
    n_asvs <- ncol(seqtab.nochim) 
    
    
    # Create the boundary of each slice
    slices = c(seq(from = 1, n_asvs, by=taxo_slice_size), n_asvs)
    #[1]  1  4  7 10 10
    
    # Remove the last slice if repeated
    slices <- unique(slices)
    # [1]  1  4  7 10
    
    for (i in 1:(length(slices)-1)){
      
      print(cat("Taxo slice = ", i, "\n"))
      
      seq_one <- seqtab.nochim[,slices[i]:slices[i+1]]
      
      taxa_one <- assignTaxonomy(seqs=seq_one,
                                 refFasta=pr2_file,
                                 taxLevels = tax_levels,
                                 minBoot = 0, outputBootstraps = TRUE,
                                 verbose = TRUE,
                                 multithread = multithread)
      boot_one <- data.frame(taxa_one$boot) %>%
        rename_all(funs(str_c(.,"_boot")))
      taxa_one <- data.frame(taxa_one$tax)  %>% 
        rownames_to_column(var = "sequence")
      
      taxa_list[[i]] <- bind_cols(taxa_one, boot_one)
    }
    
    taxa.df <- purrr::reduce(taxa_list, bind_rows)
    
    
    # saveRDS(taxa, file_dataset("_taxa.rds")) 
    
    # =================================
    #  For debugging reload the files
    # =================================
    
    #  taxa <- readRDS(file_dataset("_taxa.rds"))
    
    # =================================
    #    Create the ASV table
    # =================================  
    
    seqtab.nochim_trans <- as.data.frame(t(seqtab.nochim)) %>% 
      rownames_to_column(var = "sequence") %>%
      rowid_to_column(var = "asv_number") %>%
      mutate(asv_code = sprintf("asv_%03d_%05d", dataset_id, asv_number)) %>% 
      mutate(sequence = str_replace_all(sequence, "(-|\\.)","")) %>% 
      left_join(taxa.df) 
    
    write_tsv(seqtab.nochim_trans, file_dataset("_dada2.tsv"), na="")
    
    # =================================
    #    Create tables for import into database
    # =================================  
    
    metapr2_asv <- seqtab.nochim_trans %>% 
      mutate(gene = gene, gene_region = gene_region, organelle = organelle, dataset_id = dataset_id) %>% 
      select(asv_code,sequence, asv_code:dataset_id)
    
    metapr2_asv$sequence_hash = purrr::map_chr(metapr2_asv$sequence,digest::sha1)
    
    write_tsv(metapr2_asv, file_dataset("_metapr2_asv.txt"), na="")
    
    metapr2_asv_abundance <- seqtab.nochim_trans %>% 
      select(-asv_number, -sequence, -(kingdom:species_boot)) %>% 
      gather("file_code", "n_reads", -contains("asv_code")) %>% 
      filter(n_reads > 0 )
    
    write_tsv(metapr2_asv_abundance, file_dataset("_metapr2_asv_abundance.txt"), na="")
    
    # =================================
    #    Write fasta file with taxo
    # =================================   
    
    seq_out <- Biostrings::DNAStringSet(seqtab.nochim_trans$sequence)
    names(seq_out) <- str_c(seqtab.nochim_trans$asv_code,seqtab.nochim_trans$species, sep="|")
    Biostrings::writeXStringSet(seq_out, file_dataset("_taxo.fasta"), 
                                compress=FALSE, width = 20000)
  }
  
  
  # =================================
  #    Clean up
  # =================================  
  
  rm(list = ls())