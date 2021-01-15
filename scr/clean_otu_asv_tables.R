require(dplyr)
require(stringr)
require(here)
require(compare)

#### Mothur-processed dataset ####
OTUtabwostrdino <- read.table(here("data", "alleukV4q20N0_contigs.unique.trim.unique.good.filter.unique2.precluster.pick.wodinostramenonosingletons.an.unique_list.0.02.rep.count_table"), header=T, sep="\t", row.names=1)
dim(OTUtabwostrdino)

OTUtabstrameno <- read.table(here("data", "alleukV4q20N0_contigs.unique.trim.unique.good.filter.unique2.precluster.pick.stramenopilesnosingletons.an.unique_list.0.02.rep.count_table"), header=T, sep="\t", row.names=1)

OTUtabdino <- read.table(here("data", "alleukV4q20N0_contigs.unique.trim.unique.good.filter.unique2.precluster.pick.dinophytanosingletons.an.unique_list.0.02.rep.count_table"), header=T, sep="\t", row.names=1)

# Check that the samples are in the same order in all files
compare(colnames(OTUtabwostrdino), colnames(OTUtabstrameno), ignoreOrder = F, allowAll = F)
compare(colnames(OTUtabwostrdino), colnames(OTUtabdino), ignoreOrder = F, allowAll = F)

OTUtab <- rbind(OTUtabstrameno, OTUtabwostrdino, OTUtabdino)
colnames(OTUtab) # Original column names are very long and uninformative
OTUtab <- OTUtab %>% select(-total) # We don't need this column
newnames_oldnames <- read.table(here("data", "newnames_oldnames.txt"), header = T) #New column names with shorter sample codes OR:

# Sort samples according to order of new names:

OTUtab <- OTUtab %>% select(newnames_oldnames$old_sample_name)


# Explanation of sample codes can be found in table1_samplemeta_datapaper.xlsx
# 'ni', 'nh' = niskin bottle, net haul
# 'hi', 'mi' = hiseq, miseq
# 'ex' = replicate DNA extraction
# '60' = PCR with annealing at 60 degrees
# 'bc' = same sample with different internal barcode
# 'li' = same sample in different HiSeq library


#### Combine sequencing_events that represent the same sample_sf ####
OTUtab$MP1_15sum <- OTUtab$MP1_15miseq+OTUtab$MP1_Pr_15
OTUtab$MP1_31sum <- OTUtab$MP1_31miseq+OTUtab$MP1_Pr_31
OTUtab$MP1_3sum <- OTUtab$MP1_3miseq+OTUtab$MP1_Pr_3
OTUtab$MP1_47sum <- OTUtab$MP1_47miseq+OTUtab$MP1_Pr_47
OTUtab$MP1_51sum <- OTUtab$MP1_51miseq+OTUtab$MP1_Pr_51
OTUtab$MP1_59sum <- OTUtab$MP1_59miseq+OTUtab$MP1_Pr_59

OTUtab$MP2_17sum <- OTUtab$MP2_17miseq+OTUtab$MP2_17
OTUtab$MP2_4sum <- OTUtab$MP2_4miseq+OTUtab$MP2_4

OTUtab$MP3_NP_107sum <- OTUtab$MP3_NP_107miseq+OTUtab$MP3_NP_107
OTUtab$MP3_NP_111sum <- OTUtab$MP3_NP_111miseq+OTUtab$MP3_NP_111
OTUtab$MP3_NP_114sum <- OTUtab$MP3_NP_114miseq+OTUtab$MP3_NP_114
OTUtab$MP3_NP_118sum <- OTUtab$MP3_NP_118miseq+OTUtab$MP3_NP_118
OTUtab$MP3_NP_55sum <- OTUtab$MP3_NP_55miseq+OTUtab$MP3_NP_55
OTUtab$MP3_NP_63sum <- OTUtab$MP3_NP_63miseq+OTUtab$MP3_NP_63
OTUtab$MP3_NP_7sum <- OTUtab$MP3_NP_7miseq+OTUtab$MP3_NP_7
OTUtab$MP3_NP_87sum <- OTUtab$MP3_NP_87miseq+OTUtab$MP3_NP_87
OTUtab$MP3_NP_99sum <- OTUtab$MP3_NP_99miseq+OTUtab$MP3_NP_99

OTUtab$MP4_MI_16sum <- OTUtab$MP4_MI_16miseq+OTUtab$MP4_MI_16
OTUtab$MP4_MI_3sum <- OTUtab$MP4_MI_3miseq+OTUtab$MP4_MI_3
OTUtab$MP4_MI_48sum <- OTUtab$MP4_MI_48miseq+OTUtab$MP4_MI_48
OTUtab$MP4_MI_52sum <- OTUtab$MP4_MI_52miseq+OTUtab$MP4_MI_52
OTUtab$MP4_MI_66sum <- OTUtab$MP4_MI_66miseq+OTUtab$MP4_MI_66

OTUtab$MP4_NP_30sum <- OTUtab$MP4_NP_30miseq+OTUtab$MP4_NP_30
OTUtab$MP4_NP_34sum <- OTUtab$MP4_NP_34miseq+OTUtab$MP4_NP_34
OTUtab$MP4_NP_38sum <- OTUtab$MP4_NP_38miseq+OTUtab$MP4_NP_38
OTUtab$MP4_NP_66sum <- OTUtab$MP4_NP_66miseq+OTUtab$MP4_NP_66
OTUtab$MP4_NP_74sum <- OTUtab$MP4_NP_74miseq+OTUtab$MP4_NP_74
OTUtab$MP4_NP_86sum <- OTUtab$MP4_NP_86miseq+OTUtab$MP4_NP_86

OTUtab$MP5_MI_11sum <- OTUtab$MP5_MI_11miseq+OTUtab$MP5_MI_11
OTUtab$MP5_MI_18sum <- OTUtab$MP5_MI_18miseq+OTUtab$MP5_MI_18
OTUtab$MP5_MI_29sum <- OTUtab$MP5_MI_29miseq+OTUtab$MP5_MI_29
OTUtab$MP5_MI_36sum <- OTUtab$MP5_MI_36miseq+OTUtab$MP5_MI_36
OTUtab$MP5_NP_35sum <- OTUtab$MP5_NP_35miseq+OTUtab$MP5_NP_35
OTUtab$MP5_NP_59sum <- OTUtab$MP5_NP_59miseq+OTUtab$MP5_NP_59

#####
#Compare tests
cor.test(OTUtab$MP4_NP_30, OTUtab$MP4_NP_30_32)
cor.test(OTUtab$MP4_NP_30, OTUtab$MP4_NP_30_dl)
cor.test(OTUtab$MP4_NP_30, OTUtab$MP4_NP_30miseq)
cor.test(OTUtab$MP4_MI_22, OTUtab$MP4_MI_22_60deg)
cor.test(OTUtab$MP3_NP_87, OTUtab$MP3_NP_87miseq)
cor.test(OTUtab$MP1_15miseq, OTUtab$MP1_Pr_15)
cor.test(OTUtab$MP1_3miseq, OTUtab$MP1_Pr_3)
cor.test(OTUtab$MP1_31miseq, OTUtab$MP1_Pr_31)
cor.test(OTUtab$MP1_47miseq, OTUtab$MP1_Pr_47)
cor.test(OTUtab$MP1_51miseq, OTUtab$MP1_Pr_51)
cor.test(OTUtab$MP1_59miseq, OTUtab$MP1_Pr_59)

cor.test(OTUtab$MP4_NP_38, OTUtab$MP4_NP_38_40)

#Combine rest of duplicate samples:
OTUtab$MP4_NP_38sum2 <- OTUtab$MP4_NP_38sum + OTUtab$MP4_NP_38_40
OTUtab$MP4_MI_22sum2 <- OTUtab$MP4_MI_22 + OTUtab$MP4_MI_22_60deg
OTUtab$MP4_MI_34sum2 <- OTUtab$MP4_MI_34 + OTUtab$MP4_MI_34_ob
OTUtab$MP4_NP_30sum2 <- OTUtab$MP4_NP_30sum + OTUtab$MP4_NP_30_32 + OTUtab$MP4_NP_30_dl
OTUtab$MP4_NP_66sum2 <- OTUtab$MP4_NP_66sum +OTUtab$MP4_NP_66_60deg +OTUtab$MP4_NP_66_ob
OTUtab$MP4_NP_74sum2 <- OTUtab$MP4_NP_74sum + OTUtab$MP4_NP_74_76
OTUtab$MP4_NP_86sum2 <- OTUtab$MP4_NP_86sum + OTUtab$MP4_NP_86_88










#### ljk ####
pico <- readLines("picosamplenames2.txt")
three10 <- readLines("samplenames_3_10.txt")
three180 <- readLines("samplenames_3_180.txt")
ten50 <- readLines("samplenames_10_50.txt")
fifty200 <- readLines("samplenames_50_200.txt")
net_10_50 <- readLines("nethaul_10_50.txt")
net_50_200 <- readLines("nethaul_50_200.txt")
net_all <- c(net_10_50, net_50_200)
ten200 <- readLines("samplenames_10_200.txt")

sample_names <- read.table(file = "provenavn2.txt", header = T, row.names = 1)
View(sample_names)


pico_new <- as.character(sample_names[pico,])
three10_new <- as.character(sample_names[three10,])
three180_new <- as.character(sample_names[three180,])
ten50_new <- as.character(sample_names[ten50,])
fifty200_new <- as.character(sample_names[fifty200,])
net_all_new <- as.character(sample_names[net_all,])

#### 14.01.2021 - combine dada2 processed samples ####
#### Combine sequencing_events that represent the same sample_sf ####
OTUtab$MP1_B16_0001_0.4.3_sum <- OTUtab$MP1_15miseq+OTUtab$MP1_B16_0001_0.4.3_hiseq
OTUtab$MP1_B08_1000_0.4.3 <- OTUtab$MP1_31miseq+OTUtab$MP1_B08_1000_0.4.3_hiseq
OTUtab$MP1_B16_1000_0.4.3_sum <- OTUtab$MP1_3miseq+OTUtab$MP1_B16_1000_0.4.3_hiseq
OTUtab$MP1_B16_0001_3.180_sum <- OTUtab$MP1_47miseq+OTUtab$MP1_B16_0001_3.180_hiseq
OTUtab$MP1_B08_0001_3.180sum <- OTUtab$MP1_51miseq+OTUtab$MP1_B08_0001_3.180
OTUtab$MP1_B08_500_3.180 <- OTUtab$MP1_59miseq+OTUtab$MP1_B08_500_3.180

OTUtab$MP2_M04_0020_0.4.3sum <- OTUtab$MP2_17miseq+OTUtab$MP2_M04_0020_0.4.3
OTUtab$MP2_M04_0001_0.4.3sum <- OTUtab$MP2_4miseq+OTUtab$MP2_M04_0001_0.4.3

OTUtab$MP3_P04_500_3.10sum <- OTUtab$MP3_NP_107miseq+OTUtab$MP3_P04_500_3.10
OTUtab$MP3_P04_500_0.4.3sum <- OTUtab$MP3_NP_111miseq+OTUtab$MP3_P04_500_0.4.3
OTUtab$MP3_P04_0001_3.10sum <- OTUtab$MP3_NP_114miseq+OTUtab$MP3_P04_0001_3.10
OTUtab$MP3_P04_0001_0.4.3sum <- OTUtab$MP3_NP_118miseq+OTUtab$MP3_P04_0001_0.4.3
OTUtab$MP3_P03_0001_0.4.3sum <- OTUtab$MP3_NP_55miseq+OTUtab$MP3_P03_0001_0.4.3
OTUtab$MP3_P03_447_0.4.3sum <- OTUtab$MP3_NP_63miseq+OTUtab$MP3_P03_447_0.4.3
OTUtab$MP3_P01_0020_0.4.3sum <- OTUtab$MP3_NP_7miseq+OTUtab$MP3_P01_0020_0.4.3
OTUtab$MP3_P03_0015_0.4.3sum <- OTUtab$MP3_NP_87miseq+OTUtab$MP3_P03_0015_0.4.3
OTUtab$MP3_P04_1000_3.10sum <- OTUtab$MP3_NP_99miseq+OTUtab$MP3_P04_1000_3.10

OTUtab$MP4_P05_213_50.200sum <- OTUtab$MP4_MI_16miseq+OTUtab$MP4_P05_213_50.200
OTUtab$MP4_P05_0001_50.200sum <- OTUtab$MP4_MI_3miseq+OTUtab$MP4_P05_0001_50.200
OTUtab$MP4_P06_500_10.50sum <- OTUtab$MP4_MI_48miseq+OTUtab$MP4_P06_500_10.50
OTUtab$MP4_P06_net_50.200sum <- OTUtab$MP4_MI_52miseq+OTUtab$MP4_P06_net_50.200
OTUtab$MP4_P07_0025_10.50sum <- OTUtab$MP4_MI_66miseq+OTUtab$MP4_P07_0025_10.50

OTUtab$MP4_P06_0001_0.4.3sum <- OTUtab$MP4_NP_30miseq+OTUtab$MP4_P06_0001_0.4.3
OTUtab$MP4_P06_0024_3.10sum <- OTUtab$MP4_NP_34miseq+OTUtab$MP4_P06_0024_3.10
OTUtab$MP4_P06_0024_0.4.3sum <- OTUtab$MP4_NP_38miseq+OTUtab$MP4_P06_0024_0.4.3
OTUtab$MP4_P07_0025_3.10sum <- OTUtab$MP4_NP_66miseq+OTUtab$MP4_P07_0025_3.10
OTUtab$MP4_P07_1000_3.10sum <- OTUtab$MP4_NP_74miseq+OTUtab$MP4_P07_1000_3.10
OTUtab$MP4_P07_500_0.4.3sum <- OTUtab$MP4_NP_86miseq+OTUtab$MP4_P07_500_0.4.3

OTUtab$MP5_N02_0020_50.200sum <- OTUtab$MP5_MI_11miseq+OTUtab$MP5_N02_0020_50.200
OTUtab$MP5_N03_net_10.50sum <- OTUtab$MP5_MI_18miseq+OTUtab$MP5_N03_net_10.50
OTUtab$MP5_N03_0020_50.200sum <- OTUtab$MP5_MI_29miseq+OTUtab$MP5_N03_0020_50.200
OTUtab$MP5_N04_1000_10.50sum <- OTUtab$MP5_MI_36miseq+OTUtab$MP5_N04_1000_10.50
OTUtab$MP5_N03_0020_0.4.3sum <- OTUtab$MP5_NP_35miseq+OTUtab$MP5_N03_0020_0.4.3
OTUtab$MP5_N04_0020_0.4.3sum <- OTUtab$MP5_NP_59miseq+OTUtab$MP5_N04_0020_0.4.3

#####
#Compare tests
cor.test(OTUtab$MP4_P06_0001_0.4.3, OTUtab$MP4_NP_30_32)
cor.test(OTUtab$MP4_P06_0001_0.4.3, OTUtab$MP4_NP_30_dl)
cor.test(OTUtab$MP4_P06_0001_0.4.3, OTUtab$MP4_NP_30miseq)
cor.test(OTUtab$MP4_P05_net_50.200, OTUtab$MP4_MI_22_60deg)
cor.test(OTUtab$MP3_P03_0015_0.4.3, OTUtab$MP3_NP_87miseq)
cor.test(OTUtab$MP1_B16_0001_0.4.3, OTUtab$MP1_15miseq)
cor.test(OTUtab$MP1_B16_1000_0.4.3, OTUtab$MP1_3miseq)
cor.test(OTUtab$MP1_B08_1000_0.4.3, OTUtab$MP1_31miseq)
cor.test(OTUtab$MP1_B16_0001_3.180, OTUtab$MP1_47miseq)
cor.test(OTUtab$MP1_B08_0001_3.180, OTUtab$MP1_51miseq)
cor.test(OTUtab$MP1_B08_500_3.180, OTUtab$MP1_59miseq)

cor.test(OTUtab$MP4_NP_38, OTUtab$MP4_NP_38_40)

#Combine rest of duplicate samples:
OTUtab$MP4_P06_0024_0.4.3sum2 <- OTUtab$MP4_NP_38sum + OTUtab$MP4_NP_38_40
OTUtab$MP4_P05_net_50.200sum2 <- OTUtab$MP4_MI_22 + OTUtab$MP4_MI_22_60deg
OTUtab$MP4_P06_0024_50.200sum2 <- OTUtab$MP4_MI_34 + OTUtab$MP4_MI_34_ob
OTUtab$MP4_P06_0001_0.4.3sum2 <- OTUtab$MP4_NP_30sum + OTUtab$MP4_NP_30_32 + OTUtab$MP4_NP_30_dl
OTUtab$MP4_P07_0025_3.10sum2 <- OTUtab$MP4_NP_66sum +OTUtab$MP4_NP_66_60deg +OTUtab$MP4_NP_66_ob
OTUtab$MP4_P07_1000_3.10sum2 <- OTUtab$MP4_NP_74sum + OTUtab$MP4_NP_74_76
OTUtab$MP4_P07_500_0.4.3sum2 <- OTUtab$MP4_NP_86sum + OTUtab$MP4_NP_86_88

```