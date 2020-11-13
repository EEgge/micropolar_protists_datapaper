require(dplyr)
require(stringr)
require(here)
require(compare)

#### Mothur-processed dataset ####
OTUtabwostrdino <- read.table(here("data", "alleukV4q20N0_contigs.unique.trim.unique.good.filter.unique2.precluster.pick.wodinostramenonosingletons.an.unique_list.0.02.rep.count_table"), header=T, sep="\t", row.names=1)
dim(OTUtabwostrdino)

OTUtabstrameno <- read.table(here("data", "alleukV4q20N0_contigs.unique.trim.unique.good.filter.unique2.precluster.pick.stramenopilesnosingletons.an.unique_list.0.02.rep.count_table"), header=T, sep="\t", row.names=1)

OTUtabdino <- read.table(here("data", "alleukV4q20N0_contigs.unique.trim.unique.good.filter.unique2.precluster.pick.dinophytanosingletons.an.unique_list.0.02.rep.count_table"), header=T, sep="\t", row.names=1)

# Check that 
compare(colnames(OTUtabwostrdino), colnames(OTUtabstrameno), ignoreOrder = F, allowAll = F)
compare(colnames(OTUtabwostrdino), colnames(OTUtabdino), ignoreOrder = F, allowAll = F)

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
```