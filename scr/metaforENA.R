
getwd()
ENAmiseq <- read_delim("ENAchecklist_miseq.txt", delim = "\t")
head(ENAmiseq)
ENAmiseqsamp <- readLines("miseqsamp_old.txt")

enamiseqonly <- ENAmiseq %>% filter(old_saple_sf %in% ENAmiseqsamp)
View(enamiseqonly)
write.table(enamiseqonly, file = "ENAmiseqonly.txt")
