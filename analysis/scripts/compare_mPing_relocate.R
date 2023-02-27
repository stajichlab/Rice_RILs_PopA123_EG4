#!/usr/bin/env Rscript
library(tidyverse)
library(stringr)
library(dplyr)
library(fs)
infolder = "BED_files"

genomefai = read_tsv("genome/Nipponbare_IRGSP_1.0.fai",col_names = c("ACC","LENGTH","LENTOTAL","OFFSET","LINELEN")) %>% select(c(ACC)) %>% filter(substr(ACC,0,2) == "NC")
acc2chrom <- as_tibble(genomefai %>% add_column(CHROM = 1:length(genomefai$ACC)))
head(acc2chrom)
t<- read_tsv("BED_files/A123-2_1_relocate2_nonredundant.bed",skip=1,show_col_types = FALSE,
         col_names = c("ACC","START","END","ID","READCOUNT","STRAND"))
t
bed_files = fs::dir_ls(infolder,regexp = "_relocate2_nonredundant.bed$")
BEDfiles = bed_files %>% map_dfr(read_tsv,skip=1,show_col_types = FALSE,.id="source",
                                 col_names = c("ACC","START","END","ID","READCOUNT","STRAND")) %>% 
      mutate(RIL=str_replace(str_replace(source,"BED_files/",""),"_relocate2_nonredundant.bed","")) %>% 
      select(-c(source,READCOUNT,ID)) %>% inner_join(acc2chrom) %>% 
      mutate(LOC_ID=str_c(str_c("Chr",CHROM),START,sep="_"))

head(BEDfiles) 

write_tsv(BEDfiles %>% select(c(LOC_ID,RIL)),"RILs.mPing.tsv")

# set a value for having a "seen" location as 1
RILs.Locations = BEDfiles %>% select(c(LOC_ID,RIL)) %>% mutate(seen = 1)

# convert this to a table with pivot_wider
RIL.matrix <- RILs.Locations %>% pivot_wider(names_from = LOC_ID, values_from = seen,values_fill = 0)

write_tsv(RIL.matrix,"RILs.mPing.matrix.tsv")

# lets count A123 shared with RIL2-1
AParent = unique(RILs.Locations %>% select(RIL) %>% filter(str_detect(RIL,"^A")))
EParent = unique(RILs.Locations %>% select(RIL) %>% filter(str_detect(RIL,"^E")))
RILnames = unique(RILs.Locations %>% select(RIL) %>% filter(str_detect(RIL,"^R")))

AParent.Locations <- unique(RILs.Locations %>% inner_join(AParent) %>% select(LOC_ID))
EParent.Locations <- unique(RILs.Locations %>% inner_join(EParent) %>% select(LOC_ID))
OnlyRILS.Locations <- unique(RILs.Locations %>% inner_join(RILnames) %>% select(LOC_ID))

parentShare = AParent.Locations %>% inner_join(EParent.Locations)

# inRIL="RIL2-8"
# thisRILLoc = RILs.Locations %>% filter(RIL==inRIL) %>% select(LOC_ID)
# notthisRILLoc = RILs.Locations %>% filter(RIL != inRIL) %>% select(LOC_ID)
# ParentAShare = thisRILLoc %>% inner_join(AParent.Locations)
# print(thisRILLoc)
# ParentAShare = thisRILLoc %>% inner_join(AParent.Locations)
# ParentEShare = thisRILLoc %>% inner_join(EParent.Locations)
# ParentAllShare = unique(bind_rows(ParentAShare,ParentEShare))
# RILParentUnique = thisRILLoc %>% anti_join(ParentAllShare)
# RILUnique = thisRILLoc %>% anti_join(notthisRILLoc)
# print(RILUnique)
# print(RILParentUnique)
RIL_site_counts <- function(inRIL) {
    #print(inRIL)
    thisRILLoc = RILs.Locations %>% filter(RIL == inRIL) %>% select(LOC_ID)
    notthisRILLoc = RILs.Locations %>% filter(RIL != inRIL) %>% select(LOC_ID)
    #print(thisRILLoc)
    ParentAShare = thisRILLoc %>% inner_join(AParent.Locations)
    ParentEShare = thisRILLoc %>% inner_join(EParent.Locations)
    Parental = unique(bind_rows(ParentAShare,ParentEShare))
    NonParental = thisRILLoc %>% anti_join(Parental)
    SharedNonParental = NonParental %>% inner_join(notthisRILLoc)
    Unique = thisRILLoc %>% anti_join(notthisRILLoc) %>% anti_join(Parental)
    
    r = tibble(
      RIL=inRIL,
      Total=length(thisRILLoc$LOC_ID),
      Parental=length(Parental$LOC_ID),
      SharedAParents = length(ParentAShare$LOC_ID),
      SharedEParents = length(ParentEShare$LOC_ID), 
      NonParental = length(NonParental$LOC_ID),
      NonParentalShared = length(SharedNonParental$LOC_ID),
      Unique=length(Unique$LOC_ID),
      )
    return(r)
}

summary_table <- lapply(unique(RILs.Locations$RIL),RIL_site_counts) %>% bind_rows()
write_tsv(summary_table,"RILs_shared_pattern.tsv")
