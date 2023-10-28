rm(list = ls())
gc()
library(dplyr)
library(data.table)
library(stringr)
library(data.table)

# download Rtools
# https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html

# download SNPlocs.Hsapiens.dbSNP155.GRCh37和SNPlocs.Hsapiens.dbSNP155.GRCh38
# http://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP155.GRCh37.html
# http://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP155.GRCh38.html

# download 1000genomes
# http://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.1000genomes.hs37d5.html


SNP_LOC_DATA <- SNPlocs.Hsapiens.dbSNP155.GRCh37::SNPlocs.Hsapiens.dbSNP155.GRCh37
Genome="BSgenome.Hsapiens.1000genomes.hs37d5"


ref_match <- FALSE

# load data without rsids
new_data <- fread("GCST.tsv.gz") %>%
  rename(CHR = chromosome, # 对应染色体
         BP = base_pair_location, # 对应pos
         A1 = other_allele,  # 对应ref——allele
         A2 = effect_allele) # 对应alt——allele

# split data
if(nrow(new_data)>100000){r<-floor(nrow(new_data)/100000)}else{r<-1}
gc()

# get rsids from BSgenome
for (i in 1:r) {
  print(paste("Totol:",r,"Now:",i))
  if(i==r){sumstats_dt<-new_data[((i-1)*100000+1):nrow(new_data),]}else{sumstats_dt<-new_data[((i-1)*100000+1):(i*100000),]}
  SNP <- CHR <- i.RefSNP_id <- BP <- A1 <- A2 <- NULL
  col_headers <- names(sumstats_dt)
  gr_snp <- data.table::copy(sumstats_dt)

if(sum(c("A1", "A2") %in% col_headers) == 2){

# remove nchar(allele) > 1
num_indels <- nrow(gr_snp[(nchar(gr_snp$A1)>1 | nchar(gr_snp$A2)>1),])
if(num_indels>0){gr_snp <- gr_snp[!(nchar(A1)>1 | nchar(A2)>1),]}
}

# remove NA of CHR BP
gr_snp <- gr_snp[complete.cases(gr_snp[, c("CHR", "BP"), with = FALSE])]

# forge CHR:POS:POS class
gr_snp <- GenomicRanges::makeGRangesFromDataFrame(gr_snp,keep.extra.columns = TRUE,seqnames.field = "CHR",start.field = "BP",end.field = "BP")

# get rsids
rsids <- data.table::setDT(data.frame(BSgenome::snpsByOverlaps(SNP_LOC_DATA, ranges = gr_snp,genome=Genome,)))

# remove the alleles in alleles_as_ambig are consistent with the reference genome.
rsids <- rsids[rsids$genome_compat]
rsids <- rsids[!duplicated(paste0(rsids$seqnames,rsids$pos)),]
data.table::setnames(rsids, "seqnames", "CHR")
data.table::setnames(rsids, "pos", "BP")

# in case there is CHR8 and chr8 - but keep sex chr as upper
rsids[, CHR := tolower(as.character(CHR))]
rsids[, CHR := gsub("x|23", "X", CHR)]
rsids[, CHR := gsub("y", "Y", CHR)]
sumstats_dt[, CHR := tolower(as.character(CHR))]
sumstats_dt[, CHR := gsub("x|23", "X", CHR)]
sumstats_dt[, CHR := gsub("y", "Y", CHR)]

# ensure bp is numeric
sumstats_dt[, BP := as.numeric(BP)]

# join on CHR BP to sumstats
data.table::setkeyv(sumstats_dt, c("CHR", "BP"))
data.table::setkeyv(rsids, c("CHR", "BP"))
sumstats_dt[rsids, SNP := i.RefSNP_id]
sumstats_dt <- sumstats_dt[complete.cases(sumstats_dt[, "SNP"]), ]

if(ref_match){
  sumstats_dt <- merge(x=sumstats_dt,y=rsids,by.x="SNP",by.y="RefSNP_id")
  sumstats_dt <- sumstats_dt[sumstats_dt$A1==sumstats_dt$ref_allele,]
}
fwrite(sumstats_dt,file = paste0(i,"sumstats.csv"))
rm(list = c("sumstats_dt","rsids","gr_snp"))
gc()
}

# combine data
rm(list ="new_data")
gc()
data <- data.frame()
file_to_bind <- dir()[str_detect(dir(),"sumstats")]
for (file in file_to_bind) {
  dat <- fread(file)
  data <- rbind(data,dat)
  file.remove(file)
  rm(list ="dat")
  gc()
}
