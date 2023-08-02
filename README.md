# MR-code
MR R code
### This is a code to solve the problem of extracting data online from TwoSampleMR package or 502 error in online LD.
For example：
```
exp_clumped <- read_exposure_data(filename = exp_dir,  
                                  sep = ",",  
                                  snp_col = "rsid",  
                                  beta_col = "beta",  
                                  se_col = "se",  
                                  effect_allele_col = "eff",  
                                  other_allele_col = "ref",  
                                  chr_col = "chr",  
                                  pos_col = "pos",  
                                  pval_col = "Pvalue",  
                                  phenotype_col = "bac",  
                                  clump = TRUE)  # The following problems will occur when using clump online.
```
#### Running the above code will cause the following problems：
```
Server code: 502; Server is possibly experiencing traffic, trying again...
Server code: 502; Server is possibly experiencing traffic, trying again...
Server code: 502; Server is possibly experiencing traffic, trying again...
Server code: 502; Server is possibly experiencing traffic, trying again...
Server code: 502; Server is possibly experiencing traffic, trying again...
Server code: 502; Server is possibly experiencing traffic, trying again...
Server error: 502
Failed to retrieve results from server. See error status message in the returned object and contact the developers if the problem persists.
Removing 50 of 50 variants due to LD with other variants or absence from LD reference panel # This problem will remove all SNPs.
```
### So we can add a code loop to solve this problem.
```
while(TRUE){
  message_to_next <<- TRUE
  error_to_next <<- FALSE
  try({withCallingHandlers(exp <- expr, 
                           message = function(c) if (stringr::str_detect(as.character(c),"Failed to")) message_to_next <<- FALSE)
    error_to_next <<- TRUE})
  if(message_to_next == TRUE&error_to_next == TRUE) { break }
}
```
### Please replace 'exp <- expr' with the code above, for example：
```
while(TRUE){
  message_to_next <<- TRUE
  error_to_next <<- FALSE
  try({withCallingHandlers(exp_clumped <- read_exposure_data(filename = exp_dir,
                                  sep = ",",
                                  snp_col = "rsid",
                                  beta_col = "beta",
                                  se_col = "se",
                                  effect_allele_col = "eff",
                                  other_allele_col = "ref",
                                  chr_col = "chr",
                                  pos_col = "pos",
                                  pval_col = "Pvalue",
                                  phenotype_col = "bac",
                                  clump = TRUE), 
                           message = function(c) if (stringr::str_detect(as.character(c),"Failed to")) message_to_next <<- FALSE)
    error_to_next <<- TRUE})
  if(message_to_next == TRUE&error_to_next == TRUE) { break }
}
```


