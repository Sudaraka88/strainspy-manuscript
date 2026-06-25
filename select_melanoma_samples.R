# Melanoma meta - Thanks, Kshitij+Allyson (20260615)
meta_m = read.csv("data/melanoma_pooled/meta_melanoma_sm.csv")
meta_m2 = read.csv("data/melanoma_pooled/Clean_metadata_Aug5_Allyson_match_2025.txt", sep = '\t')
# Let's filter sequences to match with Fig. 4a

############# Matson
meta_tmp = meta_m[meta_m$Study_cohort == "Matson_2018",]
table(meta_tmp$Response_study_reported_RECIST)/nrow(meta_tmp) # looks okay 
meta_tmp$ORR[which(meta_tmp$Response_study_reported_RECIST == "SD") ] = "SD"
table(meta_tmp$ORR, useNA = 'always')/nrow(meta_tmp) # looks okay 

# drop SD
keep_matson = meta_tmp$X[meta_tmp$Response_study_reported_RECIST != "SD"]

############# McCulloch
meta_tmp = meta_m[meta_m$Study_cohort == "McCulloch_2022",]
meta_mc = readxl::read_xlsx("data/melanoma_pooled/mccloch_supp.xlsx", sheet = 1)
# date-of-collection vs. treatment-start-date
table(abs(lubridate::as_date(meta_mc$Date_of_Collection) - lubridate::as_date(meta_mc$Treatment_Start_Date)) <= 15)
# looks like these are the 37 samples used in the analysis
meta_tmp = meta_tmp %>%
  filter(Patient %in% meta_mc$Sample[which(abs(lubridate::as_date(meta_mc$Date_of_Collection) - lubridate::as_date(meta_mc$Treatment_Start_Date)) <= 15)]) # should this be 14 days? There is one sample exactly 15 days

table(meta_tmp$Response_study_reported_RECIST, useNA = 'always')/nrow(meta_tmp) # looks okay 
# CR + PR = 0.5675676
# PD = 0.21621622
# SD = 0.21621622 

meta_tmp$ORR[which(meta_tmp$Response_study_reported_RECIST == "SD") ] = "SD"
table(meta_tmp$ORR, useNA = 'always')/nrow(meta_tmp) # looks okay 

# drop SD
keep_mcculloch = meta_tmp$X[meta_tmp$Response_study_reported_RECIST != "SD"]

############# Frankel
meta_tmp = meta_m[meta_m$Study_cohort == "Frankel_2017",]
table(meta_tmp$Visit_type, useNA = 'always') # remove T1
meta_tmp = meta_tmp[meta_tmp$Visit_type == "Baseline",] 
table(meta_tmp$Response_study_reported_RECIST, useNA = 'always')/nrow(meta_tmp) # looks okay 
meta_tmp$ORR[which(meta_tmp$Response_study_reported_RECIST == "SD") ] = "SD"
table(meta_tmp$ORR, useNA = 'always')/nrow(meta_tmp) # looks okay 

# drop SD
keep_frankel = meta_tmp$X[meta_tmp$Response_study_reported_RECIST != "SD"]

############# Lee
meta_tmp = meta_m[meta_m$Study_simplified == "Lee_2022",] # merged cohorts

# doesn't look like it matches up - there is also 1 NA
table(meta_tmp$Response_study_reported_RECIST, useNA = 'always')

# let's look at the recommended package
x = curatedMetagenomicData::curatedMetagenomicData("2022-04-13.LeeKA_2022.relative_abundance", dryrun = F)
# Looks like the sample is there, but not ICB details, best to drop it

meta_tmp$ORR[which(meta_tmp$Response_study_reported_RECIST == "SD") ] = "SD"
table(meta_tmp$ORR)/nrow(meta_tmp) # looks okay 

# drop SD
keep_lee = meta_tmp$X[meta_tmp$Response_study_reported_RECIST != "SD"]
# First sample is the weird one, drop that out
keep_lee = keep_lee[-1]

keep_samples = c(keep_frankel, keep_lee, keep_matson, keep_mcculloch)
write.table(keep_samples, "data/melanoma_pooled/gunjur_samples.txt", col.names = F, row.names = F, quote = F)

