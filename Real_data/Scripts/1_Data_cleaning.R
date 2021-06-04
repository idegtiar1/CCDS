#########################################################################################################
# Title: "Data cleaning"
# Author: ""
# Date: "3/17/21"
# Description: Using latest dataset, used for GLW paper
#########################################################################################################

### Load datasets & libraries
rm(list=ls())

library(qs) # for reading in qs files
library(tidyverse) # for data manipulation
library(dplyr) # for data manipulation

### Load data
# Save dta as qs file first time
# target_sample = read_dta("glw/new/derived/glw_sample/output/analytic.dta")
# qsave(target_sample, "/data/NYS_Medicaid/Name/target_sample.qs")
target_sample = qread("/data/NYS_Medicaid/Name/target_sample.qs")

# Neighborhood-zip crosswalk
zip_to_neighborhood_crosswalk = readr::read_csv("/Helper_files/zip_to_neighborhood_crosswalk.csv", 
                                         col_types = "fff")

### Checks
nrow(target_sample)
length(unique(target_sample$recip_id))
target_sample %>% group_by(aa_sample, ac_sample) %>% dplyr::summarize(n_obs = n(),n=length(unique(recip_id)))
target_sample %>% group_by(aa_sample, ac_sample, in_sample_6) %>% dplyr::summarize(n_obs = n(),n=length(unique(recip_id)))
length(unique(target_sample$recip_id))/length(target_sample$recip_id) # percent without multiple months of observations
mean(target_sample$aa_sample==1 | target_sample$ac_sample==1) # percent auto-enrollees or active choosers = 100%
mean((target_sample$aa_sample==1 | target_sample$ac_sample==1) & target_sample$in_sample_6==1) # 28% of observations in first 6 mo: in_sample_6 flags the first 6 months of enrollment in the post-period for recipient actually in the post-period for 6 months
length(unique(target_sample$recip_id[target_sample$in_sample_6==1]))/length(unique(target_sample$recip_id)) # 81.5% of observations have 6 months of enrollment information
unique(target_sample %>% filter(in_sample_6==1) %>% select(post_assign_month)) %>% unique() # indeed subsets to post-assignment months 0-5

### Data cleaning
## Subset to first 6 months of data post-assignment and subset to auto-assignees (rand) and active choosers (obs)
target_sample_6mo = target_sample %>% filter(0 <= post_assign_month & post_assign_month <= 5) 
nrow(target_sample_6mo)
target_sample_6mo = target_sample_6mo %>% filter(aa_sample==1 | ac_sample==1)
nrow(target_sample_6mo)
target_sample_6mo %>% group_by(in_sample_6, recip_id) %>% 
  dplyr::summarise(n_months = n(), lspending = mean(lpay0_tot)) %>% 
  ungroup() %>% group_by(in_sample_6) %>% 
  dplyr::summarize(mean_months = mean(n_months),
                   mean_lspending = mean(lspending)) # mean number of months of follow-up and mean monthly log spending


## Check for duplicate benes
target_sample_6mo = target_sample_6mo %>% arrange(recip_id, year, month) # arrange chronologically
target_sample_6mo$recip_id_post_assign_month = paste0(target_sample_6mo$recip_id,
                                                             target_sample_6mo$post_assign_month)
length(unique(target_sample_6mo$recip_id))
target_sample_6mo %>% group_by(aa_sample, ac_sample) %>% dplyr::summarize(n_obs = n(),n=length(unique(recip_id)))
target_sample_6mo %>% group_by(aa_sample, ac_sample, in_sample_6) %>% dplyr::summarize(n_obs = n(),n=length(unique(recip_id)))
length(unique(target_sample_6mo$recip_id_post_assign_month))
sum(duplicated(target_sample_6mo$recip_id_post_assign_month))
sum(duplicated(target_sample_6mo[target_sample_6mo$aa_sample==1,"recip_id_post_assign_month"]))
sum(duplicated(target_sample_6mo[target_sample_6mo$ac_sample==1,"recip_id_post_assign_month"])) # half duplicates rand/half obs
sum(duplicated(target_sample_6mo$recip_id_post_assign_month))/length(unique(target_sample_6mo$recip_id_post_assign_month)) # 1.5% have multiple enrollments
duplicates = target_sample_6mo[duplicated(target_sample_6mo$recip_id_post_assign_month),] # to remove
length(unique(duplicates$recip_id_episode))

# Remove duplicates
target_sample_6mo_no_dup = target_sample_6mo %>% distinct(recip_id_post_assign_month, .keep_all = TRUE)
sum(duplicated(target_sample_6mo_no_dup$recip_id_post_assign_month)) # check that indeed no duplicates left
target_sample_6mo_no_dup %>% group_by(aa_sample, ac_sample, in_sample_6) %>% dplyr::summarize(n_obs = n(),n=length(unique(recip_id)))

## Subset to relevant variables
X_var = c("female", "race_num", "recip_age", "aid_group", "ssi", 
          "zip", "baseline_decile")# sex, race, age, aid group, SSI eligibility, neighborhood (zip code), county, baseline spending decile
relevant_var = c("recip_id", "in_sample_6", "aa_sample", "ac_sample",
                 "post_assign_month", "pay0_tot", "plan_assign", "plan",
                 X_var)
target_sample_6mo_subset = target_sample_6mo_no_dup %>% select(relevant_var)
length(unique(target_sample_6mo_subset$recip_id))
length(unique(target_sample_6mo_subset[target_sample_6mo_subset$aa_sample==1,"recip_id"]$recip_id))
length(unique(target_sample_6mo_subset[target_sample_6mo_subset$ac_sample==1,"recip_id"]$recip_id))


## Transform data long to wide and create variable for monthly spending and baseline (month 0) plan, age, aid zip, ssi
target_sample_6mo_wide = target_sample_6mo_subset %>% pivot_wider(names_from = post_assign_month, 
                                                                   values_from = c(pay0_tot, plan, plan_assign, recip_age, 
                                                                                   aid_group, zip, ssi), 
                                                                   values_fill=list(0,"0",NA,NA,NA,NA)) %>%
  select(-c(plan_1,plan_2,plan_3,plan_4,plan_5,
            plan_assign_1,plan_assign_2,plan_assign_3,plan_assign_4,plan_assign_5,
            recip_age_1,recip_age_2,recip_age_3,recip_age_4,recip_age_5,
            aid_group_1,aid_group_2,aid_group_3,aid_group_4,aid_group_5,
            zip_1,zip_2,zip_3,zip_4,zip_5,
            ssi_1,ssi_2,ssi_3,ssi_4,ssi_5)) %>%
  rename(plan_baseline = plan_0,
         plan_assign_baseline = plan_assign_0,
         age = recip_age_0,
         aid_group = aid_group_0,
         zip = zip_0,
         ssi = ssi_0)
length(unique(target_sample_6mo_wide$recip_id))
nrow(target_sample_6mo_wide)
table(target_sample_6mo_wide[,"aa_sample"], useNA = "ifany")


# Check characteristics & missingness for those with 6 months of data
table(target_sample_6mo_wide$in_sample_6, useNA = "ifany")
table(target_sample_6mo_wide[target_sample_6mo_wide$in_sample_6==1,"aa_sample"], useNA = "ifany")
table(target_sample_6mo_wide[target_sample_6mo_wide$in_sample_6==1,"ac_sample"], useNA = "ifany")
table(target_sample_6mo_wide[target_sample_6mo_wide$in_sample_6==1,"plan_baseline"], useNA = "ifany")
table(target_sample_6mo_wide[target_sample_6mo_wide$in_sample_6==1,"female"], useNA = "ifany")
table(target_sample_6mo_wide[target_sample_6mo_wide$in_sample_6==1,"race_num"], useNA = "ifany")
table(target_sample_6mo_wide[target_sample_6mo_wide$in_sample_6==1,"age"], useNA = "ifany")
table(target_sample_6mo_wide[target_sample_6mo_wide$in_sample_6==1,"aid_group"], useNA = "ifany")
table(target_sample_6mo_wide[target_sample_6mo_wide$in_sample_6==1,"ssi"], useNA = "ifany")
table(target_sample_6mo_wide[target_sample_6mo_wide$in_sample_6==1,"baseline_decile"], useNA = "ifany")

sum(is.na(target_sample_6mo_wide[target_sample_6mo_wide$in_sample_6==1,"zip"]))


# values_fill didn't work so fill NA spending values with zeros
target_sample_6mo_wide = target_sample_6mo_wide %>% mutate_at(vars(pay0_tot_0,pay0_tot_1,pay0_tot_2,pay0_tot_3,pay0_tot_4,pay0_tot_5), ~replace(., is.na(.), 0))

## Create variable for total log spending over 6 months
target_sample_6mo_wide$pay0tot_0to5 = apply(subset(target_sample_6mo_wide, select=pay0_tot_0:pay0_tot_5),1,sum) 
target_sample_6mo_wide$lpay0tot_0to5 = log(target_sample_6mo_wide$pay0tot_0to5 + 1) # matching GLW analysis

# Create missingness indicator for baseline_decile and fill in with 1's
target_sample_6mo_wide$baseline_decile_missing = ifelse(is.na(target_sample_6mo_wide$baseline_decile),1,0)
target_sample_6mo_wide = target_sample_6mo_wide %>% mutate(baseline_decile_clean = replace_na(baseline_decile, 1))

## Check other missingness
target_sample_6mo_wide %>% filter(in_sample_6==1) %>% select(-baseline_decile) %>%  filter_all(any_vars(is.na(.))) %>% nrow()
benes_missing = target_sample_6mo_wide %>% filter(in_sample_6==1) %>% select(-baseline_decile) %>% filter_all(any_vars(is.na(.))) 
length(unique(benes_missing$recip_id))

# Use first occurrence for those 30 benes missing post_assign_month 0 observations to fill in baseline values
vars_to_merge = c("recip_id","plan_baseline","plan_assign_baseline","age","aid_group","zip","ssi")
benes_to_merge = target_sample_6mo[match(unique(benes_missing$recip_id), target_sample_6mo$recip_id),] %>% 
  select(c("recip_id","plan","plan_assign","recip_age","aid_group","zip","ssi"))
colnames(benes_to_merge) = vars_to_merge

# Function taken from https://alistaire.rbind.io/blog/coalescing-joins/ for replacing missing data
coalesce_join <- function(x, y, 
                          by = NULL, suffix = c(".x", ".y"), 
                          join = dplyr::full_join, ...) {
  joined <- join(x, y, by = by, suffix = suffix, ...)
  # names of desired output
  cols <- union(names(x), names(y))
  
  to_coalesce <- names(joined)[!names(joined) %in% cols]
  suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
  # remove suffixes and deduplicate
  to_coalesce <- unique(substr(
    to_coalesce, 
    1, 
    nchar(to_coalesce) - nchar(suffix_used)
  ))
  
  coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
    joined[[paste0(.x, suffix[1])]], 
    joined[[paste0(.x, suffix[2])]]
  ))
  names(coalesced) <- to_coalesce
  
  dplyr::bind_cols(joined, coalesced)[cols]
}

target_sample_6mo_wide[,vars_to_merge] = coalesce_join(target_sample_6mo_wide[,vars_to_merge], 
                                                   benes_to_merge[,vars_to_merge], 
                                                   by = "recip_id")
# Check that there's no longer any missingness
target_sample_6mo_wide %>% filter(in_sample_6==1) %>% select(-baseline_decile) %>%  filter_all(any_vars(is.na(.))) %>% nrow()


## Merge in neighborhood
target_sample_6mo_wide$zip = as.factor(target_sample_6mo_wide$zip)
target_sample_6mo_wide = left_join(target_sample_6mo_wide, zip_to_neighborhood_crosswalk, by = "zip")

### Collapse smaller categories (n<120) for data privacy and computational efficiency
target_sample_6mo_wide$aid_group_clean120 = ifelse(target_sample_6mo_wide$aid_group %in% "SSI_AGED", "OTHER", target_sample_6mo_wide$aid_group)
target_sample_6mo_wide$age_clean120 = ifelse(target_sample_6mo_wide$age == 65, 64, target_sample_6mo_wide$age)

## Add poverty rates by neighborhood
# Poverty rates taken from https://www.cccnewyork.org/wp-content/publications/CCCReport.ConcentratedPoverty.April-2012.pdf
# When multiple neighborhoods were listed, I took the average of the the poverty rates summarized in the publication
# What was "central", "northwest", etc. was judged off a map and linked to community districts via: https://communityprofiles.planning.nyc.gov/
poverty_rates = as.data.frame(matrix(c("Borough Park", 0.322, # Borough Park = 312 Borough Park
                                       "Bronx Park and Fordham", 0.327, # Bronx Park and Fordham = 207 205 Fordham, University Heights
                                       "Bushwick and Williamsburg", 0.285, # Bushwick and Williamsburg = 304 Bushwick
                                       "Canarsie and Flatlands", 0.114, # Canarsie and Flatlands = 318 Canarsie
                                       "Central Bronx", mean(0.254,0.211), # Central Bronx = 209 211 Unionport/Soundview, Pelham Parkway
                                       "Central Brooklyn", mean(0.307,0.259,0.256), # Central Brooklyn = 303 308 309: Crown Heights North, Crown Heights South
                                       "Central Harlem", mean(0.281,0.287), # Central Harlem = 110 109 Central Harlem, Manhattanville
                                       "Central Queens", 0.137, # Central Queens =  408 Fresh Meadows/Briarwood
                                       "Chelsea and Clinton", 0.117, # Chelsea and Clinton = 104/105 Chelsea/Clinton/Midtown
                                       "East Harlem", 0.308, # East Harlem = 111 East Harlem
                                       "East New York and New Lots", mean(0.36,0.398), # East New York and New Lots = 305 316 East New York, Brownsville 
                                       "Flatbush",  mean(0.224,0.154), # Flatbush = 314 317 Flatbush/Midwood, East Flatbush
                                       "Gramercy Park and Murray Hill", 0.07, # Gramercy Park and Murray Hill = 106 Hurray Hill/Stuyvesant
                                       "Greenpoint", 0.265, # Greenpoint = 301 Williamsburg/Greenpoint
                                       "Greenwich Village and Soho", 0.099, # Manhattan 1 & 2 (Battery Park City, Greenwich Village & Soho) 101/102
                                       "High Bridge and Morrisania", mean(0.35, 0.435), # High Bridge and Morrisania = 204 203/206 Concourse/Highbridge, Morrisania/East Tremont
                                       "Hunts Point and Mott Haven", 0.411, # Hunts Point and Mott Haven = 201/202 Mott Haven/Hunts Point
                                       "Inwood and Washington Heights", 0.195, # Inwood and Washington Heights = 112 Washington Heights
                                       "Jamaica", 0.188, # Jamaica = 412 Jamaica/St. Albans
                                       "Kingsbridge and Riverdale", 0.185, # Bronx 8 (Riverdale, Fieldston & Kingsbridge) 208
                                       "Lower East Side", 0.222, # Lower East Side = 103 Lower East Side
                                       "Lower Manhattan", 0.099, # Manhattan 1 & 2 (Battery Park City, Greenwich Village & Soho) 101/102 - NOTE REPEAT B/C RATES LUMPED TOGETHER!
                                       "Mid-Island", 0.097, # Staten Island South Beach 502
                                       "North Queens", 0.143, # North Queens = 407 Flushing
                                       "Northeast Bronx", 0.212, # Northeast Bronx 212 Williamsbridge 
                                       "Northeast Queens", mean(0.073, 0.071), # Queens Bayside 411, Queens Village 413
                                       "Northwest Brooklyn", mean(0.181,0.113), # Northwest Brooklyn 302 306 Fort Greene/Brooklyn Heights, Park Slope
                                       "Northwest Queens", mean(0.19,0.224,0.192), # Northwest Queens = 401 403 404 Astoria/Long Island City, Jackson Heights, Elmhurst/Corona
                                       #"Other", mean(0.099,0.068,0.185, 0.073, 0.071, 0.097, 0.07), # Other = rich neighborhoods 101/102 108 208 411, 413, 502, 503 Battery Park/Tribeca/Greenwich Village, Upper East Side, Riverdale, Bayside, Queens Village, South Beach, Tottenville
                                       "Port Richmond", 0.179, # Staten Island 1 Willowbrook 501 - NOTE REPEAT B/C RATES LUMPED TOGETHER!
                                       "Rockaways", 0.224, # Rockaways = 414 The Rockaways
                                       "South Shore", 0.070, # Staten Island 3 Tottenville 503
                                       "Southeast Bronx", 0.164, # Southeast Bronx 210 Throngs Neck
                                       "Southeast Queens", 0.188, # Southeast Queens = 412 Jamaica/St. Albans
                                       "Southern Brooklyn", mean(0.28,0.137), # Southern Brooklyn 313 315 Coney Island, Sheepshead Bay
                                       "Southwest Brooklyn", mean(0.153,0.14), # Southwest Brooklyn 310 311 Bat Ridge, Bensonhurst
                                       "Southwest Queens", mean(0.131,0.116),# Southwest Queens = 409 410 Woodhaven, Howard Beach
                                       "Stapleton and St. George", 0.179, # Stapleton and St. George = 501 Willowbrook
                                       "Sunset Park", 0.267, # Sunset Park = 307 Sunset Park
                                       "Upper East Side", 0.068, # Manhattan Upper East Side 108
                                       "Upper West Side", 0.104, # Upper West Side = 107 Upper West Side
                                       "West Central Queens", mean(0.171,0.097), # West Central Queens = 405 406 Ridgewood/Glendale, Rego Park/Forest Hills
                                       "West Queens", 0.122), # West Queens = 402 Sunnyside/Woodside
                                     ncol=2, byrow=T)) %>% magrittr::set_colnames(c("neighborhood", "percent_poverty"))
poverty_rates$neighborhood = as.factor(as.character(poverty_rates$neighborhood))
poverty_rates$percent_poverty = as.numeric(as.character(poverty_rates$percent_poverty))


poverty_rates$poverty_cat = ifelse(poverty_rates$percent_poverty <0.15, "low",
                                   ifelse(poverty_rates$percent_poverty <0.22, "med",
                                          ifelse(poverty_rates$percent_poverty <0.3, "med-high",
                                                 "high"))) # Create categories with roughly equal neighborhoods in each category (skewing towards more neigborhoods in low-poverty category since fewer benes are from those neighborhoods)

poverty_rates$poverty_cat = factor(poverty_rates$poverty_cat, levels=c("low","med","med-high","high"))
table(poverty_rates$poverty_cat)
plot(hist(poverty_rates$percent_poverty))

target_sample_6mo_wide = left_join(target_sample_6mo_wide, poverty_rates, by = "neighborhood")
table(target_sample_6mo_wide$poverty_cat)
plot(density(target_sample_6mo_wide[target_sample_6mo_wide$in_sample_6==1,"percent_poverty"]$percent_poverty))


## Label variables for interpretability
target_sample_6mo_wide = target_sample_6mo_wide %>% mutate(plan_baseline_clean = plyr::mapvalues(plan_baseline,
                                                                                     c("1", "2", "3", "4", "5", "6", "7", "8", "9","10"),#c(" 1"," 2"," 3"," 4"," 5"," 6"," 7"," 8"," 9","10"), 
                                                                                     c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")),
                                                 race = plyr::mapvalues(race_num,
                                                                                     c("1","2","3","4","5"), 
                                                                                     c("White non-Hispanic", "Black", "Asian or Pacific Islander", "American Indian or Alaskan Native", "Other")))

## Factors
target_sample_6mo_wide$female = as.factor(target_sample_6mo_wide$female)
target_sample_6mo_wide$race = as.factor(target_sample_6mo_wide$race)
target_sample_6mo_wide$aid_group_clean120 = as.factor(target_sample_6mo_wide$aid_group_clean120)
target_sample_6mo_wide$ssi = as.factor(target_sample_6mo_wide$ssi)
target_sample_6mo_wide$neighborhood = as.factor(target_sample_6mo_wide$neighborhood)
target_sample_6mo_wide$plan_baseline_clean = as.factor(target_sample_6mo_wide$plan_baseline_clean)


## Subset to relevant variables
X_var2 = c("female", "age_clean120", "aid_group_clean120", "ssi", 
          "neighborhood", "baseline_decile_clean", "baseline_decile_missing",
          "percent_poverty")# sex, age, aid group, SSI eligibility, neighborhood, baseline spending decile
relevant_var2 = c("recip_id", "aa_sample", "in_sample_6",  "lpay0tot_0to5", "plan_baseline_clean", X_var2)
target_sample_6mo_wide_subset = target_sample_6mo_wide %>% select(relevant_var2)


## Subset to auto-assignees and active choosers with 6 months of data
target_sample_relevant = target_sample_6mo_wide_subset %>% filter(in_sample_6==1)
nrow(target_sample_relevant)/nrow(target_sample_6mo_wide_subset) # proportion with 6 months of data = 82.6%
table(target_sample_relevant[,"aa_sample"], useNA = "ifany")

## Standardize continuous covariates by rounded version of mean/SD
mean(target_sample_relevant$age_clean120)
sd(target_sample_relevant$age_clean120)
target_sample_relevant$age_scaled = (target_sample_relevant$age_clean120 - 35)/10 # Mean age is 35; unit increase = 10 years
mean(target_sample_relevant$percent_poverty)
sd(target_sample_relevant$percent_poverty)
target_sample_relevant$percent_poverty_scaled = (target_sample_relevant$percent_poverty - 0.23)/0.1 # Mean percent poverty 23%; unit increase = 10%

## Make original age and poverty variables factors
target_sample_6mo_wide$age_clean120 = as.factor(target_sample_6mo_wide$age_clean120)
target_sample_6mo_wide$percent_poverty = as.factor(target_sample_6mo_wide$percent_poverty)

## Checks
length(unique(target_sample_relevant$recip_id))
nrow(target_sample_relevant)
target_sample_relevant %>% group_by(aa_sample) %>% summarize(n=length(unique(recip_id)))

### Create random subset of the population for testing
set.seed(123)
indices = sample(1:nrow(target_sample_relevant), 1500, replace = F)
target_sample_relevant_subset = target_sample_relevant[indices,]


## Save relevant datasets
qsave(target_sample_6mo_wide, "/data/NYS_Medicaid/Name/target_sample_6mo_wide.qs") # all AC/AA subjects -- will use to compare characteristics
qsave(target_sample_relevant, "/data/NYS_Medicaid/Name/target_sample_relevant.qs") # Restricted to those with 6 mo of observations (can measure outcome) -- will use for analysis
qsave(target_sample_relevant, "/data/NYS_Medicaid/Name/target_sample_relevant1500.qs")
