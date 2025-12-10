library(dplyr)
library(lmerTest)
library(MuMIn)
library(broom)
library(broom.mixed)
library(sjPlot)
library(ggplot2)
library(gridExtra)

# ---------------------------------------------------------------------
# Load Data
# ---------------------------------------------------------------------
master.file <- read.csv(file = "./Data/master.file.csv", stringsAsFactors = FALSE)

# ---------------------------------------------------------------------
# Demographic recoding
# ---------------------------------------------------------------------
master.file <- master.file %>%
  mutate(
    interview_age_yrs = interview_age / 12,
    
    demo_sex_v2_recode = factor(
      demo_sex_v2_recode,
      levels = c("Male", "Female")
    ),
    
    race_ethnicity_recode = factor(
      race_ethnicity_recode,
      levels = c("White", "Black", "Hispanic", "Other")
    ),
    
    family_income_recode_l_binary = case_when(
      family_income_recode_l %in% c("[>=200k]", "[100-200k]") ~ ">=100k",
      family_income_recode_l %in% c("[50-100k]", "[<50k]") ~ "<100k",
      TRUE ~ NA
    ) %>% factor(levels = c(">=100k", "<100k")),
    
    parental_educ_highest_recode = case_when(
      parental_educ_highest_recode %in% c("< HS Diploma", "HS Diploma/GED") ~ "High School or Below",
      TRUE ~ parental_educ_highest_recode
    ) %>% factor(
      levels = c("Post Graduate Degree", "Bachelor",
                 "Some College", "High School or Below")
    ),
    
    parent_married_l = factor(
      parent_married_l,
      levels = c(1, 0),
      labels = c("Married", "Not Married")
    ),
    
    # Interview date
    interview_date_recode = as.Date(interview_date)
  ) %>%
  filter(interview_date_recode < as.Date("2020-03-20"))  # keep pre-COVID interviews

# ---------------------------------------------------------------------
# Create analytic groups: grp1, grp2, grp3
# ---------------------------------------------------------------------

# Y0 sample with complete variables: n=11619
dt.complete <- master.file %>%
  filter(
    eventname == 0,
    !is.na(MH_cbcl_scr_syn_totprob_r),
    !is.na(NC_nihtbx_totalcomp_uncorrected),
    !is.na(PH_sai.p_all.activity_count)
  )

# Group splits
dt.grp.tr <- dt.complete %>% filter(matched_group == 2)   # training: n=5650 
dt.grp.tt <- dt.complete %>% filter(matched_group == 1)   # large-sample validation: n=5670
dt.grp.sml <- dt.complete %>% 
  filter(!src_subject_id %in% c(dt.grp.tr$src_subject_id, dt.grp.tt$src_subject_id))  # small validation: n=299

# ---------------------------------------------------------------------
# Physical Activity annotation table
# ---------------------------------------------------------------------

pa.nm <- grep("^PH_sai\\.p_act", names(master.file), value = TRUE)[1:29]

pa.nm.label <- c(
  "Ballet/Dance", "Baseball/Softball", "Basketball", "Climbing", "Field Hockey",
  "Football", "Gymnastics", "Ice Hockey", "Horseback", "Ice or Inline Skating",
  "Martial Arts", "Lacrosse", "Rugby", "Skateboarding", "Skiing/Snowboarding",
  "Soccer", "Surfing", "Swimming/Water polo", "Tennis", "Track/Running",
  "Wrestling/MMA", "Volleyball", "Yoga/Tai Chi", "Musical Instrument",
  "Drawing/Art", "Drama/Theater", "Crafts", "Competitive Games",
  "Hobbies/Collecting"
)

anno <- data.frame(
  pa.nm = pa.nm,
  pa.nm.label = pa.nm.label,
  stringsAsFactors = FALSE
)