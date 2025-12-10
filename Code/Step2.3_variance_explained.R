#---------------------------------------------------------------------------------
# Step 2.2: variance the positive activity can explained
# 11/29/2025
#---------------------------------------------------------------------------------

setwd("/Users/xueweihan/Library/CloudStorage/Box-Box/Xuewei/ABCD/PhysicalActivity/")

library(dplyr)
library(lmerTest)
library(MuMIn)
library(broom)
library(broom.mixed)

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
dt.grp1 <- dt.complete %>% filter(matched_group == 1)   # large-sample validation: n=5670
dt.grp2 <- dt.complete %>% filter(matched_group == 2)   # training: n=5650 
dt.grp3 <- dt.complete %>% 
  filter(!src_subject_id %in% c(dt.grp1$src_subject_id, dt.grp2$src_subject_id))  # small validation: n=299

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


# ---------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------

input.dir  <- "./Data/Output_v4/"
output.dir <- "./Output_revision_Nov2025/"

index <- read.csv(
  file = paste0(input.dir, "index_lifetime_total.r_zero.cut_v2.csv"),
  stringsAsFactors = FALSE
)

outcome.var <- c("MH_cbcl_scr_syn_totprob_r", "NC_nihtbx_totalcomp_uncorrected")
label.var   <- c("cbcltotal", "totalcomp")

# ---------------------------------------------------------------------
# MH - Activity classification index
# ---------------------------------------------------------------------

outcome <- outcome.var[1]

final.index <- index %>%
  select(sai.nm, paste0(label.var[match(outcome, outcome.var)], ".tt.idx")) %>%
  setNames(c("pa.nm", "final.index"))

pa.pos.nm   <- final.index$pa.nm[final.index$final.index == "pos"]
pa.neg.nm   <- final.index$pa.nm[final.index$final.index == "neg"]
pa.mixed.nm <- final.index$pa.nm[final.index$final.index == "mixed"]

# Compute activity counts
make_p12 <- function(x) sub("PH_sai.p_act_", "PH_sai.p_act_p12_", x) #life-time and p12 discrepancy: no difference after y1; for y0, more involvement for life-time counts

master.file.sub <- master.file %>%
  filter(!is.na(PH_sai.p_all.activity_count)) %>%
  mutate(
    pa.pos.count   = rowSums(select(., all_of(pa.pos.nm)),   na.rm = TRUE),
    pa.neg.count   = rowSums(select(., all_of(pa.neg.nm)),   na.rm = TRUE),
    pa.mixed.count = rowSums(select(., all_of(pa.mixed.nm)), na.rm = TRUE),
    pa.pos.count.p12   = rowSums(select(., all_of(make_p12(pa.pos.nm))),   na.rm = TRUE),
    pa.neg.count.p12   = rowSums(select(., all_of(make_p12(pa.neg.nm))),   na.rm = TRUE),
    pa.mixed.count.p12 = rowSums(select(., all_of(make_p12(pa.mixed.nm))), na.rm = TRUE)
  )

# Create comparison index (ind = 1 vs 0)

# Identify subjects with data at all 3 waves (0–2)
yrcheck <- master.file.sub %>%
  filter(eventname < 3) %>%
  count(src_subject_id, name = "n_yr")

sub.yr3 <- yrcheck$src_subject_id[yrcheck$n_yr == 3]

# check with consistent involvement in positive activity
countcheck <- master.file.sub %>%
  filter(src_subject_id %in% sub.yr3, eventname <= 2) %>%
  group_by(src_subject_id) %>%
  summarize(nyr.pos2 = length(eventname[which(pa.pos.count>=2)]),
            nyr.pos1 = length(eventname[which(pa.pos.count>=1)]),
            nyr.pos0 = length(eventname[which(pa.pos.count==0)]),
            nyr.neg0 = length(eventname[which(pa.neg.count==0)]))

modeling <- function(sub.ind.1, sub.ind.0) {
  dt.tt <- master.file.sub %>%
    filter(
      eventname <= 2,
      src_subject_id %in% dt.grp1$src_subject_id # select subjects in large-sample validation
    ) %>%
    mutate(
      ind = case_when(
        src_subject_id %in% sub.ind.1 ~ 1,
        src_subject_id %in% sub.ind.0 ~ 0,
        TRUE ~ NA
      )
    ) %>%
    filter(!is.na(ind))
  
  table(dt.tt$ind, dt.tt$eventname)
  
  formula0 <- paste0("MH_cbcl_scr_syn_totprob_r ~ interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary + (1|src_subject_id) + (1|rel_family_id)")
  mod0 <- lmer(formula0, data = dt.tt, REML = TRUE)
  r2m.0 <- r.squaredGLMM(mod0)[1]
  
  formula1 <- paste0("MH_cbcl_scr_syn_totprob_r ~ ind + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary + (1|src_subject_id) + (1|rel_family_id)")
  mod1 <- lmer(formula1, data = dt.tt, REML = TRUE)
  r2m.1 <- r.squaredGLMM(mod1)[1]
  
  # Extract ind effects
  temp1 <- tidy(mod1, conf.int = TRUE) %>%
    filter(grepl("ind", term)) %>%
    mutate(
      n        = nobs(mod1),
      n.0      = length(unique(mod1@frame$src_subject_id[which(mod1@frame$ind==0)])),
      n.1      = length(unique(mod1@frame$src_subject_id[which(mod1@frame$ind==1)])),
      r2m.0    = as.numeric(r2m.0),
      var.expl = as.numeric(r2m.1 - r2m.0)
    )
  
  temp1
}

temp_table <- list()
### MH - pos>=2 & neg=0 for 3 yrs vs. pos=0 for 3 yrs
sub.ind.1 <- countcheck$src_subject_id[which(countcheck$nyr.pos2==3 & countcheck$nyr.neg0==3)]
sub.ind.0 <- countcheck$src_subject_id[which(countcheck$nyr.pos0==3)]
temp_table[["MH - pos>=2 & neg=0 for 3 yrs vs. pos=0 for 3 yrs"]] <- modeling(sub.ind.1, sub.ind.0)

### MH - pos>=1 & neg=0 for 3 yrs vs. pos=0 for 3 yrs
sub.ind.1 <- countcheck$src_subject_id[which(countcheck$nyr.pos1==3 & countcheck$nyr.neg0==3)]
sub.ind.0 <- countcheck$src_subject_id[which(countcheck$nyr.pos0==3)]
temp_table[["MH - pos>=1 & neg=0 for 3 yrs vs. pos=0 for 3 yrs"]] <- modeling(sub.ind.1, sub.ind.0)

### MH - pos>=2 & neg=0 for 3 yrs vs. pos>=2 for 3 yrs & neg>0 for any years
sub.ind.1 <- countcheck$src_subject_id[which(countcheck$nyr.pos2==3 & countcheck$nyr.neg0==3)]
sub.ind.0 <- countcheck$src_subject_id[which(countcheck$nyr.pos2==3 & countcheck$nyr.neg0!=3)]
temp_table[["MH - pos>=2 & neg=0 for 3 yrs vs. pos>=2 for 3 yrs & neg>0 for any years"]] <- modeling(sub.ind.1, sub.ind.0)

### MH - pos>=1 & neg=0 for 3 yrs vs. pos>=1 for 3 yrs & neg>0 for any years
sub.ind.1 <- countcheck$src_subject_id[which(countcheck$nyr.pos1==3 & countcheck$nyr.neg0==3)]
sub.ind.0 <- countcheck$src_subject_id[which(countcheck$nyr.pos1==3 & countcheck$nyr.neg0!=3)]
temp_table[["MH - pos>=1 & neg=0 for 3 yrs vs. pos>=1 for 3 yrs & neg>0 for any years"]] <- modeling(sub.ind.1, sub.ind.0)

### MH - pos>=2 & neg==0 for 3 yrs vs. pos=0 & neg=0 for 3 yrs
sub.ind.1 <- countcheck$src_subject_id[which(countcheck$nyr.pos2==3 & countcheck$nyr.neg0==3)]
sub.ind.0 <- countcheck$src_subject_id[which(countcheck$nyr.pos0==3 & countcheck$nyr.neg0==3)]
temp_table[["MH - pos>=2 & neg==0 for 3 yrs vs. pos=0 & neg=0 for 3 yrs"]] <- modeling(sub.ind.1, sub.ind.0)

### MH - pos>=2 for any 2 yrs & neg==0 for 3 yrs vs. pos=0 & neg=0 for 3 yrs
sub.ind.1 <- countcheck$src_subject_id[which(countcheck$nyr.pos2==2 & countcheck$nyr.neg0==3)]
sub.ind.0 <- countcheck$src_subject_id[which(countcheck$nyr.pos0==3 & countcheck$nyr.neg0==3)]
temp_table[["MH - pos>=2 for any 2 yrs & neg==0 for 3 yrs vs. pos=0 & neg=0 for 3 yrs"]] <- modeling(sub.ind.1, sub.ind.0)

### MH - pos>=2 for any 1 yrs & neg==0 for 3 yrs vs. pos=0 & neg=0 for 3 yrs
sub.ind.1 <- countcheck$src_subject_id[which(countcheck$nyr.pos2==1 & countcheck$nyr.neg0==3)]
sub.ind.0 <- countcheck$src_subject_id[which(countcheck$nyr.pos0==3 & countcheck$nyr.neg0==3)]
temp_table[["MH - pos>=2 for any 1 yrs & neg==0 for 3 yrs vs. pos=0 & neg=0 for 3 yrs"]] <- modeling(sub.ind.1, sub.ind.0)

# ---------------------------------------------------------------------
# Cog - Activity classification index
# ---------------------------------------------------------------------

outcome <- outcome.var[2]

final.index <- index %>%
  select(sai.nm, paste0(label.var[match(outcome, outcome.var)], ".tt.idx")) %>%
  setNames(c("pa.nm", "final.index"))

pa.pos.nm   <- final.index$pa.nm[final.index$final.index == "pos"]
pa.neg.nm   <- final.index$pa.nm[final.index$final.index == "neg"]
pa.mixed.nm <- final.index$pa.nm[final.index$final.index == "mixed"]

# Compute activity counts
make_p12 <- function(x) sub("PH_sai.p_act_", "PH_sai.p_act_p12_", x) #life-time and p12 discrepancy: no difference after y1; for y0, more involvement for life-time counts

master.file.sub <- master.file %>%
  filter(!is.na(PH_sai.p_all.activity_count)) %>%
  mutate(
    pa.pos.count   = rowSums(select(., all_of(pa.pos.nm)),   na.rm = TRUE),
    pa.neg.count   = rowSums(select(., all_of(pa.neg.nm)),   na.rm = TRUE),
    pa.mixed.count = rowSums(select(., all_of(pa.mixed.nm)), na.rm = TRUE),
    pa.pos.count.p12   = rowSums(select(., all_of(make_p12(pa.pos.nm))),   na.rm = TRUE),
    pa.neg.count.p12   = rowSums(select(., all_of(make_p12(pa.neg.nm))),   na.rm = TRUE),
    pa.mixed.count.p12 = rowSums(select(., all_of(make_p12(pa.mixed.nm))), na.rm = TRUE)
  )

# Create comparison index (ind = 1 vs 0)

# check with consistent involvement in positive activity
countcheck <- master.file.sub %>%
  filter(eventname ==0) %>%
  mutate(ind.pos2 = ifelse(pa.pos.count>=2, 1, 0),
            ind.pos1 = ifelse(pa.pos.count>=1, 1, 0),
            ind.pos0 = ifelse(pa.pos.count==0, 1, 0),
            ind.neg0 = ifelse(pa.neg.count==0, 1, 0))

modeling <- function(sub.ind.1, sub.ind.0) {
  dt.tt <- master.file.sub %>%
    filter(
      eventname ==0,
      src_subject_id %in% dt.grp1$src_subject_id # select subjects in large-sample validation
    ) %>%
    mutate(
      ind = case_when(
        src_subject_id %in% sub.ind.1 ~ 1,
        src_subject_id %in% sub.ind.0 ~ 0,
        TRUE ~ NA
      )
    ) %>%
    filter(!is.na(ind))
  
  table(dt.tt$ind)

  formula0 <- paste0("NC_nihtbx_totalcomp_uncorrected ~  interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary + (1|rel_family_id) + (1|site_id_l)")
  mod0 <- lmer(formula0, data = dt.tt, REML = TRUE)
  r2m.0 <- r.squaredGLMM(mod0)[1]
  
  formula1 <- paste0("NC_nihtbx_totalcomp_uncorrected ~ ind + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary + (1|rel_family_id) + (1|site_id_l)")
  mod1 <- lmer(formula1, data = dt.tt, REML = TRUE)
  r2m.1 <- r.squaredGLMM(mod1)[1]
  
  # Extract ind effects
  temp1 <- tidy(mod1, conf.int = TRUE) %>%
    filter(grepl("ind", term)) %>%
    mutate(
      n        = nobs(mod1),
      n.0      = length(mod1@frame$ind[which(mod1@frame$ind==0)]),
      n.1      = length(mod1@frame$ind[which(mod1@frame$ind==1)]),
      r2m.0    = as.numeric(r2m.0),
      var.expl = as.numeric(r2m.1 - r2m.0)
    )
  
  temp1
}

### Cog - pos>=2 & neg=0 at yr0 vs. pos=0 at yr0
sub.ind.1 <- countcheck$src_subject_id[which(countcheck$ind.pos2==1 & countcheck$ind.neg0==1)]
sub.ind.0 <- countcheck$src_subject_id[which(countcheck$ind.pos0==1)]
temp_table[["Cog - pos>=2 & neg=0 at yr0 vs. pos=0 at yr0"]] <- modeling(sub.ind.1, sub.ind.0)

### Cog - pos>=2 & neg=0 at yr0 vs. pos=0 & neg=0 at yr0
sub.ind.1 <- countcheck$src_subject_id[which(countcheck$ind.pos2==1 & countcheck$ind.neg0==1)]
sub.ind.0 <- countcheck$src_subject_id[which(countcheck$ind.pos0==1 & countcheck$ind.neg0==1)]
temp_table[["Cog - pos>=2 & neg=0 at yr0 vs. pos=0 & neg=0 at yr0"]] <- modeling(sub.ind.1, sub.ind.0)

### Cog - pos>=2 at yr0 vs. pos=0 at yr0
sub.ind.1 <- countcheck$src_subject_id[which(countcheck$ind.pos2==1)]
sub.ind.0 <- countcheck$src_subject_id[which(countcheck$ind.pos0==1)]
temp_table[["Cog - pos>=2 at yr0 vs. pos=0 at yr0"]] <- modeling(sub.ind.1, sub.ind.0)

### Cog - pos>=1 & neg=0 at yr0 vs. pos=0 at yr0
sub.ind.1 <- countcheck$src_subject_id[which(countcheck$ind.pos1==1 & countcheck$ind.neg0==1)]
sub.ind.0 <- countcheck$src_subject_id[which(countcheck$ind.pos0==1)]
temp_table[["Cog - pos>=1 & neg=0 at yr0 vs. pos=0 at yr0"]] <- modeling(sub.ind.1, sub.ind.0)

### Cog - pos>=2 & neg=0 at yr0 vs. pos>=2 & neg>0 at yr0
sub.ind.1 <- countcheck$src_subject_id[which(countcheck$ind.pos2==1 & countcheck$ind.neg0==1)]
sub.ind.0 <- countcheck$src_subject_id[which(countcheck$ind.pos2==1 & countcheck$ind.neg0==0)]
temp_table[["Cog - pos>=2 & neg=0 at yr0 vs. pos>=2 & neg>0 at yr0"]] <- modeling(sub.ind.1, sub.ind.0)

### Cog - pos>=1 & neg=0 at yr0 vs. pos>=1 & neg>0 at yr0
sub.ind.1 <- countcheck$src_subject_id[which(countcheck$ind.pos1==1 & countcheck$ind.neg0==1)]
sub.ind.0 <- countcheck$src_subject_id[which(countcheck$ind.pos1==1 & countcheck$ind.neg0==0)]
temp_table[["Cog - pos>=1 & neg=0 at yr0 vs. pos>=1 & neg>0 at yr0"]] <- modeling(sub.ind.1, sub.ind.0)

lmer_table <- bind_rows(temp_table, .id = "Model")
write.csv(lmer_table, file = "./Output_revision_Nov2025/Var_Exp_Table.csv", row.names = FALSE)
