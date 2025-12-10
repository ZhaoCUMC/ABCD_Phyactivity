#---------------------------------------------------------------------------------
# Hypothesis testing on time: mixed-effect regression models
#---------------------------------------------------------------------------------

# ---------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------

input.dir  <- "./Output/"
output.dir <- "./Output/"

index <- read.csv(
  file = paste0(input.dir, "index_lifetime_total_r.csv"),
  stringsAsFactors = FALSE
)

outcome.var <- c("MH_cbcl_scr_syn_totprob_r", "NC_nihtbx_totalcomp_uncorrected")
label.var   <- c("cbcltotal", "totalcomp")

### functions ###
summ.mod <- function(mod, outcome) {
  temp1 <- tidy(mod, conf.int = TRUE)
  temp2 <- temp1[grepl("pa.pos|pa.neg|pa.mixed", temp1$term),]
  temp2$n <- nobs(mod)
  temp2$Outcome <- outcome
  temp2$SHAP.outcome <- label.var[match(outcome, outcome.var)]
  temp2
}

# ---------------------------------------------------------------------
# MH ~ physical activity
# ---------------------------------------------------------------------

outcome <- outcome.var[1]

final.index <- index %>%
  select(sai.nm, paste0(label.var[match(outcome, outcome.var)], ".tt.idx")) %>%
  setNames(c("pa.nm", "final.index"))

pa.pos.nm   <- final.index$pa.nm[final.index$final.index == "pos"]
pa.neg.nm   <- final.index$pa.nm[final.index$final.index == "neg"]
pa.mixed.nm <- final.index$pa.nm[final.index$final.index == "mixed"]

# Compute activity time
make_time_y0 <- function(x) sub("PH_sai.p_act_", "PH_sai.p_act_time_hr.xyr.y0only_", x) 
make_time_p12 <- function(x) sub("PH_sai.p_act_", "PH_sai.p_act_time_hr.per.year_", x) 

# Function to compute log-sum time across columns
compute_logsum <- function(data, cols) {
  apply(data[, cols, drop = FALSE], 1, function(x) sum(log(x + 1), na.rm = TRUE))
}

master.file.sub <- master.file %>%
  filter(!is.na(PH_sai.p_all.activity_count)) %>%
  mutate(
    pa.pos.time        = compute_logsum(., make_time_y0(pa.pos.nm)),
    pa.neg.time        = compute_logsum(., make_time_y0(pa.neg.nm)),
    pa.mixed.time      = compute_logsum(., make_time_y0(pa.mixed.nm)),
    
    pa.pos.time.p12    = compute_logsum(., make_time_p12(pa.pos.nm)),
    pa.neg.time.p12    = compute_logsum(., make_time_p12(pa.neg.nm)),
    pa.mixed.time.p12  = compute_logsum(., make_time_p12(pa.mixed.nm))
  )

dt <- master.file.sub[which(master.file.sub$src_subject_id %in% dt.grp.tt$src_subject_id),]
dt$eventname <- factor(dt$eventname, levels=c(0,1,2,3))

temp_table <- list()

# ---------------------------------------------------------------------
# Baseline large-sample: life-time
formula1 <- paste0(outcome, " ~ pa.pos.time + pa.neg.time + pa.mixed.time + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary + (1|rel_family_id) + (1|site_id_l)")
dt.y0 <- dt[which(dt$eventname==0),]
mod <- lmer(formula1, data=dt.y0, REML=FALSE)
print(summary(mod))

temp_table[[paste0("Baseline large-sample: life-time")]] <- summ.mod(mod, outcome)

# ---------------------------------------------------------------------
# Longitudinal large-sample: past 12 month
formula2 <- paste0(outcome, " ~ pa.pos.time.p12 + pa.neg.time.p12 + pa.mixed.time.p12 + eventname + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary + (1|src_subject_id) + (1|rel_family_id) + (1|site_id_l)")
mod <- lmer(formula2, data=dt, REML=FALSE)
print(summary(mod))

temp_table[[paste0("Longitudinal large-sample: past 12 month")]] <- summ.mod(mod, outcome)

# ---------------------------------------------------------------------
# Baseline small-sample: life-time
formula3 <- paste0(outcome, " ~ pa.pos.time + pa.neg.time + pa.mixed.time + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary")
dt.test.y0 <- master.file.sub[which(master.file.sub$eventname==0 & (master.file.sub$src_subject_id %in% dt.grp.sml$src_subject_id)),]
mod <- lm(formula3, data=dt.test.y0)
print(summary(mod))

temp_table[[paste0("Baseline small-sample: life-time")]] <- summ.mod(mod, outcome)

# ---------------------------------------------------------------------
# Longitudinal small-sample: past 12 month
formula4 <- paste0(outcome, " ~ pa.pos.time.p12 + pa.neg.time.p12 + pa.mixed.time.p12 + eventname + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary + (1|src_subject_id)")
dt.test.long <- master.file.sub[which(master.file.sub$src_subject_id %in% dt.grp.sml$src_subject_id),]
mod <- lmer(formula4, data=dt.test.long, REML=FALSE)
print(summary(mod))

temp_table[[paste0("Longitudinal small-sample: past 12 month")]] <- summ.mod(mod, outcome)

final.table <- bind_rows(temp_table, .id = "Model")
write.csv(final.table, 
          file = paste0(output.dir, "summary_", 
                        outcome,
                        "_pa.time.csv"), 
          row.names = F)

# ---------------------------------------------------------------------
# Cog ~ physical activity
# ---------------------------------------------------------------------

outcome <- outcome.var[2]

final.index <- index %>%
  select(sai.nm, paste0(label.var[match(outcome, outcome.var)], ".tr.idx")) %>%
  setNames(c("pa.nm", "final.index"))

pa.pos.nm   <- final.index$pa.nm[final.index$final.index == "pos"]
pa.neg.nm   <- final.index$pa.nm[final.index$final.index == "neg"]
pa.mixed.nm <- final.index$pa.nm[final.index$final.index == "mixed"]

# Compute activity time
make_time_y0 <- function(x) sub("PH_sai.p_act_", "PH_sai.p_act_time_hr.xyr.y0only_", x) 

# Function to compute log-sum time across columns
compute_logsum <- function(data, cols) {
  apply(data[, cols, drop = FALSE], 1, function(x) sum(log(x + 1), na.rm = TRUE))
}

master.file.sub <- master.file %>%
  filter(!is.na(PH_sai.p_all.activity_count)) %>%
  mutate(
    pa.pos.time        = compute_logsum(., make_time_y0(pa.pos.nm)),
    pa.neg.time        = compute_logsum(., make_time_y0(pa.neg.nm)),
    pa.mixed.time      = compute_logsum(., make_time_y0(pa.mixed.nm))
  )

dt <- master.file.sub[which(master.file.sub$src_subject_id %in% dt.grp.tt$src_subject_id),]
dt$eventname <- factor(dt$eventname, levels=c(0,1,2,3))

temp_table <- list()

# ---------------------------------------------------------------------
## Baseline at large-sample validation: life-time
formula1 <- paste0(outcome, " ~ pa.pos.time*family_income_recode_l_binary + pa.neg.time + pa.mixed.time + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary + (1|rel_family_id) + (1|site_id_l)")
dt.y0 <- dt[which(dt$eventname==0),]
mod <- lmer(formula1, data=dt.y0, REML=TRUE)
summary(mod)

temp_table[[paste0("Baseline large-sample: life time")]] <- summ.mod(mod, outcome)

# ---------------------------------------------------------------------
# Baseline at small-sample validation: life-time
formula3 <- paste0(outcome, " ~ pa.pos.time*family_income_recode_l_binary + pa.neg.time + pa.mixed.time + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary")
dt.test.y0 <- master.file.sub[which(master.file.sub$eventname==0 & (master.file.sub$src_subject_id %in% dt.grp.sml$src_subject_id)),]
mod <- lm(formula3, data=dt.test.y0)
print(summary(mod))
temp_table[[paste0("Baseline small-sample: life time")]] <- summ.mod(mod, outcome)

final.table <- bind_rows(temp_table, .id = "Model")
write.csv(final.table, 
          file = paste0(output.dir, "summary_", 
                        outcome,
                        "_pa.time.csv"), 
          row.names = F)


