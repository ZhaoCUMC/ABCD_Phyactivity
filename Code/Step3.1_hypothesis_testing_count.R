#---------------------------------------------------------------------------------
# Hypothesis testing on counts: mixed-effect regression models
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

check_model_plots <- function(model, title_name) {
  
  # --- Extract residuals and fitted -------------------------
  # resid <- residuals(model)
  # fitted <- fitted(model)
  
  # used residual from fixed effects only
  fitted <- predict(model, re.form=NA)
  resid <- model@frame[,1] - fitted
  
  p1 <- ggplot(data.frame(fitted, resid), aes(fitted, resid)) +
    geom_point(alpha = 0.4) +
    geom_hline(yintercept = 0, color = "red") +
    labs(title = title_name, x = "Fitted values", y = "Residuals")
  
  p1
}

# ---------------------------------------------------------------------
# MH ~ physical activity
# ---------------------------------------------------------------------

outcome <- outcome.var[1]

final.index <- index %>%
  select(sai.nm, paste0(label.var[match(outcome, outcome.var)], ".tr.idx")) %>%
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
    pa.mixed.count.p12 = rowSums(select(., all_of(make_p12(pa.mixed.nm))), na.rm = TRUE),
    
    # categorical versions (capped at 3)
    pa.pos.count.cat        = ifelse(pa.pos.count>3, 3, pa.pos.count),
    pa.neg.count.cat        = ifelse(pa.neg.count>3, 3, pa.neg.count),
    pa.mixed.count.cat      = ifelse(pa.mixed.count>3, 3, pa.mixed.count),
    pa.pos.count.p12.cat    = ifelse(pa.pos.count.p12>3, 3, pa.pos.count.p12),
    pa.neg.count.p12.cat    = ifelse(pa.neg.count.p12>3, 3, pa.neg.count.p12),
    pa.mixed.count.p12.cat  = ifelse(pa.mixed.count.p12>3, 3, pa.mixed.count.p12)
  )

dt <- master.file.sub[which(master.file.sub$src_subject_id %in% dt.grp.tt$src_subject_id),]
dt$eventname <- factor(dt$eventname, levels=c(0,1,2,3))

temp_table <- list()
# ---------------------------------------------------------------------
# Baseline large-sample: life-time
formula1 <- paste0(outcome, " ~ pa.pos.count.cat + pa.neg.count.cat + pa.mixed.count.cat + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary + (1|rel_family_id) + (1|site_id_l)")
dt.y0 <- dt[which(dt$eventname==0),]
mod <- lmer(formula1, data=dt.y0, REML=FALSE)
print(summary(mod))

p1 <- check_model_plots(mod, "Residuals vs Fitted - CBCL Baseline Large-sample")

temp_table[[paste0("Baseline large-sample: life-time")]] <- summ.mod(mod, outcome)

# ---------------------------------------------------------------------
# Longitudinal large-sample: past 12 month
formula2 <- paste0(outcome, " ~ pa.pos.count.p12.cat + pa.neg.count.p12.cat + pa.mixed.count.p12.cat + eventname + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary + (1|src_subject_id) + (1|rel_family_id) + (1|site_id_l)")
mod <- lmer(formula2, data=dt, REML=FALSE)
print(summary(mod))

p2 <- check_model_plots(mod, "Residuals vs Fitted - CBCL Longitudinal Large-sample")

temp_table[[paste0("Longitudinal large-sample: past 12 month")]] <- summ.mod(mod, outcome)

# ---------------------------------------------------------------------
# Baseline small-sample: life-time
formula3 <- paste0(outcome, " ~ pa.pos.count.cat + pa.neg.count.cat + pa.mixed.count.cat + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary")
dt.test.y0 <- master.file.sub[which(master.file.sub$eventname==0 & (master.file.sub$src_subject_id %in% dt.grp.sml$src_subject_id)),]
mod <- lm(formula3, data=dt.test.y0)
print(summary(mod))

resid <- residuals(mod)
fitted <- fitted(mod)

p3 <- ggplot(data.frame(fitted, resid), aes(fitted, resid)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, color = "red") +
  labs(title = "Residuals vs Fitted - CBCL Baseline Small-sample", x = "Fitted values", y = "Residuals")

temp_table[[paste0("Baseline small-sample: life-time")]] <- summ.mod(mod, outcome)

# ---------------------------------------------------------------------
# Longitudinal small-sample: past 12 month
formula4 <- paste0(outcome, " ~ pa.pos.count.p12.cat + pa.neg.count.p12.cat + pa.mixed.count.p12.cat + eventname + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary + (1|src_subject_id)")
dt.test.long <- master.file.sub[which(master.file.sub$src_subject_id %in% dt.grp.sml$src_subject_id),]
mod <- lmer(formula4, data=dt.test.long, REML=FALSE)
print(summary(mod))

p4 <- check_model_plots(mod, "Residuals vs Fitted - CBCL Longitudinal Small-sample")

temp_table[[paste0("Longitudinal small-sample: past 12 month")]] <- summ.mod(mod, outcome)

final.table <- bind_rows(temp_table, .id = "Model")
write.csv(final.table, 
          file = paste0(output.dir, "summary_", 
                        outcome,
                        "_pa.count.csv"), 
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

# Compute activity counts
make_p12 <- function(x) sub("PH_sai.p_act_", "PH_sai.p_act_p12_", x) #life-time and p12 discrepancy: no difference after y1; for y0, more involvement for life-time counts

master.file.sub <- master.file %>%
  filter(!is.na(PH_sai.p_all.activity_count)) %>%
  mutate(
    pa.pos.count   = rowSums(select(., all_of(pa.pos.nm)),   na.rm = TRUE),
    pa.neg.count   = rowSums(select(., all_of(pa.neg.nm)),   na.rm = TRUE),
    pa.mixed.count = rowSums(select(., all_of(pa.mixed.nm)), na.rm = TRUE),
    
    # categorical versions (capped at 3)
    pa.pos.count.cat        = ifelse(pa.pos.count>3, 3, pa.pos.count),
    pa.neg.count.cat        = ifelse(pa.neg.count>3, 3, pa.neg.count),
    pa.mixed.count.cat      = ifelse(pa.mixed.count>3, 3, pa.mixed.count)
  )

dt <- master.file.sub[which(master.file.sub$src_subject_id %in% dt.grp.tt$src_subject_id),]

temp_table <- list()

# ---------------------------------------------------------------------
## Baseline at large-sample validation: life-time
formula1 <- paste0(outcome, " ~ pa.pos.count.cat*family_income_recode_l_binary + pa.neg.count.cat + pa.mixed.count.cat + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary + (1|rel_family_id) + (1|site_id_l)")
dt.y0 <- dt[which(dt$eventname==0),]
mod <- lmer(formula1, data=dt.y0, REML=TRUE)
summary(mod)

p5 <- check_model_plots(mod, "Residuals vs Fitted - Cog Baseline Large-sample")

temp_table[[paste0("Baseline large-sample: life time count cat")]] <- summ.mod(mod, outcome)

pdf(file.path(paste0(output.dir, "totalcomp_interaction_count_large.sample.pdf")),
    width = 6, height = 4, onefile = FALSE)
print(plot_model(mod, type = "pred", terms = c("pa.pos.count.cat", "family_income_recode_l_binary"),
                 ci.lvl = 0.95, mdrt.values = "all", #show.values = TRUE,
                 legend.title = "Family Income",
                 title = "", #show.legend = FALSE,
                 axis.title = c("Cog-Group1 Count", "Overall Cognition")) +
        theme_bw() +
        scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = seq(0, 3, by = 1)) +
        theme(legend.position = "bottom"))
dev.off()

# without interaction
formula1 <- paste0(outcome, " ~ pa.pos.count.cat + pa.neg.count.cat + pa.mixed.count.cat + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary + (1|rel_family_id) + (1|site_id_l)")
dt.y0 <- dt[which(dt$eventname==0),]
mod <- lmer(formula1, data=dt.y0, REML=TRUE)
summary(mod)
temp_table[[paste0("Baseline large-sample: life time count cat - no interaction")]] <- summ.mod(mod, outcome)

# sub-group analysis
formula1 <- paste0(outcome, " ~ pa.pos.count.cat + pa.neg.count.cat + pa.mixed.count.cat + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + (1|rel_family_id) + (1|site_id_l)")
dt.y0 <- dt[which(dt$eventname==0 & dt$family_income_recode_l_binary==">=100k"),]
mod <- lmer(formula1, data=dt.y0, REML=TRUE)
print(summary(mod))
temp_table[[paste0("Baseline large-sample: life time count cat >=100k")]] <- summ.mod(mod, outcome)

dt.y0 <- dt[which(dt$eventname==0 & dt$family_income_recode_l_binary=="<100k"),]
mod <- lmer(formula1, data=dt.y0, REML=TRUE)
print(summary(mod))
temp_table[[paste0("Baseline large-sample: life time count cat <100k")]] <- summ.mod(mod, outcome)

# ---------------------------------------------------------------------
# Baseline at small-sample validation: life-time
formula3 <- paste0(outcome, " ~ pa.pos.count.cat*family_income_recode_l_binary + pa.neg.count.cat + pa.mixed.count.cat + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_binary")
dt.test.y0 <- master.file.sub[which(master.file.sub$eventname==0 & (master.file.sub$src_subject_id %in% dt.grp.sml$src_subject_id)),]
mod <- lm(formula3, data=dt.test.y0)
print(summary(mod))
temp_table[[paste0("Baseline small-sample: life time count cat")]] <- summ.mod(mod, outcome)

resid <- residuals(mod)
fitted <- fitted(mod)

p6 <- ggplot(data.frame(fitted, resid), aes(fitted, resid)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, color = "red") +
  labs(title = "Residuals vs Fitted - Cog Baseline Small-sample", x = "Fitted values", y = "Residuals")

final_plot <- plot_grid(
  p1, p2, 
  p3, p4, 
  p5, p6,
  ncol = 2,
  align = "v",
  labels = c("A", "B", "C", "D", "E", "F")
)

final_plot

ggsave(paste0("./Output/Model_assumption_check.jpg"), 
       plot = final_plot, width = 12, height = 12, dpi = 600)

pdf(file.path(paste0(output.dir, "totalcomp_interaction_count_small.sample.pdf")),
    width = 6, height = 4, onefile = FALSE)
print(plot_model(mod, type = "pred", terms = c("pa.pos.count.cat", "family_income_recode_l_binary"),
                 ci.lvl = 0.95, mdrt.values = "all", #show.values = TRUE,
                 legend.title = "Family Income", title = "", #show.legend = FALSE,
                 axis.title = c("Cog-Group1 Count", "Overall Cognition")) +
        theme_bw() +
        scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = seq(0, 3, by = 1)) +
        theme(legend.position = "bottom"))
dev.off()

# sub-group analysis
formula3 <- paste0(outcome, " ~ pa.pos.count.cat + pa.neg.count.cat + pa.mixed.count.cat + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode")
dt.test.y0 <- master.file.sub[which(master.file.sub$eventname==0 & (master.file.sub$src_subject_id %in% dt.grp.sml$src_subject_id) 
                                    & master.file.sub$family_income_recode_l_binary==">=100k"),]
mod <- lm(formula3, data=dt.test.y0)
print(summary(mod))
temp_table[[paste0("Baseline small-sample: life time count cat >=100k")]] <- summ.mod(mod, outcome)


dt.test.y0 <- master.file.sub[which(master.file.sub$eventname==0 & (master.file.sub$src_subject_id %in% dt.grp.sml$src_subject_id) 
                                    & master.file.sub$family_income_recode_l_binary=="<100k"),]
mod <- lm(formula3, data=dt.test.y0)
print(summary(mod))
temp_table[[paste0("Baseline small-sample: life time count cat <100k")]] <- summ.mod(mod, outcome)

final.table <- bind_rows(temp_table, .id = "Model")
write.csv(final.table, 
          file = paste0(output.dir, "summary_", 
                        outcome,
                        "_pa.count.csv"), 
          row.names = F)


