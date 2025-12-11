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

# ---------------------------------------------------------------------
# configure family income
# ---------------------------------------------------------------------
# recode family income
# 6 levels: <25k, 25k-50k, 50k-75k, 75k-100k, 100k-200k, 200k+
  
# master.file <- master.file %>%
#   mutate(
#     family_income_recode_l_6lvl = case_when(
#       demo_comb_income_v2<=4 ~ "<25k",
#       demo_comb_income_v2<=6 ~ "25k-50k",
#       demo_comb_income_v2==7 ~ "50k-75k",
#       demo_comb_income_v2==8 ~ "75k-100k",
#       demo_comb_income_v2==9 ~ "100k-200k",
#       demo_comb_income_v2==10 ~ "200k+",
#       TRUE ~ NA
#     ) %>% factor(levels = c("200k+", "100k-200k", '75k-100k', "50k-75k", "25k-50k", "<25k"))
#   )

master.file <- master.file %>%
  mutate(
    family_income_recode_l_6lvl = case_when(
      demo_comb_income_v2<=4 ~ 1,
      demo_comb_income_v2<=6 ~ 2,
      demo_comb_income_v2==7 ~ 3,
      demo_comb_income_v2==8 ~ 4,
      demo_comb_income_v2==9 ~ 5,
      demo_comb_income_v2==10 ~ 6,
      TRUE ~ NA
    ) 
  )


table(master.file$family_income_recode_l_6lvl, master.file$eventname)

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

# ---------------------------------------------------------------------
## Baseline at large-sample validation: life-time
formula1 <- paste0(outcome, " ~ pa.pos.count.cat*family_income_recode_l_6lvl + pa.neg.count.cat + pa.mixed.count.cat + interview_age_yrs + demo_sex_v2_recode + race_ethnicity_recode + family_income_recode_l_6lvl + (1|rel_family_id) + (1|site_id_l)")
dt.y0 <- dt[which(dt$eventname==0),]
mod <- lmer(formula1, data=dt.y0, REML=TRUE)
summary(mod)

income_labels <- c("<25k", "25k–50k", "50k–75k", "75k–100k", "100k–200k", "200k+")
col_order <- c("#bdd7e7", "#6baed6", "#3182bd", "#9e9ac8", "#756bb1", "#54278f")
col_order <- c(
  "#6EC5FF",  # light blue
  "#40A3D8",  # blue
  "#34BFA1",  # teal
  "#8BC34A",  # green
  "#F9A825",  # yellow-orange
  "#E64A19"   # red-orange
)

p <- plot_model(mod, type = "pred", terms = c("pa.pos.count.cat", "family_income_recode_l_6lvl [1, 2, 3, 4, 5, 6]"),
           ci.lvl = 0.95, mdrt.values = "all", #show.values = TRUE,
           legend.title = "Family Income",
           title = "Recreational Activity x Family Income Interaction on Cognitive Function", #show.legend = FALSE,
           axis.title = c("Cog-Group1 Count", "Overall Cognition")) +
  theme_bw() +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = seq(0, 3, by = 1)) +
  scale_color_manual(
    values = col_order,
    labels = income_labels
  ) +
  scale_fill_manual(
    values = col_order,
    labels = income_labels
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.width = unit(1, "cm"))+ guides(color = guide_legend(nrow = 1))
ggsave(paste0("./Figures/Interaction_cog_famincome6lvl_v2.jpg"), 
       plot = p, width = 8, height = 5, dpi = 600)


