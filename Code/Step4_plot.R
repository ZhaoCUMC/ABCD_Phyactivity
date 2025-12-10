#---------------------------------------------------------------------------------
# Step 4: Visualization
# 12/06/2025
#---------------------------------------------------------------------------------

# ---------------------------------------------------------------------
# Figure xx: forestplot 
# ---------------------------------------------------------------------

# large sample holdout
plot.dt <- read.csv("./Output/table_for_forestplot_largesample.csv")
plot.dt$estimate.v2 <- paste0(sprintf("%.2f", plot.dt$estimate), " [", sprintf("%.2f", plot.dt$conf.low), ", ", sprintf("%.2f", plot.dt$conf.high), "]")
plot.dt$p.value.v2 <- ifelse(plot.dt$p.value<0.001, "<0.001", sprintf("%.3f", plot.dt$p.value))

tt <- plot.dt
tt$estimate.v2[which(is.na(tt$estimate))] <- ""
tt$p.value.v2[which(is.na(tt$p.value.v2))] <- ""
tabletext <- data.frame(Model=tt$Name,
                        Exposure=tt$Name.2,
                        Est=tt$estimate.v2,
                        p.value=tt$p.value.v2)
### full
pdf(file.path("./Output/Forestplot_CBCL_largesample.pdf"),
    width = 10, height = 5, onefile = FALSE)
print(tt |>
        forestplot(labeltext=tabletext, graph.pos=5,
                   mean=estimate,
                   lower=conf.low, upper=conf.high,
                   title="", xlab="Est",
                   col=fpColors(box="black", lines="black", zero = "gray50"),
                   zero=0, cex=0.9, lineheight = "auto", boxsize=0.3, colgap=unit(6,"mm"),
                   lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.2, bg = "white"
        ) |>
        fp_add_header(Model = c("Model"),
                      Exposure = c("Exposure"),
                      Est = c("Est (95%CI)"),
                      p.value = c("P Value")
        ) |>
        fp_set_zebra_style("#EFEFEF") |>
        fp_add_lines() |>
        fp_set_style(align = "llrr",
                     txt_gp = fpTxtGp(ticks = gpar(cex = 0.8),
                                      xlab  = gpar(cex = 0.8))
        ))
dev.off()

# small sample holdout
plot.dt <- read.csv("./Output/table_for_forestplot_smallsample.csv")
plot.dt$estimate.v2 <- paste0(sprintf("%.2f", plot.dt$estimate), " [", sprintf("%.2f", plot.dt$conf.low), ", ", sprintf("%.2f", plot.dt$conf.high), "]")
plot.dt$p.value.v2 <- ifelse(plot.dt$p.value<0.001, "<0.001", sprintf("%.3f", plot.dt$p.value))

tt <- plot.dt
tt$estimate.v2[which(is.na(tt$estimate))] <- ""
tt$p.value.v2[which(is.na(tt$p.value.v2))] <- ""
tabletext <- data.frame(Model=tt$Name,
                        Exposure=tt$Name.2,
                        Est=tt$estimate.v2,
                        p.value=tt$p.value.v2)
### full
pdf(file.path("./Output/Forestplot_CBCL_smallsample.pdf"),
    width = 10, height = 5, onefile = FALSE)
print(tt |>
        forestplot(labeltext=tabletext, graph.pos=5,
                   mean=estimate,
                   lower=conf.low, upper=conf.high,
                   title="", xlab="Est",
                   col=fpColors(box="black", lines="black", zero = "gray50"),
                   zero=0, cex=0.9, lineheight = "auto", boxsize=0.3, colgap=unit(6,"mm"),
                   lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.2, bg = "white"
        ) |>
        fp_add_header(Model = c("Model"),
                      Exposure = c("Exposure"),
                      Est = c("Est (95%CI)"),
                      p.value = c("P Value")
        ) |>
        fp_set_zebra_style("#EFEFEF") |>
        fp_add_lines() |>
        fp_set_style(align = "llrr",
                     txt_gp = fpTxtGp(ticks = gpar(cex = 0.8),
                                      xlab  = gpar(cex = 0.8))
        ))
dev.off()

#--------------------------------------------------------------------------------
# Reviewer 3 #7 calculate SE of SHAP values per participant-feature pair
#--------------------------------------------------------------------------------
# for cbcl
shap.cbcl<-read.csv("./SHAP_Result/cbcl_discovery/shap_mat_rf_training.n.runs.csv")

head(shap.cbcl)

feature_cols <- grep("^PH_sai", names(shap.cbcl), value = TRUE)

shap_long <- shap.cbcl %>%
  select(src_subject_id, run, all_of(feature_cols)) %>%
  pivot_longer(
    cols = all_of(feature_cols),
    names_to = "feature",
    values_to = "shap_value"
  )

shap_summary <- shap_long %>%
  group_by(src_subject_id, feature) %>%
  summarise(
    shap_mean = mean(shap_value, na.rm = TRUE),
    # shap_se   = sd(shap_value, na.rm = TRUE) / sqrt(n()),  # change back to SD
    shap_se   = sd(shap_value, na.rm = TRUE), # use this one in the final version to align with CV definition
    .groups = "drop"
  )

threshold <- 5e-2

shap_rel <- shap_summary %>%
  mutate(relative_se = ifelse(abs(shap_mean) < threshold,
                              NA,  # drop or mark as missing
                              shap_se / shap_mean))

shap_rel <- merge(shap_rel, anno %>% setNames(c("feature", "feature_anno")), all.x = TRUE)

# report this one 12/03/2025
p1 <- ggplot(shap_rel,
            aes(x = reorder(feature_anno, relative_se, median),
                y = relative_se)) +
  geom_boxplot(outlier.alpha = 0.3) +
  coord_flip() +
  labs(
    x = "",
    y = "Coefficient of Variation",
    title = "Stability of SHAP Values Across Activities (Behavioral Outcome)"
  ) +
  theme_bw(base_size = 14)

# for cognitive
shap.cog<-read.csv("./SHAP_Result/cog_discovery/shap_mat_rf_training.n.runs.csv")

head(shap.cog)

feature_cols <- grep("^PH_sai", names(shap.cog), value = TRUE)

shap_long <- shap.cog %>%
  select(src_subject_id, run, all_of(feature_cols)) %>%
  pivot_longer(
    cols = all_of(feature_cols),
    names_to = "feature",
    values_to = "shap_value"
  )

shap_summary <- shap_long %>%
  group_by(src_subject_id, feature) %>%
  summarise(
    shap_mean = mean(shap_value, na.rm = TRUE),
    # shap_se   = sd(shap_value, na.rm = TRUE) / sqrt(n()),  # SE across runs
    shap_se   = sd(shap_value, na.rm = TRUE), # report this
    .groups = "drop"
  )

threshold <- 5e-2

shap_rel <- shap_summary %>%
  mutate(relative_se = ifelse(abs(shap_mean) < threshold,
                              NA,  # drop or mark as missing
                              shap_se / shap_mean))

shap_rel <- merge(shap_rel, anno %>% setNames(c("feature", "feature_anno")), all.x = TRUE)

# report this one 12/03/2025
p2 <- ggplot(shap_rel,
       aes(x = reorder(feature_anno, relative_se, median),
           y = relative_se)) +
  geom_boxplot(outlier.alpha = 0.3) +
  coord_flip() +
  labs(
    x = "",
    y = "Coefficient of Variation",
    title = "Stability of SHAP Values Across Activities (Cognitive Outcome)"
  ) +
  theme_bw(base_size = 14)

final_plot <- plot_grid(
  p1, p2,
  ncol = 1,
  align = "v",
  labels = c("A", "B")
)

final_plot

ggsave(paste0("./Output/Relative_Uncertainty_combine.jpg"), 
       plot = final_plot, width = 12, height = 15, dpi = 600)

# ---------------------------------------------------------------------
# eFigure xx: comparison of feature importance and SHAP values 
# ---------------------------------------------------------------------

df <- readxl::read_excel("./Output/table_for_varimp.xlsx")
df_long <- df %>%
  pivot_longer(
    cols = -c(Predictors, Annotation),
    names_to = "Metric",
    values_to = "Importance"
  )

plot_metric <- function(metric_name) {
  df_long %>%
    filter(Metric == metric_name) %>%
    arrange(Importance) %>%
    mutate(Annotation = factor(Annotation, levels = Annotation)) %>%
    ggplot(aes(x = Importance, y = Annotation)) +
    geom_col(fill = "#4C72B0") +
    theme_bw() +
    labs(title = metric_name, x = "Importance", y = "") +
    theme(
      plot.title = element_text(size = 15, face = "bold"),
      axis.text.y = element_text(size = 10)
    )
}

p1 <- plot_metric("CBCL RF-based Feature Importance (Training)")
p2 <- plot_metric("CBCL SHAP-based Feature Importance (Training)")
p3 <- plot_metric("CBCL SHAP-based Feature Importance (Holdout)")
p4 <- plot_metric("Cog RF-based Feature Importance (Training)")
p5 <- plot_metric("Cog SHAP-based Feature Importance (Training)")
p6 <- plot_metric("Cog SHAP-based Feature Importance (Holdout)")

final_plot <- plot_grid(
  p1, p2, p3,
  p4, p5, p6,
  ncol = 3,
  align = "v",
  labels = c("A", "B", "C", "D", "E", "F")
)

final_plot

ggsave(paste0("./Output/Var_Imp.jpg"), 
       plot = final_plot, width = 20, height = 15, dpi = 600)

# ---------------------------------------------------------------------
# quadrant plot between MH and cognition
# ---------------------------------------------------------------------

index <- read.csv("./Output/index_lifetime_total_r.csv")
index$pa.nm <- index$sai.nm
index <- merge(anno, index)

tt <- index
tt$cbcltotal.tt <- -tt$cbcltotal.tt
#tt$pa.nm.label[which(tt$cbcltotal.tt.idx=="mixed" & tt$totalcomp.tt.idx=="mixed")] <- ""
tt$cbcltotal.tt.idx <- factor(tt$cbcltotal.tt.idx, levels = c("pos", "neg", "mixed"))
tt$totalcomp.tt.idx <- factor(tt$totalcomp.tt.idx, levels = c("pos", "neg", "mixed"))
tt$size.value <- abs(tt$cbcltotal.tt) + abs(tt$totalcomp.tt)
tt$color.value <- paste0(tt$cbcltotal.tt.idx, "_", tt$totalcomp.tt.idx)
tt$color.value <- factor(ifelse(tt$color.value=="pos_pos", 1, 
                                ifelse(tt$color.value=="neg_neg", 2, 
                                       ifelse(tt$color.value=="neg_pos", 3, 4))))


p <- ggplot(tt, aes(x=cbcltotal.tt, y=totalcomp.tt, size=size.value, color=color.value)) +
  geom_point(alpha=0.8) +
  scale_size(range = c(.5, 10), name="Mean SHAP Values") + 
  geom_label_repel(aes(label = pa.nm.label), color = "black", #angle=30,
                   #position = "jitter", hjust=0.5, vjust=1.5, 
                   size = 3, max.overlaps = 40) +
  labs(title = "", x = "CBCL Total Problems", y = "Overall Cognition") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white")) +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ylim(-0.75, 2) + xlim(-1.25, 1.5) + 
  # annotate("text", x = 0, y = Inf, label = "Positive", hjust = 0.5, vjust = 1, color = "red") +
  # annotate("text", x = 0, y = -Inf, label = "Negative", hjust = 0.5, vjust = -1, color = "blue") +
  # annotate("text", x = Inf, y = 0, label = "Positive", hjust = 1, vjust = 0, color = "red") +
  # annotate("text", x = -Inf, y = 0, label = "Negative", hjust = 0, vjust = 0, color = "blue") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
ggsave(paste0("./Output/quadrant_plot.pdf"), 
       plot = p, width = 9, height = 9)
