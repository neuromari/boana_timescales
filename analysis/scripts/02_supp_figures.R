################## B. pulchellus Timescales of Call Variability #################
#### author: M. Rodriguez-Santiago
#### last update: 7.8.2026
####
#### this script creates plots and tables in the supplemental materials
#### tables were subsequently modified in word for aesthetics and ease of reading 
##
##
#
#### packages needed to run this code ####
#
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(flextable)
library(officer)
#
#### Supp Fig 1: histogram distributions of ICOs ####
#
#within males
output_dir_within_calls <- "./analysis/graphs/histograms/males"
# Fixed histogram breaks for comparability
breaks_fixed <- seq(0, 10, length.out = 101)  # 100 bins
bin_width <- 0.1
# Vertical line values per frogID
vline_values <- data.frame(
  frogID = c(1281, 1303, 1314, 1325, 1346, 1552, 1573, 1594, 1605, 1626, 
             1682, 1741, 1752, 1763, 1774, 1785, 1806, 1817),
  vline_x = c(0.41, 0.57, 0.57, 0.32, 0.44, 0.6, 0.64, 0.7, 0.38, 0.76,
              0.46, 0.36, 0.41, 0.4, 0.65, 0.41, 0.41, 0.36))
# Split data by frogID
all_groups <- all_calls_bin %>%
  filter(!is.na(ICI), ICI < 3) %>%
  group_by(frogID) %>%
  group_split()

# Compute maximum count across all males for consistent y-axis
max_count <- max(sapply(all_groups, function(df) {
  hist(df$call_period, breaks = breaks_fixed, plot = FALSE)$counts
}))

summary_rows <- list()

for (df in all_groups) {
  frog_id <- df$frogID[1]
  vline_x <- vline_values$vline_x[vline_values$frogID == frog_id]
  
  h  <- hist(df$call_period, breaks = breaks_fixed, plot = FALSE)
  tr <- find_trough(df$call_period)                         
  bounds <- find_boundaries(h$counts, h$breaks, trough = tr$trough)
  
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    frogID    = frog_id,
    trough    = ifelse(is.na(tr$trough),    NA, round(tr$trough, 2)),
    bimodal   = tr$bimodal,
    sep_depth = ifelse(is.na(tr$sep_depth), NA, round(tr$sep_depth, 2)),
    dist1_min = round(bounds$dist1_min, 2),
    dist1_max = round(bounds$dist1_max, 2),
    dist2_min = ifelse(is.na(bounds$dist2_min), NA, round(bounds$dist2_min, 2)),
    dist2_max = ifelse(is.na(bounds$dist2_max), NA, round(bounds$dist2_max, 2)),
    peak      = vline_x)
  
  ann_y <- max_count * 0.95
  
  if (!is.na(tr$trough)) {
    df$dist_group <- ifelse(df$call_period < tr$trough, "dist1", "dist2")
  } else {
    df$dist_group <- "dist1"
  }
  df$dist_group <- factor(df$dist_group, levels = c("dist1", "dist2"))
  
  line_col <- "gray30"
  p_hist <- ggplot(df, aes(x = call_period, fill = dist_group)) +
    geom_histogram(breaks = breaks_fixed, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c(dist1 = "#e4cee6", dist2 = "#9f839d")) +
    geom_vline(xintercept = bounds$dist1_min, color = line_col,
               linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = bounds$dist1_max, color = line_col,
               linetype = "dashed", linewidth = 0.5) +
    {if (!is.na(bounds$dist2_min))
      geom_vline(xintercept = bounds$dist2_min, color = line_col,
                 linetype = "dashed", linewidth = 0.5)} +
    {if (!is.na(bounds$dist2_max))
      geom_vline(xintercept = bounds$dist2_max, color = line_col,
                 linetype = "dashed", linewidth = 0.5)} +
    annotate("text", x = bounds$dist1_min, y = ann_y,
             label = round(bounds$dist1_min, 2),
             hjust = 1.1, vjust = 1, size = 5, color = "gray30") +
    annotate("text", x = bounds$dist1_max, y = ann_y,
             label = round(bounds$dist1_max, 2),
             hjust = 1.1, vjust = 1, size = 5, color = "red") +
    {if (!is.na(bounds$dist2_min))
      annotate("text", x = bounds$dist2_min, y = ann_y,
               label = round(bounds$dist2_min, 2),
               hjust = 1.1, vjust = 1, size = 5, color = "gray30")} +
    {if (!is.na(bounds$dist2_max))
      annotate("text", x = bounds$dist2_max, y = ann_y,
               label = round(bounds$dist2_max, 2),
               hjust = 1.1, vjust = 1, size = 5, color = "gray30")} +
    ylab("count") +
    scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1)) +
    scale_y_continuous(limits = c(0, max_count)) +
    theme(panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.background  = element_blank(),
      axis.line         = element_line(colour = "black", linewidth = 0.5),
      axis.title.y      = element_text(size = 20),
      axis.title.x      = element_blank(),
      axis.text.x       = element_blank(),
      axis.text.y       = element_text(size = 16),
      legend.position   = "none")
  
  p_box <- ggplot(df, aes(x = call_period, y = "")) +
    geom_boxplot(outlier.alpha = 0, width = 0.5) +
    geom_point(position = position_jitter(height = 0.1), alpha = 0.3) +
    xlab("ICO (s)") +
    scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1)) +
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line        = element_line(colour = "black", linewidth = 0.5),
      axis.title.y     = element_blank(),
      axis.title.x     = element_text(size = 20),
      axis.text.x      = element_text(size = 16),
      axis.text.y      = element_blank(),
      legend.position  = "none")
  
  combo_plot <- plot_grid(p_hist, p_box, ncol = 1, align = "v", rel_heights = c(2, 1))
  
  pdf_file <- paste0(output_dir_within_calls, "hist_ico_frog_", frog_id, ".pdf")
  ggsave(pdf_file, plot = combo_plot, width = 10, height = 5, dpi = 300, limitsize = FALSE)
}
#
#### Supp Table 3: CV spectral features summary ####
#
cv_vars <- c("meandom", "maxdom","meanpeakf", "dfrange", "modindx",
             "sfm","skew", "entropy",
             "time.median", "time.IQR")
# 
# within bouts
spectral_cvw <- all_notes_final %>%
  group_by(frogID, note_type_relabel) %>%
  filter(INI < 10,
         note_type_relabel %in% c("1", "2")) %>%
  summarise(n = n(),
            across(all_of(cv_vars),
              list(mean = ~ mean(.x, na.rm = TRUE),
                   sd   = ~ sd(.x,   na.rm = TRUE),
                   sem  = ~ sem(.x), 
                   CVw  = ~ (sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE)) * 100),.names = "{.col}.{.fn}"),
            .groups = "drop")
spectral_cvw_long <- spectral_cvw %>%
  pivot_longer(cols = -c(frogID, note_type_relabel, n),
               names_to  = c("feature", ".value"),
               names_pattern = "(.+)\\.(mean|sd|sem|CVw)$") %>%
  mutate(dynamics = case_when(CVw >= 12 ~ "dynamic",
                              CVw >= 5  ~ "intermediate",
                              CVw >= 0  ~ "static"))
spectral_cvb <- spectral_cvw_long %>%
  group_by(feature, note_type_relabel) %>%
  summarise(n= n(),
            grand_mean = mean(mean, na.rm = TRUE),
            grand_sd   = sd(mean,   na.rm = TRUE),
            grand_sem  = sem(mean),
            grand_CVw  = mean(CVw,  na.rm = TRUE),
            CVb= (grand_sd / grand_mean) * 100,
            .groups = "drop")

spectral_cvb_label <- spectral_cvb %>%
  mutate(dynamics_w = case_when(grand_CVw >= 12 ~ "dynamic",
                                grand_CVw >= 5  ~ "intermediate",
                                grand_CVw >= 0  ~ "static"),
         dynamics_b = case_when(CVb >= 12 ~ "dynamic",
                                CVb >= 5  ~ "intermediate",
                                CVb >= 0  ~ "static"),
         CV_ratio = CVb / grand_CVw)
# table
cvb_table_data_spectral <- spectral_cvb_label %>%
  group_by(note_type_relabel) %>%
  select(feature, n, grand_mean, grand_sd, grand_CVw, CVb, CV_ratio)

cvb_ft_spectral <- cvb_table_data_spectral %>%
  flextable() %>%
  set_header_labels(
    note_type_relabel = "note",
    feature           = "feature",
    n                 = "N",
    grand_mean        = "mean",
    grand_sd          = "sd",
    grand_CVw         = "mean CVw",
    CVb               = "CVb (%)",
    CV_ratio          = "CVb/CVw Ratio") %>%
  colformat_double(j= c("grand_mean", "grand_sd", "grand_CVw", "CVb", "CV_ratio"),
                   digits = 2) %>%
  bg(i = ~ CVb >= 0 & CVb < 5,  j = "CVb", bg = CVColors[1]) %>%
  bg(i = ~ CVb >= 5 & CVb < 12, j = "CVb", bg = CVColors[2]) %>%
  bg(i = ~ CVb >= 12,            j = "CVb", bg = CVColors[3]) %>%
  bold(part = "header") %>%
  add_header_lines("Variation in Spectral Note Features in Bouts") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvb_ft_spectral
save_as_docx(cvb_ft_spectral, path = "./manuscript/figures/supp/supp_CVb_notes_spec_bouts.docx")
#
#
# within series
spectral_cvw_series <- all_notes_final %>%
  group_by(frogID, note_type_relabel) %>%
  filter(INI < 10,
         series_type %in% c("three", "four", "five"),
         note_type_relabel %in% c("1", "2")) %>%
  summarise(n = n(),
            across(all_of(cv_vars),
                   list(mean = ~ mean(.x, na.rm = TRUE),
                        sd   = ~ sd(.x,   na.rm = TRUE),
                        sem  = ~ sem(.x), 
                        CVw  = ~ (sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE)) * 100),.names = "{.col}.{.fn}"),
            .groups = "drop")
spectral_cvw_long_series <- spectral_cvw_series %>%
  pivot_longer(cols = -c(frogID, note_type_relabel, n),
               names_to  = c("feature", ".value"),
               names_pattern = "(.+)\\.(mean|sd|sem|CVw)$") %>%
  mutate(dynamics = case_when(CVw >= 12 ~ "dynamic",
                              CVw >= 5  ~ "intermediate",
                              CVw >= 0  ~ "static"))
spectral_cvb_series <- spectral_cvw_long_series %>%
  group_by(feature, note_type_relabel) %>%
  summarise(n= n(),
            grand_mean = mean(mean, na.rm = TRUE),
            grand_sd   = sd(mean,   na.rm = TRUE),
            grand_sem  = sem(mean),
            grand_CVw  = mean(CVw,  na.rm = TRUE),
            CVb= (grand_sd / grand_mean) * 100,
            .groups = "drop")

spectral_cvb_label_series <- spectral_cvb_series %>%
  mutate(dynamics_w = case_when(grand_CVw >= 12 ~ "dynamic",
                                grand_CVw >= 5  ~ "intermediate",
                                grand_CVw >= 0  ~ "static"),
         dynamics_b = case_when(CVb >= 12 ~ "dynamic",
                                CVb >= 5  ~ "intermediate",
                                CVb >= 0  ~ "static"),
         CV_ratio = CVb / grand_CVw)
# table
cvb_table_data_spectral_series <- spectral_cvb_label_series %>%
  group_by(note_type_relabel) %>%
  select(feature, n, grand_mean, grand_sd, grand_CVw, CVb, CV_ratio)

cvb_ft_spectral_series <- cvb_table_data_spectral_series %>%
  flextable() %>%
  set_header_labels(
    note_type_relabel = "note",
    feature           = "feature",
    n                 = "N",
    grand_mean        = "mean",
    grand_sd          = "sd",
    grand_CVw         = "mean CVw",
    CVb               = "CVb (%)",
    CV_ratio          = "CVb/CVw Ratio") %>%
  colformat_double(j= c("grand_mean", "grand_sd", "grand_CVw", "CVb", "CV_ratio"),
                   digits = 2) %>%
  bg(i = ~ CVb >= 0 & CVb < 5,  j = "CVb", bg = CVColors[1]) %>%
  bg(i = ~ CVb >= 5 & CVb < 12, j = "CVb", bg = CVColors[2]) %>%
  bg(i = ~ CVb >= 12,            j = "CVb", bg = CVColors[3]) %>%
  bold(part = "header") %>%
  add_header_lines("Variation in Spectral Note Features in Series") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvb_ft_spectral_series
save_as_docx(cvb_ft_spectral_series, path = "./manuscript/figures/supp/supp_CVb_notes_spec_series.docx")
#
#
# within series by sequence
spectral_cvw_seq <- all_notes_final %>%
  group_by(frogID, note_type_relabel, series_seq_cat) %>%
  filter(INI < 10,
         series_type %in% c("three", "four", "five"),
         series_seq_cat %in% c("1","50"),
         note_type_relabel %in% c("1", "2")) %>%
  summarise(n = n(),
            across(all_of(cv_vars),
                   list(mean = ~ mean(.x, na.rm = TRUE),
                        sd   = ~ sd(.x,   na.rm = TRUE),
                        sem  = ~ sem(.x), 
                        CVw  = ~ (sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE)) * 100),.names = "{.col}.{.fn}"),
            .groups = "drop")
spectral_cvw_long_seq <- spectral_cvw_seq %>%
  pivot_longer(cols = -c(frogID, note_type_relabel,series_seq_cat, n),
               names_to  = c("feature", ".value"),
               names_pattern = "(.+)\\.(mean|sd|sem|CVw)$") %>%
  mutate(dynamics = case_when(CVw >= 12 ~ "dynamic",
                              CVw >= 5  ~ "intermediate",
                              CVw >= 0  ~ "static"))
spectral_cvb_seq <- spectral_cvw_long_seq %>%
  group_by(feature, note_type_relabel, series_seq_cat) %>%
  summarise(n= n(),
            grand_mean = mean(mean, na.rm = TRUE),
            grand_sd   = sd(mean,   na.rm = TRUE),
            grand_sem  = sem(mean),
            grand_CVw  = mean(CVw,  na.rm = TRUE),
            CVb= (grand_sd / grand_mean) * 100,
            .groups = "drop")

spectral_cvb_label_seq <- spectral_cvb_seq %>%
  mutate(dynamics_w = case_when(grand_CVw >= 12 ~ "dynamic",
                                grand_CVw >= 5  ~ "intermediate",
                                grand_CVw >= 0  ~ "static"),
         dynamics_b = case_when(CVb >= 12 ~ "dynamic",
                                CVb >= 5  ~ "intermediate",
                                CVb >= 0  ~ "static"),
         CV_ratio = CVb / grand_CVw)
# table
cvb_table_data_spectral_seq <- spectral_cvb_label_seq %>%
  group_by(series_seq_cat, note_type_relabel) %>%
  select(feature, n, grand_mean, grand_sd, grand_CVw, CVb, CV_ratio)

cvb_ft_spectral_seq <- cvb_table_data_spectral_seq %>%
  flextable() %>%
  set_header_labels(series_seq_cat = "seq order",
                    note_type_relabel = "note",
                    feature= "feature",
                    n= "N",
                    grand_mean= "mean",
                    grand_sd= "sd",
                    grand_CVw= "mean CVw",
                    CVb= "CVb (%)",
                    CV_ratio= "CVb/CVw Ratio") %>%
  colformat_double(j= c("grand_mean", "grand_sd", "grand_CVw", "CVb", "CV_ratio"),
                   digits = 2) %>%
  bg(i = ~ CVb >= 0 & CVb < 5,  j = "CVb", bg = CVColors[1]) %>%
  bg(i = ~ CVb >= 5 & CVb < 12, j = "CVb", bg = CVColors[2]) %>%
  bg(i = ~ CVb >= 12,            j = "CVb", bg = CVColors[3]) %>%
  bold(part = "header") %>%
  add_header_lines("Variation in Spectral Note Features in Series") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvb_ft_spectral_seq
save_as_docx(cvb_ft_spectral_seq, path = "./manuscript/figures/supp/supp_CVb_notes_spec_seq.docx")
