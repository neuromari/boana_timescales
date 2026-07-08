################## B. pulchellus Timescales of Call Variability #################
#### author: M. Rodriguez-Santiago
#### last update: 7.8.2026
####
#### this script contains all calculations for call-level analyses
#### need only be called once and stored in environment
####
#### variables:
##   call_type_relabel: doublet (notes 1 and 2), squeak (note 3), solo (a singular note 1 and 2 not in doublet) 
##   call_period is difference in onsets between consec calls (regardless call type)
##   call_type: calls in calls (solo if only 1, two for 2 calls in calls, etc)
##   call_seq: sequence of call within the calls (0 = solo, 1 = first in calls, 2=second in calls, etc)
##   series_bin = # of calls within a series
##   series_seq_cat = within-series category where 
##     0: solo call
##     1: part of first call in series
##     50: part of mid call in series
##     100: part of last call in series
####
####
#
#### all packages needed for this script ####
##
library(dplyr)
library(tidyr)
library(flextable)
library(ggplot2)
library(tibble)
#
#
#### data ####
#
all_calls <- read.csv("./data/all_df_calls_final.csv")
colnames(all_calls) 

all_calls_bin <- all_calls %>%
  arrange(frogID, begin.time.s) %>%  # ensure calls are in time order
  group_by(frogID) %>%
  mutate(series_bin = case_when(series_type == "solo" ~ "solo", 
                                series_type == "two" ~ "two", 
                                series_type == "three" ~ "three",
                                series_type == "four" ~ "four+",
                                series_type == "five" ~ "four+"),
         series_seq_bin = case_when(series_seq =="0" ~ "0",
                                  series_seq =="1" ~ "1",
                                  series_seq == "2" ~ "2",
                                  series_seq == "3" ~ "3",
                                  series_seq == "4" ~ "4",
                                  series_seq == "5" ~ "4",
                                  series_seq == "6" ~ "4"),
         series_seq_cat = case_when(series_type == "solo" ~ "0", #solo call
                                  series_type == "two" & series_seq == "1" ~ "1", 
                                  series_type == "two" & series_seq == "2" ~ "100", 
                                  series_type == "three" & series_seq == "1" ~ "1", 
                                  series_type == "three" & series_seq == "2" ~ "50",
                                  series_type == "three" & series_seq == "3" ~ "100",
                                  series_type == "four" & series_seq == "1" ~ "1",
                                  series_type == "four" & series_seq == "2" ~ "50",
                                  series_type == "four" & series_seq == "3" ~ "50",
                                  series_type == "four" & series_seq == "4" ~ "100",
                                  series_type == "five" & series_seq == "1" ~ "1",
                                  series_type == "five" & series_seq == "2" ~ "50",
                                  series_type == "five" & series_seq == "3" ~ "50",
                                  series_type == "five" & series_seq == "4" ~ "50",
                                  series_type == "five" & series_seq == "5" ~ "100",
                                  series_type == "six" & series_seq == "1" ~ "1",
                                  series_type == "six" & series_seq == "2" ~ "50",
                                  series_type == "six" & series_seq == "3" ~ "50",
                                  series_type == "six" & series_seq == "4" ~ "50",
                                  series_type == "six" & series_seq == "5" ~ "50",
                                  series_type == "six" & series_seq == "6" ~ "100"),
         bout_mid = cumsum(lag(ICI, default = 0) > 5) + 1) %>%
  ungroup()

all_calls_bin$frogID<-as.factor(all_calls_bin$frogID)

call_ir <- all_calls_bin %>%
  group_by(frogID) %>%
  arrange(begin.time.s, .by_group = TRUE) %>% 
  mutate(int_ratio_ici = ICI / (ICI + lead(ICI)),
         int_ratio_ioi = call_period / (call_period + lead(call_period))) %>%
  filter(ICI<10) %>%
  ungroup()
#
#
#### call feature summary stats ####
#
# summary stats within bouts
avg_all_frog_bout <- all_calls_bin %>%
  filter(ICI<10) %>%
  group_by(frogID) %>%
  arrange(begin.time.s) %>%
  summarise(n_calls = n(),
            mean_dur = mean(call_type_dur, na.rm = TRUE),
            sd_dur = sd(call_type_dur, na.rm = TRUE),
            CVw_dur = sd_dur / mean_dur,
            mean_ici = mean(ICI, na.rm = TRUE),
            sd_ici = sd(ICI, na.rm = TRUE),
            CVw_ici = sd_ici/mean_ici,
            mean_ioi = mean(call_period, na.rm = TRUE),
            sd_ioi = sd(call_period, na.rm = TRUE),
            CVw_ioi = sd_ioi/mean_ioi) %>%
  mutate(frogID = as.factor(frogID))
# summary for between males
avg_all_bout <- all_calls_bin %>%
  filter(ICI<10) %>%
  summarise(n_calls = n(),
            mean_dur = mean(call_type_dur, na.rm = TRUE),
            sd_dur = sd(call_type_dur, na.rm = TRUE),
            CVb_dur = sd_dur / mean_dur,
            mean_ici = mean(ICI, na.rm = TRUE),
            sd_ici = sd(ICI, na.rm = TRUE),
            CVb_ici = sd_ici/mean_ici,
            mean_ioi = mean(call_period, na.rm = TRUE),
            sd_ioi = sd(call_period, na.rm = TRUE),
            CVb_ioi = sd_ioi/mean_ioi)
#
#### call interval ratio observed and null ####
#
## observed
# IR histogram
call_ir %>%
  filter(ICI<10) %>%
  ggplot(aes(x = int_ratio_ioi)) +
  geom_histogram(position = "identity", alpha = 1, bins = 100, color = "black", fill="#cc99cc", linewidth=0.3) +
  xlab("observed ioi ratio") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .5),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x =  element_text(size = 16),
        axis.text.y = element_text(size = 16))
ggsave("./analysis/graphs/IRs/int_ratio_ioi_all_seq_cat_1_50.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./analysis/graphs/IRs/int_ratio_ioi_all_seq_cat_1_50.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
#
#
## null ioi ratio 
#
ioi_pool <- split(call_ir$call_period, call_ir$frogID)
ioi_pool <- lapply(ioi_pool, function(v) v[is.finite(v) & !is.na(v) & v > 0])
# Calculate how many observed values each frog has
n_emp <- split(call_ir$int_ratio_ioi, call_ir$frogID)
n_emp_clean <- lapply(n_emp, function(x) x[is.finite(x) & !is.na(x)])
# For plotting: generate null matching observed size per frog
n_frogs <- length(ioi_pool)
obs_per_frog <- sapply(n_emp_clean, length)
total_obs <- sum(obs_per_frog)
set.seed(1)
# Generate null ratios matching the exact count per frog for PLOTTING
null_ioi_plot <- lapply(names(ioi_pool), function(fid) {
  n_samples <- obs_per_frog[fid]
  replicate(n_samples, { 
    x <- sample(ioi_pool[[fid]], 2, replace = TRUE)
    x 
  }, simplify = FALSE)  # Keep as list of vectors
})
names(null_ioi_plot) <- names(ioi_pool)

# Flatten properly
nulldist_IOI_plot <- tibble(
  frogID = rep(names(null_ioi_plot), times = sapply(null_ioi_plot, function(x) length(unlist(x)))),
  call_period = as.numeric(unlist(null_ioi_plot))) %>%
  filter(is.finite(call_period), !is.na(call_period)) %>%
  dplyr::group_by(frogID) %>%
  mutate(null_int_ratio = call_period / (call_period + lead(call_period))) %>%
  filter(row_number() %% 2 == 1) %>%
  dplyr::ungroup() %>%
  filter(is.finite(null_int_ratio), !is.na(null_int_ratio))
# Generate larger null pool for resampling tests (100,000 samples like the seal paper)
set.seed(1)
null_ioi_large <- lapply(ioi_pool, function(v) replicate(100000, { x <- sample(v, 2, replace = TRUE); x }))

nulldist_IOI_large <- tibble(
  frogID = rep(names(null_ioi_large), each = 2 * 100000), 
  call_period = as.numeric(unlist(null_ioi_large))) %>%
  filter(is.finite(call_period), !is.na(call_period)) %>%
  dplyr::group_by(frogID) %>%
  mutate(null_int_ratio = call_period / (call_period + lead(call_period))) %>%
  filter(row_number() %% 2 == 1) %>%
  dplyr::ungroup() %>%
  filter(is.finite(null_int_ratio), !is.na(null_int_ratio))

# 1,000 resampled KS tests (matching seal paper method)
set.seed(123)
observed_ratios_ioi <- call_ir %>%
  filter(is.finite(int_ratio_ioi), !is.na(int_ratio_ioi)) %>%
  pull(int_ratio_ioi)

n_obs <- length(observed_ratios_ioi)

ks_resampled_results <- replicate(1000, {
  resampled_null <- sample(nulldist_IOI_large$null_int_ratio, size = n_obs, replace = FALSE)
  ks_result <- suppressWarnings(ks.test(resampled_null, observed_ratios_ioi))
  c(D = ks_result$statistic, p = ks_result$p.value)
})
D_values <- ks_resampled_results[1,]
p_values <- ks_resampled_results[2,]
# Final statistics 
final_stats_ioi_resampled <- tibble(
  D_min = min(D_values),
  D_max = max(D_values),
  p_mean = mean(p_values),
  p_median = median(p_values))
print(final_stats_ioi_resampled)
#
# Plot null distribution (using matched-size null for equal y-axis)
nulldist_IOI_plot %>%
  ggplot(aes(x = null_int_ratio)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 100, color = "black", linewidth=0.3, fill="gray") +
  xlab("null ioi interval ratio") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .5),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x =  element_text(size = 16),
        axis.text.y = element_text(size = 16))
ggsave("./analysis/graphs/IRs/int_ioi_null.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./analysis/graphs/IRs/int_ioi_null.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
##
#
#### call interval ratio within male ####
#
# ICI ratios
output_dir <- "./analysis/graphs/histograms/ratios/"
# Fixed histogram breaks for comparability
breaks_fixed <- seq(0, 10, length.out = 101)  # 100 bins
# Vertical line values per frogID
vline_values <- data.frame(
  frogID = c(1281, 1303, 1314, 1325, 1346, 1552, 1573, 1594, 1605, 1626, 
             1682, 1741, 1752, 1763, 1774, 1785, 1806, 1817),
  vline_x = c(0.41, 0.57, 0.57, 0.32, 0.44, 0.6, 0.64, 0.7, 0.38, 0.76,
              0.46, 0.36, 0.41, 0.4, 0.65, 0.41, 0.41, 0.36))
# Split data by frogID
all_groups <- call_ir %>%
  filter(!is.na(ICI), ICI < 10) %>%
  group_by(frogID) %>%
  group_split()
#
# IOI
# Compute maximum count across all males for consistent y-axis
max_count <- max(sapply(all_groups, function(df) {
  hist(df$int_ratio_ioi, breaks = breaks_fixed, plot = FALSE)$counts
}))

# Loop over each male
for (df in all_groups) {
  frog_id <- df$frogID[1]
  
  # Lookup the vertical line x-value
  vline_x <- vline_values$vline_x[vline_values$frogID == frog_id]
  
  p_hist <- ggplot(df, aes(x = int_ratio_ioi)) +
    geom_histogram(position = "identity", alpha = 0.6, bins = 100) +
    xlab("observed ioi ratio") +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black", size = .5),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.text.x =  element_text(size = 16),
          axis.text.y = element_text(size = 16))
  
  pdf_file <- paste0(output_dir, "hist_ratio_ioi_frog_", frog_id, ".pdf")
  ggsave(pdf_file, plot = p_hist, width = 7, height = 3.55, dpi = 300, limitsize = FALSE)
}
#
#### CV calcs: within and between males in bouts ####
#
# CVw in bouts
calls_cvw <- all_calls_bin  %>%
  group_by(frogID) %>%
  filter(ICI<10) %>%
  summarise(n = n(),
            callsdur.mean = mean(call_type_dur, na.rm = TRUE),
            callsdur.sd = sd(call_type_dur, na.rm = TRUE),
            callsdur.sem = sem(call_type_dur),
            callsdur.CVw = (callsdur.sd/callsdur.mean)*100,
            callsici.mean = mean(ICI, na.rm = TRUE),
            callsici.sd = sd(ICI, na.rm = TRUE),
            callsici.sem = sem(ICI),
            callsici.CVw = (callsici.sd/callsici.mean)*100,
            callsioi.mean = mean(call_period, na.rm = TRUE),
            callsioi.sd = sd(call_period, na.rm = TRUE),
            callsioi.sem = sem(call_period),
            callsioi.CVw = (callsioi.sd/callsioi.mean)*100)
calls_cvw_long <- calls_cvw %>%
  pivot_longer(cols = -c(frogID, n),
               names_to = c("feature", ".value"),
               names_pattern = "(.+)\\.(mean|sd|sem|CVw)$") %>%
  mutate(dynamics = case_when(CVw >= 12 ~ 'dynamic',
                              CVw >= 5 ~ 'intermediate',
                              CVw >= 0 ~ 'static'))

calls_cvb <- calls_cvw_long %>%
  group_by(feature) %>%
  summarize(n = n(),
            grand_mean = mean(mean, na.rm = TRUE),
            grand_sd = sd(mean, na.rm = TRUE),
            grand_sem = sem(mean),
            grand_CVw = mean(CVw, na.rm = TRUE),
            CVb = (grand_sd / grand_mean) * 100,
            .groups = "drop")
calls_cvb_label <- calls_cvb %>%
  mutate(dynamics_w = case_when(grand_CVw >= 12 ~ 'dynamic',
                                grand_CVw >= 5 ~ 'intermediate',
                                grand_CVw >= 0 ~ 'static'),
         dynamics_b = case_when(CVb >= 12 ~ 'dynamic',
                                CVb >= 5 ~ 'intermediate',
                                CVb >= 0 ~ 'static'),
         CV_ratio = CVb / grand_CVw)
#
#### CV tables: between and within males in bouts ####
#
# within
cvw_table_data_calls <- calls_cvw_long %>%
  filter(!is.na(CVw)) %>%
  arrange(frogID, CVw) %>%
  select(frogID, feature, n, mean, sd, sem, CVw)
cvw_ft_calls <- cvw_table_data_calls %>%
  flextable() %>%
  set_header_labels(frogID = "frog",
                    feature = "feature",
                    n = "N",
                    mean = "mean",
                    sd = "sd",
                    sem = "sem",
                    CVw = "CVw (%)") %>%
  colformat_double(j = c("mean", "sd", "CVw"),
                   digits = 2) %>%
  bg(i = ~ CVw >= 0 & CVw < 5,
     j = "CVw",
     bg = CVColors[1]) %>%
  bg(i = ~ CVw >= 5 & CVw < 12,
     j = "CVw",
     bg = CVColors[2]) %>%
  bg(i = ~ CVw >= 12,
     j = "CVw",
     bg = CVColors[3]) %>%
  bold(part = "header") %>%
  add_header_lines("Variation in Call Features in Bouts by Male") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvw_ft_calls
save_as_docx(cvw_ft_calls, path = "./analysis/graphs/tables/CVw_table_calls.docx")
#
#
# between
cvb_table_data_calls <- calls_cvb_label %>%
  select(feature, n, grand_mean, grand_sd, grand_CVw, CVb, CV_ratio)
cvb_ft_calls <- cvb_table_data_calls %>%
  flextable() %>%
  set_header_labels(feature = "feature",
                    n = "N",
                    grand_mean = "mean",
                    grand_sd = "sd",
                    grand_CVw = "mean CVw",
                    CVb = "CVb (%)",
                    CV_ratio = "CVb/CVw Ratio") %>%
  colformat_double(j = c("grand_mean", "grand_sd", "grand_CVw", "CVb", "CV_ratio"),
                   digits = 2) %>%
  bg(i = ~ CVb >= 0 & CVb < 5,
     j = "CVb",
     bg = CVColors[1]) %>%
  bg(i = ~ CVb >= 5 & CVb < 12,
     j = "CVb",
     bg = CVColors[2]) %>%
  bg(i = ~ CVb >= 12,
     j = "CVb",
     bg = CVColors[3]) %>%
  bold(part = "header") %>%
  add_header_lines("Variation in Call Features in Bouts") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvb_ft_calls
save_as_docx(cvb_ft_calls, path = "./analysis/graphs/tables/CVb_table_calls.docx")
#
#### CV calcs: within and between males in series ####
#
# within
calls_cvw_series <- all_calls_bin %>%
  group_by(frogID) %>%
  filter(ICI < 10, 
         series_seq_cat %in% c("1","50"),
         series_type %in% c("three", "four", "five")) %>%
  summarise(n = n(),
            callsdur.mean = mean(call_type_dur, na.rm = TRUE),
            callsdur.sd = sd(call_type_dur, na.rm = TRUE),
            callsdur.sem = sem(call_type_dur),
            callsdur.CVw = (callsdur.sd/callsdur.mean)*100,
            callsici.mean = mean(ICI, na.rm = TRUE),
            callsici.sd = sd(ICI, na.rm = TRUE),
            callsici.sem = sem(ICI),
            callsici.CVw = (callsici.sd/callsici.mean)*100,
            callsioi.mean = mean(call_period, na.rm = TRUE),
            callsioi.sd = sd(call_period, na.rm = TRUE),
            callsioi.sem = sem(call_period),
            callsioi.CVw = (callsioi.sd/callsioi.mean)*100,
            .groups = "drop")
calls_cvw_long_series <- calls_cvw_series %>%
  pivot_longer(cols = -c(frogID, n),
               names_to = c("feature", ".value"),
               names_pattern = "(.+)\\.(mean|sd|sem|CVw)$") %>%
  mutate(dynamics = case_when(CVw >= 12 ~ 'dynamic',
                              CVw >= 5 ~ 'intermediate',
                              CVw >= 0 ~ 'static'))
# between
calls_cvb_series <- calls_cvw_long_series %>%
  group_by(feature) %>%
  summarize(n = n(),
            grand_mean = mean(mean, na.rm = TRUE),
            grand_sd = sd(mean, na.rm = TRUE),
            grand_sem = sem(mean),
            grand_CVw = mean(CVw, na.rm = TRUE),
            CVb = (grand_sd / grand_mean) * 100,
            .groups = "drop")
calls_cvb_series_label <- calls_cvb_series %>%
  mutate(dynamics_w = case_when(grand_CVw >= 12 ~ 'dynamic',
                                grand_CVw >= 5 ~ 'intermediate',
                                grand_CVw >= 0 ~ 'static'),
         dynamics_b = case_when(CVb >= 12 ~ 'dynamic',
                                CVb >= 5 ~ 'intermediate',
                                CVb >= 0 ~ 'static'),
         CV_ratio = CVb / grand_CVw)
#
#### CV tables: between and within males in series ####
#
cvw_table_data_calls_series <- calls_cvw_long_series %>%
  filter(!is.na(CVw)) %>%
  arrange(frogID, CVw) %>%
  select(frogID, feature, n, mean, sd, sem, CVw)
cvw_ft_calls_series <- cvw_table_data_calls_series %>%
  flextable() %>%
  set_header_labels(frogID = "frog",
                    feature = "feature",
                    n = "N",
                    mean = "mean",
                    sd = "sd",
                    sem = "sem",
                    CVw = "CVw (%)") %>%
  colformat_double(j = c("mean", "sd", "CVw"),
                   digits = 2) %>%
  bg(i = ~ CVw >= 0 & CVw < 5,
     j = "CVw",
     bg = CVColors[1]) %>%
  bg(i = ~ CVw >= 5 & CVw < 12,
     j = "CVw",
     bg = CVColors[2]) %>%
  bg(i = ~ CVw >= 12,
     j = "CVw",
     bg = CVColors[3]) %>%
  bold(part = "header") %>%
  add_header_lines("Variation in Call Features by Male Within Series") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvw_ft_calls_series
save_as_docx(cvw_ft_calls_series, path = "./analysis/graphs/tables/CVw_table_calls_in_series.docx")

# between
cvb_table_data_calls_series <- calls_cvb_series_label %>%
  select(feature, n, grand_mean, grand_sd, grand_CVw, CVb, CV_ratio)
cvb_ft_calls_series <- cvb_table_data_calls_series %>%
  flextable() %>%
  set_header_labels(feature = "feature",
                    n = "N",
                    grand_mean = "mean",
                    grand_sd = "sd",
                    grand_CVw = "mean CVw",
                    CVb = "CVb (%)",
                    CV_ratio = "CVb/CVw Ratio") %>%
  colformat_double(j = c("grand_mean", "grand_sd", "grand_CVw", "CVb", "CV_ratio"),
                   digits = 2) %>%
  bg(i = ~ CVb >= 0 & CVb < 5,
     j = "CVb",
     bg = CVColors[1]) %>%
  bg(i = ~ CVb >= 5 & CVb < 12,
     j = "CVb",
     bg = CVColors[2]) %>%
  bg(i = ~ CVb >= 12,
     j = "CVb",
     bg = CVColors[3]) %>%
  bold(part = "header") %>%
  add_header_lines("Variation in Call Timing Features") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvb_ft_calls_series
save_as_docx(cvb_ft_calls_series, path = "./analysis/graphs/tables/CVb_table_calls_in_series.docx")
#
#### CV calcs: within and between males in series by seq order ####
#
calls_cvw_seq <- all_calls_bin %>%
  group_by(frogID, series_seq_cat) %>% #change group here to plot by seq_order or series duration
  filter(ICI < 10, 
         series_seq_cat %in% c("1","50"),
         series_type %in% c("three", "four", "five")) %>%
  summarise(n = n(),
            callsdur.mean = mean(call_type_dur, na.rm = TRUE),
            callsdur.sd = sd(call_type_dur, na.rm = TRUE),
            callsdur.sem = sem(call_type_dur),
            callsdur.CVw = (callsdur.sd/callsdur.mean)*100,
            callsici.mean = mean(ICI, na.rm = TRUE),
            callsici.sd = sd(ICI, na.rm = TRUE),
            callsici.sem = sem(ICI),
            callsici.CVw = (callsici.sd/callsici.mean)*100,
            callsioi.mean = mean(call_period, na.rm = TRUE),
            callsioi.sd = sd(call_period, na.rm = TRUE),
            callsioi.sem = sem(call_period),
            callsioi.CVw = (callsioi.sd/callsioi.mean)*100,
            .groups = "drop")
calls_cvw_seq_long <- calls_cvw_seq %>%
  pivot_longer(cols = -c(frogID, series_seq_cat, n),
               names_to = c("feature", ".value"),
               names_pattern = "(.+)\\.(mean|sd|sem|CVw)$") %>%
  mutate(dynamics = case_when(CVw >= 12 ~ 'dynamic',
                              CVw >= 5 ~ 'intermediate',
                              CVw >= 0 ~ 'static'))
# between
calls_cvb_seq <- calls_cvw_seq_long %>%
  group_by(series_seq_cat, feature) %>%
  summarize(n = n(),
            grand_mean = mean(mean, na.rm = TRUE),
            grand_sd = sd(mean, na.rm = TRUE),
            grand_sem = sem(mean),
            grand_CVw = mean(CVw, na.rm = TRUE),
            CVb = (grand_sd / grand_mean) * 100,
            .groups = "drop")
calls_cvb_seq_label <- calls_cvb_seq %>%
  mutate(dynamics_w = case_when(grand_CVw >= 12 ~ 'dynamic',
                                grand_CVw >= 5 ~ 'intermediate',
                                grand_CVw >= 0 ~ 'static'),
         dynamics_b = case_when(CVb >= 12 ~ 'dynamic',
                                CVb >= 5 ~ 'intermediate',
                                CVb >= 0 ~ 'static'),
         CV_ratio = CVb / grand_CVw)
#
#### CV tables: between and within males in series by seq order ####
#
cvw_table_data_calls_seq <- calls_cvw_seq_long %>%
  filter(!is.na(CVw)) %>%
  group_by(series_seq_cat) %>%
  arrange(frogID, CVw) %>%
  select(frogID, feature, n, mean, sd, sem, CVw)
cvw_ft_calls_seq <- cvw_table_data_calls_seq %>%
  flextable() %>%
  set_header_labels(series_seq_cat = "seq order",
                    frogID = "frog",
                    feature = "feature",
                    n = "N",
                    mean = "mean",
                    sd = "sd",
                    sem = "sem",
                    CVw = "CVw (%)") %>%
  colformat_double(j = c("mean", "sd", "CVw"),
                   digits = 2) %>%
  bg(i = ~ CVw >= 0 & CVw < 5,
     j = "CVw",
     bg = CVColors[1]) %>%
  bg(i = ~ CVw >= 5 & CVw < 12,
     j = "CVw",
     bg = CVColors[2]) %>%
  bg(i = ~ CVw >= 12,
     j = "CVw",
     bg = CVColors[3]) %>%
  bold(part = "header") %>%
  add_header_lines("Variation in Call Features by Male by Seq Within Series") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvw_ft_calls_seq
save_as_docx(cvw_ft_calls_seq, path = "./analysis/graphs/tables/CVw_table_calls_in_series_seq.docx")

# between
cvb_table_data_calls_seq <- calls_cvb_seq_label %>%
  group_by(series_seq_cat) %>%
  select(feature, n, grand_mean, grand_sd, grand_CVw, CVb, CV_ratio)
cvb_ft_calls_seq <- cvb_table_data_calls_seq %>%
  flextable() %>%
  set_header_labels(series_seq_cat = "seq order",
                    feature = "feature",
                    n = "N",
                    grand_mean = "mean",
                    grand_sd = "sd",
                    grand_CVw = "mean CVw",
                    CVb = "CVb (%)",
                    CV_ratio = "CVb/CVw Ratio") %>%
  colformat_double(j = c("grand_mean", "grand_sd", "grand_CVw", "CVb", "CV_ratio"),
                   digits = 2) %>%
  bg(i = ~ CVb >= 0 & CVb < 5,
     j = "CVb",
     bg = CVColors[1]) %>%
  bg(i = ~ CVb >= 5 & CVb < 12,
     j = "CVb",
     bg = CVColors[2]) %>%
  bg(i = ~ CVb >= 12,
     j = "CVb",
     bg = CVColors[3]) %>%
  bold(part = "header") %>%
  add_header_lines("Variation in Call Features by Seq Within Series") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvb_ft_calls_seq
save_as_docx(cvb_ft_calls_seq, path = "./analysis/graphs/tables/CVb_table_calls_in_series_seq.docx")
#