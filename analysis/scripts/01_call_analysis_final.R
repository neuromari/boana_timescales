################# B. pulchellus Timescales of Call Variability #################
#### author: M. Rodriguez-Santiago
####
####
#### this script contains all calculations for call-level analyses
#### need only be called once and stored in environment for entire project
####
#### last update: 3.21.2026
##
##
#
#### all packages needed for this script ####
##
library(tidyverse)
library(tidyr)
library(cowplot)
library(ggforce)
library(moments)
library(ggpmisc)
library(lmerTest)
library(car)
library(emmeans)
library(performance)
library(flextable)
library(ggpubr)
library(reshape2)
library(gridExtra)
library(grid)
library(ggpmisc)
library(officer)
#
### ALL CODE FOR CALL TYPE-LEVEL ANALYSES IN SEQUENCE OF PRESENTATION IN PAPER ####
#
#### data and variable explanations ####
#
all_calls <- read.csv("all_df_calls_final.csv")
colnames(all_calls) 
# call_period is difference in onsets between consec calls (regardless call type)
# call_type: calls in calls (solo if only 1, two for 2 calls in calls, etc)
# call_seq: sequence of call within the calls (0 = solo, 1 = first in calls, 2=second in calls, etc)
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
ggsave("./v15_graphs/final/calls/int_ratio_ioi_all_seq_cat_1_50.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./v15_graphs/final/calls/int_ratio_ioi_all_seq_cat_1_50.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
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
ggsave("./v15_graphs/final/calls/int_ioi_null.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./v15_graphs/final/calls/int_ioi_null.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
##
#
#### call interval ratio within series observed and null ####
#
call_ir %>%
  filter(ICI < 10,
         series_seq_cat %in% c("1", "50"),
         series_type %in% c("three", "four", "five")) %>%
  count(series_type, series_seq_cat) %>%
  pivot_wider(names_from  = series_seq_cat,
              values_from = n,
              names_prefix = "bin_",
              values_fill  = 0) %>%
  arrange(series_type)
## observed
# IR histogram
call_ir %>%
  filter(ICI < 10, 
         series_seq_cat %in% c("1","50"),
         series_type %in% c("three", "four", "five")) %>%
  ggplot(aes(x = int_ratio_ioi, fill = series_seq_cat)) +
  geom_histogram(position = "identity", alpha = 1, bins = 100, 
                 color = "black",linewidth=0.3) +
  xlab("observed ioi ratio") +
  scale_fill_manual(name = "  within\n  series\n position",
                     values = c("1" = "#e4cee6", "50" = "#cc99cc"),
                     labels = c("first", "mid")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .5),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x =  element_text(size = 16),
        axis.text.y = element_text(size = 16))
ggsave("./v15_graphs/final/calls/int_ratio_ico_within_seq.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./v15_graphs/final/calls/int_ratio_ico_within_seq.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
#
##
call_ir %>%
  filter(ICI < 10,
         (series_type == "three" & series_seq == 1) |
           (series_type == "four"  & series_seq %in% c(1, 2)) |
           (series_type == "five"  & series_seq %in% c(1, 2, 3))) %>%
  mutate(series_seq = factor(series_seq)) %>%
  ggplot(aes(x = int_ratio_ioi, fill = series_seq)) +
  geom_histogram(position = "identity", alpha = 1, bins = 100,
                 color = "black", linewidth = 0.3) +
  xlab("observed ioi ratio") +
  scale_fill_manual(name= "  within\n  series\n position",
    values = c("1" = "#e4cee6", "2" = "#cc99cc", "3" = "#9966aa"),
    labels = c("1" = "first", "2" = "second", "3" = "third")) +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line        = element_line(colour = "black", size = .5),
    axis.title.y     = element_text(size = 20),
    axis.title.x     = element_text(size = 20),
    axis.text.x      = element_text(size = 16),
    axis.text.y      = element_text(size = 16))
ggsave("./v15_graphs/final/calls/int_ratio_ico_within_seq_nopenul.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./v15_graphs/final/calls/int_ratio_ico_within_seq_nopenul.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
#
#
#### call interval ratio within male ####
#
# ICI ratios
output_dir <- "./v15_graphs/final/calls/histograms/ratios/"
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

# Compute maximum count across all males for consistent y-axis
max_count <- max(sapply(all_groups, function(df) {
  hist(df$int_ratio_ici, breaks = breaks_fixed, plot = FALSE)$counts
}))

# Loop over each male
for (df in all_groups) {
  frog_id <- df$frogID[1]
  
  # Lookup the vertical line x-value
  vline_x <- vline_values$vline_x[vline_values$frogID == frog_id]
  
  p_hist <- ggplot(df, aes(x = int_ratio_ici)) +
    geom_histogram(position = "identity", alpha = 0.6, bins = 100) +
    xlab("observed ici ratio") +
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
  
  pdf_file <- paste0(output_dir, "hist_ratio_ici_frog_", frog_id, ".pdf")
  ggsave(pdf_file, plot = p_hist, width = 7, height = 3.55, dpi = 300, limitsize = FALSE)
}
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
#
## graph IR histogram colored by calls order PER MALE
#
frog_ids <- unique(call_ir$frogID)
for(frog in frog_ids) {
  
  frog_data <- call_ir %>%
    filter(frogID == frog)
  
  p <- frog_data %>%
    ggplot(aes(x = int_ratio_ioi, fill = series_seq_cat)) +
    geom_histogram(position = "identity", alpha = 0.6, bins = 100) +
    xlab("observed ioi ratio") +
    ggtitle(paste("Frog ID:", frog)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    scale_fill_manual(values = callsSeqCatMid,
                      name = "call order",
                      breaks = c("1", "50", "100", "0"),
                      labels = c("1st", "mid", "last", "solo")) +
    scale_color_manual(values = callsSeqCatMid,
                       name = "call order",
                       breaks = c("1", "50", "100", "0"),
                       labels = c("1st", "mid", "last", "solo")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black", size = .5),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          plot.title = element_text(size = 18, hjust = 0.5))
  
  ggsave(paste0("./v15_graphs/calls/histograms/int_ratio_ioi_call_seq_", frog, ".pdf"), 
         plot = p, width = 7, height = 3.5, dpi = 300, limitsize = FALSE)
}
#
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
save_as_docx(cvw_ft_calls, path = "./v15_graphs/final/tables/CVw_table_calls.docx")
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
save_as_docx(cvb_ft_calls, path = "./v15_graphs/final/tables/CVb_table_calls.docx")
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
save_as_docx(cvw_ft_calls_series, path = "./v15_graphs/final/tables/CVw_table_calls_in_series.docx")

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
save_as_docx(cvb_ft_calls_series, path = "./v15_graphs/final/tables/CVb_table_calls_in_series.docx")
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
save_as_docx(cvw_ft_calls_seq, path = "./v15_graphs/final/tables/CVw_table_calls_in_series_seq.docx")

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
save_as_docx(cvb_ft_calls_seq, path = "./v15_graphs/final/tables/CVb_table_calls_in_series_seq.docx")
#






notes_cvw_series <- notes_ir %>%
  group_by(frogID, note_type_relabel, series_seq_cat) %>%
  filter(INI < 10,
         series_seq_cat %in% c("1","50"),
         series_type %in% c("three", "four", "five"),
         note_type_relabel %in% c("1", "2")) %>%
  summarise(n = n(),
            notesdur.mean = mean(note.dur.s, na.rm = TRUE),
            notesdur.sd = sd(note.dur.s, na.rm = TRUE),
            notesdur.sem = sem(note.dur.s),
            notesdur.CVw = (notesdur.sd / notesdur.mean) * 100,
            notesini.mean = mean(INI, na.rm = TRUE),
            notesini.sd = sd(INI, na.rm = TRUE),
            notesini.sem = sem(INI),
            notesini.CVw = (notesini.sd / notesini.mean) * 100,
            notesino.mean = mean(note.period, na.rm = TRUE),
            notesino.sd = sd(note.period, na.rm = TRUE),
            notesino.sem = sem(note.period),
            notesino.CVw = (notesino.sd / notesino.mean) * 100,
            notesmeandom.mean = mean(meandom, na.rm = TRUE),
            notesmeandom.sd = sd(meandom, na.rm = TRUE),
            notesmeandom.sem = sem(meandom),
            notesmeandom.CVw = (notesmeandom.sd / notesmeandom.mean) * 100,
            .groups = "drop")
notes_cvw_long_series <- notes_cvw_series %>%
  pivot_longer(cols = -c(frogID, note_type_relabel, series_seq_cat, n),
               names_to = c("feature", ".value"),
               names_pattern = "(.+)\\.(mean|sd|sem|CVw)$") %>%
  mutate(dynamics = case_when(CVw >= 12 ~ 'dynamic',
                              CVw >= 5 ~ 'intermediate',
                              CVw >= 0 ~ 'static'))
notes_cvb_series <- notes_cvw_long_series %>%
  group_by(series_seq_cat, feature, note_type_relabel) %>%
  summarize(n = n(),
            grand_mean = mean(mean, na.rm = TRUE),
            grand_sd = sd(mean, na.rm = TRUE),
            grand_sem = sem(mean),
            grand_CVw = mean(CVw, na.rm = TRUE),
            CVb = (grand_sd / grand_mean) * 100,
            .groups = "drop")
notes_cvb_series_label <- notes_cvb_series %>%
  mutate(dynamics_w = case_when(grand_CVw >= 12 ~ 'dynamic',
                                grand_CVw >= 5 ~ 'intermediate',
                                grand_CVw >= 0 ~ 'static'),
         dynamics_b = case_when(CVb >= 12 ~ 'dynamic',
                                CVb >= 5 ~ 'intermediate',
                                CVb >= 0 ~ 'static'),
         CV_ratio = CVb / grand_CVw)
calls_plot_seq <- calls_cvb_series_label %>%
  mutate(feature_label = case_when(
    feature == "callsdur" ~ "call dur",
    feature == "callsici" ~ "ICI",
    feature == "callsioi" ~ "ICO"))
notes_plot_seq <- notes_cvb_series_label %>%
  filter((feature %in% c("notesdur", "notesmeandom")) |
           (feature %in% c("notesini", "notesino") & note_type_relabel == "1")) %>%
  mutate(feature_label = case_when(
    feature == "notesdur" & note_type_relabel == "1" ~ "dur N1",
    feature == "notesdur" & note_type_relabel == "2" ~ "dur N2",
    feature == "notesini" & note_type_relabel == "1" ~ "INI",
    feature == "notesino" & note_type_relabel == "1" ~ "INO",
    feature == "notesmeandom" & note_type_relabel == "1" ~ "DF N1",
    feature == "notesmeandom" & note_type_relabel == "2" ~ "DF N2"))
combined_seq <- bind_rows(calls_plot_seq, notes_plot_seq) %>%
  mutate(series_label = case_when(
    series_seq_cat == "1" ~ "first in series",
    series_seq_cat == "50" ~ "mid series"),
    feature_label = factor(feature_label, levels = c("call dur", "ICI", "ICO", "dur N1", "dur N2", "INI", "INO", "DF N1", "DF N2")))
label_colors <- c("call dur" = "#9f839d",
                  "ICI" = "#cfa9d1",
                  "ICO" = "#e4cee6",
                  "dur N1" = "#839c7e",
                  "dur N2" = "#839c7e",
                  "INI" = "#99cc99",
                  "INO" = "#cde4cc",
                  "DF N1" = "#a8c4a0",
                  "DF N2" = "#a8c4a0")
label_shapes <- c("call dur" = 16,
                  "ICI" = 16,
                  "ICO" = 16,
                  "dur N1" = 17,
                  "dur N2" = 15,
                  "INI" = 17,
                  "INO" = 17,
                  "DF N1" = 17,
                  "DF N2" = 15)
p <- ggplot(combined_seq, aes(x = feature_label, y = CV_ratio, color = feature_label, shape = feature_label)) +
  geom_point(size = 7, alpha = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.7) +
  facet_wrap(~series_label, scales = "free_x") +
  xlab("") +
  ylab("CVb / CVw") +
  scale_color_manual(values = label_colors) +
  scale_shape_manual(values = label_shapes) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = .3),
        legend.position = "none",
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.placement = "outside",
        strip.background = element_blank())

g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-b', g$layout$name))

for(i in seq_along(strips)) {
  strip_grob <- g$grobs[[strips[i]]]
  if(i == 1) {
    strip_grob$grobs[[1]]$children[[1]]$gp$fill <- "#cccccc"
  } else {
    strip_grob$grobs[[1]]$children[[1]]$gp$fill <- "#666666"
  }
  g$grobs[[strips[i]]] <- strip_grob
}

grid::grid.draw(g)
ggsave("./v15_graphs/CV/CVratio_by_series_position.pdf", width=10, height=5, dpi = 300, limitsize = FALSE)
ggsave("./v15_graphs/CV/CVratio_by_series_position.svg", width=10, height=5, dpi = 300, limitsize = FALSE)
#
#