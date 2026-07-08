################## B. pulchellus Timescales of Call Variability #################
#### author: M. Rodriguez-Santiago
####
####
#### this script contains all calculations for series-level analyses
#### need only be called once and stored in environment for entire project
####
#### last update: 3.21.2026
##
##
#
#### packages ####
#
library(tidyverse)
library(cowplot)
library(ggforce)
library(moments)
library(ggpmisc)
library(lmerTest)
library(car)
library(emmeans)
library(flextable)
library(ggpubr)
library(ggeffects)
#
#### data import and tidy ####
#
all_series <- read.csv("all_df_series_final.csv")
colnames(all_series)
all_series_bin <- all_series %>%
  arrange(frogID, begin.time.s) %>%
  group_by(frogID) %>%
  mutate(series_bin = case_when(series_type == "solo" ~ "1", 
                                series_type == "two" ~ "2", 
                                series_type == "three" ~ "3",
                                series_type == "four" ~ "4+",
                                series_type == "five" ~ "4+"),
         series_seq_bin = case_when(series_seq =="0" ~ "0", #grouping calls
                                    series_seq =="1" ~ "1",
                                    series_seq == "2" ~ "2",
                                    series_seq == "3" ~ "3",
                                    series_seq == "4" ~ "4",
                                    series_seq == "5" ~ "4"),
         series_seq_cat = case_when(series_type == "solo" ~ "0", #solo call
                                    series_type == "two" & series_seq == "1" ~ "1", #first call
                                    series_type == "two" & series_seq == "2" ~ "100", #last call
                                    series_type == "three" & series_seq == "1" ~ "1", #first call
                                    series_type == "three" & series_seq == "2" ~ "50", #second 
                                    series_type == "three" & series_seq == "3" ~ "100", #last
                                    series_type == "four" & series_seq == "1" ~ "1",
                                    series_type == "four" & series_seq == "2" ~ "50",
                                    series_type == "four" & series_seq == "3" ~ "50",
                                    series_type == "four" & series_seq == "4" ~ "100",
                                    series_type == "five" & series_seq == "1" ~ "1",
                                    series_type == "five" & series_seq == "2" ~ "50",
                                    series_type == "five" & series_seq == "3" ~ "50",
                                    series_type == "five" & series_seq == "4" ~ "50",
                                    series_type == "five" & series_seq == "5" ~ "100")) %>%
  ungroup()
all_series_bin_10 <- all_series_bin %>%
  filter(ISI<10) %>%
  mutate(series_type = factor(series_bin, 
                              levels = c("1", "2", "3", "4+"),
                              labels = c("solo", "two", "three", "4+")))
#
# eco variables
eco_var <- read.csv("eco_variables_2020.csv")
colnames(eco_var)
all_series_eco <- left_join(all_series_bin_10, eco_var, by = "frogID")
colnames(all_series_eco)
#
#### series feature summary stats within bouts ####
#
avg_all_frog_bout_ser <- all_series_bin %>%
  filter(ISI<10) %>%
  group_by(frogID) %>%
  arrange(begin.time.s) %>%
  summarise(n_series = n(),
            mean_dur = mean(series_dur, na.rm = TRUE),
            sd_dur = sd(series_dur, na.rm = TRUE),
            CVw_dur = sd_dur / mean_dur,
            mean_ici = mean(ISI, na.rm = TRUE),
            sd_ici = sd(ISI, na.rm = TRUE),
            CVw_ici = sd_ici/mean_ici,
            mean_ioi = mean(ISO, na.rm = TRUE),
            sd_ioi = sd(ISO, na.rm = TRUE),
            CVw_ioi = sd_ioi/mean_ioi) %>%
  mutate(frogID = as.factor(frogID))
# summary for between males
avg_all_bout_ser <- all_series_bin %>%
  filter(ISI<10) %>%
  summarise(n_series = n(),
            mean_dur = mean(series_dur, na.rm = TRUE),
            sd_dur = sd(series_dur, na.rm = TRUE),
            CVb_dur = sd_dur / mean_dur,
            mean_ici = mean(ISI, na.rm = TRUE),
            sd_ici = sd(ISI, na.rm = TRUE),
            CVb_ici = sd_ici/mean_ici,
            mean_ioi = mean(ISO, na.rm = TRUE),
            sd_ioi = sd(ISO, na.rm = TRUE),
            CVb_ioi = sd_ioi/mean_ioi)
#
#### series interval ratio - observed and null ####
#
## observed
series_ir <- all_series_eco %>%
  group_by(frogID) %>%
  arrange(begin.time.s, .by_group = TRUE) %>% 
  mutate(int_ratio_isi = ISI / (ISI + lead(ISI)),
         int_ratio_iso = ISO / (ISO + lead(ISO))) %>%
  filter(ISI < 10) %>%
  ungroup()
# graph  IR histogram
series_ir %>%
  ggplot(aes(x = int_ratio_iso)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 100, color = "black", linewidth=0.3, fill="#66cccc") +
  xlab("observed iso ratio") +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.25)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .3),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x =  element_text(size = 16),
        axis.text.y = element_text(size = 16))
ggsave("./v15_graphs/final/series/int_ratio_iso_all.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./v15_graphs/final/series/int_ratio_iso_all.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
#
#
## null isi ratio 
iso_pool <- split(series_ir$ISO, series_ir$frogID)
iso_pool <- lapply(iso_pool, function(v) v[is.finite(v) & !is.na(v) & v > 0])
# Calculate how many observed values each frog has
n_emp <- split(series_ir$int_ratio_iso, series_ir$frogID)
n_emp_clean <- lapply(n_emp, function(x) x[is.finite(x) & !is.na(x)])
# For plotting: generate null matching observed size per frog
n_frogs <- length(iso_pool)
obs_per_frog <- sapply(n_emp_clean, length)
total_obs <- sum(obs_per_frog)
set.seed(1)
# Generate null ratios matching the exact count per frog for PLOTTING
null_iso_plot <- lapply(names(iso_pool), function(fid) {
  n_samples <- obs_per_frog[fid]
  replicate(n_samples, { x <- sample(iso_pool[[fid]], 2, replace = TRUE); x })
})
names(null_iso_plot) <- names(iso_pool)

nulldist_ISO_plot <- tibble(
  frogID = rep(names(null_iso_plot), times = sapply(null_iso_plot, function(x) ncol(x)) * 2),
  ISO = as.numeric(unlist(null_iso_plot))) %>%
  filter(is.finite(ISO), !is.na(ISO)) %>%
  dplyr::group_by(frogID) %>%
  mutate(null_int_ratio = ISO / (ISO + lead(ISO))) %>%
  filter(row_number() %% 2 == 1) %>%
  dplyr::ungroup() %>%
  filter(is.finite(null_int_ratio), !is.na(null_int_ratio))

# Generate larger null pool for resampling tests (100,000 samples like the seal paper)
set.seed(1)
null_iso_large <- lapply(iso_pool, function(v) replicate(100000, { x <- sample(v, 2, replace = TRUE); x }))

nulldist_ISO_large <- tibble(
  frogID = rep(names(null_iso_large), each = 2 * 100000), 
  ISO = as.numeric(unlist(null_iso_large))) %>%
  filter(is.finite(ISO), !is.na(ISO)) %>%
  dplyr::group_by(frogID) %>%
  mutate(null_int_ratio = ISO / (ISO + lead(ISO))) %>%
  filter(row_number() %% 2 == 1) %>%
  dplyr::ungroup() %>%
  filter(is.finite(null_int_ratio), !is.na(null_int_ratio))

# 1,000 resampled KS tests (matching seal paper method)
set.seed(123)
observed_ratios_iso <- series_ir %>%
  filter(is.finite(int_ratio_iso), !is.na(int_ratio_iso)) %>%
  pull(int_ratio_iso)
n_obs <- length(observed_ratios_iso)
ks_resampled_results <- replicate(1000, {
  # Resample null distribution (without replacement) to match observed size
  resampled_null <- sample(nulldist_ISO_large$null_int_ratio, size = n_obs, replace = FALSE)
  ks_result <- ks.test(resampled_null, observed_ratios_iso)
  c(D = ks_result$statistic, p = ks_result$p.value)
})

D_values <- ks_resampled_results[1,]
p_values <- ks_resampled_results[2,]

final_stats_iso_resampled <- tibble(
  D_min = min(D_values),
  D_max = max(D_values),
  p_mean = mean(p_values),
  p_median = median(p_values))
print(final_stats_iso_resampled)
#
nulldist_ISO_plot %>%
  ggplot(aes(x = null_int_ratio)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 100, color = "black", linewidth=0.3, fill="gray") +
  xlab("null onset interval ratio") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .3),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x =  element_text(size = 16),
        axis.text.y = element_text(size = 16))
ggsave("./v15_graphs/final/series/int_iso_null.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./v15_graphs/final/series/int_iso_null.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
#
#
#### calculation of series feature CVs between and within males ####
#
# within
series_cvw <- series_ir  %>%
  group_by(frogID) %>%
  filter(ISI<10) %>%
  summarise(n = n(),
            seriesdur.mean = mean(series_dur, na.rm = TRUE),
            seriesdur.sd = sd(series_dur, na.rm = TRUE),
            seriesdur.sem = sem(series_dur),
            seriesdur.CVw = (seriesdur.sd/seriesdur.mean)*100,
            seriesisi.mean = mean(ISI, na.rm = TRUE),
            seriesisi.sd = sd(ISI, na.rm = TRUE),
            seriesisi.sem = sem(ISI),
            seriesisi.CVw = (seriesisi.sd/seriesisi.mean)*100,
            seriesiso.mean = mean(ISO, na.rm = TRUE),
            seriesiso.sd = sd(ISO, na.rm = TRUE),
            seriesiso.sem = sem(ISO),
            seriesiso.CVw = (seriesiso.sd/seriesiso.mean)*100)
series_cvw_long <- series_cvw %>%
  pivot_longer(cols = -c(frogID, n),
               names_to = c("feature", ".value"),
               names_pattern = "(.+)\\.(mean|sd|sem|CVw)$") %>%
  mutate(dynamics = case_when(CVw >= 12 ~ 'dynamic',
                              CVw >= 5 ~ 'intermediate',
                              CVw >= 0 ~ 'static'))
# between
series_cvb <- series_cvw_long %>%
  group_by(feature) %>%
  summarize(n = n(),
            grand_mean = mean(mean, na.rm = TRUE),
            grand_sd = sd(mean, na.rm = TRUE),
            grand_sem = sem(mean),
            grand_CVw = mean(CVw, na.rm = TRUE),
            CVb = (grand_sd / grand_mean) * 100,
            .groups = "drop")
series_cvb_label <- series_cvb %>%
  mutate(dynamics_w = case_when(grand_CVw >= 12 ~ 'dynamic',
                                grand_CVw >= 5 ~ 'intermediate',
                                grand_CVw >= 0 ~ 'static'),
         dynamics_b = case_when(CVb >= 12 ~ 'dynamic',
                                CVb >= 5 ~ 'intermediate',
                                CVb >= 0 ~ 'static'),
         CV_ratio = CVb / grand_CVw)
#
#### CV tables: series ####
#
cvw_table_data_series <- series_cvw_long %>%
  filter(!is.na(CVw)) %>%
  arrange(frogID, CVw) %>%
  select(frogID, feature, n, mean, sd, sem, CVw)
cvw_ft_series <- cvw_table_data_series %>%
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
  add_header_lines("Variation in Series Features in Calling Bouts By Male") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvw_ft_series
save_as_docx(cvw_ft_series, path = "./v15_graphs/final/tables/CVw_table_series.docx")

cvb_table_data_series <- series_cvb_label %>%
  select(feature, n, grand_mean, grand_sd, grand_CVw, CVb, CV_ratio)
cvb_ft_series <- cvb_table_data_series %>%
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
  add_header_lines("Variation in Series Features in Calling Bouts") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvb_ft_series
save_as_docx(cvb_ft_series, path = "./v15_graphs/final/tables/CVb_table_series.docx")
#
