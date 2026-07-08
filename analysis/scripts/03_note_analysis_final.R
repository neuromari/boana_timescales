################## B. pulchellus Timescales of Call Variability #################
#### author: M. Rodriguez-Santiago
#### last update: 3.21.2026
####
#### this script contains all calculations for note-level analyses
#### need only be called once and stored in environment
####
#### variables:
##   note_type_uc = categories of automated note selections where '1=1, 2=2'uc' is notes previously not described in literature (these numbers were originally coded in Raven as part of manual annotation)
##   note_type_relabel = final note categories described in manuscript
##     1, 2 — the two note types described in the literature, which co-occur as a 1–2 sequence
##     11 — a note whose oscillogram visually matches note 1 but is not followed by a note 2
##     12 — a note whose oscillogram matches note 2 but is not preceded by a note 1
##   call_type_relabel: doublet (notes 1 and 2), squeak (note 3), solo (a singular note 1 and 2 not in doublet) 
##   series_bin = # of calls within a series
##   series_seq_cat = within-series category where 
##     0: solo call
##     1: part of first call in series
##     50: part of mid call in series
##     100: part of last call in series
####
#### input: all_df_notes.csv generated in 00_DataImport.R 
#### output: all_df_notes_final.csv with series variables
####
####
#
#
#### packages needed to run this code ####
#
library(dplyr)
library(flextable)
library(ggpubr)
library(lmerTest)
library(car)
library(emmeans)
library(performance)
library(reshape2)
library(gridExtra)
library(grid)
library(ggpmisc)
#
#
#### data import & tidy ####
#
all_notes <- read.csv("./data/all_df_notes.csv")
unique(all_notes$Note_type) #1 2 and 3 11 12 are solo notes
colnames(all_notes)
all_notes <- all_notes %>%
  select(FileName, frogID, begin.time.s, end.time.s, 
         Note_type, note.dur.s, INI, note.period, note_num_series,
         call_type, series_type, series_seq) %>%
  mutate(note_type_uc = case_when(Note_type=="1" ~ "1",
                                  Note_type == "2" ~ "2",
                                  Note_type == "3" ~ "uc", 
                                  Note_type == "11" ~ "uc",
                                  Note_type == "12" ~ "uc",
                                  Note_type == "uc" ~ "uc"),
         note_type_relabel = case_when(Note_type == "1" ~ "1", 
                                       Note_type == "2" ~ "2", 
                                       Note_type == "3" ~ "3",
                                       Note_type == "11" ~ "uc",
                                       Note_type == "12" ~ "uc",
                                       Note_type == "uc" ~ "uc"),
         call_type_relabel = case_when(call_type == "doublet" ~ "doublet",
                                       call_type == "solo" & Note_type=="3" ~ "squeak",
                                       call_type == "solo" ~ "solo"),
         series_bin = case_when(series_type == "solo" ~ "1", 
                                series_type == "two" ~ "2", 
                                series_type == "three" ~ "3",
                                series_type == "four" ~ "4+",
                                series_type == "five" ~ "4+"),
         note_num_series_bin = case_when(note_num_series =="1" ~ "1",
                                         note_num_series == "2" ~ "2",
                                         note_num_series == "3" ~ "3",
                                         note_num_series == "4" ~ "4",
                                         note_num_series == "5" ~ "5",
                                         note_num_series =="6" ~ "6",
                                         note_num_series =="7" ~ "8+",
                                         note_num_series =="8" ~ "8+",
                                         note_num_series =="9" ~ "8+",
                                         note_num_series =="10" ~ "8+"),
         series_seq_cat = case_when(
           note_type_relabel == "uc" ~ NA_character_,
           series_type == "solo" & note_num_series == "1" ~ "0",
           series_type == "solo" & note_num_series == "2" ~ "0",
           series_type == "two" & note_num_series == "1" ~ "1", 
           series_type == "two" & note_num_series == "2" ~ "1", 
           series_type == "two" & note_num_series == "3" ~ "100", 
           series_type == "two" & note_num_series == "4" ~ "100",
           series_type == "three" & note_num_series == "1" ~ "1", 
           series_type == "three" & note_num_series == "2" ~ "1",
           series_type == "three" & note_num_series == "3" ~ "50",
           series_type == "three" & note_num_series == "4" ~ "50",
           series_type == "three" & note_num_series == "5" ~ "100",
           series_type == "three" & note_num_series == "6" ~ "100",
           series_type == "four" & note_num_series == "1" ~ "1", 
           series_type == "four" & note_num_series == "2" ~ "1",
           series_type == "four" & note_num_series == "3" ~ "50",
           series_type == "four" & note_num_series == "4" ~ "50",
           series_type == "four" & note_num_series == "5" ~ "50",
           series_type == "four" & note_num_series == "6" ~ "50",
           series_type == "four" & note_num_series == "7" ~ "100",
           series_type == "four" & note_num_series == "8" ~ "100",
           series_type == "five" & note_num_series == "1" ~ "1", 
           series_type == "five" & note_num_series == "2" ~ "1",
           series_type == "five" & note_num_series == "3" ~ "50",
           series_type == "five" & note_num_series == "4" ~ "50",
           series_type == "five" & note_num_series == "5" ~ "50",
           series_type == "five" & note_num_series == "6" ~ "50",
           series_type == "five" & note_num_series == "7" ~ "50",
           series_type == "five" & note_num_series == "8" ~ "50",
           series_type == "five" & note_num_series == "9" ~ "100",
           series_type == "five" & note_num_series == "10" ~ "100",
           TRUE ~ NA_character_))
eco_var <- read.csv("./data/eco_variables_2020.csv")
colnames(eco_var)
all_notes_eco <- left_join(all_notes, eco_var, by = "frogID")
unique(all_notes$note_type_relabel) #1 2 and 3 11 12 are solo notes

# Kim code output
colnaspecsK <- read.csv("SpectrotemporalOutputs_boanaKH.csv") 
colnames(specsK)
# X1, X2 columns are probability assignment to group 1 and group 2 using random forest
specsK_subset <- specsK %>%
  select(frogID, begin.time.s, 
         meanfreq, meandom, maxdom, enddom, dfrange,meanpeakf, freq.IQR, 
         time.median, time.IQR, 
         skew, entropy, sfm, modindx)
all_notes_final <- all_notes_eco %>%
  left_join(specsK_subset, by = c("frogID", "begin.time.s")) %>%
  filter(maxdom<5)
colnames(all_notes_final)
write.csv(all_notes, "./data/all_notes_final.csv") #run only once

#
notes_ir <- all_notes_final %>%
  group_by(frogID) %>%
  arrange(begin.time.s, .by_group = TRUE) %>% 
  mutate(int_ratio_ino = note.period / (note.period + lead(note.period)),
         int_ratio_ini = INI/(INI+lead(INI))) %>%
  filter(note.period<10) %>%
  ungroup()
#
#### note interval ratios ####
#
## observed
#
#ino
notes_ir %>%
  ggplot(aes(x = int_ratio_ino)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 100) +
  xlab("observed ino ratio") +
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
ggsave("./v15_graphs/int_ratio_ino.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./v15_graphs/int_ratio_ino.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
#
##
# null ino
##

ino_pool <- split(notes_ir$note.period, notes_ir$frogID)
ino_pool <- lapply(ino_pool, function(v) v[is.finite(v) & !is.na(v) & v > 0])
# Calculate how many observed values each frog has
n_emp <- split(notes_ir$int_ratio_ino, notes_ir$frogID)
n_emp_clean <- lapply(n_emp, function(x) x[is.finite(x) & !is.na(x)])
max_obs_per_frog <- max(sapply(n_emp_clean, length))
# Generate enough null values per frog (at least as many as max observed)
n_frogs <- length(ino_pool)
samples_per_frog <- max(ceiling(6000 / n_frogs), max_obs_per_frog + 100)
set.seed(1)
null_ino  <- lapply(ino_pool, function(v) replicate(samples_per_frog, { x <- sample(v, 2, replace = TRUE); x }))
nulldist_ino <- tibble(frogID = rep(names(null_ino), each = 2 * samples_per_frog), 
                       note.period = as.numeric(unlist(null_ino))) %>%
  filter(is.finite(note.period), !is.na(note.period)) %>%
  dplyr::group_by(frogID) %>%
  mutate(null_int_ratio = note.period / (note.period + lead(note.period))) %>%
  filter(row_number() %% 2 == 1) %>%      # keep only first of each sampled pair (within frog)
  dplyr::ungroup() %>%
  filter(is.finite(null_int_ratio), !is.na(null_int_ratio))

#resample 1,000 null distributions (WITHOUT replacement)
ks_results <- do.call(rbind, Map(function(fid, emp) {
  emp <- emp[is.finite(emp) & !is.na(emp)]
  res <- replicate(1000, {
    samp <- sample(nulldist_ino$null_int_ratio[nulldist_ino$frogID == fid], size = length(emp), replace = FALSE)
    out  <- suppressWarnings(stats::ks.test(samp, emp)) # 2-sample KS (empirical vs resampled null)
    c(D = unname(out$statistic), p = out$p.value)
  })
  cbind(as.data.frame(t(res)), frogID = fid)
}, names(n_emp), n_emp))

ks_results_ino <- ks_results

final_stats_ino <- tibble(frogID = unique(ks_results_ino$frogID),
  D_min  = tapply(ks_results_ino$D, ks_results_ino$frogID, min, na.rm = TRUE),
  D_max  = tapply(ks_results_ino$D, ks_results_ino$frogID, max, na.rm = TRUE),
  p_mean = tapply(ks_results_ino$p, ks_results_ino$frogID, mean, na.rm = TRUE))
#
nulldist_ino %>%
  ggplot(aes(x = null_int_ratio)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 100) +
  xlab("null ino ratio") +
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
ggsave("./v15_graphs/int_ratio_null_ino.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./v15_graphs/int_ratio_null_ino.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
#
#### note interval ratios by seq order within series####
#
## observed
#
#ino
notes_ir %>%
  filter(INI < 10, 
         series_seq_cat %in% c("1","50"),
         series_type %in% c("three", "four", "five"),
         note_type_relabel %in% c("1")) %>%
  ggplot(aes(x = int_ratio_ino, fill = series_seq_cat)) +
  geom_histogram(position = "identity", alpha = 1, bins = 100, 
                 color = "black",linewidth=0.3) +
  xlab("observed ino ratio") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_fill_manual(name = "  within\n  series\n position",
                     values = c("1" = "#bcdab7", "50" = "#336633"),
                     labels = c("1st", "mid")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .5),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x =  element_text(size = 16),
        axis.text.y = element_text(size = 16))
ggsave("./v15_graphs/final/notes/int_ratio_ino_within_seq.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./v15_graphs/final/notes/int_ratio_ino_within_seq.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
#
#### note interval ratios within male ####
#
##
breaks_fixed <- seq(0, 1, length.out = 101)  # 100 bins
all_groups <- notes_ir %>%
  filter(is.finite(int_ratio_ino), !is.na(int_ratio_ino)) %>%
  group_by(frogID) %>%
  group_split()
# Compute maximum count across all males for consistent y-axis
max_count <- max(sapply(all_groups, function(df) {
  hist(df$int_ratio_ino, breaks = breaks_fixed, plot = FALSE)$counts
}))

for (df in all_groups) {
  frog_id <- df$frogID[1]
  
  p_hist <- ggplot(df, aes(x = int_ratio_ino)) +
    geom_histogram(breaks = breaks_fixed,
                   fill = "gray", color = "black", alpha = 0.6) +
    ylab("count") +
    scale_x_continuous(limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    scale_y_continuous(limits = c(0, max_count)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", size = 0.5),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
  
  pdf_file <- paste0(output_dir_notes, "hist_ino_ratio_frog_", frog_id, ".pdf")
  ggsave(pdf_file, plot = p_hist, width = 10, height = 5, dpi = 300, limitsize = FALSE)
}
##
#### within male INO histograms ####
#
## this code identifies trophs in the first distribution bound per male within bouts
## plots all males with their distinct first distribution min (troph) as vertical lines
output_dir_notes <- "./v15_graphs/final/notes/histograms/"
# Fixed histogram breaks for comparability
breaks_fixed <- seq(0, 10, length.out = 101)  # 100 bins
# Vertical line values per frogID
vline_values <- data.frame(
  frogID = c(1281, 1303, 1314, 1325, 1346, 1552, 1573, 1594, 1605, 1626, 
             1682, 1741, 1752, 1763, 1774, 1785, 1806, 1817),
  vline_x = c(0.41, 0.57, 0.57, 0.32, 0.44, 0.6, 0.64, 0.7, 0.38, 0.76,
              0.46, 0.36, 0.41, 0.4, 0.65, 0.41, 0.41, 0.36))
all_groups <- all_notes_eco %>%
  filter(!is.na(note.period), note.period < 5) %>%
  group_by(frogID) %>%
  group_split()

# Compute maximum count across all males for consistent y-axis
max_count <- max(sapply(all_groups, function(df) {
  hist(df$note.period, breaks = breaks_fixed, plot = FALSE)$counts
}))

# Loop over each male
for (df in all_groups) {
  frog_id <- df$frogID[1]
  
  # Lookup the vertical line x-value
  vline_x <- vline_values$vline_x[vline_values$frogID == frog_id]
  
  p_hist <- ggplot(df, aes(x = note.period)) +
    geom_histogram(breaks = breaks_fixed,
                   fill = "gray", color = "black", alpha = 0.6) +
    geom_vline(xintercept = vline_x, color = "#A16928", linetype = "dashed", size = 0.3) +
    ylab("count") +
    scale_x_continuous(limits = c(0, 5),
                       breaks = seq(0, 5, by = 1)) +
    scale_y_continuous(limits = c(0, max_count)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", size = 0.3),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 12))
  
  p_box <- ggplot(df, aes(x = note.period, y = "")) +
    geom_boxplot(outlier.alpha = 0, width = 0.5) +
    geom_point(position = position_jitter(height = 0.1), alpha = 0.3) +
    geom_vline(xintercept = vline_x, color = "#A16928", 
               linetype = "dashed", size = 0.5) +
    xlab("ino (s)") +
    scale_x_continuous(limits = c(0, 5),
                       breaks = seq(0, 5, by = 1)) +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black", size = .3),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_blank(),
      legend.position = "none")
  
  combo_plot <- plot_grid(p_hist, p_box, ncol = 1, align = "v", rel_heights = c(2,1))
  
  pdf_file <- paste0(output_dir_notes, "hist_ino_frog_", frog_id, ".pdf")
  ggsave(pdf_file, plot = combo_plot, width = 10, height = 5, dpi = 300, limitsize = FALSE)
}
#
#### CV calcs: within and between males in bouts ####
#
# CVw in bouts
notes_cvw <- all_notes_final %>%
  group_by(frogID, note_type_relabel) %>%
  filter(INI < 10,
         note_type_relabel %in% c("1", "2")) %>%
  summarise(n = n(),
            notesdur.mean    = mean(note.dur.s, na.rm = TRUE),
            notesdur.sd      = sd(note.dur.s, na.rm = TRUE),
            notesdur.sem     = sem(note.dur.s),
            notesdur.CVw     = (notesdur.sd / notesdur.mean) * 100,
            notesini.mean    = mean(INI, na.rm = TRUE),
            notesini.sd      = sd(INI, na.rm = TRUE),
            notesini.sem     = sem(INI),
            notesini.CVw     = (notesini.sd / notesini.mean) * 100,
            notesino.mean    = mean(note.period, na.rm = TRUE),
            notesino.sd      = sd(note.period, na.rm = TRUE),
            notesino.sem     = sem(note.period),
            notesino.CVw     = (notesino.sd / notesino.mean) * 100,
            notesmeandom.mean = mean(meandom, na.rm = TRUE),
            notesmeandom.sd   = sd(meandom, na.rm = TRUE),
            notesmeandom.sem  = sem(meandom),
            notesmeandom.CVw  = (notesmeandom.sd / notesmeandom.mean) * 100,
            .groups = "drop")

notes_cvw_long <- notes_cvw %>%
  pivot_longer(cols = -c(frogID, note_type_relabel, n),
               names_to = c("feature", ".value"),
               names_pattern = "(.+)\\.(mean|sd|sem|CVw)$") %>%
  mutate(dynamics = case_when(CVw >= 12 ~ 'dynamic',
                              CVw >= 5  ~ 'intermediate',
                              CVw >= 0  ~ 'static'))
notes_cvb <- notes_cvw_long %>%
  group_by(feature, note_type_relabel) %>%
  summarize(n          = n(),
            grand_mean = mean(mean, na.rm = TRUE),
            grand_sd   = sd(mean, na.rm = TRUE),
            grand_sem  = sem(mean),
            grand_CVw  = mean(CVw, na.rm = TRUE),
            CVb        = (grand_sd / grand_mean) * 100,
            .groups    = "drop")
notes_cvb_label <- notes_cvb %>%
  mutate(dynamics_w = case_when(grand_CVw >= 12 ~ 'dynamic',
                                grand_CVw >= 5  ~ 'intermediate',
                                grand_CVw >= 0  ~ 'static'),
         dynamics_b = case_when(CVb >= 12 ~ 'dynamic',
                                CVb >= 5  ~ 'intermediate',
                                CVb >= 0  ~ 'static'),
         CV_ratio = CVb / grand_CVw)
#
#### CV tables: between and within males in bouts ####
#
# within
cvw_table_data_notes <- notes_cvw_long %>%
  filter(!is.na(CVw)) %>%
  group_by(frogID, note_type_relabel) %>%
  arrange(frogID, note_type_relabel, CVw) %>%
  select(frogID, feature, n, mean, sd, sem, CVw)
cvw_ft_notes <- cvw_table_data_notes %>%
  flextable() %>%
  set_header_labels(note_type_relabel = "note",
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
  add_header_lines("Variation in Note Features in Bouts by Male") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvw_ft_notes
save_as_docx(cvw_ft_notes, path = "./v15_graphs/final/tables/CVw_table_notes.docx")
#
# between
cvb_table_data_notes <- notes_cvb_label %>%
  group_by(note_type_relabel) %>%
  select(feature, n, grand_mean, grand_sd, grand_CVw, CVb, CV_ratio)
cvb_ft_notes <- cvb_table_data_notes %>%
  flextable() %>%
  set_header_labels(note_type_relabel = "note",
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
  add_header_lines("Variation in Note Features in Bouts") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvb_ft_notes
save_as_docx(cvb_ft_notes, path = "./v15_graphs/final/tables/CVb_table_notes.docx")
#
#### CV calcs: within and between males in series ####
# within
notes_cvw_series <- all_notes_final %>%
  group_by(frogID, note_type_relabel) %>%
  filter(INI < 10,
         series_type %in% c("three", "four", "five"),
         note_type_relabel %in% c("1", "2")) %>%
  summarise(n = n(),
            notesdur.mean    = mean(note.dur.s, na.rm = TRUE),
            notesdur.sd      = sd(note.dur.s, na.rm = TRUE),
            notesdur.sem     = sem(note.dur.s),
            notesdur.CVw     = (notesdur.sd / notesdur.mean) * 100,
            notesini.mean    = mean(INI, na.rm = TRUE),
            notesini.sd      = sd(INI, na.rm = TRUE),
            notesini.sem     = sem(INI),
            notesini.CVw     = (notesini.sd / notesini.mean) * 100,
            notesino.mean    = mean(note.period, na.rm = TRUE),
            notesino.sd      = sd(note.period, na.rm = TRUE),
            notesino.sem     = sem(note.period),
            notesino.CVw     = (notesino.sd / notesino.mean) * 100,
            notesmeandom.mean = mean(meandom, na.rm = TRUE),
            notesmeandom.sd   = sd(meandom, na.rm = TRUE),
            notesmeandom.sem  = sem(meandom),
            notesmeandom.CVw  = (notesmeandom.sd / notesmeandom.mean) * 100,
            .groups = "drop")

notes_cvw_long_series <- notes_cvw_series %>%
  pivot_longer(cols = -c(frogID, note_type_relabel, n),
               names_to = c("feature", ".value"),
               names_pattern = "(.+)\\.(mean|sd|sem|CVw)$") %>%
  mutate(dynamics = case_when(CVw >= 12 ~ 'dynamic',
                              CVw >= 5  ~ 'intermediate',
                              CVw >= 0  ~ 'static'))
# between
notes_cvb_series <- notes_cvw_long_series %>%
  group_by(feature, note_type_relabel) %>%
  summarize(n          = n(),
            grand_mean = mean(mean, na.rm = TRUE),
            grand_sd   = sd(mean, na.rm = TRUE),
            grand_sem  = sem(mean),
            grand_CVw  = mean(CVw, na.rm = TRUE),
            CVb        = (grand_sd / grand_mean) * 100,
            .groups    = "drop")

notes_cvb_label_series <- notes_cvb_series %>%
  mutate(dynamics_w = case_when(grand_CVw >= 12 ~ 'dynamic',
                                grand_CVw >= 5  ~ 'intermediate',
                                grand_CVw >= 0  ~ 'static'),
         dynamics_b = case_when(CVb >= 12 ~ 'dynamic',
                                CVb >= 5  ~ 'intermediate',
                                CVb >= 0  ~ 'static'),
         CV_ratio = CVb / grand_CVw)
#
#### CV tables: within and between males in series ####
#
# within
cvw_table_data_notes_series <- notes_cvw_long_series %>%
  filter(!is.na(CVw)) %>%
  group_by(frogID, note_type_relabel) %>%
  arrange(frogID, note_type_relabel, CVw) %>%
  select(frogID, feature, n, mean, sd, sem, CVw)
cvw_ft_notes_series <- cvw_table_data_notes_series %>%
  flextable() %>%
  set_header_labels(note_type_relabel = "note",
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
  add_header_lines("Variation in Note Features in Series by Male") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvw_ft_notes_series
save_as_docx(cvw_ft_notes_series, path = "./v15_graphs/final/tables/CVw_table_notes_in_series.docx")
#
# between
cvb_table_data_notes_series <- notes_cvb_label_series %>%
  group_by(note_type_relabel) %>%
  select(feature, n, grand_mean, grand_sd, grand_CVw, CVb, CV_ratio)
cvb_ft_notes_series <- cvb_table_data_notes_series %>%
  flextable() %>%
  set_header_labels(note_type_relabel = "note",
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
  add_header_lines("Variation in Note Features in Series") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvb_ft_notes_series
save_as_docx(cvb_ft_notes_series, path = "./v15_graphs/final/tables/CVb_table_notes_in_series.docx")
#
#### CV calcs: within and between males in series by seq order ####
# within
notes_cvw_series_seq <- all_notes_final %>%
  group_by(frogID, note_type_relabel, series_seq_cat) %>%
  filter(INI < 10,
         series_type %in% c("three", "four", "five"),
         series_seq_cat %in% c("1","50"),
         note_type_relabel %in% c("1", "2")) %>%
  summarise(n = n(),
            notesdur.mean    = mean(note.dur.s, na.rm = TRUE),
            notesdur.sd      = sd(note.dur.s, na.rm = TRUE),
            notesdur.sem     = sem(note.dur.s),
            notesdur.CVw     = (notesdur.sd / notesdur.mean) * 100,
            notesini.mean    = mean(INI, na.rm = TRUE),
            notesini.sd      = sd(INI, na.rm = TRUE),
            notesini.sem     = sem(INI),
            notesini.CVw     = (notesini.sd / notesini.mean) * 100,
            notesino.mean    = mean(note.period, na.rm = TRUE),
            notesino.sd      = sd(note.period, na.rm = TRUE),
            notesino.sem     = sem(note.period),
            notesino.CVw     = (notesino.sd / notesino.mean) * 100,
            notesmeandom.mean = mean(meandom, na.rm = TRUE),
            notesmeandom.sd   = sd(meandom, na.rm = TRUE),
            notesmeandom.sem  = sem(meandom),
            notesmeandom.CVw  = (notesmeandom.sd / notesmeandom.mean) * 100,
            .groups = "drop")

notes_cvw_long_series_seq <- notes_cvw_series_seq %>%
  pivot_longer(cols = -c(series_seq_cat, frogID, note_type_relabel, n),
               names_to = c("feature", ".value"),
               names_pattern = "(.+)\\.(mean|sd|sem|CVw)$") %>%
  mutate(dynamics = case_when(CVw >= 12 ~ 'dynamic',
                              CVw >= 5  ~ 'intermediate',
                              CVw >= 0  ~ 'static'))
# between
notes_cvb_series_seq <- notes_cvw_long_series_seq %>%
  group_by(series_seq_cat, feature, note_type_relabel) %>%
  summarize(n          = n(),
            grand_mean = mean(mean, na.rm = TRUE),
            grand_sd   = sd(mean, na.rm = TRUE),
            grand_sem  = sem(mean),
            grand_CVw  = mean(CVw, na.rm = TRUE),
            CVb        = (grand_sd / grand_mean) * 100,
            .groups    = "drop")

notes_cvb_label_series_seq <- notes_cvb_series_seq %>%
  mutate(dynamics_w = case_when(grand_CVw >= 12 ~ 'dynamic',
                                grand_CVw >= 5  ~ 'intermediate',
                                grand_CVw >= 0  ~ 'static'),
         dynamics_b = case_when(CVb >= 12 ~ 'dynamic',
                                CVb >= 5  ~ 'intermediate',
                                CVb >= 0  ~ 'static'),
         CV_ratio = CVb / grand_CVw)
#
#### CV tables: between and within males in series by seq order ####
#
# within
cvw_table_data_notes_series_seq <- notes_cvw_long_series_seq %>%
  filter(!is.na(CVw)) %>%
  group_by(series_seq_cat, frogID, note_type_relabel) %>%
  arrange(frogID, note_type_relabel, CVw) %>%
  select(frogID, feature, n, mean, sd, sem, CVw)
cvw_ft_notes_series_seq <- cvw_table_data_notes_series_seq %>%
  flextable() %>%
  set_header_labels(series_seq_cat = "seq order",
                    note_type_relabel = "note",
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
  add_header_lines("Variation in Note Features in Series by Seq Order by Male") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvw_ft_notes_series_seq
save_as_docx(cvw_ft_notes_series_seq, path = "./v15_graphs/final/tables/CVw_table_notes_in_series_seq.docx")
#
# between
cvb_table_data_notes_series_seq <- notes_cvb_label_series_seq %>%
  group_by(series_seq_cat,note_type_relabel) %>%
  select(feature, n, grand_mean, grand_sd, grand_CVw, CVb, CV_ratio)
cvb_ft_notes_series_seq <- cvb_table_data_notes_series_seq %>%
  flextable() %>%
  set_header_labels(series_seq_cat="seq order",
                    note_type_relabel="note",
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
  add_header_lines("Variation in Note Features in Series by Seq Order") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvb_ft_notes_series_seq
save_as_docx(cvb_ft_notes_series_seq, path = "./v15_graphs/final/tables/CVb_table_notes_series_seq.docx")
#















# graph within series
all_colors <- c(TimeCalls, TimeNotes)
calls_plot_A <- cvb_table_data_calls_series %>%
  filter(feature %in% c("callsdur", "callsici", "callsioi")) %>%
  mutate(note_type_relabel = NA,
         feature_note = feature)
notes_plot_A <- notes_cvb_label %>%
  filter((feature %in% c("notesdur", "notesmeandom")) |
      (feature %in% c("notesini" ,"notesino") & note_type_relabel == "1")) %>%
  mutate(feature_note = paste0(feature, "_note", note_type_relabel))

combined_A <- bind_rows(calls_plot_A, notes_plot_A)
shape_values <- c("callsdur"            = 16,
                  "callsici"            = 16,
                  "callsioi"            = 16,
                  "notesdur_note1"      = 17,
                  "notesdur_note2"      = 15,
                  "notesini_note1"      = 17,
                  "notesino_note1"      = 17,
                  "notesmeandom_note1"  = 17,
                  "notesmeandom_note2"  = 15)
ggscatter(combined_A, x = "grand_CVw", y = "CV_ratio",
          color = "feature_note", shape = "feature_note", size = 7, alpha = 1) +
  xlab("CVw") +
  ylab("CVb / CVw") +
  geom_vline(xintercept = 12, linetype = "dashed", alpha = 0.7) +
  geom_hline(yintercept = 1,  linetype = "dashed", alpha = 0.7) +
  scale_color_manual(values = all_colors,
                     labels = c("callsdur"           = "call dur",
                                "callsici"           = "ICI",
                                "callsioi"           = "ICO",
                                "notesdur_note1"     = "dur N1",
                                "notesdur_note2"     = "dur N2",
                                "notesini_note1"     = "INI",
                                "notesino_note1"     = "INO",
                                "notesmeandom_note1" = "DF N1",
                                "notesmeandom_note2" = "DF N2")) +
  scale_shape_manual(values = shape_values,
                     labels = c("callsdur"           = "call dur",
                                "callsici"           = "ICI",
                                "callsioi"           = "ICO",
                                "notesdur_note1"     = "dur N1",
                                "notesdur_note2"     = "dur N2",
                                "notesini_note1"     = "INI",
                                "notesino_note1"     = "INO",
                                "notesmeandom_note1" = "DF N1",
                                "notesmeandom_note2" = "DF N2")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line        = element_line(colour = "black", size = .3),
        legend.position  = "right",
        axis.title.x     = element_text(size = 20),
        axis.title.y     = element_text(size = 20),
        axis.text.x      = element_text(size = 16),
        axis.text.y      = element_text(size = 16))
ggsave("./v15_graphs/cvw_cvratio_panelA.pdf", width = 7, height = 5, dpi = 300, limitsize = FALSE)
ggsave("./v15_graphs/cvw_cvratio_panelA.svg", width = 7, height = 5, dpi = 300, limitsize = FALSE)
#
##
# CV ratio of call and note features by sequence order
#
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
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16, angle = 30, hjust = 0.7,vjust = 0.9),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 20),
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
#### CV tables ####
#
## ALL DATA
cvw_table_data_notes <- notes_cvw_long %>%
  filter(!is.na(CVw)) %>%
  arrange(frogID, CVw) %>%
  select(frogID, feature, n, mean, sd, sem, CVw)
cvw_ft_notes <- cvw_table_data_notes %>%
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
  add_header_lines("Variation in Note Features by Male") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvw_ft_notes
save_as_docx(cvw_ft_notes, path = "./v15_graphs/tables/notes_CVw_table.docx")

cvb_table_data_notes <- notes_cvb_label %>%
  select(feature, n, grand_mean, grand_sd, grand_CVw, CVb, CV_ratio)
cvb_ft_notes <- cvb_table_data_notes %>%
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
  add_header_lines("Variation in Note Timing Features") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvb_ft_notes
save_as_docx(cvb_ft_notes, path = "./v15_graphs/tables/notes_cvb_table.docx")
#
#
## within series, JUST NOTE 1
cvw_table_data_notes_series <- notes_cvw_long_series %>%
  filter(!is.na(CVw)) %>%
  arrange(frogID, CVw) %>%
  select(frogID, feature, n, mean, sd, sem, CVw)
cvw_ft_notes_series <- cvw_table_data_notes_series %>%
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
  add_header_lines("Variation in Note Features by Male") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvw_ft_notes_series
save_as_docx(cvw_ft_notes_series, path = "./v15_graphs/tables/notes_CVw_table_within_series.docx")

cvb_table_data_notes_series <- notes_cvb_label_series %>%
  select(feature, n, grand_mean, grand_sd, grand_CVw, CVb, CV_ratio)
cvb_ft_notes_series <- cvb_table_data_notes_series %>%
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
  add_header_lines("Variation in Note Timing Features") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvb_ft_notes_series
save_as_docx(cvb_ft_notes, path = "./v15_graphs/tables/notes_cvb_table_within_series.docx")
#
#### within male note feature summary ####
#
table2_summary <- all_notes_final %>%
  group_by(frogID, note_type_relabel) %>%
  filter(INI < 10) %>%
  filter(note_type_relabel != 'uc') %>% 
  dplyr::summarise(n = n(),
    across(all_of(cv_vars), list(
      mean = ~mean(., na.rm = TRUE),
      sd = ~sd(., na.rm = TRUE),
      min = ~min(., na.rm = TRUE),
      max = ~max(., na.rm = TRUE)), 
      .names = "{.col}.{.fn}"),
    .groups = "drop")
table2_list <- list()

for (frog in unique(table2_summary$frogID)) {
  frog_data <- table2_summary %>% filter(frogID == frog)
  
  temp_table <- data.frame(
    frogID = frog,
    Feature = cv_vars,
    stringsAsFactors = FALSE)
  
  for (note_type in note_types) {
    col_name_base <- gsub(" ", "_", note_type)
    note_data <- frog_data %>% filter(note_type_relabel == note_type)
    
    if (nrow(note_data) > 0) {
      mean_sd_col <- sapply(cv_vars, function(var) {
        mean_val <- note_data[[paste0(var, ".mean")]]
        sd_val <- note_data[[paste0(var, ".sd")]]
        if (length(mean_val) > 0 && !is.na(mean_val)) {
          sprintf("%.3f ± %.3f", mean_val, sd_val) 
        } else {
          "NA"
        }
      })
      
      range_col <- sapply(cv_vars, function(var) {
        min_val <- note_data[[paste0(var, ".min")]]
        max_val <- note_data[[paste0(var, ".max")]]
        if (length(min_val) > 0 && !is.na(min_val)) {
          sprintf("(%.3f - %.3f)", min_val, max_val) 
        } else {
          "NA"
        }
      })
      
      temp_table[[paste0(col_name_base, "_mean_sd")]] <- mean_sd_col
      temp_table[[paste0(col_name_base, "_range")]] <- range_col
    } else {
      temp_table[[paste0(col_name_base, "_mean_sd")]] <- "NA"
      temp_table[[paste0(col_name_base, "_range")]] <- "NA"
    }
  }
  
  table2_list[[length(table2_list) + 1]] <- temp_table
}

table2_display <- bind_rows(table2_list)
ft2 <- flextable(table2_display)
header_labels2 <- list(frogID = "Frog ID", Feature = "Note Feature")
header_row_values2 <- c("Frog ID", "Note Feature")
header_row_colwidths2 <- c(1, 1)
col_positions_for_italic2 <- c()
current_col2 <- 3  
for (i in seq_along(note_types)) {
  note_type <- note_types[i]
  col_name_base <- gsub(" ", "_", note_type)
  
  header_labels2[[paste0(col_name_base, "_mean_sd")]] <- "mean ± sd"
  header_labels2[[paste0(col_name_base, "_range")]] <- "range (min - max)"
  
  header_row_values2 <- c(header_row_values2, note_type)
  header_row_colwidths2 <- c(header_row_colwidths2, 2)
  
  col_positions_for_italic2 <- c(col_positions_for_italic2, current_col2)
  current_col2 <- current_col2 + 2
}

if (sum(header_row_colwidths2) != ncol(table2_display)) {
  stop("ERROR: Sum of colwidths does not match number of columns for Table 2!")
}

ft2 <- set_header_labels(ft2, frogID = "Frog ID", Feature = "Note Feature")

for (i in seq_along(note_types)) {
  note_type <- note_types[i]
  col_name_base <- gsub(" ", "_", note_type)
  
  labels_to_set <- list()
  labels_to_set[[paste0(col_name_base, "_mean_sd")]] <- "mean ± sd"
  labels_to_set[[paste0(col_name_base, "_range")]] <- "range (min - max)"
  
  ft2 <- do.note(set_header_labels, c(list(x = ft2), labels_to_set))
}

ft2 <- ft2 %>%
  add_header_row(values = header_row_values2,
                 colwidths = header_row_colwidths2) %>%
  italic(i = 2, j = col_positions_for_italic2, part = "header") %>%
  align(align = "center", part = "header") %>%
  align(j = 1:2, align = "left", part = "body") %>%
  align(j = 3:ncol(table2_display), align = "center", part = "body") %>%
  border_inner() %>%
  border_outer(border = fp_border(width = 2)) %>%
  merge_v(j = 1) %>%
  autofit()
save_as_docx(ft2, path = "./v15_graphs/tables/table2_note_features_by_frog.docx")
##
#
#### supp: note categorization (no spectral features) ####
#
# graphs produced in R, figure aesthetics modified in Adobe Illustrator
#

#
#note type proportions
notes_prop <- all_notes %>%
  mutate(Note_type = case_when(
    Note_type %in% c(11, 12) ~ "uc",
    TRUE ~ as.character(Note_type)  # convert to character to mix with "uc"
  )) %>%
  group_by(Note_type) %>%
  summarise(n = n()) %>%
  mutate(
    proportion = n / sum(n),
    percent = proportion * 100
  )
#pie chart
ggplot(notes_prop, aes(x="", y=proportion, fill=Note_type)) +
  geom_bar(stat = "identity", width=1, alpha=0.7) +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c(NoteColors)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  guides(fill=guide_legend(title="note type"))
# note_type    n  proportion percent
# 1          2890    0.507    50.7  
# 2          2747    0.482    48.2   
# 3           48     0.00842   0.842
# uc          16    0.00281   0.281
ggsave("./v15_graphs/notes_pie.pdf", width=6, height=6, dpi = 300, limitsize = FALSE) 
ggsave("./v15_graphs/notes_pie.svg", width=6, height=6, dpi = 300, limitsize = FALSE) 
#
# duration boxplot comparison
all_notes_final_10 %>%
  filter(note_type_relabel != 'uc') %>%
  ggplot(aes(x=note_type_relabel, y=note.dur.s, color=note_type_relabel)) +
  geom_boxplot(outlier.alpha = 0, width = 0.5) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.2) +
  xlab("note type") +
  ylab("duration (s)") + 
  scale_color_manual(values=c(NoteColors)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .5),
        axis.title.y=element_text(size=20),  # X axis title
        axis.title.x=element_text(size=20),  # Y axis title
        axis.text.x=element_text(size=16),  # X axis text
        axis.text.y=element_text(size=16),
        legend.position = "none")
ggsave("./v15_graphs/note_dur_type_boxplot.pdf", width=4, height=4, dpi = 300) 
ggsave("./v15_graphs/note_dur_type_boxplot.svg", width=4, height=4, dpi = 300) 
#
# DF boxplot comparison
all_notes_final_10 %>%
  filter(note_type_relabel != 'uc') %>%
  ggplot(aes(x=note_type_relabel, y=meandom, color=note_type_relabel)) +
  geom_boxplot(outlier.alpha = 0, width = 0.5) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.2) +
  xlab("note type") +
  ylab("dominant freq (kHz)") + 
  scale_color_manual(values=c(NoteColors)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .5),
        axis.title.y=element_text(size=20),  # X axis title
        axis.title.x=element_text(size=20),  # Y axis title
        axis.text.x=element_text(size=16),  # X axis text
        axis.text.y=element_text(size=16),
        legend.position = "none")
ggsave("./v15_graphs/note_domfreq_type_boxplot.pdf", width=4, height=4, dpi = 300) 
ggsave("./v15_graphs/note_domfreq_boxplot.svg", width=4, height=4, dpi = 300) 
#
# ini boxplot comparison
all_notes_final_10 %>%
  filter(note_type_relabel != 'uc') %>%
  ggplot(aes(x=note_type_relabel, y=INI, color=note_type_relabel)) +
  geom_boxplot(outlier.alpha = 0, width = 0.5) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.2) +
  xlab("note type") +
  ylab("ini (s)") + 
  ylim(0,5) +
  scale_color_manual(values=c(NoteColors)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .5),
        axis.title.y=element_text(size=20),  # X axis title
        axis.title.x=element_text(size=20),  # Y axis title
        axis.text.x=element_text(size=16),  # X axis text
        axis.text.y=element_text(size=16),
        legend.position = "none")
ggsave("./v15_graphs/note_ini_type_boxplot.pdf", width=4, height=4, dpi = 300) 
ggsave("./v15_graphs/note_ini_type_boxplot.svg", width=4, height=4, dpi = 300) 
#
##
#### CVw and CVb note tables ####
#
cvw_table_data <- notes_cvw_long %>%
  filter(!is.na(series_seq_cat)) %>%
  filter(!is.na(CVw)) %>%
  mutate(note_type_relabel = case_when(note_type_relabel == "1" ~ "note 1",
                                       note_type_relabel == "2" ~ "note 2",
                                       note_type_relabel == "3" ~ "note 3",
                                       TRUE ~ note_type_relabel)) %>%
  mutate(note_type_relabel = factor(note_type_relabel, levels = c('note 1', 'note 2', 'note 3'))) %>%
  mutate(series_seq_cat = factor(series_seq_cat, levels = c('0', '1', '50', '100'))) %>%
  arrange(series_seq_cat, note_type_relabel, frogID, CVw) %>%
  select(series_seq_cat, note_type_relabel, frogID, feature, n, mean, sd, sem, CVw)
cvb_table_data <- notes_cvb_label %>%
  filter(!is.na(series_seq_cat)) %>%
  filter(note_type_relabel %in% c("1", "2")) %>%
  mutate(note_type_relabel = case_when(note_type_relabel == "1" ~ "note 1",
                                       note_type_relabel == "2" ~ "note 2",
                                       TRUE ~ note_type_relabel)) %>%
  mutate(note_type_relabel = factor(note_type_relabel, levels = c('note 1', 'note 2'))) %>%
  mutate(series_seq_cat = factor(series_seq_cat, levels = c('0', '1', '50', '100'))) %>%
  arrange(series_seq_cat, note_type_relabel, CVb) %>%
  select(series_seq_cat, note_type_relabel, feature, n, grand_mean, grand_sd, grand_CVw, CVb, CV_ratio)

cvb_ft <- cvb_table_data %>%
  flextable() %>%
  set_header_labels(series_seq_cat = "position",
                    note_type_relabel = "note",
                    feature = "feature",
                    n = "N",
                    grand_mean = "mean",
                    grand_sd = "sd",
                    grand_CVw = "mean CVw",
                    CVb = "CVb (%)",
                    CV_ratio = "CVb/CVw Ratio") %>%
  colformat_double(j = c("grand_mean", "grand_sd", "grand_CVw", "CVb", "CV_ratio"),
                   digits = 2) %>%
  bg(i = ~ note_type_relabel == "note 1",
     j = "note_type_relabel",
     bg = NoteColors3["note 1"]) %>%
  bg(i = ~ note_type_relabel == "note 2",
     j = "note_type_relabel",
     bg = NoteColors3["note 2"]) %>%
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
  add_header_lines("Variation in Note Features by Series Position") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()
cvb_ft
save_as_docx(cvb_ft, path = "./v15_graphs/tables/notes_doublets_cvb_table.docx")
#

#### feature pearson corr by note ####
##
#
# note 1 
#
note1 <- all_notes_final %>%
  filter(note_type_relabel=="1") %>%
  select(note.dur.s,INI,note.period,dfrange,modindx,skew,entropy,sfm,
  meanfreq,meanpeakf,meandom,maxdom,enddom,freq.IQR,time.median,time.IQR) %>%
  filter(INI<10)
#
cor_matrix1 <- cor(note1, use = "complete.obs", method = "pearson")
cor_melted1 <- melt(cor_matrix1)
names(cor_melted1) <- c("Var1", "Var2", "Correlation")

p1 <- ggplot(cor_melted1, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "#778868", high = "#de8a5a", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), 
                       name = "pearson\ncorrelation") +
  ggtitle("note 1") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        legend.position = "right") +
  coord_fixed()

high_cor_indices1 <- which(abs(cor_matrix1) >= 0.85 & cor_matrix1 != 1, arr.ind = TRUE)
high_cor_df1 <- data.frame(Variable1 = rownames(cor_matrix1)[high_cor_indices1[, 1]],
                          Variable2 = colnames(cor_matrix1)[high_cor_indices1[, 2]],
                          Correlation = cor_matrix1[high_cor_indices1])

high_cor_df1 <- high_cor_df1[!duplicated(t(apply(high_cor_df1[, 1:2], 1, sort))), ]
high_cor_df1 <- high_cor_df1[order(abs(high_cor_df1$Correlation), decreasing = TRUE), ]

high_cor_df1$Correlation <- round(high_cor_df1$Correlation, 3)

table_grob1 <- tableGrob(high_cor_df1, rows = NULL, 
                        theme = ttheme_default(base_size = 10))

title_grob1 <- textGrob("correlations r >= 0.85)", 
                       gp = gpar(fontsize = 14, fontface = "bold"))
combined_plot1 <- grid.arrange(p1,
                              arrangeGrob(title_grob1, table_grob1, heights = c(0.1, 1)),
                              ncol = 1,
                              heights = c(1, 0.4))
ggsave("./v15_graphs/note1_spec_corr.pdf", combined_plot1, 
       width = 8, height = 8, dpi = 300, limitsize = FALSE)
##
#
# note 2
#
note2 <- all_notes_final %>%
  filter(note_type_relabel=="2") %>%
  select(note.dur.s,INI,note.period,dfrange,modindx,skew,entropy,sfm,
         meanfreq,meanpeakf,meandom,maxdom,enddom,freq.IQR,time.median,time.IQR) %>%
  filter(INI<10)
#
cor_matrix2 <- cor(note2, use = "complete.obs", method = "pearson")
cor_melted2 <- melt(cor_matrix2)
names(cor_melted2) <- c("Var1", "Var2", "Correlation")

p2 <- ggplot(cor_melted2, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "#778868", high = "#de8a5a", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), 
                       name = "pearson\ncorrelation") +
  ggtitle("note 2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        legend.position = "right") +
  coord_fixed()

high_cor_indices2 <- which(abs(cor_matrix2) >= 0.85 & cor_matrix2 != 1, arr.ind = TRUE)
high_cor_df2 <- data.frame(Variable1 = rownames(cor_matrix2)[high_cor_indices2[, 1]],
                           Variable2 = colnames(cor_matrix2)[high_cor_indices2[, 2]],
                           Correlation = cor_matrix2[high_cor_indices2])

high_cor_df2 <- high_cor_df2[!duplicated(t(apply(high_cor_df2[, 1:2], 1, sort))), ]
high_cor_df2 <- high_cor_df2[order(abs(high_cor_df2$Correlation), decreasing = TRUE), ]
high_cor_df2$Correlation <- round(high_cor_df2$Correlation, 3)

table_grob2 <- tableGrob(high_cor_df2, rows = NULL, 
                         theme = ttheme_default(base_size = 10))

title_grob2 <- textGrob("correlations r >= 0.85)", 
                        gp = gpar(fontsize = 14, fontface = "bold"))
combined_plot2 <- grid.arrange(p2,
                               arrangeGrob(title_grob2, table_grob2, heights = c(0.1, 1)),
                               ncol = 1,
                               heights = c(1, 0.4))
ggsave("./v15_graphs/note2_spec_corr.pdf", combined_plot2, 
       width = 8, height = 8, dpi = 300, limitsize = FALSE)
##
#
# note 3
#
note3 <- all_notes_final %>%
  filter(note_type_relabel=="3") %>%
  select(note.dur.s,INI,note.period,dfrange,modindx,skew,entropy,sfm,
         meanfreq,meanpeakf,meandom,maxdom,enddom,freq.IQR,time.median,time.IQR) %>%
  filter(INI<10)
#
cor_matrix3 <- cor(note3, use = "complete.obs", method = "pearson")
cor_melted3 <- melt(cor_matrix3)
names(cor_melted3) <- c("Var1", "Var2", "Correlation")

p3 <- ggplot(cor_melted3, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "#778868", high = "#de8a5a", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), 
                       name = "pearson\ncorrelation") +
  ggtitle("note 3") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        legend.position = "right") +
  coord_fixed()

high_cor_indices3 <- which(abs(cor_matrix3) >= 0.85 & cor_matrix3 != 1, arr.ind = TRUE)
high_cor_df3 <- data.frame(Variable1 = rownames(cor_matrix3)[high_cor_indices3[, 1]],
                           Variable2 = colnames(cor_matrix3)[high_cor_indices3[, 2]],
                           Correlation = cor_matrix3[high_cor_indices3])

high_cor_df3 <- high_cor_df3[!duplicated(t(apply(high_cor_df3[, 1:2], 1, sort))), ]
high_cor_df3 <- high_cor_df3[order(abs(high_cor_df3$Correlation), decreasing = TRUE), ]
high_cor_df3$Correlation <- round(high_cor_df3$Correlation, 3)

table_grob3 <- tableGrob(high_cor_df3, rows = NULL, 
                         theme = ttheme_default(base_size = 10))
title_grob3 <- textGrob("correlations r >= 0.85)", 
                        gp = gpar(fontsize = 14, fontface = "bold"))
combined_plot3 <- grid.arrange(p3,
                               arrangeGrob(title_grob3, table_grob3, heights = c(0.1, 1)),
                               ncol = 1,
                               heights = c(1, 0.4))
ggsave("./v15_graphs/note3_spec_corr.pdf", combined_plot3, 
       width = 8, height = 8, dpi = 300, limitsize = FALSE)
#
#
#### LMMs: note features based on series patterning ####
#- WITHIN SERIES ONLY #
all_notes_within <- all_notes_final_10 %>%
  filter(series_seq_cat %in% c("1", "50"))
all_notes_within$series_type <- factor(all_notes_within$series_type, 
                                       levels = c("three", "four", "five"))
all_notes_within$series_seq_cat <- factor(all_notes_within$series_seq_cat, 
                                          levels = c("1", "50"))
all_notes_within$note_type_relabel <- relevel(as.factor(all_notes_within$note_type_relabel), 
                                              ref = "1")

# NOTE DURATION MODEL (WITHIN-SERIES) - BY NOTE TYPE
note_types <- c("1", "2")
for (nt in note_types) {
  
  sink(paste0("./v15_graphs/stats/raw_output_note.dur.s_note", nt, "_within_series.txt"))
  
  tryCatch({
    
    note_data <- all_notes_within %>%
      filter(note_type_relabel == nt) %>%
      droplevels()
    
    cat("Sample size:", nrow(note_data), "\n")
    cat("series_type levels:", levels(note_data$series_type), "\n")
    cat("series_seq_cat levels:", levels(note_data$series_seq_cat), "\n")
    print(table(note_data$series_type, note_data$series_seq_cat))
    
    final_model <- lmer(note.dur.s ~ series_type * series_seq_cat + (1|frogID),
                        data = note_data,
                        REML = TRUE)
    
    print(summary(final_model))
    
    anova_full <- Anova(final_model, type = "III", test.statistic = "Chisq")
    print(anova_full)
    
    anova_df <- as.data.frame(anova_full)
    anova_df$term <- rownames(anova_df)
    interaction_p <- anova_df[anova_df$term == "series_type:series_seq_cat", "Pr(>Chisq)"]
    
    if (interaction_p >= 0.05) {
      
      print("Interaction not significant, removing and refitting additive model.")
      
      final_model <- lmer(note.dur.s ~ series_type + series_seq_cat + (1|frogID),
                          data = note_data,
                          REML = TRUE)
      
      print(summary(final_model))
      
      anova_final <- Anova(final_model, type = "III", test.statistic = "Chisq")
      print(anova_final)
      
      interaction_significant <- FALSE
      
    } else {
      
      print("Interaction is significant, keeping full model.")
      interaction_significant <- TRUE
    }
    
    if (interaction_significant) {
      
      emm_series_pw <- emmeans(final_model, pairwise ~ series_type | series_seq_cat, adjust = "bonferroni")
      emm_seq_pw <- emmeans(final_model, pairwise ~ series_seq_cat | series_type, adjust = "bonferroni")
      
      print(emm_series_pw)
      print(emm_seq_pw)
      
      tryCatch({
        print(eff_size(emm_series_pw$contrasts, sigma = sigma(final_model), edf = df.residual(final_model), method = "identity"))
      }, error = function(e) {
        print("Effect size calculation failed for series_type | series_seq_cat")
      })
      
      tryCatch({
        print(eff_size(emm_seq_pw$contrasts, sigma = sigma(final_model), edf = df.residual(final_model), method = "identity"))
      }, error = function(e) {
        print("Effect size calculation failed for series_seq_cat | series_type")
      })
    }
    
    r2_values <- r2(final_model)
    icc_value <- icc(final_model)
    
    print(paste("R2 marginal:", round(r2_values$R2_marginal, 3)))
    print(paste("R2 conditional:", round(r2_values$R2_conditional, 3)))
    print(paste("ICC:", round(icc_value$ICC_adjusted, 3)))
    
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
  })
  
  sink()
}
#
# INI MODEL (WITHIN-SERIES) - NOTE TYPE 1 ONLY
sink("./v15_graphs/stats/raw_output_INI_within_series.txt")
tryCatch({
  
  ini_data <- all_notes_within %>%
    filter(note_type_relabel == "1") %>%
    droplevels()
  
  cat("Sample size:", nrow(ini_data), "\n")
  cat("series_type levels:", levels(ini_data$series_type), "\n")
  cat("series_seq_cat levels:", levels(ini_data$series_seq_cat), "\n")
  print(table(ini_data$series_type, ini_data$series_seq_cat))
  
  final_model <- lmer(INI ~ series_type * series_seq_cat + (1|frogID),
                      data = ini_data,
                      REML = TRUE)
  
  print(summary(final_model))
  
  anova_full <- Anova(final_model, type = "III", test.statistic = "Chisq")
  print(anova_full)
  
  anova_df <- as.data.frame(anova_full)
  anova_df$term <- rownames(anova_df)
  interaction_p <- anova_df[anova_df$term == "series_type:series_seq_cat", "Pr(>Chisq)"]
  
  if (interaction_p >= 0.05) {
    
    print("Interaction not significant, removing and refitting additive model.")
    
    final_model <- lmer(INI ~ series_type + series_seq_cat + (1|frogID),
                        data = ini_data,
                        REML = TRUE)
    
    print(summary(final_model))
    
    anova_final <- Anova(final_model, type = "III", test.statistic = "Chisq")
    print(anova_final)
    
    interaction_significant <- FALSE
    
  } else {
    
    print("Interaction is significant, keeping full model.")
    interaction_significant <- TRUE
  }
  
  if (interaction_significant) {
    
    emm_series_pw <- emmeans(final_model, pairwise ~ series_type | series_seq_cat, adjust = "bonferroni")
    emm_seq_pw <- emmeans(final_model, pairwise ~ series_seq_cat | series_type, adjust = "bonferroni")
    
    print(emm_series_pw)
    print(emm_seq_pw)
    
    tryCatch({
      print(eff_size(emm_series_pw$contrasts, sigma = sigma(final_model), edf = df.residual(final_model), method = "identity"))
    }, error = function(e) {
      print("Effect size calculation failed for series_type | series_seq_cat")
    })
    
    tryCatch({
      print(eff_size(emm_seq_pw$contrasts, sigma = sigma(final_model), edf = df.residual(final_model), method = "identity"))
    }, error = function(e) {
      print("Effect size calculation failed for series_seq_cat | series_type")
    })
  }
  
  r2_values <- r2(final_model)
  icc_value <- icc(final_model)
  
  print(paste("R2 marginal:", round(r2_values$R2_marginal, 3)))
  print(paste("R2 conditional:", round(r2_values$R2_conditional, 3)))
  print(paste("ICC:", round(icc_value$ICC_adjusted, 3)))
  
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
})

sink()
#
##### supp: note features by males and with eco variables ####
#
eco_vars <- list(
  list(var = "humidity_temp", label = "humidity"),
  list(var = "air_temp_c", label = "temperature (C)"),
  list(var = "weight_g", label = "weight (g)"),
  list(var = "svl_mm", label = "svl (mm)"),
  list(var = "femur_length_mm", label = "femur length (mm)"))

note_types <- list(list(note = 1, color = "#f6edbd"),
                   list(note = 2, color = "#b5b991"),
                   list(note = 3, color = "#de8a5a"))

for (note_info in note_types) {
  note_num <- note_info$note
  note_color <- note_info$color
  
  note_data <- all_notes_final_10 %>%
    filter(note_type_relabel == note_num)
  
  for (eco_info in eco_vars) {
    eco_var <- eco_info$var
    eco_label <- eco_info$label
    
    eco_filename <- gsub(" ", "_", gsub("\\(|\\)", "", eco_label))
    eco_filename <- gsub("_+", "_", eco_filename)  # Remove multiple underscores
    eco_filename <- tolower(eco_filename)
    
    for (y_var in cv_vars) {
      
      plot_data <- note_data %>%
        filter(!is.na(.data[[eco_var]])) %>%
        mutate(!!eco_var := as.numeric(.data[[eco_var]])) %>%
        filter(!is.na(.data[[y_var]]))
      
      if (nrow(plot_data) == 0) next
      
      p <- ggplot(plot_data, aes(x = .data[[eco_var]], y = .data[[y_var]])) +
        geom_point(size = 3, alpha = 0.5, color = note_color) +
        geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
        stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
          formula = y ~ x,
          parse = TRUE,
          color = "black",
          size = 5,
          label.x.npc = "right",
          label.y.npc = "top",
          hjust = 0,
          vjust = 1) +
        xlab(eco_label) +
        ylab(y_var) +
        theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size = .7),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          legend.position = "none")
      
      if (y_var %in% c("INI", "note.period") && note_num == 2) {
        p <- p + ylim(0, 5)
      }
      
      filename_base <- paste0(eco_filename, "_", y_var, "_note", note_num)
      
      ggsave(paste0("./v15_graphs/eco/", filename_base, ".pdf"),
        plot = p,
        width = 6,
        height = 4,
        dpi = 300,
        limitsize = FALSE)
      
      ggsave(paste0("./v15_graphs/eco/", filename_base, ".svg"),
        plot = p,
        width = 6,
        height = 4,
        dpi = 300,
        limitsize = FALSE)
      
      cat(paste0("Saved: ", filename_base, "\n"))
    }
  }
}

#### covariation between note 1 and 2 temporal features ####
#
#
## within series INI
#
all_notes_final %>%
  filter(series_type=="three")
  
  filter(note_type == "doublet") %>%
  filter(INI < 5) %>%
  select(frogID, note_num, Note_type, INI) %>%
  pivot_wider(id_cols = c(frogID, note_num),
    names_from = note_type_relabel,
    values_from = INI,
    names_prefix = "INI_note_") %>%
  filter(!is.na(INI_note_1) & !is.na(INI_note_2)) %>%
  filter(frogID=="1281") %>%
  ggplot(aes(x = INI_note_1, y = INI_note_2)) +
  geom_point(size = 3, alpha = 0.5) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    color = "black",
    size = 5,
    label.x = "right",
    label.y = "top") 
  labs(x = "INI Note 1",
       y = "INI Note 2") +
  theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "none")
ggsave("./v15_graphs/series/isi_iso.pdf", width=6, height=4, dpi = 300, limitsize = FALSE)
#
#
#### note analysis: squeaks only ####
#
#
## CV of squeaks
cvb_table_data_note3 <- notes_cvb_label %>%
  filter(!is.na(series_seq_cat)) %>%
  filter(note_type_relabel == "3") %>%  # Filter for only note 3
  mutate(note_type_relabel = case_when(note_type_relabel == "3" ~ "note 3",
                                       TRUE ~ note_type_relabel)) %>%
  mutate(note_type_relabel = factor(note_type_relabel, levels = c('note 3'))) %>%
  mutate(series_seq_cat = factor(series_seq_cat, levels = c('0', '1', '50', '100'))) %>%
  arrange(series_seq_cat, note_type_relabel, CVb) %>%
  select(series_seq_cat, note_type_relabel, feature, n, grand_mean, grand_sd, grand_CVw, CVb, CV_ratio)

cvb_ft_note3 <- cvb_table_data_note3 %>%
  flextable() %>%
  set_header_labels(series_seq_cat = "position",
                    note_type_relabel = "note",
                    feature = "feature",
                    n = "N",
                    grand_mean = "mean",
                    grand_sd = "sd",
                    grand_CVw = "mean CVw",
                    CVb = "CVb (%)",
                    CV_ratio = "CVb/CVw Ratio") %>%
  colformat_double(j = c("grand_mean", "grand_sd", "grand_CVw", "CVb", "CV_ratio"),
                   digits = 2) %>%
  bg(i = ~ note_type_relabel == "note 3",
     j = "note_type_relabel",
     bg = NoteColors3["note 3"]) %>%
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
  add_header_lines("Variation in Note Features by Series Position") %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "feature", align = "left", part = "body") %>%
  autofit()