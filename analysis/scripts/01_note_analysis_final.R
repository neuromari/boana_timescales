################## B. pulchellus Timescales of Call Variability #################
#### author: M. Rodriguez-Santiago
#### last update: 7.8.2026
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
library(tidyr)
library(flextable)
library(ggplot2)
library(tibble)
#
#
#### data import & tidy ####
#
all_notes <- read.csv("./data/all_df_notes_final.csv")
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
specsK <- read.csv("./data/SpectrotemporalOutputs_boanaKH.csv") 
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
ggsave("./analysis/graphs/IRs/int_ratio_ino.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./analysis/graphs/IRs/int_ratio_ino.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
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
tibble(frogID = unique(ks_results_ino$frogID),
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
ggsave("./analysis/graphs/IRs/int_ratio_null_ino.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./analysis/graphs/IRs/int_ratio_null_ino.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
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
save_as_docx(cvw_ft_notes, path = "./analysis/graphs/tables/CVw_table_notes.docx")
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
save_as_docx(cvb_ft_notes, path = "./analysis/graphs/tables/CVb_table_notes.docx")
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
save_as_docx(cvw_ft_notes_series, path = "./analysis/graphs/tables/CVw_table_notes_in_series.docx")
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
save_as_docx(cvb_ft_notes_series, path = "./analysis/graphs/tables/CVb_table_notes_in_series.docx")
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
save_as_docx(cvw_ft_notes_series_seq, path = "./analysis/graphs/tables/CVw_table_notes_in_series_seq.docx")
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
save_as_docx(cvb_ft_notes_series_seq, path = "./analysis/graphs/tables/CVb_table_notes_series_seq.docx")
#
##