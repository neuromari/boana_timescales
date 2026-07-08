################## B. pulchellus Timescales of Call Variability #################
#### author: M. Rodriguez-Santiago
#### last update: 7.8.2026
####
#### this script plots all graphs in final manuscript
####
#### last update: 7.08.2026
##
##
##
#
#### packages needed to run this code ####
#
library(dplyr)
library(ggplot2)
library(ggforce)
library(lmerTest)
library(emmeans)
library(performance)
library(cowplot)
library(patchwork)
library(gridExtra)
library(diptest)
#
#### Figure 2 IR graphs ####
#
# panel a
observed_data_iso <- series_ir %>%
  filter(is.finite(int_ratio_iso), !is.na(int_ratio_iso)) %>%
  mutate(ratio = int_ratio_iso, type = "observed") %>%
  select(ratio, type)
null_data_iso <- nulldist_ISO_plot %>%
  filter(is.finite(null_int_ratio), !is.na(null_int_ratio)) %>%
  mutate(ratio = null_int_ratio, type = "null") %>%
  select(ratio, type)
combined_data_iso <- bind_rows(observed_data_iso, null_data_iso)
ggplot() +
  geom_histogram(data = null_data_iso, aes(x = ratio),
                 color = "black", linewidth=0.3, fill="gray", bins = 100, alpha = 0.8) +
  geom_histogram(data = observed_data_iso, aes(x = ratio), 
                 color = "black", linewidth=0.3, fill="#66cccc", bins = 100, alpha=0.8) +
  xlab("series interval ratio") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .3),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))
ggsave("./manuscript/figures/fig2_a_int_ratio_iso.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./manuscript/figures/fig2_a_int_ratio_iso.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE)
## stats
observed_ratios_iso <- series_ir %>%
  filter(is.finite(int_ratio_iso), !is.na(int_ratio_iso)) %>%
  pull(int_ratio_iso)
null_ratios_iso <- nulldist_ISO_plot %>%
  filter(is.finite(null_int_ratio), !is.na(null_int_ratio)) %>%
  pull(null_int_ratio)
ks_test_iso <- ks.test(observed_ratios_iso, null_ratios_iso)
print(ks_test_iso)
## extract mean and sd of distributions
mean(observed_ratios_iso, na.rm = TRUE)
sd(observed_ratios_iso, na.rm = TRUE)
mean(null_ratios_iso, na.rm = TRUE)
sd(null_ratios_iso, na.rm = TRUE)
#
#
# panel b
all_series_bin %>%
  filter(ISI<10) %>%
  ggplot(aes(x = ISO)) +
  geom_histogram(position = "identity", alpha = 1, bins = 100, color = "black", linewidth=0.3, fill="#66cccc") +
  xlab("inter-series onset (s)") +
  scale_x_continuous(limits = c(0, 5),
                     breaks = seq(0, 5, by = 1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .3),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x =  element_text(size = 16),
        axis.text.y = element_text(size = 16))
ggsave("./manuscript/figures/fig2_b_hist_iso_bout.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./manuscript/figures/fig2_b_hist_iso_bout.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE)
#
#
# panel c
observed_data_ioi <- call_ir %>%
  filter(is.finite(int_ratio_ioi), !is.na(int_ratio_ioi)) %>%
  mutate(ratio = int_ratio_ioi, type = "observed") %>%
  select(ratio, type)
null_data_ioi <- nulldist_IOI_plot %>%
  filter(is.finite(null_int_ratio), !is.na(null_int_ratio)) %>%
  mutate(ratio = null_int_ratio, type = "null") %>%
  select(ratio, type)
combined_data_ioi <- bind_rows(observed_data_ioi, null_data_ioi)
ggplot() +
  geom_histogram(data = null_data_ioi, aes(x = ratio), 
                 color = "black", linewidth=0.3, fill="gray", bins = 100, alpha = 1) +
  geom_histogram(data = observed_data_ioi, aes(x = ratio), 
                 color = "black", linewidth=0.3, fill="#cc99cc", bins = 100, alpha = 0.6) +
  xlab("call interval ratio") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .3),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))
ggsave("./manuscript/figures/fig2_c_int_ratio_ioi_both.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./manuscript/figures/fig2_c_int_ratio_ioi_both.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE)
## stats
observed_ratios_ioi <- call_ir %>%
  filter(is.finite(int_ratio_ioi), !is.na(int_ratio_ioi)) %>%
  pull(int_ratio_ioi)
null_ratios_ioi <- nulldist_IOI_plot %>%
  filter(is.finite(null_int_ratio), !is.na(null_int_ratio)) %>%
  pull(null_int_ratio)
ks_test_ioi <- ks.test(observed_ratios_ioi, null_ratios_ioi)
print(ks_test_ioi)
## extract mean and sd of distributions
mean(null_ratios_ioi, na.rm = TRUE)
sd(null_ratios_ioi, na.rm = TRUE)
mean(observed_ratios_ioi, na.rm = TRUE)
sd(observed_ratios_ioi, na.rm = TRUE)
## diptest for multimodality (new for revision1)
dip_test_ioi <- dip.test(observed_ratios_ioi)
print(dip_test_ioi)
#
#
# panel d
all_calls_bin %>%
  filter(ICI < 10) %>%
  ggplot(aes(x = call_period)) +
  xlab("inter-call onset (s)") +
  geom_histogram(position = "identity", 
                 alpha = 1, 
                 bins = 100, color = "black", linewidth=0.3, fill="#cc99cc") +
  scale_x_continuous(limits = c(0, 5),
                     breaks = seq(0, 5, by = 1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", linewidth = 0.3),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none")
ggsave("./manuscript/figures/fig2_d_hist_ioi_bout_all.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./manuscript/figures/fig2_d_hist_ioi_bout_all.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
# inset: color by two call rhythms
all_calls_bin %>%
  filter(ICI < 1) %>%
  ggplot(aes(x = call_period, fill=series_seq_cat)) +
  xlab("inter-call onset (s)") +
  geom_histogram(position = "identity", 
                 alpha = 1, 
                 bins = 100, color = "black", linewidth=0.3) +
  scale_fill_manual(values = c("0" = "#9f839d",
                               "1" = "#e4cee6",
                               "50" = "#e4cee6",
                               "100" = "#9f839d")) +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_blank(),
        axis.ticks=element_blank(),
        legend.position = "none")
ggsave("./manuscript/figures/fig2_d_inset_hist_ici_bout_zoom.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./manuscript/figures/fig2_d_inset_hist_ici_bout_zoom.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
#
#
# panel e
observed_data <- notes_ir %>%
  filter( is.finite(int_ratio_ino), !is.na(int_ratio_ino)) %>%
  mutate(ratio = int_ratio_ino, type = "Observed") %>%
  select(ratio, type)
null_data <- nulldist_ino %>%
  filter(is.finite(null_int_ratio), !is.na(null_int_ratio)) %>%
  mutate(ratio = null_int_ratio, type = "Null") %>%
  select(ratio, type)
combined_data <- bind_rows(observed_data, null_data)

ggplot() +
  geom_histogram(data = null_data, aes(x = ratio), 
                 color = "black", linewidth=0.3, fill="gray", bins = 100, alpha = 0.8) +
  geom_histogram(data = observed_data, aes(x = ratio), 
                 color = "black", linewidth=0.3, fill="#99cc99", bins = 100, alpha = 1) +
  xlab("note interval ratio") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .3),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))
ggsave("./manuscript/figures/fig2_e_int_ratio_ino.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./manuscript/figures/fig2_e_int_ratio_ino.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
## stats
observed_ratios <- notes_ir %>%
  filter(is.finite(int_ratio_ino), !is.na(int_ratio_ino)) %>%
  pull(int_ratio_ino)
null_ratios <- nulldist_ino %>%
  filter(is.finite(null_int_ratio), !is.na(null_int_ratio)) %>%
  pull(null_int_ratio)
# Two-sample KS test
ks_result <- ks.test(observed_ratios, null_ratios)
print(ks_result)
## extract mean and sd of distributions
mean(null_ratios, na.rm = TRUE)
sd(null_ratios, na.rm = TRUE)
mean(observed_ratios, na.rm = TRUE)
sd(observed_ratios, na.rm = TRUE)
## diptest for multimodality (new for revision1)
dip_test_ino <- dip.test(observed_ratios)
print(dip_test_ino)
#
#
# panel f
all_notes_final %>%
  filter(note.period<10) %>%
  ggplot(aes(x = note.period)) +
  geom_histogram(position = "identity", alpha = 1, bins = 100, color = "black", linewidth=0.3, fill="#99cc99") +
  scale_x_continuous(limits = c(0, 5),
                     breaks = seq(0, 5, by = 1)) +
  xlab("inter-note onset (s)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .3),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x =  element_text(size = 16),
        axis.text.y = element_text(size = 16))
ggsave("./manuscript/figures/fig2_f_ino.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./manuscript/figures/fig2_f_ino.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
#
# inset: color by two note types
all_notes_final %>%
  filter(note.period<1) %>%
  filter(note_type_relabel=="1" | note_type_relabel=="2") %>%
  ggplot(aes(x = note.period, fill=note_type_relabel)) +
  geom_histogram(position = "identity", alpha = 1, bins = 100, color = "black", linewidth=0.3) +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1)) +
  scale_fill_manual(values = c("1" = "#bcdab7",
                               "2" = "#336633")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_blank(),
        axis.ticks=element_blank(),
        legend.position = "none")
ggsave("./manuscript/figures/fig2_f_inset_ino_note_types_zoom.pdf", width=6, height=3.5, dpi = 300, limitsize = FALSE) 
ggsave("./manuscript/figures/fig2_f_inset_ino_note_types_zoom.svg", width=6, height=3.5, dpi = 300, limitsize = FALSE)
#
#### Figure 3: ISI by series length ####
#
# ISI by series bin
all_series_bin %>%
  mutate(series_type = factor(series_bin, 
                              levels = c("1", "2", "3", "4+"),
                              labels = c("solo", "two", "three", "4+"))) %>%
  ggplot(aes(x = series_type, y = ISI, color = series_type)) +
  geom_sina(alpha = 0.5, size = 2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, 
               color = "black", linewidth = 0.6, 
               show.legend = FALSE) +
  xlab("series length") +
  ylab("ISI (s)") +
  ylim(0, 5) +
  scale_color_manual(values = SeriesTypeColors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", linewidth = 0.3),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none")
ggsave("./manuscript/figures/fig3_series_dur_isi.pdf", width=6, height=4, dpi = 300, limitsize = FALSE)
ggsave("./manuscript/figures/fig3_series_dur_isi.svg", width=6, height=4, dpi = 300, limitsize = FALSE) 
#stats
series_type_10 <- all_series_bin %>%
  filter(ISI < 10) %>%
  mutate(frogID = as.factor(frogID),
         log_ISI = log(ISI),
         series_type_bin = factor(series_bin,
                                  levels = c("1", "2", "3", "4+"),
                                  labels = c("solo", "two", "three", "4+")))
mod_type_isi <- lmer(ISI ~ series_type_bin + (1|frogID), data = series_type_10)
summary(mod_type_isi)
emm_isi <- emmeans(mod_type_isi, ~ series_type_bin)
pairs(emm_isi, adjust = "holm")
#
#### Figure 4: temporal variables of CALLS within series ####
#
calls_within <- all_calls_bin %>%
  filter(ICI < 10, 
         series_seq_cat %in% c("1", "50"),
         series_type %in% c("three", "four", "five")) %>%
  mutate(series_seq_cat = factor(series_seq_cat, levels = c("1", "50")),
         series_type = factor(series_type, levels = c("three", "four", "five")))
#
#
# panel a
ggplot(calls_within, aes(x = series_type, y = call_period, color = series_seq_cat)) +
  geom_sina(alpha = 0.7, size = 4,
            position = position_dodge(width = 0.7)) +
  stat_summary(aes(group = series_seq_cat),
               fun = mean, geom = "crossbar", width = 0.3,
               color = "black", linewidth = 0.2,
               show.legend = FALSE,
               position = position_dodge(width = 0.7)) +
  xlab("series length") +
  ylab("inter-call onset (s)") +
  scale_color_manual(name = "  within\n  series\n position",
                     values = c("1" = "#e4cee6", "50" = "#cc99cc"),
                     labels = c("first","mid")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = .3),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14))
ggsave("./manuscript/figures/fig4_a_within_series_ico_type.pdf", width=6, height=4, dpi = 300, limitsize = FALSE)
ggsave("./manuscript/figures/fig4_a_within_series_ico_type.svg", width=6, height=4, dpi = 300, limitsize = FALSE) 
# stats
ico_lm<-lmer(call_period ~ series_type * series_seq_cat + (1|frogID), data = calls_within, REML = TRUE)
summary(ico_lm)
emm_series_ico <- emmeans(ico_lm, pairwise ~ series_seq_cat | series_type,
                          adjust = "bonferroni")
summary(emm_series_ico)
r2(ico_lm)
#
#
# panel b
ggplot(calls_within, aes(x = series_type, y = call_type_dur, color = series_seq_cat)) +
  geom_sina(alpha = 0.7, size = 4,
            position = position_dodge(width = 0.7)) +
  stat_summary(aes(group = series_seq_cat),
               fun = mean, geom = "crossbar", width = 0.3,
               color = "black", linewidth = 0.2,
               show.legend = FALSE,
               position = position_dodge(width = 0.7)) +
  xlab("series length") +
  ylab("call duration (s)") +
  scale_color_manual(name = "  within\n  series\n position",
                     values = c("1" = "#e4cee6", "50" = "#cc99cc"),
                     labels = c("first", "mid")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.3),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
ggsave("./manuscript/figures/fig4_b_within_series_dur_type.pdf", width=6, height=4, dpi = 300, limitsize = FALSE)
ggsave("./manuscript/figures/fig4_b_within_series_dur_type.svg", width=6, height=4, dpi = 300, limitsize = FALSE) 
# stats
calldur_lm<-lmer(call_type_dur ~ series_type * series_seq_cat + (1|frogID), data = calls_within, REML = TRUE)
summary(calldur_lm)
emm_series_pw <- emmeans(calldur_lm, pairwise ~ series_seq_cat | series_type,
                         adjust = "bonferroni")
summary(emm_series_pw)
r2(calldur_lm)
icc(calldur_lm)
#
#
# panel c
ggplot(calls_within, aes(x = series_type, y = ICI, color = series_seq_cat)) +
  geom_sina(alpha = 0.7, size = 4,
            position = position_dodge(width = 0.7)) +
  stat_summary(aes(group = series_seq_cat),
               fun = mean, geom = "crossbar", width = 0.3,
               color = "black", linewidth = 0.2,
               show.legend = FALSE,
               position = position_dodge(width = 0.7)) +
  xlab("series length") +
  ylab("inter-call interval (s)") +
  scale_color_manual(name = "  within\n  series\n position",
                     values = c("1" = "#e4cee6", "50" = "#cc99cc"),
                     labels = c("first","mid")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = .3),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14))
ggsave("./manuscript/figures/fig4_c_within_series_ici_type.pdf", width=6, height=4, dpi = 300, limitsize = FALSE)
ggsave("./manuscript/figures/fig4_c_within_series_ici_type.svg", width=6, height=4, dpi = 300, limitsize = FALSE) 
# stats
ici_lm<-lmer(ICI ~ series_type * series_seq_cat + (1|frogID), data = calls_within, REML = TRUE)
summary(ici_lm)
emm_series_ici <- emmeans(ici_lm, pairwise ~ series_seq_cat | series_type,
                         adjust = "bonferroni")
summary(emm_series_ici)
r2(ici_lm)
#
#
# note subsets
notes_data <- all_notes_final %>%
  mutate(series_seq_cat = factor(series_seq_cat, levels = c("1", "50")),
         series_type = factor(series_type, levels = c("three", "four", "five")),
         note_type_relabel = factor(note_type_relabel)) %>%
  filter(INI < 10,
         series_seq_cat %in% c("1", "50"),
         series_type %in% c("three", "four", "five"))
note1_within <- notes_data %>%
  filter(note_type_relabel =="1") %>%
  filter(INI<10) %>%
  mutate(series_seq_cat = factor(series_seq_cat, levels = c("1", "50")),
         series_type = factor(series_type, levels = c("three", "four", "five")))
note2_within <- notes_data %>%
  filter(note_type_relabel =="2") %>%
  filter(INI<10) %>%
  mutate(series_seq_cat = factor(series_seq_cat, levels = c("1", "50")),
         series_type = factor(series_type, levels = c("three", "four", "five")))
#
#
# panel d
note1_within %>%
  filter(INI < 10) %>%
  ggplot(aes(x = series_type, y = note.period, color = series_seq_cat)) +
  geom_sina(alpha = 0.7, size = 4,
            position = position_dodge(width = 0.7)) +
  stat_summary(aes(group = series_seq_cat),
               fun = mean, geom = "crossbar", width = 0.3,
               color = "black", linewidth = 0.2,
               show.legend = FALSE,
               position = position_dodge(width = 0.7)) +
  scale_y_continuous(limits = c(0, 0.1),
                     labels = scales::number_format(accuracy = 0.01)) +
  xlab("series length") +
  ylab("inter-note onset (s)") +
  scale_color_manual(name = "within series sequence",
                     values = c("1" = "#bcdab7", "50" = "#336633"),
                     labels = c("1st","mid")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = .3),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14))
ggsave("./manuscript/figures/fig4_d_within_series_ino_type.pdf", width=5, height=4, dpi = 300, limitsize = FALSE)
ggsave("./manuscript/figures/fig4_d_within_series_ino_type.svg", width=5, height=4, dpi = 300, limitsize = FALSE)
# stats
note1ino_lm<-lmer(note.period ~ series_type * series_seq_cat + (1|frogID), data = note1_within, REML = TRUE)
summary(note1ino_lm)
emm_series_note1ino <- emmeans(note1ino_lm, pairwise ~ series_seq_cat | series_type,
                               adjust = "bonferroni")
summary(emm_series_note1ino)
r2(note1ino_lm)
#
#
# panel e
base_theme <- theme(panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black", size = 0.3),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  legend.position = "none",
  strip.text = element_text(size = 14))
# set y-axis scale
y_scale <- scale_y_continuous(limits = c(0, 0.1),
                              labels = scales::number_format(accuracy = 0.01))

p1 <- notes_data %>%
  filter(note_type_relabel == "1") %>%
  ggplot(aes(x = series_type, y = note.dur.s, color = series_seq_cat)) +
  geom_sina(alpha = 0.7, size = 4,
            position = position_dodge(width = 0.7)) +
  stat_summary(aes(group = series_seq_cat),
               fun = mean, geom = "crossbar", width = 0.3,
               color = "black", linewidth = 0.2,
               show.legend = FALSE,
               position = position_dodge(width = 0.7)) +
  xlab("series length") +
  ylab("note duration (s)") +
  y_scale +
  scale_color_manual(name = "within series sequence",
                     values = c("1" = "#bcdab7", "50" = "#336633"),
                     labels = c("1st", "mid")) +
  base_theme
p2 <- notes_data %>%
  filter(note_type_relabel == "2") %>%
  ggplot(aes(x = series_type, y = note.dur.s, color = series_seq_cat)) +
  geom_sina(alpha = 0.7, size = 4,
            position = position_dodge(width = 0.7)) +
  stat_summary(aes(group = series_seq_cat),
               fun = mean, geom = "crossbar", width = 0.3,
               color = "black", linewidth = 0.2,
               show.legend = FALSE,
               position = position_dodge(width = 0.7)) +
  xlab("series length") +
  ylab(NULL) +
  y_scale +
  scale_color_manual(name = "within series sequence",
                     values = c("1" = "#bcdab7", "50" = "#336633"),
                     labels = c("1st", "mid")) +
  base_theme +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
legend <- get_legend(
  p1 + theme(legend.position = "top",
             legend.text = element_text(size = 12),
             legend.title = element_text(size = 14)))
plot_row <- plot_grid(p1, p2, ncol = 2, align = "h", axis = "bt",
                      rel_widths = c(1, 0.85))
plot_grid(legend, plot_row, ncol = 1, rel_heights = c(0.1, 0.9))
ggsave("./manuscript/figures/fig4_e_within_series_note_dur_type.pdf", width=7, height=4, dpi = 300, limitsize = FALSE)
ggsave("./manuscript/figures/fig4_e_within_series_note_dur_type.svg", width=7, height=4, dpi = 300, limitsize = FALSE) 
# stats
note1dur_lm<-lmer(note.dur.s ~ series_type * series_seq_cat + (1|frogID), data = note1_within, REML = TRUE)
summary(note1dur_lm)
emm_series_note1dur <- emmeans(note1dur_lm, pairwise ~ series_seq_cat | series_type,
                          adjust = "bonferroni")
summary(emm_series_note1dur)
r2(note1dur_lm)
note2dur_lm<-lmer(note.dur.s ~ series_type * series_seq_cat + (1|frogID), data = note2_within, REML = TRUE)
summary(note2dur_lm)
emm_series_note2dur <- emmeans(note2dur_lm, pairwise ~ series_seq_cat | series_type,
                               adjust = "bonferroni")
summary(emm_series_note2dur)
r2(note2dur_lm)
#
#
# panel f
note1_within %>%
  filter(INI < 10) %>%
  ggplot(aes(x = series_type, y = INI, color = series_seq_cat)) +
  geom_sina(alpha = 0.7, size = 4,
            position = position_dodge(width = 0.7)) +
  stat_summary(aes(group = series_seq_cat),
               fun = mean, geom = "crossbar", width = 0.3,
               color = "black", linewidth = 0.2,
               show.legend = FALSE,
               position = position_dodge(width = 0.7)) +
  scale_y_continuous(limits = c(0, 0.1),
                     labels = scales::number_format(accuracy = 0.01)) +
  xlab("series length") +
  ylab("inter-note interval (s)") +
  scale_color_manual(name = "within series sequence",
                     values = c("1" = "#bcdab7", "50" = "#336633"),
                     labels = c("1st","mid")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = .3),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14))
ggsave("./manuscript/figures/fig4_f_within_series_ini_type.pdf", width=5, height=4, dpi = 300, limitsize = FALSE)
ggsave("./manuscript/figures/fig4_f_within_series_ini_type.svg", width=5, height=4, dpi = 300, limitsize = FALSE)
# stats
note1ini_lm<-lmer(INI ~ series_type * series_seq_cat + (1|frogID), data = note1_within, REML = TRUE)
summary(note1ini_lm)
emm_series_note1ini <- emmeans(note1ini_lm, pairwise ~ series_seq_cat | series_type,
                               adjust = "bonferroni")
summary(emm_series_note1ini)
r2(note1ini_lm)
#
#### Figure 5: new panels ####
#
## data frames:
# series_cvb_label
# calls_cvb_label
# notes_cvb_label
#
series_cvb_fig <- series_cvb_label %>%
  mutate(timescale="series")
calls_cvb_fig <- calls_cvb_label %>%
  mutate(timescale="calls")
notes_cvb_fig <- notes_cvb_label %>%
  mutate(timescale="notes")

notes_filtered <- notes_cvb_fig %>%
  filter((feature %in% c("notesdur", "notesmeandom")) | 
           (feature %in% c("notesini", "notesino") & note_type_relabel == "1"))
combined_cvb <- bind_rows(series_cvb_fig %>% 
                            mutate(note_type_relabel = NA_character_),
                          calls_cvb_fig  %>% 
                            mutate(note_type_relabel = NA_character_),
                          notes_filtered)
combined_cvb <- combined_cvb %>%
  mutate(feature_label = case_when(
    feature == "seriesdur"~ "duration (s)",
    feature == "seriesisi"~ "ISI",
    feature == "seriesiso"~ "ISO",
    feature == "callsdur" ~ "duration (c)",
    feature == "callsici" ~ "ICI",
    feature == "callsioi" ~ "ICO",
    feature == "notesdur" & note_type_relabel == "1"~ "dur N1",
    feature == "notesdur" & note_type_relabel == "2"~ "dur N2",
    feature == "notesini" & note_type_relabel == "1"~ "INI",
    feature == "notesino" & note_type_relabel == "1"~ "INO",
    feature == "notesmeandom" & note_type_relabel == "1"~ "DF N1",
    feature == "notesmeandom" & note_type_relabel == "2"~ "DF N2"))
feature_label_levels <- c("duration (s)", "ISI", "ISO",
                          "duration (c)", "ICI", "ICO",
                          "dur N1", "dur N2", "INI", "INO", "DF N1", "DF N2")
combined_cvb <- combined_cvb %>%
  mutate(feature_label = factor(feature_label, levels = feature_label_levels))
label_colors <- c("duration (s)" = "#2a5674",
                  "ISI"          = "#2a5674",
                  "ISO"          = "#2a5674",
                  "duration (c)" = "#9f839d",
                  "ICI"          = "#9f839d",
                  "ICO"          = "#9f839d",
                  "dur N1"       = "#cde4cc",
                  "dur N2"       = "#99cc99",
                  "INI"          = "#839c7e",
                  "INO"          = "#839c7e",
                  "DF N1"        = "#cde4cc",
                  "DF N2"        = "#99cc99")
label_shapes <- c("duration (s)" = 16,
                  "ISI"          = 17,
                  "ISO"          = 15,
                  "duration (c)" = 16,
                  "ICI"          = 17,
                  "ICO"          = 15,
                  "dur N1"       = 16,
                  "dur N2"       = 16,
                  "INI"          = 17,
                  "INO"          = 15,
                  "DF N1"        = 8,
                  "DF N2"        = 8)

timescale_x <- combined_cvb %>%
  mutate(timescale = case_when(
    feature_label %in% c("duration (s)", "ISI", "ISO")         ~ "series",
    feature_label %in% c("duration (c)", "ICI", "ICO")         ~ "calls",
    feature_label %in% c("dur N1","dur N2","INI","INO","DF N1","DF N2") ~ "notes")) %>%
  group_by(timescale) %>%
  summarise(x = mean(grand_CVw, na.rm = TRUE), .groups = "drop")

#
## panel a: CVw x ratio between series within bout
#
p_a <- ggplot(combined_cvb,
              aes(x = grand_CVw, y = CV_ratio,
                  color = feature_label, shape = feature_label)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "gray60", linewidth = 0.4) +
  geom_point(size = 8, alpha = 1) +
  scale_color_manual(values = label_colors) +
  scale_shape_manual(values = label_shapes) +
  xlim(0, 75) +
  ylim(0,5) +
  xlab("CVw (%)") +
  ylab("CVb / CVw") +
  coord_cartesian(clip = "off") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line        = element_line(colour = "black", size = 0.3),
        legend.position  = "none",
        axis.title.y     = element_text(size = 26),
        axis.title.x     = element_text(size = 26),
        axis.text.x      = element_text(size = 20),
        axis.text.y      = element_text(size = 20),
        plot.margin      = margin(t = 10, r = 20, b = 10, l = 10))
p_a
ggsave("./manuscript/figures/fig5_a_CVw_ratio.pdf",
       width = 6, height = 5, dpi = 300, limitsize = FALSE)
ggsave("./manuscript/figures/fig5_a_CVw_ratio.svg",
       width = 8, height = 5, dpi = 300, limitsize = FALSE)
#
#
## panel b: CVw x ratio within series
calls_cvb_fig <- calls_cvb_series_label %>%
  mutate(timescale="calls")
notes_cvb_fig <- notes_cvb_label_series %>%
  mutate(timescale="notes")
notes_filtered <- notes_cvb_label_series %>%
  filter((feature %in% c("notesdur", "notesmeandom")) | 
           (feature %in% c("notesini", "notesino") & note_type_relabel == "1")) %>%
  mutate(timescale="notes")
combined_cvb_ws <- bind_rows(calls_cvb_series_label %>%
                               mutate(note_type_relabel = NA_character_),
                             notes_filtered)
combined_cvb_ws <- combined_cvb_ws %>%
  mutate(feature_label = case_when(
    feature == "callsdur" ~ "duration (c)",
    feature == "callsici" ~ "ICI",
    feature == "callsioi" ~ "ICO",
    feature == "notesdur" & note_type_relabel == "1" ~ "dur N1",
    feature == "notesdur" & note_type_relabel == "2" ~ "dur N2",
    feature == "notesini" & note_type_relabel == "1" ~ "INI",
    feature == "notesino" & note_type_relabel == "1" ~ "INO",
    feature == "notesmeandom" & note_type_relabel == "1" ~ "DF N1",
    feature == "notesmeandom" & note_type_relabel == "2" ~ "DF N2"))
feature_label_levels_ws <- c("duration (c)", "ICI", "ICO",
                             "dur N1", "dur N2", "INI", "INO", "DF N1", "DF N2")
combined_cvb_ws <- combined_cvb_ws %>%
  mutate(feature_label = factor(feature_label, levels = feature_label_levels_ws))
label_colors_ws <- c("duration (c)" = "#9f839d",
                     "ICI"          = "#9f839d",
                     "ICO"          = "#9f839d",
                     "dur N1"       = "#cde4cc",
                     "dur N2"       = "#99cc99",
                     "INI"          = "#839c7e",
                     "INO"          = "#839c7e",
                     "DF N1"        = "#cde4cc",
                     "DF N2"        = "#99cc99")
label_shapes_ws <- c("duration (c)" = 16,
                     "ICI"          = 17,
                     "ICO"          = 15,
                     "dur N1"       = 16,
                     "dur N2"       = 16,
                     "INI"          = 17,
                     "INO"          = 15,
                     "DF N1"        = 8,
                     "DF N2"        = 8)

p_ws <- ggplot(combined_cvb_ws,
               aes(x = grand_CVw, y = CV_ratio,
                   color = feature_label, shape = feature_label)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "gray60", linewidth = 0.4) +
  geom_point(size = 8, alpha = 1) +
  scale_color_manual(values = label_colors_ws) +
  scale_shape_manual(values = label_shapes_ws) +
  xlim(0, 75) +
  ylim(0,5) +
  xlab("CVw (%)") +
  ylab("CVb / CVw") +
  coord_cartesian(clip = "off") +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line        = element_line(colour = "black", size = 0.3),
    legend.position  = "none",
    axis.title.y     = element_text(size = 26),
    axis.title.x     = element_text(size = 26),
    axis.text.x      = element_text(size = 20),
    axis.text.y      = element_text(size = 20),
    plot.margin      = margin(t = 10, r = 20, b = 10, l = 10))
p_ws
ggsave("./manuscript/figures/fig5_b_CVw_ratio_within_series.pdf", width = 6, height = 5, dpi = 300, limitsize = FALSE)
ggsave("./manuscript/figures/fig5_b_CVw_ratio_within_series.svg", width = 6, height = 5, dpi = 300, limitsize = FALSE)
#
#
# CVb by sequence order
#
# data frames:
# calls_cvb_seq_label
# notes_cvb_label_series
#
calls_cvb_fig <- calls_cvb_seq_label %>%
  mutate(timescale="calls")
notes_filtered_seq <- notes_cvb_label_series_seq %>%
  filter((feature %in% c("notesdur", "notesmeandom")) | 
           (feature %in% c("notesini", "notesino") & note_type_relabel == "1"))
combined_cvb_seq <- bind_rows(calls_cvb_seq_label %>% mutate(note_type_relabel = NA_character_),
  notes_filtered_seq)
combined_cvb_seq <- combined_cvb_seq %>%
  mutate(feature_label = case_when(
    feature == "callsdur"  & series_seq_cat == "1"  ~ "dur_c1",
    feature == "callsici"  & series_seq_cat == "1"  ~ "ICI_1",
    feature == "callsioi"  & series_seq_cat == "1"  ~ "ICO_1",
    feature == "callsdur"  & series_seq_cat == "50" ~ "dur_c50",
    feature == "callsici"  & series_seq_cat == "50" ~ "ICI_50",
    feature == "callsioi"  & series_seq_cat == "50" ~ "ICO_50",
    feature == "notesdur"  & note_type_relabel == "1" & series_seq_cat == "1"  ~ "dur_n1_1",
    feature == "notesdur"  & note_type_relabel == "2" & series_seq_cat == "1"  ~ "dur_n2_1",
    feature == "notesini"  & note_type_relabel == "1" & series_seq_cat == "1"  ~ "INI_1",
    feature == "notesino"  & note_type_relabel == "1" & series_seq_cat == "1"  ~ "INO_1",
    feature == "notesmeandom" & note_type_relabel == "1" & series_seq_cat == "1" ~ "DF_n1_1",
    feature == "notesmeandom" & note_type_relabel == "2" & series_seq_cat == "1" ~ "DF_n2_1",
    feature == "notesdur"  & note_type_relabel == "1" & series_seq_cat == "50" ~ "dur_n1_50",
    feature == "notesdur"  & note_type_relabel == "2" & series_seq_cat == "50" ~ "dur_n2_50",
    feature == "notesini"  & note_type_relabel == "1" & series_seq_cat == "50" ~ "INI_50",
    feature == "notesino"  & note_type_relabel == "1" & series_seq_cat == "50" ~ "INO_50",
    feature == "notesmeandom" & note_type_relabel == "1" & series_seq_cat == "50" ~ "DF_n1_50",
    feature == "notesmeandom" & note_type_relabel == "2" & series_seq_cat == "50" ~ "DF_n2_50"))
feature_label_levels_seq <- c(
  "dur_c1",   "ICI_1",   "ICO_1",   "dur_n1_1", "dur_n2_1", "INI_1",  "INO_1", "DF_n1_1", "DF_n2_1",
  "dur_c50",  "ICI_50",  "ICO_50",  "dur_n1_50","dur_n2_50","INI_50", "INO_50", "DF_n1_50", "DF_n2_50")
combined_cvb_seq <- combined_cvb_seq %>%
  mutate(feature_label = factor(feature_label, levels = feature_label_levels_seq))
label_colors_seq <- c("dur_c1"    = "#9f839d", "dur_c50"   = "#9f839d",
                      "ICI_1"     = "#9f839d", "ICI_50"    = "#9f839d",
                      "ICO_1"     = "#9f839d", "ICO_50"    = "#9f839d",
                      "dur_n1_1"  = "#cde4cc", "dur_n1_50" = "#cde4cc",
                      "dur_n2_1"  = "#99cc99", "dur_n2_50" = "#99cc99",
                      "INI_1"     = "#839c7e", "INI_50"    = "#839c7e",
                      "INO_1"     = "#839c7e", "INO_50"    = "#839c7e",
                      "DF_n1_1"   = "#cde4cc", "DF_n1_50"  = "#cde4cc",
                      "DF_n2_1"   = "#99cc99", "DF_n2_50"  = "#99cc99")
label_shapes_seq <- c("dur_c1"    = 16, "dur_c50"   = 16,
                      "ICI_1"     = 17, "ICI_50"    = 17,
                      "ICO_1"     = 15, "ICO_50"    = 15,
                      "dur_n1_1"  = 16, "dur_n1_50" = 16,
                      "dur_n2_1"  = 16, "dur_n2_50" = 16,
                      "INI_1"     = 17, "INI_50"    = 17,
                      "INO_1"     = 15, "INO_50"    = 15,
                      "DF_n1_1"   = 8, "DF_n1_50"  = 8,
                      "DF_n2_1"   = 8, "DF_n2_50"  = 8)

p_seq <- ggplot(combined_cvb_seq,
                aes(x = grand_CVw, y = CV_ratio,
                    color = feature_label, shape = feature_label)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "gray60", linewidth = 0.4) +
  geom_point(size = 8, alpha = 1) +
  scale_color_manual(values = label_colors_seq) +
  scale_shape_manual(values = label_shapes_seq) +
  xlim(0, 75) +
  ylim(0,5) +
  xlab("CVw (%)") +
  ylab("CVb / CVw") +
  coord_cartesian(clip = "off") +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line        = element_line(colour = "black", size = 0.3),
    legend.position  = "none",
    axis.title.y     = element_text(size = 26),
    axis.title.x     = element_text(size = 26),
    axis.text.x      = element_text(size = 20),
    axis.text.y      = element_text(size = 20),
    plot.margin      = margin(t = 10, r = 20, b = 10, l = 10))
p_seq
#
#
## seq position split
#
combined_cvb_seq_1  <- combined_cvb_seq %>%
  filter(series_seq_cat == "1")
combined_cvb_seq_50 <- combined_cvb_seq %>%
  filter(series_seq_cat == "50")
p_seq_1  <- plot_seq(combined_cvb_seq_1,  show_y_axis = TRUE,  bg_color = "#ccffff") 
p_seq_50 <- plot_seq(combined_cvb_seq_50, show_y_axis = FALSE, bg_color = "#99cccc") 
p_seq_combined <- p_seq_1 + p_seq_50 +
  plot_layout(widths = c(1, 1)) &
  theme(plot.margin = margin(t = 10, r = 5, b = 10, l = 5))
p_seq_combined_final <- wrap_elements(full = {
  gt <- patchwork::patchworkGrob(p_seq_combined)
  arrangeGrob(gt, ncol = 1)
})

p_divider <- ggplot() +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             color = "gray60", linewidth = 0.6) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))
p_seq_final <- p_seq_1 + p_divider + p_seq_50 +
  plot_layout(widths = c(1, 0.03, 1))
p_seq_final
ggsave("./manuscript/figures/fig5_b_seq_split.pdf",
       p_seq_final, width = 12, height = 5, dpi = 300, limitsize = FALSE)
ggsave("./manuscript/figures/fig5_b_seq_split.svg",
       p_seq_final, width = 12, height = 5, dpi = 300, limitsize = FALSE)
#
#
## combined panels with seq position split
#
p_seq_1_combined <- plot_seq(combined_cvb_seq_1,  show_y_axis = TRUE,  bg_color = "#ccffff") 
p_seq_50_combined <- plot_seq(combined_cvb_seq_50, show_y_axis = FALSE, bg_color = "#99cccc")
p_c_in_fig <- p_seq_1_combined + p_divider + p_seq_50_combined +
  plot_layout(widths = c(1, 0.03, 1))

p_b_combined <- p_ws + theme(axis.title.y = element_blank())

fig5 <- (p_a | p_b_combined | p_c_in_fig) +
  plot_layout(widths = c(1, 1, 2))
fig5

ggsave("./manuscript/figures/fig5_combined_seq_split.pdf",
       fig5, width = 24, height = 6, dpi = 300, limitsize = FALSE)
ggsave("./manuscript/figures/fig5_combined_seq_split.svg",
       fig5, width = 24, height = 6, dpi = 300, limitsize = FALSE)
#
#### Figure 6: note features ####
#
## dur x INI
#
ggplot(all_notes_final, aes(x=note.dur.s, y=INI, color=note_type_relabel)) +
  geom_point(size=3.5, alpha=0.7) +
  xlab("duration (s)") +
  ylab("ini (s)") + 
  labs(color=" note\n type") +
  ylim(0,5) +
  xlim(0,0.15) +
  scale_color_manual(values = NoteColors,
                     labels = c("uc" = "solo")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = .3),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.position = "none")
ggsave("./manuscript/figures/fig6_b_note_dur_ino.pdf", width=6, height=4, dpi = 300) 
ggsave("./manuscript/figures/fig6_b_note_dur_ino.svg", width=6, height=4, dpi = 300)