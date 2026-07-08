################## B. pulchellus Timescales of Call Variability #################
#### author: M. Rodriguez-Santiago
#### sound analysis software: Raven
####
####
#### this script imports and combines individual .txt files exported from Raven
####
#### variables:
##   Delta Time = note duration (from Raven, waveform selection)
##   INI = inter-note interval (computed per individual)
##   note.period = inter-note onset (computed per individual)
####
#### input: automated detector .txt files from Raven
#### output: all_df_notes.csv with note feature calculations
####
####

library(dplyr)
library(ggplot2)
library(data.table)
library(readr)


## Read all note text files
list_of_files <- list.files(path = "./data/raw_data", recursive = TRUE,
                            pattern = "\\.txt$",
                            full.names = TRUE)
# Read all the files into a data table and create a FileName column to store filenames
DT <- rbindlist(sapply(list_of_files, fread, simplify = FALSE),
                use.names = TRUE, idcol = "FileName", fill=TRUE)


#### add frog ID column,etc. & subset by waveform and spectrogram ####
#
#frog 1
frog1_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_1_two_call_bouts.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01281",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 2
frog2_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_2_two_call_bouts.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01303",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 3
frog3_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_3_six_call_bouts.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01314",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 4
frog4_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_4_two_call_bouts.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01325",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 5
frog5_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_5_four_call_bouts.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01346",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 6
frog6_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_6_two_call_bouts.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01552",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 7
frog7_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_7_three_call_bouts.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01573",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 8
frog8_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_8_three_call_bouts.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01594",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog9
frog9_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_9_two_call_bouts.Amplitude.Detector.selections.txt",
         View == "Waveform 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01605",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 10
frog10_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_10_four_call_bouts.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01626",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 11
frog11_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_11_one_call_bout.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01682",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 12
frog12_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_12_five_call_bouts.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01741",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 13
frog13_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_13_seven_call_bouts.Amplitude.Detector.selections.txt",
         View == "Waveform 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01752",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())

#frog 14
frog14_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_14_two_call_bouts.Amplitude.Detector.selections.txt",
         View == "Waveform 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01763",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 15
frog15_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_15_four_call_bouts.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01774",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 16
frog16_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_16_three_call_bouts.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01785",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 17
frog17_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_17_three_call_bouts.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01806",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())
#frog 18
frog18_s <- DT %>%
  filter(FileName == "./amplitude_detector_tables/male_18_three_call_bouts.Amplitude.Detector.selections.txt",
         View == "Spectrogram 1") %>%
  rename("begin.time.s" = "Begin Time (s)",
         "end.time.s" = "End Time (s)",
         "low.freq.hz"= "Low Freq (Hz)",
         "high.freq.hz" = "High Freq (Hz)",
         "note.dur.s" = "Delta Time (s)",
         "max.freq" = "Max Freq (Hz)",
         "peak.amp" = "Peak Amp (U)",
         "freq.5" = "Freq 5% (Hz)",
         "freq.25" = "Freq 25% (Hz)",
         "freq.75" = "Freq 75% (Hz)",
         "freq.95" = "Freq 95% (Hz)") %>%
  arrange(begin.time.s) %>%
  mutate(frogID = "01817",
         INI = lead(begin.time.s)-end.time.s,
         note.period = lead(begin.time.s)-begin.time.s,
         total_note_num = row_number())

#### COMBINE ALL DATAFRAMES INTO CSV -------------------------------------------

all_male_s <- rbind(frog1_s, frog2_s, frog3_s, frog4_s, frog5_s, frog6_s, frog7_s, frog8_s, frog9_s, frog10_s, frog11_s, frog12_s, frog13_s, frog14_s, frog15_s, frog16_s, frog17_s, frog18_s)

# write.csv(all_male_s, "all_df_notes.csv", row.names = FALSE)