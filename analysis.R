source("scripts/sounds.R")
source("scripts/geolocation.R")
source("scripts/plotting.R")

print("### Reading sound files ###")
raws <- read_sounds(c(file1, file2, file3), time_expansion = TRUE)
raws$bit <- 12 # only 12 bits out of 16 are used
names(raws$sounds) <- c("masterL", "masterR", "slave1L", "slave1R",
                        "slave2L", "slave2R")

print("### Aligning sounds ###")
print("# filtering synchronization signals")
SYNC_SIGNAL_FREQ <- 8e3
records_filtered <- clean_freq(raws,
                               freq_min = SYNC_SIGNAL_FREQ / 2,
                               freq_max = SYNC_SIGNAL_FREQ * 2)
print("# finding synchronization signals")
SYNC_SIGNAL_LENGTH <- 6.5e-4
# sensitivity is best between 40 and 110 (10 to 40 is OK but slower,
# above 110 no longer works)
SYNC_SIGNALS_SENSITIVITY <- 20
sync_signals <- find_calls(records_filtered,
                           duration = SYNC_SIGNAL_LENGTH,
                           sensitivity = SYNC_SIGNALS_SENSITIVITY)
NB_SECONDS_BETWEEN_SIGNALS <- 3
sync_signals.nb <- ceiling(length(raws$sounds[[1]]) / raws$samp.rate /
                           NB_SECONDS_BETWEEN_SIGNALS) * length(raws$sounds)
sync_signals <- sync_signals %>%
  arrange(desc(power)) %>%
  head(sync_signals.nb)
## Using a ratio of 2 is enough to differentiate with
## non-synchronization signals
#  filter(power > max(power) / 2 & max_signal > max(max_signal) / 2)
first_signals <- sync_signals %>%
  group_by(recorder) %>%
  arrange(time) %>%
  slice(1)
stopifnot(all(first_signals$recorder == names(raws$sounds)))
print(paste("# lags are", toString(format(first_signals$time))))
records <- align_sounds(raws, syncs = first_signals$time,
                        sync_length = SYNC_SIGNAL_LENGTH,
                        signal_freq = SYNC_SIGNAL_FREQ)

print("### Filter unrelated frequencies (takes a minute) ###")
records_filtered <- clean_freq(records, freq_min = CALL_FREQ_MIN,
                               freq_max = CALL_FREQ_MAX)

print("### Find calls (takes a minute) ###")
# No need to have a too small sensitivity as it produces too much
# low-power calls that will anyway be discarded afterwards when
# computing the TDOAs (around 10 calls per second seems sensible)
calls <- find_calls(records_filtered,
                    duration = CALL_DURATION,
                    sensitivity = CALL_SENSITIVITY)
calls.per_second <- nrow(calls) /
  (length(records_filtered$sounds[[1]]) / records_filtered$samp.rate) /
  length(records_filtered$sounds)
print(paste("# there are", nrow(calls), "detected calls, which makes",
            format(calls.per_second), "calls per second on each recorder"))

# Compute maximum delay between each recorder
sound_speed_ms <- 331.3 + 0.606 * TEMPERATURE_C
library(readr)
recorders <- tibble(id = names(records$sounds),
                    x = c(3, 3, 0, 0, 0, 0),
                    y = c(0, 0, 0, 0, 3, 3),
                    z = c(2.20, 4.15, 1.31, 2.38, 1.56, 2.62))
distance <- position2distance(recorders)
max_delay <- max(distance$distance) / sound_speed_ms

print("### Compute exact TDOA from recordings (take a long time) ###")
TDOA <- extract_TDOA(calls, records_filtered, max_delay,
                     CALL_DURATION)
print(paste("# there are", nrow(TDOA), "TDOAs in total"))

# Clean data based on ACF to save time when computing positions
library(dplyr)
MIN_SIMULTANEOUS_OBS <- 5 # should be 5 to get meaningful positions
TDOA_selection <- TDOA %>%
  filter(acf >= MIN_ACF) %>%
  group_by(id) %>%
  filter(n() >= MIN_SIMULTANEOUS_OBS) %>%
  ungroup()
print(paste("# there are", nrow(TDOA_selection), "selected TDOAs"))

print("### Compute position of each call (take a few minutes) ###")
pos_all <- TDOA2position(TDOA_selection, recorders)
print(paste("# there are", nrow(pos_all), "positions in total"))

# Clean data to remove inconsistent data
speeds <- speed_locations(pos_all)
positions <- clean_locations(pos_all, recorders, speeds,
                             MAX_DISTANCE_M, MAX_SPEED_KMH)
print(paste("# there are", nrow(positions), "selected positions"))

# Plot the trajectory
MIN_DISTANCE_M <- 10
MIN_SPEED_KMH <- 10
plot_trajectory(positions, recorders, speeds, MIN_DISTANCE_M,
                MIN_SPEED_KMH, filename = OUTPUT_HTML)
