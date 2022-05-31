library(dplyr)
library(tuneR)

# Each file is assumed to contain a left and right track with a unique
# sample rate.
read_sounds <- function(files, time_expansion = TRUE) {
  require(tuneR)
  sounds <- list()
  for (file in files) {
    sndObj <- readWave(file)
    sounds[[length(sounds) + 1]] <- sndObj@left
    sounds[[length(sounds) + 1]] <- sndObj@right
    samp.rate <- sndObj@samp.rate * ifelse(time_expansion, 10, 1)
    bit <- sndObj@bit
  }
  list(sounds = sounds, samp.rate = samp.rate, bit = bit)
}

FFT_filter <- function(snd, samp.rate, freq_min = 0, freq_max = Inf) {
  require(tuneR)
  # Fast Fourier Transform algorithm to get the frequencies
  p <- fft(snd)
  # select just the first half since the second half is a mirror image
  # of the first
  n <- length(snd)
  nUniquePts <- ceiling((n + 1) / 2)
  # create the frequencies array
  freq <- (0:(nUniquePts - 1)) * (samp.rate / n)
  p[freq < freq_min | freq_max < freq] <- 0
  p[n:(nUniquePts + 1)] <- Conj(p[2:(n - nUniquePts + 1)])
  Re(fft(p, inverse = TRUE) / n)
}

padding_snd <- function(snd, samp.rate, lag) {
  if (lag > 0)
    return(tail(c(snd, rep(0, lag * samp.rate), length(snd))))
  else
    return(head(c(rep(0, -lag * samp.rate), snd), length(snd)))
}

# Function to compute lags and align sounds.
# - records: must contain a list `sounds` (a list of sound vectors)
#   and the sample rate `samp.rate`.
# - syncs: manual approximate times of sync signals, must be precise
#   up to sync_length.
# - sync_length: sync signal must appear during at most sync_length.
align_sounds <- function(records, syncs, sync_length,
                         signal_freq = NULL) {
  syncs <- syncs - sync_length
  sync_length <- 3 * sync_length
  for (i in seq_along(records$sounds)) {
    lag <- compute_lag(records$sounds[[1]], records$sounds[[i]],
                       records$samp.rate, syncs[1], syncs[i],
                       sync_length, signal_freq = signal_freq)
    aligned_sound <- padding_snd(records$sounds[[i]],
                                 records$samp.rate, lag$lag)
    records$sounds[[i]] <- aligned_sound
  }
  records
}

# Filter lower frequencies (takes a minute)
clean_freq <- function(records, freq_min = 0, freq_max = Inf) {
  pb <- txtProgressBar(min = 0, max = length(records$sounds),
                       style = 3)
  for (i in seq_along(records$sounds)) {
    setTxtProgressBar(pb, i - 1)
    records$sounds[[i]] <- FFT_filter(records$sounds[[i]],
                                      records$samp.rate,
                                      freq_min, freq_max)
  }
  setTxtProgressBar(pb, length(records$sounds))
  close(pb)
  records
}

# Find calls.
# - records: must contain a list `sounds` (a list of sound vectors)
#   and the sample rate `samp.rate`.
# - duration: duration of a call.
# - sensitivity: sensitivity of the method (start with large value and
#   decrease until enough call are found)
find_calls <- function(records, duration, sensitivity) {
  require(dplyr)
  samp.rate <- records$samp.rate
  calls <- list()
  pb <- txtProgressBar(min = 0, max = length(records$sounds),
                       style = 3)
  for (i in seq_along(records$sounds)) {
    setTxtProgressBar(pb, i - 1)
    snd <- abs(records$sounds[[i]])
    cs <- cumsum(snd)
    index <- head(seq_along(snd), - duration * samp.rate)
    # Compute the average signal over a duration of duration_min
    snd_avg <- (cs[index + duration * samp.rate] - cs[index]) /
      (duration * samp.rate)
    upper <- rep(TRUE, length(snd_avg))
    threshold <- sensitivity * median(snd)
    upper <- snd_avg > threshold
    while (any(upper)) {
      # TODO optimize this loop
      # Find max power
      power <- max(snd_avg[upper])
      j <- which(snd_avg == power)[1]
      j_min <- max(0, j - duration * samp.rate)
      j_max <- min(length(upper), j + duration * samp.rate)
      upper[j_min:j_max] <- FALSE
      # Center time of max power (because moving average is shifted)
      calls[[length(calls) + 1]] <- list(recorder = names(records$sounds)[i],
                                         time = j / samp.rate,
                                         power = power,
                                         max_signal = max(snd[j_min:j_max]))
    }
  }
  setTxtProgressBar(pb, length(records$sounds))
  close(pb)
  bind_rows(calls)
}

# Compute exact lag between two signals using cross-correlation
# coefficient on resampled signals with spline method.
# - snd1,snd2: signals
# - start1, start2: starting time
# - duration of the signals
# Failed idea: smoothed absolute signal to prevent the ccf from
# optimizing the lag with best phase, rather than the lag with best
# start and maximum power. Rolling mean on 4 values because frequency
# is around 45 kHz and 384 / 45 is around 8 (since we take the
# absolute, the signal repeat itself about every 4 samples)
#part1 <- rollmax(abs(snd1_s$y), k = 4 * RESAMPLING, na.pad = TRUE)
#TODO compute min and max lags with a given confidence level
#TODO optimize this function by iterating on resampling (with
#aggregated signals) to determine an efficient lag.max
compute_lag <- function(snd1, snd2, samp.rate, start1, start2,
                        duration, resampling = 4, signal_freq = NULL) {
  i_length <- duration * samp.rate
  i1_s <- max(start1 * samp.rate, 0)
  i1_f <- min(i1_s + i_length, length(snd1))
  i2_s <- max(start2 * samp.rate, 0)
  i2_f <- min(i2_s + i_length, length(snd2))
  stopifnot(i1_s <= i1_f && i2_s <= i2_f)
  snd1_s <- spline(i1_s:i1_f, snd1[i1_s:i1_f],
                   resampling * (i1_f - i1_s + 1))
  snd2_s <- spline(i2_s:i2_f, snd2[i2_s:i2_f],
                   resampling * (i2_f - i2_s + 1))
  if (!is.null(signal_freq)) {
    snd1_s$y <- FFT_filter(snd1_s$y, resampling * samp.rate,
                           signal_freq / 2, signal_freq * 2)
    snd2_s$y <- FFT_filter(snd2_s$y, resampling * samp.rate,
                           signal_freq / 2, signal_freq * 2)
  }
  cc <- ccf(snd1_s$y, snd2_s$y, plot = FALSE,
            lag.max = resampling * i_length, na.action = na.pass)
  lag <- cc$lag[which.max(cc$acf)] / (resampling * samp.rate)
  list(lag = start2 - start1 - lag, acf = max(cc$acf))
}

# Keep only one lag per recorder given all possible pairwise lags and
# confidence value (ACF) during a previously defined time window
best_TDOA <- function(TDOA) {
  require(dplyr)
  #TODO replace with igraph.mst algorithm
  TDOA_symmetric <- TDOA %>%
    mutate(rec = recorder2, recorder2 = recorder1, recorder1 = rec,
           TDOA = -TDOA) %>%
    select(-rec) %>%
    rbind(TDOA)
  most_central <- TDOA_symmetric %>%
    group_by(recorder1) %>%
    summarise(sum_acf = sum(acf)) %>%
    arrange(desc(sum_acf)) %>%
    slice(1) %>%
    select(recorder1) %>%
    unlist()
  TDOA_symmetric %>%
    filter(recorder1 == most_central) %>%
    mutate(recorder = recorder2) %>%
    select(-recorder1, -recorder2) %>%
    rbind(list(recorder = most_central, TDOA = 0, acf = 1))
}

# Extract time differences of arrival by matching calls one all
# recorders:
# - calls: must contain the columns "recorder", "time", "power" and
#   "max_signal".
# - records: must contain a list `sounds` (a list of sound vectors)
#   and the sample rate `samp.rate`.
# Return the columns "id" (time of main call), "acf" (confidence),
# "recorder", "power", "max_signal" and "TDOA".
extract_TDOA <- function(calls, records, max_delay, duration) {
  require(dplyr)
  TDOA <- list()
  calls <- calls %>% arrange(desc(power))
  recorders <- unique(calls$recorder)
  pb.max <- nrow(calls)
  pb <- txtProgressBar(min = 0, max = pb.max, style = 3)
  # Start with calls with maximum energy
  while (nrow(calls) != 0) {
    setTxtProgressBar(pb, pb.max - nrow(calls))
    most_power <- calls %>% slice(1)
    # get the best call for each recorder in the admissible time window
    calls_close <- calls %>%
      filter(abs(time - most_power$time) <= max_delay) %>%
      group_by(recorder) %>%
      slice(1) %>%
      ungroup()
    if (nrow(calls_close) != 1) {
      # Add each missing recorder for which a trace signal could be
      # found by chance
      #TODO only keep time with max power for each missing recorder
      missing_recorders <- recorders[!recorders %in% calls_close$recorder]
      if (length(missing_recorders) >= 1)
        calls_close <- crossing(recorder = missing_recorders,
                                time = unique(calls_close$time)) %>%
          mutate(time_max = sapply(records$sounds[recorder], length) /
                   records$samp.rate) %>%
          filter(time + duration < time_max) %>%
          select(-time_max) %>%
          mutate(power = 0, max_signal = 0) %>%
          rbind(calls_close)
      # If there are other calls, compute the best possible lag for
      # each pair of calls
      TDOA_local <- list()
      for (i in 1:(nrow(calls_close) - 1))
        for (j in (i + 1):nrow(calls_close)) {
          if (calls_close[j,"power"] == 0)
            next
          recorder1 <- calls_close[i,"recorder"] %>% unlist()
          recorder2 <- calls_close[j,"recorder"] %>% unlist()
          lag <- compute_lag(records$sounds[[recorder1]],
                             records$sounds[[recorder2]],
                             records$samp.rate,
                             calls_close[i,"time"] %>% unlist(),
                             calls_close[j,"time"] %>% unlist(),
                             duration)
          TDOA_local[[length(TDOA_local) + 1]] <- list(recorder1 = recorder1,
                                                       recorder2 = recorder2,
                                                       TDOA = lag$lag,
                                                       acf = lag$acf)
        }
      TDOA_local <- bind_rows(TDOA_local)
      # Remove redundancies caused by missing recorder
      TDOA_local <- TDOA_local %>%
        group_by(recorder1, recorder2) %>%
        arrange(desc(acf)) %>%
        slice(1) %>%
        ungroup()
      # Keep only one call per recorder such that the confidence
      # values are the best
      lags <- best_TDOA(TDOA_local)
      TDOA <- lags %>%
        mutate(id = most_power$time) %>%
        right_join(calls_close %>%
                   select(-time) %>%
                 distinct(), by = c("recorder")) %>%
        rbind(TDOA)
    }
    # Remove all calls from the given time window
    calls <- calls %>%
      filter(abs(time - most_power$time) > max_delay)
  }
  setTxtProgressBar(pb, pb.max)
  close(pb)
  TDOA
}
