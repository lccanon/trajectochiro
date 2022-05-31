## Parameter spectific to the specie
CALL_FREQ_MIN <- 40e3
CALL_FREQ_MAX <- Inf
# Calls last 6 ms for PIPPIP and 12 ms for NYCNOC
CALL_DURATION <- 0.006

## Parameters specific to the time of the recording
CALL_SENSITIVITY <- 10
TEMPERATURE_C <- 15

## General parameters
MIN_ACF <- 0.5
# Bat positions cannot be identified at more than 20 meters
MAX_DISTANCE_M <- 20
# Bats do not travel/hunt at more than 50 km/s
MAX_SPEED_KMH <- 50
