
################## FP analysis for NPM recordings ################################################
# Written for demo purposes by Arthur Godino
# From raw fluorescence time series files created by running the v3_dualcolor_FP Bonsai script on NPM FP3002 rig (read and organized in the first part of this script)
# Analysis method (second part of this script) adapted from Martianova et al, J Vis Exp, 2019 (DOI: 10.3791/60278) with minor modifications ###########################

require(zoo)
require(MASS)
require(MESS)
require(signal)
require(tidyverse)  # /!\ there are name conflicts between signal and tidyverse packages - for select and filter functions, always call package::function

################# INPUTS #########################################################################

wd.rig = '/Users/arthur/Desktop/rawNPM/' ; setwd(wd.rig)  # path to the raw photometry files (Bonsai output) - should be one per LED used (for this example, only use the L415 and L470)
sampf = 26         # effective sampling frequency for each channel /!\ = actual sampling frequency divided by number of LED states (default for us is 78Hz/3=26Hz)
animal = "T1-3"    # should match the subject name used in naming the Bonsai output files

################# PART 1 : DATA IMPORT, CLEANING AND TIDYING ###########################################################
# Read input files and add a column specifying which LED was ON
rig.test <- read.table(list.files(pattern = paste0("470_", animal, "_")), sep = ",", dec = ".", fill = T, header = T) %>%
            mutate(Channel = 470)  ;  head(rig.test)
rig.iso <- read.table(list.files(pattern = paste0("415_", animal)), sep = ",", dec = ".", fill = T, header = T)  %>%
            mutate(Channel = 415)  ;  head(rig.iso)

# Combine chronologically
comb <- merge(rig.iso, rig.test, all = T)  ;  head(comb)
# Delete unused red 'photodetector' zone data & create two columns, one for each channel (470 and 415) - this introduces NAs because LED states are interleaved in time
comb <- mutate(comb, Region0R = NULL) %>% 
        pivot_wider(names_from = Channel, values_from = Region1G, names_prefix = "raw.C")  ;  head(comb)
# Extrapolate raw C415 isosbestic values at C470 timestamps using linear interpolation
comb <- mutate(comb, raw.C415 = as.numeric(na.approx(zoo(raw.C415, Timestamp), na.rm = F)))  ;  head(comb)
# Keep only the timestamps where C470 was on, and the associated extrapolated C415 values
comb <- dplyr::filter(comb, is.na(raw.C470)==F & is.na(raw.C415)==F)  ;  head(comb)  # Note dplyr::filter, to avoid conflict with signal::filter
# Create a re-zeroed time column, note how the difference between each row/timestamp is 1/26 seconds, ie sampling period is 38.46 ms
comb <- mutate(comb, Relative.time = Timestamp - Timestamp[1])  ;  head(comb)

# Quick plot
demo_plot <- function(c) {
c %>% pivot_longer(matches("\\.C|zdff"), names_to = "data", values_to = "value") %>% 
      mutate(Channel = str_extract(data, "C[[:digit:]]+"), Channel = ifelse(data == "fitz.C415" | data == "zdff", "C470", Channel)) %>%
      mutate(data = factor(data, levels = c("raw.C415", "raw.C470", "lpf.C415", "lpf.C470", "bl.C415", "bl.C470", "cor.C415", "cor.C470", "z.C415", "z.C470", "fitz.C415", "zdff")),
             Channel = factor(Channel, levels = c("C415", "C470"))) %>%
      ggplot(aes(x=Relative.time, y=value)) + 
             facet_wrap(~Channel, ncol = 1, scales = "free", drop = F) +
             geom_line(aes(color=data)) +
             theme_classic()
}
demo_plot(comb)  # Note the bleaching effect esp. in the first 2-3 mins of recording

################# PART 2 : SIGNAL ANALYSIS PER SE ######################################################################

# Useful custom low-pass 4th-order Butterworth filter function /!\ Can not filter to higher frequencies than the Nyquist frequency (sampling frequency/2)
low_pass_filter <- function(lpf, sampf, v) { # lpf: frequency of the filter, sampf: original sampling frequency, v: time series to filter
  bf <- butter(4, (lpf/(sampf/2)), type="low")  # 4 for 4th-order
  output <- c(rep(v[1], round(0.05*length(v))), v, rep(v[length(v)], round(0.05*length(v)))) %>% # filtfilt function struggles with extremities, so add "fake" data on both sides to absorb the weird behavior
  filtfilt(bf, .) %>%
  .[(round(0.05*length(v))+1):(length(.)-round(0.05*length(v)))] # remove "fake" data
  return(output)
} 

# Step 1: low-pass filter both 470 and 415 signals at 5Hz to remove some meaningless noise (electrical...)
comb <- mutate(comb, lpf.C470 = low_pass_filter(5, sampf, comb$raw.C470), 
                     lpf.C415 = low_pass_filter(5, sampf, comb$raw.C415))  ;  head(comb)
demo_plot(comb) # filtering/smoothing not really visible at that level (as it should not have removed biologically-relevant signal), but if we zoom in a lot...
demo_plot(comb) + scale_x_continuous(limits = c(200,202)) # ... some smoothing can be seen. This step used to be more important for older photodetection systems where sampling was done at much higher frequencies (300-400Hz)

# Step 2: detect baseline (using LOWESS smoother) & subtract baseline from lpf signals.
# This corrects for bleaching mainly, and does so independently for each channel (esp. relevant when doing dual-color imaging, and for longer recordings)
# Alternatives for baseline detection are airPLS or window moving average (not coded in here but suggested in original manuscript)
comb <- mutate(comb, bl.C470 = lowess(comb$lpf.C470)$y,
                     bl.C415 = lowess(comb$lpf.C415)$y) %>%
        mutate(cor.C470 = lpf.C470 - bl.C470,
               cor.C415 = lpf.C415 - bl.C415)  ;  head(comb)
demo_plot(dplyr::select(comb, c(Relative.time, matches("lpf|bl"))))
demo_plot(dplyr::select(comb, c(Relative.time, matches("cor"))))

# Step 3: standardize each baseline-corrected signal using robust z-scoring 
# takes into account both signal intensity and signal-to-noise ratio
comb <- mutate(comb, z.C470 = (cor.C470 - median(comb$cor.C470))/mad(comb$cor.C470), 
                     z.C415 = (cor.C415 - median(comb$cor.C415))/mad(comb$cor.C415))  ;  head(comb)
demo_plot(dplyr::select(comb, c(Relative.time, matches("z"))))

# Step 4: fit (using a robust linear model, more resistant to outliers) standardized C2 to standardized C1 
#         & extract fitted z.C415 values (model predicted values using model linear coefficients)
model <- rlm(z.C470 ~ z.C415, data = comb)
model$model %>% ggplot(aes(x=z.C415, y=z.C470)) + # this is a visual explanation or what rlm does
  geom_point(size = .1, alpha = .5) + 
  geom_abline(slope = model$coefficients[2], intercept = model$coefficients[1], color = "blue") +
  coord_fixed() + scale_x_continuous(limits = c(-5, 5)) +
  theme_classic()
  
comb <- mutate(comb, fitz.C415 = as.numeric(predict(model)))  ;  head(comb)
demo_plot(dplyr::select(comb, c(Relative.time, matches("z"))))  # Note the weird artefact at the very beginning - happens sometimes, good practice is usually to exclude the first 5-10sec of recording at the very beginning before anything else

# Step 5: subtract fitted standardized C415 to standardized C470 to obtain final z(dF/F)
comb <- mutate(comb, zdff = z.C470 - fitz.C415)  ;  head(comb)
demo_plot(dplyr::select(comb, c(Relative.time, matches("zdff")))) 

# Zoom in to see transient shape
demo_plot(dplyr::select(comb, c(Relative.time, matches("zdff")))) +
  scale_x_continuous(limits = c(535, 595))

### THE END #####################################################################################

# All data analysis can be performed as one function:
signal_correction <- function(rig, sampf) {
  low_pass_filter <- function(lpf, sampf, v) { # lpf: frequency of the filter, sampf: original sampling frequency, v: time series to filter
    bf <- butter(4, (lpf/(sampf/2)), type="low")  # 4 for 4th-order
    output <- c(rep(v[1], round(0.05*length(v))), v, rep(v[length(v)], round(0.05*length(v)))) %>% # filtfilt function struggles with extremities, so add "fake" data on both sides to absorb the weird behavior
      filtfilt(bf, .) %>%
      .[(round(0.05*length(v))+1):(length(.)-round(0.05*length(v)))] # remove "fake" data
    return(output)
  } 
  
  rig <- rig %>% mutate(lpf.C470 = low_pass_filter(5, sampf, .$raw.C470), lpf.C415 = low_pass_filter(5, sampf, .$raw.C415)) %>%          
    mutate(bl.C470 = lowess(.$lpf.C470)$y, bl.C415 = lowess(.$lpf.C415)$y) %>%                                      
    mutate(cor.C470 = lpf.C470 - bl.C470, cor.C415 = lpf.C415 - bl.C415) %>%                                 
    mutate(z.C470 = (cor.C470-median(.$cor.C470))/mad(.$cor.C470), z.C415 = (cor.C415-median(.$cor.C415))/mad(.$cor.C415)) %>%   
    mutate(fitz.C415 = as.numeric(predict(rlm(z.C470 ~ z.C415, data = .)))) %>%  
    mutate(zdff = z.C470 - fitz.C415) 
  return(rig)
}

