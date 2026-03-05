library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(zoo)

setwd("/data2/lackey_lab/randazza/microexons/SHAPE/")
## looking at reactivity levels in human data we have yUnAy
annotate_windows <- function(df,
                             window_size = 5,
                             branchpoint_positions = NULL,
                             five_ss_positions = NULL,
                             three_ss_positions = NULL) {
  
  # Ensure sorted by nucleotide
  df <- df %>% arrange(Nucleotide) %>% 
    mutate(
      Norm_profile = pmin(pmax(Norm_profile, 0), 4)
    )
  
  # Rolling mean of Norm_profile
  df$mean <- zoo::rollapply(
    df$Norm_profile,
    width = window_size,
    FUN = mean,
    align = "center",
    fill = NA
  )
  
  # Compute window ranges
  half_win <- floor(window_size / 2)
  
  df <- df %>%
    mutate(
      window_start = Nucleotide - half_win,
      window_end   = Nucleotide + half_win
    )
  
  # Add default element column
  df$element <- "background"
  
  # Helper function to check overlap
  contains_pos <- function(start, end, positions) {
    sapply(seq_along(start), function(i) any(positions >= start[i] & positions <= end[i]))
  }
  
  # Add annotations
  if (!is.null(branchpoint_positions)) {
    df$element[contains_pos(df$window_start, df$window_end, branchpoint_positions)] <- "branchpoint"
  }
  
  if (!is.null(five_ss_positions)) {
    df$element[contains_pos(df$window_start, df$window_end, five_ss_positions)] <- "5ss"
  }
  
  if (!is.null(three_ss_positions)) {
    df$element[contains_pos(df$window_start, df$window_end, three_ss_positions)] <- "3ss"
  }
  
  return(df)
}
#microexon  
agap1 <- read.delim("combined_5NIA/AGAP1_A_HsaEX0003019/AGAP1_A_HsaEX0003019_AGAP1_A_HsaEX0003019::chr2:235717931-235718831_+__profile.txt", sep = "\t")
agap1 <- agap1 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
agap1 <- agap1[, c("Nucleotide", "Norm_profile")]
agap1_anno <- annotate_windows(
  df = agap1,
  window_size = 5,
  branchpoint_positions = c(446),
  five_ss_positions = c(504),
  three_ss_positions = c(493)
)
agap1_anno <- agap1_anno %>%
  mutate(gene = "AGAP1")
apbb2 <- read.delim("combined_5NIA/APBB2_A_HsaEX0005072/APBB2_A_HsaEX0005072_APBB2_A_HsaEX0005072::chr4:40841925-40842825_-__profile.txt", sep = "\t")
apbb2 <- apbb2 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
apbb2 <- apbb2[, c("Nucleotide", "Norm_profile")]
apbb2_anno <- annotate_windows(
  df = apbb2,
  window_size = 5,
  branchpoint_positions = c(450),
  five_ss_positions = c(499),
  three_ss_positions = c(494)
)
apbb2_anno <- apbb2_anno %>%
  mutate(gene = "APBB2")
asap2 <- read.delim("combined_5NIA/short/ASAP2_A_HsaEX0006174/ASAP2_A_HsaEX0006174_ASAP2_A_HsaEX0006174_short_profile.txt", sep = "\t")
asap2 <- asap2 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 443)
asap2 <- asap2[, c("Nucleotide", "Norm_profile")]
asap2_anno <- annotate_windows(
  df = asap2,
  window_size = 5,
  branchpoint_positions = c(193),
  five_ss_positions = c(252),
  three_ss_positions = c(244)
)
asap2_anno <- asap2_anno %>%
  mutate(gene = "ASAP2")
robo1 <- read.delim("combined_5NIA/ROBO1_A_HsaEX0054203/ROBO1_A_HsaEX0054203_ROBO1_A_HsaEX0054203::chr3:78692906-78693806_-__profile.txt", sep = "\t")
robo1 <- robo1 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
robo1 <- robo1[, c("Nucleotide", "Norm_profile")]
robo1_anno <- annotate_windows(
  df = robo1,
  window_size = 5,
  branchpoint_positions = c(410),
  five_ss_positions = c(502),
  three_ss_positions = c(494)
)
robo1_anno <- robo1_anno %>%
  mutate(gene = "ROBO1")
cpeb4 <- read.delim("combined_5NIA/short/CPEB4_A_HsaEX0016976/CPEB4_A_HsaEX0016976_CPEB4_A_HsaEX0016976_short_profile.txt", sep = "\t")
cpeb4 <- cpeb4 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 443)
cpeb4 <- cpeb4[, c("Nucleotide", "Norm_profile")]
cpeb4_anno <- annotate_windows(
  df = cpeb4,
  window_size = 5,
  branchpoint_positions = c(175),
  five_ss_positions = c(267),
  three_ss_positions = c(244)
)
cpeb4_anno <- cpeb4_anno %>%
  mutate(gene = "CPEB4")
emc1 <- read.delim("combined_5NIA/EMC1_A_HsaEX0022280/EMC1_A_HsaEX0022280_EMC1_A_HsaEX0022280::chr1:19233738-19234638_-__profile.txt", sep = "\t")
emc1 <- emc1 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
emc1 <- emc1[, c("Nucleotide", "Norm_profile")]
emc1_anno <- annotate_windows(
  df = emc1,
  window_size = 5,
  branchpoint_positions = c(461),
  five_ss_positions = c(502),
  three_ss_positions = c(494)
)
emc1_anno <- emc1_anno %>%
  mutate(gene = "EMC1")
fry <- read.delim("combined_5NIA/FRY_A_HsaEX0026448/FRY_A_HsaEX0026448_FRY_A_HsaEX0026448::chr13:32257514-32258414_+__profile.txt", sep = "\t")
fry <- fry %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
fry <- fry[, c("Nucleotide", "Norm_profile")]
fry_anno <- annotate_windows(
  df = fry,
  window_size = 5,
  branchpoint_positions = c(402),
  five_ss_positions = c(501),
  three_ss_positions = c(493)
)
fry_anno <- fry_anno %>%
  mutate(gene = "FRY")
mef2a <- read.delim("combined_5NIA/short/MEF2A_A_HsaEX0038690/MEF2A_A_HsaEX0038690_MEF2A_A_HsaEX0038690_short_profile.txt", sep = "\t")
mef2a <- mef2a %>%
  filter(Nucleotide >= 20 & Nucleotide <= 443)
mef2a <- mef2a[, c("Nucleotide", "Norm_profile")]
mef2a_anno <- annotate_windows(
  df = mef2a,
  window_size = 5,
  branchpoint_positions = c(212),
  five_ss_positions = c(267),
  three_ss_positions = c(244)
)
mef2a_anno <- mef2a_anno %>%
  mutate(gene = "MEF2A")
dock7 <- read.delim("combined_5NIA/short/DOCK7_A_HsaEX0020638/DOCK7_A_HsaEX0020638_DOCK7_A_HsaEX0020638_short_profile.txt", sep = "\t")
dock7 <- dock7 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 443)
dock7 <- dock7[, c("Nucleotide", "Norm_profile")]
dock7_anno <- annotate_windows(
  df = dock7,
  window_size = 5,
  branchpoint_positions = c(214),
  five_ss_positions = c(258),
  three_ss_positions = c(244)
)
dock7_anno <- dock7_anno %>%
  mutate(gene = "DOCK7")
itsn1 <- read.delim("combined_5NIA/ITSN1_A_HsaEX0032608/ITSN1_A_HsaEX0032608_ITSN1_A_HsaEX0032608::chr21:33801980-33802880_+__profile.txt", sep = "\t")
itsn1 <- itsn1 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
itsn1 <- itsn1[, c("Nucleotide", "Norm_profile")]
itsn1_anno <- annotate_windows(
  df = itsn1,
  window_size = 5,
  branchpoint_positions = c(447),
  five_ss_positions = c(507),
  three_ss_positions = c(493)
)
itsn1_anno <- itsn1_anno %>%
  mutate(gene = "ITSN1")
ptprk <- read.delim("combined_5NIA/PTPRK_A_HsaEX0051167/PTPRK_A_HsaEX0051167_PTPRK_A_HsaEX0051167::chr6:128000757-128001657_-__profile.txt", sep = "\t")
ptprk <- ptprk %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
ptprk <- ptprk[, c("Nucleotide", "Norm_profile")]
ptprk_anno <- annotate_windows(
  df = ptprk,
  window_size = 5,
  branchpoint_positions = c(448),
  five_ss_positions = c(505),
  three_ss_positions = c(494)
)
ptprk_anno <- ptprk_anno %>%
  mutate(gene = "PTPRK")
clec16a <- read.delim("combined_5NIA/CLEC16A_A_HsaEX0015596/CLEC16A_A_HsaEX0015596_CLEC16A_A_HsaEX0015596::chr16:10972104-10973004_+__profile.txt", sep = "\t")
clec16a <- clec16a %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
clec16a <- clec16a[, c("Nucleotide", "Norm_profile")]
clec16a_anno <- annotate_windows(
  df = clec16a,
  window_size = 5,
  branchpoint_positions = c(427),
  five_ss_positions = c(498),
  three_ss_positions = c(493)
)
clec16a_anno <- clec16a_anno %>%
  mutate(gene = "CLEC16A")
dctn1 <- read.delim("combined_5NIA/short/DCTN1_A_HsaEX0018764/DCTN1_A_HsaEX0018764_DCTN1_A_HsaEX0018764_short_profile.txt", sep = "\t")
dctn1 <- dctn1 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 443)
dctn1 <- dctn1[, c("Nucleotide", "Norm_profile")]
dctn1_anno <- annotate_windows(
  df = dctn1,
  window_size = 5,
  branchpoint_positions = c(195),
  five_ss_positions = c(258),
  three_ss_positions = c(244)
)
dctn1_anno <- dctn1_anno %>%
  mutate(gene = "DCTN1")

combo <- rbind(agap1_anno, apbb2_anno, asap2_anno, clec16a_anno, cpeb4_anno, dctn1_anno, dock7_anno,
               emc1_anno, fry_anno, itsn1_anno, mef2a_anno, ptprk_anno, robo1_anno)
combo$dataset <- "Microexon"
combo %>%
  ggplot(aes(x = element, y = mean, color = gene)) +
  geom_point()
write.table(combo, "/data2/lackey_lab/randazza/microexons/SHAPE/combined_5NIA/human_element_structures/my_micros_5window.txt", sep = "\t")

##remove windows that overlap with window centered around branchpoint and splice sites
agap1_less <- agap1_anno %>% 
  filter(!Nucleotide %in% c(491, 492, 494, 495, 502, 503, 505, 506, 444, 445, 447, 448))
apbb2_less <- apbb2_anno %>% 
  filter(!Nucleotide %in% c(492, 493, 495, 496, 497, 498, 500, 501, 448, 449, 451, 452))
asap2_less <- asap2_anno %>% 
  filter(!Nucleotide %in% c(242, 243, 245, 246, 250, 251, 253, 254, 191, 192, 194, 195))
clec_less <- clec16a_anno %>% 
  filter(!Nucleotide %in% c(491, 492, 494, 495, 496, 497, 499, 500, 425, 426, 428, 429))
cpeb4_less <- cpeb4_anno %>% 
  filter(!Nucleotide %in% c(242, 243, 245, 246, 265, 266, 268, 269, 173, 174, 176, 177))
dctn1_less <- dctn1_anno %>% 
  filter(!Nucleotide %in% c(242, 243, 245, 246, 256, 257, 259, 260, 193, 194, 196, 197))
dock7_less <- dock7_anno %>% 
  filter(!Nucleotide %in% c(242, 243, 245, 246, 256, 257, 259, 260, 202, 203, 205, 206))
emc1_less <- emc1_anno %>% 
  filter(!Nucleotide %in% c(492, 493, 495, 496, 500, 501, 503, 504, 459, 460, 462, 463))
fry_less <- fry_anno %>% 
  filter(!Nucleotide %in% c(491, 492, 494, 495, 499, 500, 502, 503, 400, 401, 403, 404))
itsn1_less <- itsn1_anno %>% 
  filter(!Nucleotide %in% c(491, 492, 494, 495, 505, 506, 508, 509, 445, 446, 448, 449))
mef2a_less <- mef2a_anno %>% 
  filter(!Nucleotide %in% c(242, 243, 245, 246, 265, 266, 268, 269, 210, 211, 213, 214))
ptprk_less <- ptprk_anno %>% 
  filter(!Nucleotide %in% c(492, 493, 495, 496, 503, 504, 506, 507, 446, 447, 449, 450))
robo1_less <- robo1_anno %>% 
  filter(!Nucleotide %in% c(492, 493, 495, 496, 500, 501, 503, 504, 408, 409, 411, 412))

combo <- rbind(agap1_less, apbb2_less, asap2_less, clec_less, cpeb4_less, dctn1_less, dock7_less,
               emc1_less, fry_less, itsn1_less, mef2a_less, ptprk_less, robo1_less)
combo$dataset <- "Microexon"
combo %>%
  ggplot(aes(x = element, y = mean, color = gene)) +
  geom_point()
write.table(combo, "/data2/lackey_lab/randazza/microexons/SHAPE/combined_5NIA/human_element_structures/my_micros_5window_single_bpss.txt", sep = "\t")

#midexon 
tnrc6a <- read.delim("combined_5NIA/TNRC6A_A_HsaEX0066372/TNRC6A_A_HsaEX0066372_TNRC6A_A_HsaEX0066372::chr16:24804255-24805155_+__profile.txt", sep = "\t")
tnrc6a <- tnrc6a %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
tnrc6a <- tnrc6a[, c("Nucleotide", "Norm_profile")]
tnrc6a_anno <- annotate_windows(
  df = tnrc6a,
  window_size = 5,
  branchpoint_positions = c(409),
  five_ss_positions = c(639),
  three_ss_positions = c(493)
)
tnrc6a_anno <- tnrc6a_anno %>%
  mutate(gene = "TNRC6A")
arhgef7 <- read.delim("combined_5NIA/ARHGEF7_A_HsaEX0005734/ARHGEF7_A_HsaEX0005734_ARHGEF7_A_HsaEX0005734::chr13:111291668-111292568_+__profile.txt", sep = "\t")
arhgef7 <- arhgef7 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
arhgef7 <- arhgef7[, c("Nucleotide", "Norm_profile")]
arhgef7_anno <- annotate_windows(
  df = arhgef7,
  window_size = 5,
  branchpoint_positions = c(450),
  five_ss_positions = c(669),
  three_ss_positions = c(493)
)
arhgef7_anno <- arhgef7_anno %>%
  mutate(gene = "ARHGEF7")
rps6kb2 <- read.delim("combined_5NIA/RPS6KB2_A_HsaEX0055640/RPS6KB2_A_HsaEX0055640_RPS6KB2_A_HsaEX0055640::chr11:67431906-67432806_+__profile.txt", sep = "\t")
rps6kb2 <- rps6kb2 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
rps6kb2 <- rps6kb2[, c("Nucleotide", "Norm_profile")]
rps6kb2_anno <- annotate_windows(
  df = rps6kb2,
  window_size = 5,
  branchpoint_positions = c(459),
  five_ss_positions = c(629),
  three_ss_positions = c(493)
)
rps6kb2_anno <- rps6kb2_anno %>%
  mutate(gene = "RPS6KB2")
cers5 <- read.delim("combined_5NIA/CERS5_A_HsaEX0035443/CERS5_A_HsaEX0035443_CERS5_A_HsaEX0035443::chr12:50132614-50133514_-__profile.txt", sep = "\t")
cers5 <- cers5 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
cers5 <- cers5[, c("Nucleotide", "Norm_profile")]
cers5_anno <- annotate_windows(
  df = cers5,
  window_size = 5,
  branchpoint_positions = c(425),
  five_ss_positions = c(596),
  three_ss_positions = c(494)
)
cers5_anno <- cers5_anno %>%
  mutate(gene = "CERS5")
lrrfip2 <- read.delim("combined_5NIA/LRRFIP2_A_HsaEX0036760/LRRFIP2_A_HsaEX0036760_LRRFIP2_A_HsaEX0036760::chr3:37065875-37066775_-__profile.txt", sep = "\t")
lrrfip2 <- lrrfip2 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
lrrfip2 <- lrrfip2[, c("Nucleotide", "Norm_profile")]
lrrfip2_anno <- annotate_windows(
  df = lrrfip2,
  window_size = 5,
  branchpoint_positions = c(442),
  five_ss_positions = c(595),
  three_ss_positions = c(494)
)
lrrfip2_anno <- lrrfip2_anno %>%
  mutate(gene = "LRRFIP2")
gabbr1 <- read.delim("combined_5NIA/GABBR1_A_HsaEX0026760/GABBR1_A_HsaEX0026760_GABBR1_A_HsaEX0026760::chr6:29608929-29609829_-__profile.txt", sep = "\t")
gabbr1 <- gabbr1 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
gabbr1 <- gabbr1[, c("Nucleotide", "Norm_profile")]
gabbr1_anno <- annotate_windows(
  df = gabbr1,
  window_size = 5,
  branchpoint_positions = c(433),
  five_ss_positions = c(644),
  three_ss_positions = c(494)
)
gabbr1_anno <- gabbr1_anno %>%
  mutate(gene = "GABBR1")
phf21a <- read.delim("combined_5NIA/PHF21A_A_HsaEX0047036/PHF21A_A_HsaEX0047036_PHF21A_A_HsaEX0047036::chr11:45945553-45946453_-__profile.txt", sep = "\t")
phf21a <- phf21a %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
phf21a <- phf21a[, c("Nucleotide", "Norm_profile")]
phf21a_anno <- annotate_windows(
  df = phf21a,
  window_size = 5,
  branchpoint_positions = c(429),
  five_ss_positions = c(657),
  three_ss_positions = c(494)
)
phf21a_anno <- phf21a_anno %>%
  mutate(gene = "PHF21A")

combo <- rbind(arhgef7_anno, cers5_anno, gabbr1_anno, lrrfip2_anno, phf21a_anno, rps6kb2_anno, tnrc6a_anno)
combo$dataset <- "Midexon"
combo %>%
  ggplot(aes(x = element, y = mean, color = gene)) +
  geom_point()
write.table(combo, "/data2/lackey_lab/randazza/microexons/SHAPE/combined_5NIA/human_element_structures/my_avg_5window.txt", sep = "\t")

##remove windows that overlap with window centered around branchpoint and splice sites
arhgef7_less <- arhgef7_anno %>% 
  filter(!Nucleotide %in% c(491, 492, 494, 495, 627, 628, 670, 671, 448, 449, 451, 452))
cers5_less <- cers5_anno %>% 
  filter(!Nucleotide %in% c(492, 493, 495, 496, 594, 595, 597, 598, 423, 424, 426, 427))
gabbr1_less <- gabbr1_anno %>% 
  filter(!Nucleotide %in% c(492, 493, 495, 496, 642, 643, 645, 646, 431, 432, 434, 435))
lrrfip2_less <- lrrfip2_anno %>% 
  filter(!Nucleotide %in% c(492, 493, 495, 496, 593, 594, 596, 597, 440, 441, 443, 444))
phf21a_less <- phf21a_anno %>% 
  filter(!Nucleotide %in% c(492, 493, 495, 496, 655, 656, 658, 659, 427, 428, 430, 431))
rps6kb2_less <- rps6kb2_anno %>% 
  filter(!Nucleotide %in% c(491, 492, 494, 495, 627, 628, 630, 631, 457, 458, 460, 461))
tnrc6a_less <- tnrc6a_anno %>% 
  filter(!Nucleotide %in% c(491, 492, 494, 495, 627, 628, 630, 631, 407, 408, 410, 411))

combo <- rbind(arhgef7_less, cers5_less, gabbr1_less, lrrfip2_less, phf21a_less, rps6kb2_less, tnrc6a_less)
combo$dataset <- "Midexon"
combo %>%
  ggplot(aes(x = element, y = mean, color = gene)) +
  geom_point()
write.table(combo, "/data2/lackey_lab/randazza/microexons/SHAPE/combined_5NIA/human_element_structures/my_avg_5window_single_bpss.txt", sep = "\t")

mymicro <- read.delim("my_micros_5window_single_bpss.txt", sep ="\t")
myavg <- read.delim("my_avg_5window_single_bpss.txt", sep = "\t")
combo <- rbind(mymicro, myavg)
combo$element <- factor(combo$element, levels = c("branchpoint", "3ss", "5ss", "background"))
combo %>%
  ggplot(aes(x = element, y = mean, fill = dataset)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Element", y = "Mean Reactivity")
test <- wilcox.test(combo[combo$element == "branchpoint" & combo$dataset == "Midexon", ]$mean,
                    combo[combo$element == "branchpoint" & combo$dataset == "Microexon", ]$mean) #0.887
test <- wilcox.test(combo[combo$element == "3ss" & combo$dataset == "Midexon", ]$mean,
                    combo[combo$element == "3ss" & combo$dataset == "Microexon", ]$mean) #0.86
test <- wilcox.test(combo[combo$element == "5ss" & combo$dataset == "Midexon", ]$mean,
                    combo[combo$element == "5ss" & combo$dataset == "Microexon", ]$mean) #0.611
test <- wilcox.test(combo[combo$element == "background" & combo$dataset == "Midexon", ]$mean,
                    combo[combo$element == "background" & combo$dataset == "Microexon", ]$mean) #0.0612