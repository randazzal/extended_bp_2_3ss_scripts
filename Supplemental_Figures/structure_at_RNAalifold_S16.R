library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(zoo)
library(ggbeeswarm)

setwd("/data2/lackey_lab/randazza/microexons/SHAPE/")
agap1_h <- read.delim("combined_5NIA/AGAP1_A_HsaEX0003019/AGAP1_A_HsaEX0003019_AGAP1_A_HsaEX0003019::chr2:235717931-235718831_+__profile.txt", sep = "\t")
apbb2_h <- read.delim("combined_5NIA/APBB2_A_HsaEX0005072/APBB2_A_HsaEX0005072_APBB2_A_HsaEX0005072::chr4:40841925-40842825_-__profile.txt", sep = "\t")
asap2_h <- read.delim("combined_5NIA/short/ASAP2_A_HsaEX0006174/ASAP2_A_HsaEX0006174_ASAP2_A_HsaEX0006174_short_profile.txt", sep = "\t")
robo1_h <- read.delim("combined_5NIA/ROBO1_A_HsaEX0054203/ROBO1_A_HsaEX0054203_ROBO1_A_HsaEX0054203::chr3:78692906-78693806_-__profile.txt", sep = "\t")
cpeb4_h <- read.delim("combined_5NIA/short/CPEB4_A_HsaEX0016976/CPEB4_A_HsaEX0016976_CPEB4_A_HsaEX0016976_short_profile.txt", sep = "\t")
emc1_h <- read.delim("combined_5NIA/EMC1_A_HsaEX0022280/EMC1_A_HsaEX0022280_EMC1_A_HsaEX0022280::chr1:19233738-19234638_-__profile.txt", sep = "\t")
fry_h <- read.delim("combined_5NIA/FRY_A_HsaEX0026448/FRY_A_HsaEX0026448_FRY_A_HsaEX0026448::chr13:32257514-32258414_+__profile.txt", sep = "\t")
mef2a_h <- read.delim("combined_5NIA/short/MEF2A_A_HsaEX0038690/MEF2A_A_HsaEX0038690_MEF2A_A_HsaEX0038690_short_profile.txt", sep = "\t")
dock7_h <- read.delim("combined_5NIA/short/DOCK7_A_HsaEX0020638/DOCK7_A_HsaEX0020638_DOCK7_A_HsaEX0020638_short_profile.txt", sep = "\t")
itsn1_h <- read.delim("combined_5NIA/ITSN1_A_HsaEX0032608/ITSN1_A_HsaEX0032608_ITSN1_A_HsaEX0032608::chr21:33801980-33802880_+__profile.txt", sep = "\t")
ptprk_h <- read.delim("combined_5NIA/PTPRK_A_HsaEX0051167/PTPRK_A_HsaEX0051167_PTPRK_A_HsaEX0051167::chr6:128000757-128001657_-__profile.txt", sep = "\t")
agap1_c <- read.delim("combined_5NIA/short/AGAP1_A_GgaEX0005391/AGAP1_A_GgaEX0005391_AGAP1_A_GgaEX0005391_short_profile.txt", sep = "\t")
apbb2_c <- read.delim("combined_5NIA/short/APBB2_A_GgaEX1008412/APBB2_A_GgaEX1008412_APBB2_A_GgaEX1008412_short_profile.txt", sep = "\t")
asap2_c <- read.delim("combined_5NIA/short/ASAP2_A_GgaEX0029181/ASAP2_A_GgaEX0029181_ASAP2_A_GgaEX0029181_short_profile.txt", sep = "\t")
robo1_c <- read.delim("combined_5NIA/short/ROBO1_A_GgaEX0022749/ROBO1_A_GgaEX0022749_ROBO1_A_GgaEX0022749_short_profile.txt", sep = "\t")
cpeb4_c <- read.delim("combined_5NIA/short/CPEB4_A_HsaEX0016976/CPEB4_A_HsaEX0016976_CPEB4_A_HsaEX0016976_short_profile.txt", sep = "\t")
emc1_c <- read.delim("combined_5NIA/short/EMC1_A_GgaEX0005294/EMC1_A_GgaEX0005294_EMC1_A_GgaEX0005294_short_profile.txt", sep = "\t")
fry_c <- read.delim("combined_5NIA/short/FRY_A_GgaEX0029190/FRY_A_GgaEX0029190_FRY_A_GgaEX0029190_short_profile.txt", sep = "\t")
mef2a_c <- read.delim("combined_5NIA/short/MEF2A_A_GgaEX0010134/MEF2A_A_GgaEX0010134_MEF2A_A_GgaEX0010134_short_profile.txt", sep = "\t")
dock7_c <- read.delim("combined_5NIA/short/DOCK7_A_GgaEX0016434/DOCK7_A_GgaEX0016434_DOCK7_A_GgaEX0016434_short_profile.txt", sep = "\t")
itsn1_c <- read.delim("combined_5NIA/short/ITSN1_A_GgaEX0023569/ITSN1_A_GgaEX0023569_ITSN1_A_GgaEX0023569_short_profile.txt", sep = "\t")
ptprk_c <- read.delim("combined_5NIA/short/PTPRK_A_GgaEX0021629/PTPRK_A_GgaEX0021629_PTPRK_A_GgaEX0021629_short_profile.txt", sep = "\t")

clamp_reactivity <- function(df) {
  df %>%
    dplyr::select(Nucleotide, Sequence, Norm_profile) %>%
    dplyr::mutate(
      React1 = pmin(pmax(Norm_profile, 0), 4),
    )
}
#hairpin 23-28 & 34-39; 29-33 loop
agap1_h <- clamp_reactivity(agap1_h)
agap1_h <- agap1_h %>%
  filter(Nucleotide %in% c(23:39, 472:549))
apbb2_h <- clamp_reactivity(apbb2_h)
apbb2_h <- apbb2_h %>%
  filter(Nucleotide %in% c(23:39, 311:523))
asap2_h <- clamp_reactivity(asap2_h)
asap2_h <- asap2_h %>%
  filter(Nucleotide %in% c(23:39, 223:258))
robo1_h <- clamp_reactivity(robo1_h)
robo1_h <- robo1_h %>%
  filter(Nucleotide %in% c(23:39, 497:674))
cpeb4_h <- clamp_reactivity(cpeb4_h)
cpeb4_h <- cpeb4_h %>%
  filter(Nucleotide %in% c(23:39, 157:267))
emc1_h <- clamp_reactivity(emc1_h)
emc1_h <- emc1_h %>%
  filter(Nucleotide %in% c(23:39, 463:555))
fry_h <- clamp_reactivity(fry_h)
fry_h <- fry_h %>%
  filter(Nucleotide %in% c(23:39, 455:551))
mef2a_h <- clamp_reactivity(mef2a_h)
mef2a_h <- mef2a_h %>%
  filter(Nucleotide %in% c(23:39, 207:301))
dock7_h <- clamp_reactivity(dock7_h)
dock7_h <- dock7_h %>%
  filter(Nucleotide %in% c(23:39, 247:352))
itsn1_h <- clamp_reactivity(itsn1_h)
itsn1_h <- itsn1_h %>%
  filter(Nucleotide %in% c(23:39, 459:515))
ptprk_h <- clamp_reactivity(ptprk_h)
ptprk_h <- ptprk_h %>%
  filter(Nucleotide %in% c(23:39, 369:644))
agap1_c <- clamp_reactivity(agap1_c)
agap1_c <- agap1_c %>%
  filter(Nucleotide %in% c(23:39, 223:293))
apbb2_c <- clamp_reactivity(apbb2_c)
apbb2_c <- apbb2_c %>%
  filter(Nucleotide %in% c(23:39, 75:268))
asap2_c <- clamp_reactivity(asap2_c)
asap2_c <- asap2_c %>%
  filter(Nucleotide %in% c(23:39, 223:258))
robo1_c <- clamp_reactivity(robo1_c)
robo1_c <- robo1_c %>%
  filter(Nucleotide %in% c(23:39, 248:417))
cpeb4_c <- clamp_reactivity(cpeb4_c)
cpeb4_c <- cpeb4_c %>%
  filter(Nucleotide %in% c(23:39, 173:267))
emc1_c <- clamp_reactivity(emc1_c)
emc1_c <- emc1_c %>%
  filter(Nucleotide %in% c(23:39, 211:303))
fry_c <- clamp_reactivity(fry_c)
fry_c <- fry_c %>%
  filter(Nucleotide %in% c(23:39, 205:297))
mef2a_c <- clamp_reactivity(mef2a_c)
mef2a_c <- mef2a_c %>%
  filter(Nucleotide %in% c(23:39, 203:298))
dock7_c <- clamp_reactivity(dock7_c)
dock7_c <- dock7_c %>%
  filter(Nucleotide %in% c(23:39, 247:348))
itsn1_c <- clamp_reactivity(itsn1_c)
itsn1_c <- itsn1_c %>%
  filter(Nucleotide %in% c(23:39, 211:266))
ptprk_c <- clamp_reactivity(ptprk_c)
ptprk_c <- ptprk_c %>%
  filter(Nucleotide %in% c(23:39, 132:382))

assign_class <- function(df, paired_range, unpaired_range, hairpin_stem_range, hairpin_loop_range) {
  df$class <- "outside"
  df$class[df$Nucleotide %in% paired_range]   <- "paired"
  df$class[df$Nucleotide %in% unpaired_range] <- "unpaired"
  df$class[df$Nucleotide %in% hairpin_stem_range] <- "hairpin_stem"
  df$class[df$Nucleotide %in% hairpin_loop_range] <- "hairpin_loop"
  
  df$class <- factor(df$class,
                     levels = c("hairpin_stem", "hairpin_loop", "paired", "unpaired", "outside"))
  return(df)
}

assign_class_shuffled <- function(df, paired_range, unpaired_range, hairpin_stem_range, hairpin_loop_range, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Shuffle Norm_profile and create React1/React2
  shuffled <- sample(df$React1)
  df$React1 <- shuffled
  
  # Assign class same as above
  df$class <- "outside"
  df$class[df$Nucleotide %in% paired_range]   <- "paired"
  df$class[df$Nucleotide %in% unpaired_range] <- "unpaired"
  df$class[df$Nucleotide %in% hairpin_stem_range] <- "hairpin_stem"
  df$class[df$Nucleotide %in% hairpin_loop_range] <- "hairpin_loop"
  
  df$class <- factor(df$class,
                     levels = c("hairpin_stem", "hairpin_loop", "paired", "unpaired", "outside"))
  return(df)
}

hsa_agap1_class <- assign_class(agap1_h,
                                paired_range = c(472, 473,475:482, 484:486, 491:493, 496:499, 501:504, 506, 507,
                                                 520:524, 545:549),
                                unpaired_range = c(474, 483, 487:490, 494, 495, 500, 505,
                                                   525:544),
                                hairpin_stem_range = c(23:28, 34:39),
                                hairpin_loop_range = c(29:33))
hsa_agap1_class$dataset <- "experimental"
hsa_agap1_shuf <- assign_class_shuffled(agap1_h,
                                        paired_range = c(472, 473,475:482, 484:486, 491:493, 496:499, 501:504, 506, 507,
                                                         520:524, 545:549),
                                        unpaired_range = c(474, 483, 487:490, 494, 495, 500, 505,
                                                           525:544),
                                        hairpin_stem_range = c(23:28, 34:39),
                                        hairpin_loop_range = c(29:33), seed = 123)
hsa_agap1_shuf$dataset <- "shuffled"
combo <- rbind(hsa_agap1_class, hsa_agap1_shuf) %>%
  mutate(gene = "AGAP1", species = "human")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/agap1_exp_shuffle.txt", sep = "\t")

hsa_class <- assign_class(apbb2_h,
                          paired_range = c(311:315,317, 318, 340, 341, 343:347,
                                           382:384, 389:391,
                                           429:433, 438:443, 474:479, 482:486,
                                           495:501, 517:523),
                          unpaired_range = c(316, 319:339, 342,
                                             385:388,
                                             434:437, 444:473, 480, 481,
                                             502:516),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
hsa_class$dataset <- "experimental"
hsa_shuf <- assign_class_shuffled(apbb2_h,
                                  paired_range = c(311:315,317, 318, 340, 341, 343:347,
                                                   382:384, 389:391,
                                                   429:433, 438:443, 474:479, 482:486,
                                                   495:501, 517:523),
                                  unpaired_range = c(316, 319:339, 342,
                                                     385:388,
                                                     434:437, 444:473, 480, 481,
                                                     502:516),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
hsa_shuf$dataset <- "shuffled"
combo <- rbind(hsa_class, hsa_shuf) %>%
  mutate(gene = "APBB2", species = "human")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/apbb2_exp_shuffle.txt", sep = "\t")

hsa_class <- assign_class(asap2_h,
                          paired_range = c(223:234, 236:238, 243:251, 253:258),
                          unpaired_range = c(235, 239:242, 252),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
hsa_class$dataset <- "experimental"
hsa_shuf <- assign_class_shuffled(asap2_h,
                                  paired_range = c(223:234, 236:238, 243:251, 253:258),
                                  unpaired_range = c(235, 239:242, 252),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
hsa_shuf$dataset <- "shuffled"
combo <- rbind(hsa_class, hsa_shuf) %>%
  mutate(gene = "ASAP2", species = "human")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/asap2_exp_shuffle.txt", sep = "\t")

hsa_class <- assign_class(cpeb4_h,
                          paired_range = c(157:160, 162:164, 166:168, 173:178, 180:183,
                                           228:230, 232:235, 237:240, 252:255, 258:261, 265:267),
                          unpaired_range = c(161, 165, 169:172, 179,
                                             231, 236, 241:251, 256, 257, 262:264),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
hsa_class$dataset <- "experimental"
hsa_shuf <- assign_class_shuffled(cpeb4_h,
                                  paired_range = c(157:160, 162:164, 166:168, 173:178, 180:183,
                                                   228:230, 232:235, 237:240, 252:255, 258:261, 265:267),
                                  unpaired_range = c(161, 165, 169:172, 179,
                                                     231, 236, 241:251, 256, 257, 262:264),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
hsa_shuf$dataset <- "shuffled"
combo <- rbind(hsa_class, hsa_shuf) %>%
  mutate(gene = "CPEB4", species = "human")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/cpeb4_exp_shuffle.txt", sep = "\t")

hsa_class <- assign_class(dock7_h,
                          paired_range = c(247:256, 270:279,
                                           328:335, 337, 338, 343:352),
                          unpaired_range = c(257:269,
                                             336, 339:342),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
hsa_class$dataset <- "experimental"
hsa_shuf <- assign_class_shuffled(dock7_h,
                                  paired_range = c(247:256, 270:279,
                                                   328:335, 337, 338, 343:352),
                                  unpaired_range = c(257:269,
                                                     336, 339:342),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
hsa_shuf$dataset <- "shuffled"
combo <- rbind(hsa_class, hsa_shuf) %>%
  mutate(gene = "DOCK7", species = "human")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/dock7_exp_shuffle.txt", sep = "\t")

hsa_class <- assign_class(emc1_h,
                          paired_range = c(463:468, 479:482, 491:494, 502:507,
                                           537, 538, 541:543, 549:551, 554, 555),
                          unpaired_range = c(469:478, 483:490, 495:501,
                                             539, 540, 544:548, 552, 553),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
hsa_class$dataset <- "experimental"
hsa_shuf <- assign_class_shuffled(emc1_h,
                                  paired_range = c(463:468, 479:482, 491:494, 502:507,
                                                   537, 538, 541:543, 549:551, 554, 555),
                                  unpaired_range = c(469:478, 483:490, 495:501,
                                                     539, 540, 544:548, 552, 553),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
hsa_shuf$dataset <- "shuffled"
combo <- rbind(hsa_class, hsa_shuf) %>%
  mutate(gene = "EMC1", species = "human")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/emc1_exp_shuffle.txt", sep = "\t")

hsa_class <- assign_class(fry_h,
                          paired_range = c(455:463, 466:471, 490:495, 499:507,
                                           539:543, 547:551),
                          unpaired_range = c(464, 465, 472:489, 496:498,
                                             544:546),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
hsa_class$dataset <- "experimental"
hsa_shuf <- assign_class_shuffled(fry_h,
                                  paired_range = c(455:463, 466:471, 490:495, 499:507,
                                                   539:543, 547:551),
                                  unpaired_range = c(464, 465, 472:489, 496:498,
                                                     544:546),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
hsa_shuf$dataset <- "shuffled"
combo <- rbind(hsa_class, hsa_shuf) %>%
  mutate(gene = "FRY", species = "human")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/fry_exp_shuffle.txt", sep = "\t")

hsa_class <- assign_class(itsn1_h,
                          paired_range = c(459:467, 473:478, 480:484, 489:493, 495:500, 506:512, 514, 515),
                          unpaired_range = c(468:472, 479, 485:488, 494, 501:505, 513),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
hsa_class$dataset <- "experimental"
hsa_shuf <- assign_class_shuffled(itsn1_h,
                                  paired_range = c(459:467, 473:478, 480:484, 489:493, 495:500, 506:512, 514, 515),
                                  unpaired_range = c(468:472, 479, 485:488, 494, 501:505, 513),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
hsa_shuf$dataset <- "shuffled"
combo <- rbind(hsa_class, hsa_shuf) %>%
  mutate(gene = "ITSN1", species = "human")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/itsn1_exp_shuffle.txt", sep = "\t")

hsa_class <- assign_class(mef2a_h,
                          paired_range = c(207, 208, 210:219, 223:227, 230:234, 238:242, 246:250, 254:257, 262:267, 269, 270,
                                           279:282, 284:287, 293:296, 298:301),
                          unpaired_range = c(209, 220:222, 228, 229, 235:237, 243:245, 251:253, 258:261, 268,
                                             283, 288:292, 297),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
hsa_class$dataset <- "experimental"
hsa_shuf <- assign_class_shuffled(mef2a_h,
                                  paired_range = c(207, 208, 210:219, 223:227, 230:234, 238:242, 246:250, 254:257, 262:267, 269, 270,
                                                   279:282, 284:287, 293:296, 298:301),
                                  unpaired_range = c(209, 220:222, 228, 229, 235:237, 243:245, 251:253, 258:261, 268,
                                                     283, 288:292, 297),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
hsa_shuf$dataset <- "shuffled"
combo <- rbind(hsa_class, hsa_shuf) %>%
  mutate(gene = "MEF2A", species = "human")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/mef2a_exp_shuffle.txt", sep = "\t")

hsa_class <- assign_class(ptprk_h,
                          paired_range = c(369:373, 377:380, 386:389, 393:397,
                                           462:469, 473:479, 484:487, 496:499, 505:507, 509:512, 516:523,
                                           607:616, 630:635, 641:644),
                          unpaired_range = c(374:376, 381:385, 390:392,
                                             470:472, 480:483, 488:495, 500:504, 508, 513:515,
                                             617:629, 636:640),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
hsa_class$dataset <- "experimental"
hsa_shuf <- assign_class_shuffled(ptprk_h,
                                  paired_range = c(369:373, 377:380, 386:389, 393:397,
                                                   462:469, 473:479, 484:487, 496:499, 505:507, 509:512, 516:523,
                                                   607:616, 630:635, 641:644),
                                  unpaired_range = c(374:376, 381:385, 390:392,
                                                     470:472, 480:483, 488:495, 500:504, 508, 513:515,
                                                     617:629, 636:640),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
hsa_shuf$dataset <- "shuffled"
combo <- rbind(hsa_class, hsa_shuf) %>%
  mutate(gene = "PTPRK", species = "human")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/ptprk_exp_shuffle.txt", sep = "\t")

hsa_class <- assign_class(robo1_h,
                          paired_range = c(497, 498, 500:507, 515:524,
                                           528:531, 534:537, 539, 540, 544, 545, 547:550, 553:556,
                                           610:618, 621:625, 627:629, 634:636, 638:642, 649:652, 659:662, 666:674),
                          unpaired_range = c(499, 508:514, 525:527,
                                             532, 533, 538, 541:543, 546, 551, 552,
                                             619, 620, 626, 630:633, 637, 643:648, 653:658, 663:665),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
hsa_class$dataset <- "experimental"
hsa_shuf <- assign_class_shuffled(robo1_h,
                                  paired_range = c(497, 498, 500:507, 515:524,
                                                   528:531, 534:537, 539, 540, 544, 545, 547:550, 553:556,
                                                   610:618, 621:625, 627:629, 634:636, 638:642, 649:652, 659:662, 666:674),
                                  unpaired_range = c(499, 508:514, 525:527,
                                                     532, 533, 538, 541:543, 546, 551, 552,
                                                     619, 620, 626, 630:633, 637, 643:648, 653:658, 663:665),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
hsa_shuf$dataset <- "shuffled"
combo <- rbind(hsa_class, hsa_shuf) %>%
  mutate(gene = "ROBO1", species = "human")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/robo1_exp_shuffle.txt", sep = "\t")

###chicken
gga_agap1_class <- assign_class(agap1_c,
                                paired_range = c(223, 224,226:233, 235:237, 242:244, 247:250, 252:255, 257, 258,
                                                 271:275, 291:293),
                                unpaired_range = c(225, 234, 238:241, 245, 246, 251, 256,
                                                   276:290),
                                hairpin_stem_range = c(23:28, 34:39),
                                hairpin_loop_range = c(29:33))
gga_agap1_class$dataset <- "experimental"
gga_agap1_shuf <- assign_class_shuffled(agap1_c,
                                        paired_range = c(223, 224,226:233, 235:237, 242:244, 247:250, 252:255, 257, 258,
                                                         271:275, 291:293),
                                        unpaired_range = c(225, 234, 238:241, 245, 246, 251, 256,
                                                           276:290),
                                        hairpin_stem_range = c(23:28, 34:39),
                                        hairpin_loop_range = c(29:33), seed = 123)
gga_agap1_shuf$dataset <- "shuffled"
combo <- rbind(gga_agap1_class, gga_agap1_shuf) %>%
  mutate(gene = "AGAP1", species = "chicken")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/gga_agap1_exp_shuffle.txt", sep = "\t")

gga_class <- assign_class(apbb2_c,
                          paired_range = c(75:76, 80, 85, 110, 113, 116:118,
                                           142:143, 151:153,
                                           175:179, 189:198, 223:229, 231:234,
                                           245:251, 263:268),
                          unpaired_range = c(81:84, 79, 86:109, 115,
                                             144:150,
                                             182:188, 195:222, 229, 230,
                                             252:262),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
gga_class$dataset <- "experimental"
gga_shuf <- assign_class_shuffled(apbb2_c,
                                  paired_range = c(75:76, 80, 85, 110, 113, 116:118,
                                                   142:143, 151:153,
                                                   175:179, 189:198, 223:229, 231:234,
                                                   245:251, 263:268),
                                  unpaired_range = c(81:84, 79, 86:109, 115,
                                                     144:150,
                                                     182:188, 195:222, 229, 230,
                                                     252:262),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
gga_shuf$dataset <- "shuffled"
combo <- rbind(gga_class, gga_shuf) %>%
  mutate(gene = "APBB2", species = "chicken")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/gga_apbb2_exp_shuffle.txt", sep = "\t")

gga_class <- assign_class(asap2_c,
                          paired_range = c(223:234, 236:238, 243:251, 253:258),
                          unpaired_range = c(235, 239:242, 252),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
gga_class$dataset <- "experimental"
gga_shuf <- assign_class_shuffled(asap2_c,
                                  paired_range = c(223:234, 236:238, 243:251, 253:258),
                                  unpaired_range = c(235, 239:242, 252),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
gga_shuf$dataset <- "shuffled"
combo <- rbind(gga_class, gga_shuf) %>%
  mutate(gene = "ASAP2", species = "chicken")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/gga_asap2_exp_shuffle.txt", sep = "\t")

gga_class <- assign_class(cpeb4_c,
                          paired_range = c(173:175, 177:178, 180:182, 185:187, 189:192,
                                           228:230, 232:235, 237:240, 252:255, 258:261, 265:267),
                          unpaired_range = c(176, 179, 183:184, 188,
                                             231, 236, 241:251, 256, 257, 262:264),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
gga_class$dataset <- "experimental"
gga_shuf <- assign_class_shuffled(cpeb4_c,
                                  paired_range = c(173:175, 177:178, 180:182, 185:187, 189:192,
                                                   228:230, 232:235, 237:240, 252:255, 258:261, 265:267),
                                  unpaired_range = c(176, 179, 183:184, 188,
                                                     231, 236, 241:251, 256, 257, 262:264),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
gga_shuf$dataset <- "shuffled"
combo <- rbind(gga_class, gga_shuf) %>%
  mutate(gene = "CPEB4", species = "chicken")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/gga_cpeb4_exp_shuffle.txt", sep = "\t")

gga_class <- assign_class(dock7_c,
                          paired_range = c(247:256, 273:281,
                                           325:331, 333, 334, 339:348),
                          unpaired_range = c(257:272,
                                             332, 335:338),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
gga_class$dataset <- "experimental"
gga_shuf <- assign_class_shuffled(dock7_c,
                                  paired_range = c(247:256, 273:281,
                                                   325:331, 333, 334, 339:348),
                                  unpaired_range = c(257:272,
                                                     332, 335:338),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
gga_shuf$dataset <- "shuffled"
combo <- rbind(gga_class, gga_shuf) %>%
  mutate(gene = "DOCK7", species = "chicken")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/gga_dock7_exp_shuffle.txt", sep = "\t")

gga_class <- assign_class(emc1_c,
                          paired_range = c(211:213, 228:231, 241:244, 252:256,
                                           286, 287, 289, 298:299, 303),
                          unpaired_range = c(214:226, 233:241, 246:251,
                                             288, 290:297, 300:302),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
gga_class$dataset <- "experimental"
gga_shuf <- assign_class_shuffled(emc1_c,
                                  paired_range = c(211:213, 228:231, 241:244, 252:256,
                                                   286, 287, 289, 298:299, 303),
                                  unpaired_range = c(214:226, 233:241, 246:251,
                                                     288, 290:297, 300:302),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
gga_shuf$dataset <- "shuffled"
combo <- rbind(gga_class, gga_shuf) %>%
  mutate(gene = "EMC1", species = "chicken")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/gga_emc1_exp_shuffle.txt", sep = "\t")

gga_class <- assign_class(fry_c,
                          paired_range = c(205:214, 217:222, 241:246, 250:258,
                                           287:291, 295:297),
                          unpaired_range = c(215, 216, 223:240, 247:249,
                                             292:294),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
gga_class$dataset <- "experimental"
gga_shuf <- assign_class_shuffled(fry_c,
                                  paired_range = c(205:214, 217:222, 241:246, 250:258,
                                                   287:291, 295:297),
                                  unpaired_range = c(215, 216, 223:240, 247:249,
                                                     292:294),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
gga_shuf$dataset <- "shuffled"
combo <- rbind(gga_class, gga_shuf) %>%
  mutate(gene = "FRY", species = "chicken")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/gga_fry_exp_shuffle.txt", sep = "\t")

gga_class <- assign_class(itsn1_c,
                          paired_range = c(211:219, 225:229, 231:235, 240:244, 246:252, 257:263, 265, 266),
                          unpaired_range = c(220:223, 230, 236:239, 245, 253:256, 264),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
gga_class$dataset <- "experimental"
gga_shuf <- assign_class_shuffled(itsn1_c,
                                  paired_range = c(211:219, 225:229, 231:235, 240:244, 246:252, 257:263, 265, 266),
                                  unpaired_range = c(220:223, 230, 236:239, 245, 253:256, 264),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
gga_shuf$dataset <- "shuffled"
combo <- rbind(gga_class, gga_shuf) %>%
  mutate(gene = "ITSN1", species = "chicken")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/gga_itsn1_exp_shuffle.txt", sep = "\t")

gga_class <- assign_class(mef2a_c,
                          paired_range = c(203, 208, 210:219, 223:227, 229:234, 238:242, 246:250, 254:257, 262:267, 269, 270,
                                           280:281, 284:286, 292:293, 294:298),
                          unpaired_range = c(209, 220:222, 228, 235:237, 243:245, 251:252, 258:261, 268,
                                             283, 288:291),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
gga_class$dataset <- "experimental"
gga_shuf <- assign_class_shuffled(mef2a_c,
                                  paired_range = c(203, 208, 210:219, 223:227, 229:234, 238:242, 246:250, 254:257, 262:267, 269, 270,
                                                   280:281, 284:286, 292:293, 294:298),
                                  unpaired_range = c(209, 220:222, 228, 235:237, 243:245, 251:252, 258:261, 268,
                                                     283, 288:291),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
gga_shuf$dataset <- "shuffled"
combo <- rbind(gga_class, gga_shuf) %>%
  mutate(gene = "MEF2A", species = "chicken")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/gga_mef2a_exp_shuffle.txt", sep = "\t")

gga_class <- assign_class(ptprk_c,
                          paired_range = c(132:134, 137:138, 143:144, 146,
                                           207:219, 223:229, 235:237, 246:249, 255:257, 259:262, 266:280,
                                           346:354, 370:375, 379:382),
                          unpaired_range = c(135:136, 139:142, 145,
                                             220:222, 230:234, 238:245, 250:254, 258, 263:265,
                                             356:369, 376:3780),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
gga_class$dataset <- "experimental"
gga_shuf <- assign_class_shuffled(ptprk_c,
                                  paired_range = c(132:134, 137:138, 143:144, 146,
                                                   207:219, 223:229, 235:237, 246:249, 255:257, 259:262, 266:280,
                                                   346:354, 370:375, 379:382),
                                  unpaired_range = c(135:136, 139:142, 145,
                                                     220:222, 230:234, 238:245, 250:254, 258, 263:265,
                                                     356:369, 376:3780),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
gga_shuf$dataset <- "shuffled"
combo <- rbind(gga_class, gga_shuf) %>%
  mutate(gene = "PTPRK", species = "chicken")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/gga_ptprk_exp_shuffle.txt", sep = "\t")

gga_class <- assign_class(robo1_c,
                          paired_range = c(248, 250:257, 265:274,
                                           278:281, 284:288, 289, 290, 294, 295, 297:300, 302:306,
                                           356:365, 368:372, 374:376, 378:380, 381:383, 392:395, 401:405, 410:417),
                          unpaired_range = c(249, 258:264, 275:277,
                                             282, 283, 291:293, 296, 301,
                                             366, 367, 373, 377, 384:391, 396:400, 406:409),
                          hairpin_stem_range = c(23:28, 34:39),
                          hairpin_loop_range = c(29:33))
gga_class$dataset <- "experimental"
gga_shuf <- assign_class_shuffled(robo1_c,
                                  paired_range = c(248, 250:257, 265:274,
                                                   278:281, 284:288, 289, 290, 294, 295, 297:300, 302:306,
                                                   356:365, 368:372, 374:376, 378:380, 381:383, 392:395, 401:405, 410:417),
                                  unpaired_range = c(249, 258:264, 275:277,
                                                     282, 283, 291:293, 296, 301,
                                                     366, 367, 373, 377, 384:391, 396:400, 406:409),
                                  hairpin_stem_range = c(23:28, 34:39),
                                  hairpin_loop_range = c(29:33), seed = 123)
gga_shuf$dataset <- "shuffled"
combo <- rbind(gga_class, gga_shuf) %>%
  mutate(gene = "ROBO1", species = "chicken")
write.table(combo, "combined_5NIA/human_v_chicken_images/RNAalifold_structure/gga_robo1_exp_shuffle.txt", sep = "\t")

## plotting
agap1_hsa <- read.delim("agap1_exp_shuffle.txt", sep = "\t")
apbb2_hsa <- read.delim("apbb2_exp_shuffle.txt", sep = "\t")
asap2_hsa <- read.delim("asap2_exp_shuffle.txt", sep = "\t")
cpeb4_hsa <- read.delim("cpeb4_exp_shuffle.txt", sep = "\t")
dock7_hsa <- read.delim("dock7_exp_shuffle.txt", sep = "\t")
emc1_hsa <- read.delim("emc1_exp_shuffle.txt", sep = "\t")
fry_hsa <- read.delim("fry_exp_shuffle.txt", sep = "\t")
itsn1_hsa <- read.delim("itsn1_exp_shuffle.txt", sep = "\t")
mef2a_hsa <- read.delim("mef2a_exp_shuffle.txt", sep = "\t")
ptprk_hsa <- read.delim("ptprk_exp_shuffle.txt", sep = "\t")
robo1_hsa <- read.delim("robo1_exp_shuffle.txt", sep = "\t")

agap1_gga <- read.delim("gga_agap1_exp_shuffle.txt", sep = "\t")
apbb2_gga <- read.delim("gga_apbb2_exp_shuffle.txt", sep = "\t")
asap2_gga <- read.delim("gga_asap2_exp_shuffle.txt", sep = "\t")
cpeb4_gga <- read.delim("gga_cpeb4_exp_shuffle.txt", sep = "\t")
dock7_gga <- read.delim("gga_dock7_exp_shuffle.txt", sep = "\t")
emc1_gga <- read.delim("gga_emc1_exp_shuffle.txt", sep = "\t")
fry_gga <- read.delim("gga_fry_exp_shuffle.txt", sep = "\t")
itsn1_gga <- read.delim("gga_itsn1_exp_shuffle.txt", sep = "\t")
mef2a_gga <- read.delim("gga_mef2a_exp_shuffle.txt", sep = "\t")
ptprk_gga <- read.delim("gga_ptprk_exp_shuffle.txt", sep = "\t")
robo1_gga <- read.delim("gga_robo1_exp_shuffle.txt", sep = "\t")

chick <- rbind(agap1_gga, apbb2_gga, asap2_gga, cpeb4_gga, dock7_gga, emc1_gga, fry_gga,
               itsn1_gga, mef2a_gga, ptprk_gga, robo1_gga)
human <- rbind(agap1_hsa, apbb2_hsa, asap2_hsa, cpeb4_hsa, dock7_hsa, emc1_hsa, fry_hsa,
               itsn1_hsa, mef2a_hsa, ptprk_hsa, robo1_hsa)

combo <- rbind(chick, human)
combo <- combo %>%
  filter(!class == "outside")
combo %>%
  ggplot(aes(x = class, y = React1, fill = dataset)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  theme_bw() +
  labs(x = "Class",
       y = "React1",
       fill = "Species")

chick <- chick %>%
  filter(!class == "outside")
chick %>%
  ggplot(aes(x = class, y = React1, fill = dataset)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  theme_minimal() +
  labs(x = "Class",
       y = "Reactivity",
       fill = "Dataset",
       color = "Gene")
test <- wilcox.test(chick[chick$dataset == "experimental" & chick$class == "hairpin_loop", ]$React1,
                    chick[chick$dataset == "shuffled" & chick$class == "hairpin_loop", ]$React1) #5.25e-12
test <- wilcox.test(chick[chick$dataset == "experimental" & chick$class == "hairpin_stem", ]$React1,
                    chick[chick$dataset == "shuffled" & chick$class == "hairpin_stem", ]$React1) #3.32e-18
test <- wilcox.test(chick[chick$dataset == "experimental" & chick$class == "paired", ]$React1,
                    chick[chick$dataset == "shuffled" & chick$class == "paired", ]$React1) #0.624
test <- wilcox.test(chick[chick$dataset == "experimental" & chick$class == "unpaired", ]$React1,
                    chick[chick$dataset == "shuffled" & chick$class == "unpaired", ]$React1) #0.0235
test <- wilcox.test(chick[chick$dataset == "experimental" & chick$class == "paired", ]$React1,
                    chick[chick$dataset == "experimental" & chick$class == "unpaired", ]$React1) #0.00172
test <- wilcox.test(chick[chick$dataset == "shuffled" & chick$class == "paired", ]$React1,
                    chick[chick$dataset == "shuffled" & chick$class == "unpaired", ]$React1) #0.882

human <- human %>%
  filter(!class == "outside")
human %>%
  ggplot(aes(x = class, y = React1, fill = dataset)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  theme_minimal() +
  labs(x = "Class",
       y = "Reactivity",
       fill = "Dataset",
       color = "Gene")
test <- wilcox.test(human[human$dataset == "experimental" & human$class == "hairpin_loop", ]$React1,
                    human[human$dataset == "shuffled" & human$class == "hairpin_loop", ]$React1) #0.18
test <- wilcox.test(human[human$dataset == "experimental" & human$class == "hairpin_stem", ]$React1,
                    human[human$dataset == "shuffled" & human$class == "hairpin_stem", ]$React1) #4.05e-25
test <- wilcox.test(human[human$dataset == "experimental" & human$class == "paired", ]$React1,
                    human[human$dataset == "shuffled" & human$class == "paired", ]$React1) #0.814
test <- wilcox.test(human[human$dataset == "experimental" & human$class == "unpaired", ]$React1,
                    human[human$dataset == "shuffled" & human$class == "unpaired", ]$React1) #0.0107
test <- wilcox.test(human[human$dataset == "experimental" & human$class == "paired", ]$React1,
                    human[human$dataset == "experimental" & human$class == "unpaired", ]$React1) #0.00349
test <- wilcox.test(human[human$dataset == "shuffled" & human$class == "paired", ]$React1,
                    human[human$dataset == "shuffled" & human$class == "unpaired", ]$React1) #0.915
