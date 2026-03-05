library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(zoo)
library(ggbeeswarm)

setwd("/data2/lackey_lab/randazza/microexons/SHAPE/")
#microexon
agap1 <- read.delim("combined_5NIA/AGAP1_A_HsaEX0003019/AGAP1_A_HsaEX0003019_AGAP1_A_HsaEX0003019::chr2:235717931-235718831_+__profile.txt", sep = "\t")
agap1 <- agap1 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
agap1 <- agap1[, c("Nucleotide", "Norm_profile")]
apbb2 <- read.delim("combined_5NIA/APBB2_A_HsaEX0005072/APBB2_A_HsaEX0005072_APBB2_A_HsaEX0005072::chr4:40841925-40842825_-__profile.txt", sep = "\t")
apbb2 <- apbb2 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
apbb2 <- apbb2[, c("Nucleotide", "Norm_profile")]
asap2 <- read.delim("combined_5NIA/short/ASAP2_A_HsaEX0006174/ASAP2_A_HsaEX0006174_ASAP2_A_HsaEX0006174_short_profile.txt", sep = "\t")
asap2 <- asap2 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 443)
asap2 <- asap2[, c("Nucleotide", "Norm_profile")]
robo1 <- read.delim("combined_5NIA/ROBO1_A_HsaEX0054203/ROBO1_A_HsaEX0054203_ROBO1_A_HsaEX0054203::chr3:78692906-78693806_-__profile.txt", sep = "\t")
robo1 <- robo1 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
robo1 <- robo1[, c("Nucleotide", "Norm_profile")]
cpeb4 <- read.delim("combined_5NIA/short/CPEB4_A_HsaEX0016976/CPEB4_A_HsaEX0016976_CPEB4_A_HsaEX0016976_short_profile.txt", sep = "\t")
cpeb4 <- cpeb4 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 443)
cpeb4 <- cpeb4[, c("Nucleotide", "Norm_profile")]
emc1 <- read.delim("combined_5NIA/EMC1_A_HsaEX0022280/EMC1_A_HsaEX0022280_EMC1_A_HsaEX0022280::chr1:19233738-19234638_-__profile.txt", sep = "\t")
emc1 <- emc1 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
emc1 <- emc1[, c("Nucleotide", "Norm_profile")]
fry <- read.delim("combined_5NIA/FRY_A_HsaEX0026448/FRY_A_HsaEX0026448_FRY_A_HsaEX0026448::chr13:32257514-32258414_+__profile.txt", sep = "\t")
fry <- fry %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
fry <- fry[, c("Nucleotide", "Norm_profile")]
mef2a <- read.delim("combined_5NIA/short/MEF2A_A_HsaEX0038690/MEF2A_A_HsaEX0038690_MEF2A_A_HsaEX0038690_short_profile.txt", sep = "\t")
mef2a <- mef2a %>%
  filter(Nucleotide >= 20 & Nucleotide <= 443)
mef2a <- mef2a[, c("Nucleotide", "Norm_profile")]
dock7 <- read.delim("combined_5NIA/short/DOCK7_A_HsaEX0020638/DOCK7_A_HsaEX0020638_DOCK7_A_HsaEX0020638_short_profile.txt", sep = "\t")
dock7 <- dock7 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 443)
dock7 <- dock7[, c("Nucleotide", "Norm_profile")]
itsn1 <- read.delim("combined_5NIA/ITSN1_A_HsaEX0032608/ITSN1_A_HsaEX0032608_ITSN1_A_HsaEX0032608::chr21:33801980-33802880_+__profile.txt", sep = "\t")
itsn1 <- itsn1 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
itsn1 <- itsn1[, c("Nucleotide", "Norm_profile")]
ptprk <- read.delim("combined_5NIA/PTPRK_A_HsaEX0051167/PTPRK_A_HsaEX0051167_PTPRK_A_HsaEX0051167::chr6:128000757-128001657_-__profile.txt", sep = "\t")
ptprk <- ptprk %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
ptprk <- ptprk[, c("Nucleotide", "Norm_profile")]
clec16a <- read.delim("combined_5NIA/CLEC16A_A_HsaEX0015596/CLEC16A_A_HsaEX0015596_CLEC16A_A_HsaEX0015596::chr16:10972104-10973004_+__profile.txt", sep = "\t")
clec16a <- clec16a %>%
  filter(Nucleotide >= 20 & Nucleotide <= 943)
clec16a <- clec16a[, c("Nucleotide", "Norm_profile")]
dctn1 <- read.delim("combined_5NIA/short/DCTN1_A_HsaEX0018764/DCTN1_A_HsaEX0018764_DCTN1_A_HsaEX0018764_short_profile.txt", sep = "\t")
dctn1 <- dctn1 %>%
  filter(Nucleotide >= 20 & Nucleotide <= 443)
dctn1 <- dctn1[, c("Nucleotide", "Norm_profile")]

#want to compare the reactivities between branchpoint and 3'ss 
#(lets try boxplots into 4 groups; upstream 100, bp-3'ss, exon, downstream 100; all genes combined)
clamp_reactivity <- function(df) {
  df %>%
    dplyr::mutate(
      Norm_profile = pmin(pmax(Norm_profile, 0), 4),
    )
}
agap1 <- clamp_reactivity(agap1)
apbb2 <- clamp_reactivity(apbb2)
asap2 <- clamp_reactivity(asap2)
clec16a <- clamp_reactivity(clec16a)
cpeb4 <- clamp_reactivity(cpeb4)
dctn1 <- clamp_reactivity(dctn1)
dock7 <- clamp_reactivity(dock7)
emc1 <- clamp_reactivity(emc1)
fry <- clamp_reactivity(fry)
itsn1 <- clamp_reactivity(itsn1)
mef2a <- clamp_reactivity(mef2a)
ptprk <- clamp_reactivity(ptprk)
robo1 <- clamp_reactivity(robo1)

assign_class <- function(df, up_range, bp3ss_range, microexon_range, down_range) {
  df$class <- "outside"
  df$class[df$Nucleotide %in% up_range]   <- "upstream"
  df$class[df$Nucleotide %in% bp3ss_range] <- "Branchpoint_2_3ss"
  df$class[df$Nucleotide %in% microexon_range] <- "Microexon"
  df$class[df$Nucleotide %in% down_range] <- "downstream"
  
  df$class <- factor(df$class,
                     levels = c("upstream", "Branchpoint_2_3ss", "Microexon", "downstream", "outside"))
  return(df)
}
agap1_class <- assign_class(agap1,
                            up_range = c(346:445),
                            bp3ss_range = c(446:492),
                            microexon_range = c(493:504),
                            down_range = c(505:604))
agap1_class$gene <- "AGAP1"
apbb2_class <- assign_class(apbb2,
                            up_range = c(530:449),
                            bp3ss_range = c(450:493),
                            microexon_range = c(494:499),
                            down_range = c(500:599))
apbb2_class$gene <- "APBB2"
asap2_class <- assign_class(asap2,
                            up_range = c(93:192),
                            bp3ss_range = c(193:243),
                            microexon_range = c(244:252),
                            down_range = c(253:352))
asap2_class$gene <- "ASAP2"
clec_class <- assign_class(clec16a,
                           up_range = c(327:426),
                           bp3ss_range = c(427:492),
                           microexon_range = c(493:498),
                           down_range = c(499:598))
clec_class$gene <- "CLEC16A"
cpeb4_class_shared <- assign_class(cpeb4,
                                   up_range = c(94:193),
                                   bp3ss_range = c(194:243),
                                   microexon_range = c(244:267),
                                   down_range = c(268:367))
cpeb4_class_shared$gene <- "CPEB4"
cpeb4_class_hum <- assign_class(cpeb4,
                                up_range = c(75:174),
                                bp3ss_range = c(175:243),
                                microexon_range = c(244:267),
                                down_range = c(268:367))
cpeb4_class_hum$gene <- "CPEB4"
dctn1_class <- assign_class(dctn1,
                            up_range = c(95:194),
                            bp3ss_range = c(195:243),
                            microexon_range = c(244:258),
                            down_range = c(259:358))
dctn1_class$gene <- "DCTN1"
dock7_class <- assign_class(dock7,
                            up_range = c(114:213),
                            bp3ss_range = c(214:243),
                            microexon_range = c(244:258),
                            down_range = c(259:358))
dock7_class$gene <- "DOCK7"
emc1_class_shared <- assign_class(emc1,
                                  up_range = c(295:394),
                                  bp3ss_range = c(395:493),
                                  microexon_range = c(494:502),
                                  down_range = c(503:602))
emc1_class_shared$gene <- "EMC1"
emc1_class_hum <- assign_class(emc1,
                               up_range = c(361:460),
                               bp3ss_range = c(461:493),
                               microexon_range = c(494:502),
                               down_range = c(503:602))
emc1_class_hum$gene <- "EMC1"
fry_class <- assign_class(fry,
                          up_range = c(302:401),
                          bp3ss_range = c(402:492),
                          microexon_range = c(493:501),
                          down_range = c(502:601))
fry_class$gene <- "FRY"
itsn1_class_shared <- assign_class(itsn1,
                                   up_range = c(370:469),
                                   bp3ss_range = c(470:492),
                                   microexon_range = c(493:507),
                                   down_range = c(508:607))
itsn1_class_shared$gene <- "ITSN1"
itsn1_class_hum <- assign_class(itsn1,
                                up_range = c(347:446),
                                bp3ss_range = c(447:492),
                                microexon_range = c(493:507),
                                down_range = c(508:607))
itsn1_class_hum$gene <- "ITSN1"
mef2a_class <- assign_class(mef2a,
                            up_range = c(112:211),
                            bp3ss_range = c(212:243),
                            microexon_range = c(244:267),
                            down_range = c(268:367))
mef2a_class$gene <- "MEF2A"
ptprk_class_shared <- assign_class(ptprk,
                                   up_range = c(317:416),
                                   bp3ss_range = c(417:493),
                                   microexon_range = c(494:505),
                                   down_range = c(506:605))
ptprk_class_shared$gene <- "PTPRK"
ptprk_class_hum <- assign_class(ptprk,
                                up_range = c(348:447),
                                bp3ss_range = c(448:493),
                                microexon_range = c(494:505),
                                down_range = c(506:605))
ptprk_class_hum$gene <- "PTPRK"
robo1_class_shared <- assign_class(robo1,
                                   up_range = c(324:423),
                                   bp3ss_range = c(424:493),
                                   microexon_range = c(494:502),
                                   down_range = c(503:602))
robo1_class_shared$gene <- "ROBO1"
robo1_class_hum <- assign_class(robo1,
                                up_range = c(310:409),
                                bp3ss_range = c(410:493),
                                microexon_range = c(494:502),
                                down_range = c(503:602))
robo1_class_hum$gene <- "ROBO1"

#microA  
#splice sites
#agap1 493-504; apbb2 494-499; asap2 244-252; clec16a 493-498; cpeb4 244-267; dctn1 244-258; dock7 244-258; 
#emc1 494-502; fry 493-501; itsn1 493-507; mef2a 244-267; ptprk 494-505; robo1 494-502
#branch point (highest svm_scr in human but also exist in chick sequence)
#shared best in agap1 (446), best apbb2 (450); best asap2 (193); clec16a (427); cpeb4 (shared 194; best 175)
#best dock7 (214); dctn1 (195); emc1 (best 461, shared 395)
#fry (402); itsn1 (best 447, shared 470); mef2a (212); ptprk (best 448; shared 417); robo1 (best 410; shared 424)
combo <- rbind(agap1_class, apbb2_class, asap2_class, clec_class, cpeb4_class_hum, dctn1_class, dock7_class,
               emc1_class_hum, fry_class, itsn1_class_hum, mef2a_class, ptprk_class_hum, robo1_class_hum)
combo %>%
  ggplot(aes(x = class, y = Norm_profile, color = gene)) +
  geom_boxplot() +
  #geom_quasirandom(width = 0.15, size = 1.5, alpha = 0.7) +
  theme_minimal() 

#midexon 
clamp_reactivity <- function(df) {
  df %>%
    dplyr::mutate(
      Norm_profile = pmin(pmax(Norm_profile, 0), 4),
    )
}

tnrc6a <- read.delim("combined_5NIA/TNRC6A_A_HsaEX0066372/TNRC6A_A_HsaEX0066372_TNRC6A_A_HsaEX0066372::chr16:24804255-24805155_+__profile.txt", sep = "\t")
arhgef7 <- read.delim("combined_5NIA/ARHGEF7_A_HsaEX0005734/ARHGEF7_A_HsaEX0005734_ARHGEF7_A_HsaEX0005734::chr13:111291668-111292568_+__profile.txt", sep = "\t")
rps6kb2 <- read.delim("combined_5NIA/RPS6KB2_A_HsaEX0055640/RPS6KB2_A_HsaEX0055640_RPS6KB2_A_HsaEX0055640::chr11:67431906-67432806_+__profile.txt", sep = "\t")
cers5 <- read.delim("combined_5NIA/CERS5_A_HsaEX0035443/CERS5_A_HsaEX0035443_CERS5_A_HsaEX0035443::chr12:50132614-50133514_-__profile.txt", sep = "\t")
lrrfip2 <- read.delim("combined_5NIA/LRRFIP2_A_HsaEX0036760/LRRFIP2_A_HsaEX0036760_LRRFIP2_A_HsaEX0036760::chr3:37065875-37066775_-__profile.txt", sep = "\t")
gabbr1 <- read.delim("combined_5NIA/GABBR1_A_HsaEX0026760/GABBR1_A_HsaEX0026760_GABBR1_A_HsaEX0026760::chr6:29608929-29609829_-__profile.txt", sep = "\t")
phf21a <- read.delim("combined_5NIA/PHF21A_A_HsaEX0047036/PHF21A_A_HsaEX0047036_PHF21A_A_HsaEX0047036::chr11:45945553-45946453_-__profile.txt", sep = "\t")

tnrc6a <- clamp_reactivity(tnrc6a)
arhgef7 <- clamp_reactivity(arhgef7)
rps6kb2 <- clamp_reactivity(rps6kb2)
cers5 <- clamp_reactivity(cers5)
lrrfip2 <- clamp_reactivity(lrrfip2)
gabbr1 <- clamp_reactivity(gabbr1)
phf21a <- clamp_reactivity(phf21a)

tnrc6a_class <- assign_class(tnrc6a,
                             up_range = c(309:408),
                             bp3ss_range = c(409:492),
                             microexon_range = c(493:639),
                             down_range = c(640:739))
tnrc6a_class$gene <- "TNRC6A"
arhgef7_class <- assign_class(arhgef7,
                              up_range = c(350:449),
                              bp3ss_range = c(450:492),
                              microexon_range = c(493:669),
                              down_range = c(670:769))
arhgef7_class$gene <- "ARHGEF7"
rps6kb2_class <- assign_class(rps6kb2,
                              up_range = c(359:458),
                              bp3ss_range = c(459:492),
                              microexon_range = c(493:629),
                              down_range = c(630:729))
rps6kb2_class$gene <- "RPS6KB2"
cers5_class <- assign_class(cers5,
                            up_range = c(325:424),
                            bp3ss_range = c(425:493),
                            microexon_range = c(494:596),
                            down_range = c(597:696))
cers5_class$gene <- "CERS5"
lrrfip2_class <- assign_class(lrrfip2,
                              up_range = c(342:441),
                              bp3ss_range = c(442:493),
                              microexon_range = c(494:595),
                              down_range = c(596:695))
lrrfip2_class$gene <- "LRRFIP2"
gabbr1_class <- assign_class(gabbr1,
                             up_range = c(333:432),
                             bp3ss_range = c(433:493),
                             microexon_range = c(494:644),
                             down_range = c(645:744))
gabbr1_class$gene <- "GABBR1"
phf21a_class <- assign_class(phf21a,
                             up_range = c(329:428),
                             bp3ss_range = c(429:493),
                             microexon_range = c(494:657),
                             down_range = c(658:757))
phf21a_class$gene <- "PHF21A"
#splice sites
#tnrc6a 493-639, arhgef7 493-669, rps6kb2 493-629, cers5 494-596, lrrfip2 494-595, gabbr1 494-644, phf21a 494-657
#branchpoint (best)
#tnrc6a 409, arhgef7 450, rps6kb2 459, cers5 425, lrrfip2 442, gabbr1 433, phf21a 429

more <- rbind(tnrc6a_class, arhgef7_class, rps6kb2_class, cers5_class, lrrfip2_class, gabbr1_class, phf21a_class)
more %>%
  ggplot(aes(x = class, y = Norm_profile, color = gene)) +
  geom_boxplot() +
  #geom_quasirandom(width = 0.15, size = 1.5, alpha = 0.7) +
  theme_minimal() 

combo$dataset <- "Microexons"
combo <- combo[, c("Nucleotide", "Norm_profile", "gene", "dataset", "class")]
more$dataset <- "Midexons"
more <- more[, c("Nucleotide", "Norm_profile", "gene", "dataset", "class")]
mega <- rbind(combo, more)

mega %>%
  ggplot(aes(x = class, y = Norm_profile, fill = dataset)) +
  geom_boxplot() +
  #geom_quasirandom(width = 0.15, size = 1.5, alpha = 0.7) +
  theme_minimal()

test <- wilcox.test(more[more$class == "upstream", ]$Norm_profile, combo[combo$class == "upstream", ]$Norm_profile) #0.286
test <- wilcox.test(more[more$class == "Branchpoint_2_3ss", ]$Norm_profile, combo[combo$class == "Branchpoint_2_3ss", ]$Norm_profile) #0.000224
test <- wilcox.test(more[more$class == "Microexon", ]$Norm_profile, combo[combo$class == "Microexon", ]$Norm_profile) #0.384
test <- wilcox.test(more[more$class == "downstream", ]$Norm_profile, combo[combo$class == "downstream", ]$Norm_profile) #0.311
