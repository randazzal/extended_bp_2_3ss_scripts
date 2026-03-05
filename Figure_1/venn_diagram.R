library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(ggeasy)
library(VennDiagram)

## psi changes and UGC sequence comparisons
psi <- read.delim("/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/all_psi_mandf.txt", sep = "\t")
avg_noUGC <- read.delim("/data2/lackey_lab/randazza/microexons/motifs/UGC_motif_chick/DIY/gga_avg_psi_noUGC.txt", sep = "\t")
avg_noUGC[, 4] <- NULL
avg_UGC <- read.delim("/data2/lackey_lab/randazza/microexons/motifs/UGC_motif_chick/DIY/gga_avg_psi_1plusUGC.txt", sep = "\t")
avg_UGC[, 4] <- NULL
micro_noUGC <- read.delim("/data2/lackey_lab/randazza/microexons/motifs/UGC_motif_chick/DIY/gga_micro_psi_noUGC.txt", sep = "\t")
micro_noUGC[, 4] <- NULL
micro_UGC <- read.delim("/data2/lackey_lab/randazza/microexons/motifs/UGC_motif_chick/DIY/gga_micro_psi_1plusUGC.txt", sep = "\t")
micro_UGC[, 4] <- NULL

micro_noUGC <- micro_noUGC %>%
  separate(EVENT, into = c("GeneID", "exonStart_0base", "exonEnd", "length"), sep = "_", remove = FALSE)
micro_UGC <- micro_UGC %>%
  separate(EVENT, into = c("GeneID", "exonStart_0base", "exonEnd", "length"), sep = "_", remove = FALSE)
avg_noUGC <- avg_noUGC %>%
  separate(EVENT, into = c("GeneID", "exonStart_0base", "exonEnd", "length"), sep = "_", remove = FALSE)
avg_UGC <- avg_UGC %>%
  separate(EVENT, into = c("GeneID", "exonStart_0base", "exonEnd", "length"), sep = "_", remove = FALSE)

## start with 0.3 difference
little <- read.delim("/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/HH15_HH27_changes/change_less03_lowexp_removed.txt", sep = "\t")
little$length <- little$exonEnd - little$exonStart_0base
little$name <- paste(little$GeneID, little$exonStart_0base, little$exonEnd, little$length, sep = "_")
m_little <- little[little$length <= 27, ]
a_little <- little[little$length > 100 & little$length <= 200, ]

nomicro_l1 <- merge(micro_noUGC, m_little, by.x = "EVENT", by.y = "name")
yepmicro_l1 <- merge(micro_UGC, m_little, by.x = "EVENT", by.y = "name")
noavg_l1 <- merge(avg_noUGC, a_little, by.x = "EVENT", by.y = "name")
yepavg_l1 <- merge(avg_UGC, a_little, by.x = "EVENT", by.y = "name")

drastic <- read.delim("/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/HH15_HH27_changes/change_more03_lowexp_removes.txt", sep = "\t")
drastic$length <- drastic$exonEnd - drastic$exonStart_0base
drastic$name <- paste(drastic$GeneID, drastic$exonStart_0base, drastic$exonEnd, drastic$length, sep = "_")
m_drastic <- drastic[drastic$length <= 27, ]
a_drastic <- drastic[drastic$length > 100 & drastic$length <= 200, ]

nomicro_d1 <- merge(micro_noUGC, m_drastic, by.x = "EVENT", by.y = "name")
yepmicro_d1 <- merge(micro_UGC, m_drastic, by.x = "EVENT", by.y = "name")
noavg_d1 <- merge(avg_noUGC, a_drastic, by.x = "EVENT", by.y = "name")
yepavg_d1 <- merge(avg_UGC, a_drastic, by.x = "EVENT", by.y = "name")

temp <- rbind(micro_noUGC, micro_UGC)
change_only <- anti_join(m_drastic, temp, by = c("name" = "EVENT"))
m_drastic <- m_drastic %>%
  filter(!name %in% change_only$name)
change_only <- anti_join(m_little, temp, by = c("name" = "EVENT"))
m_little <- m_little %>%
  filter(!name %in% change_only$name)
change_only <- anti_join(m_nochange, temp, by = c("name" = "EVENT"))
m_nochange <- m_nochange %>%
  filter(!name %in% change_only$name)

list_events <- list(UGC = micro_UGC$EVENT,
                    Drastic= m_drastic$name)
venn.plot <- venn.diagram(
  x = list_events,
  filename = NULL, # Use NULL to draw the diagram directly on R's plotting device
  fill = c("#f1a340", "#998ec3"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  na = "none",
  main = "Microexons with UGC and Big PSI Change"
)
grid::grid.draw(venn.plot)

grid::grid.newpage()

list_events <- list(UGC = micro_noUGC$EVENT,
                    Drastic= m_drastic$name)
venn.plot <- venn.diagram(
  x = list_events,
  filename = NULL, # Use NULL to draw the diagram directly on R's plotting device
  fill = c("#f1a340", "#998ec3"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  na = "none",
  main = "Microexons without UGC and Big PSI Change"
)
grid::grid.draw(venn.plot)

grid::grid.newpage()
