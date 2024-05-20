#generate .dat file for SAINT score

library(tidyverse)


#first import txt proteinGroup file

#inter.dat
a <- proteinGroups |>
  filter(Peptides > 1,
         Unique.peptides > 0,
         !grepl("CON", Protein.IDs),
         !grepl("REV", Protein.IDs)) |>
  mutate(Gene = vapply(strsplit(Gene.names, ";"), "[", 1, FUN.VALUE = character(1)),
         Protein = vapply(strsplit(Majority.protein.IDs, ";"), "[", 1, FUN.VALUE = character(1)),
         ID = paste0(Gene, "_", Protein)) |>
  dplyr::select(-Gene.names, -Majority.protein.IDs) |>
  dplyr::select(contains("LFQ"), ID)

a <- a |>
  mutate(Protein = vapply(strsplit(ID, "_"), "[", 2, FUN.VALUE = character(1))) |>
  dplyr::select(contains("LFQ"), Protein)

colnames(a) <- gsub("LFQ.intensity.", "", colnames(a))

a <- a |> pivot_longer(!Protein) |> 
  mutate(condition = case_when(
    grepl("P1", name) ~ "P1",
    grepl("P2", name) ~ "P2",
    TRUE ~ "CTRL"))

int <- data.frame(
  a$name, a$condition, a$Protein, a$value
)

int |> write.table(file = "inter.txt", 
                   col.names = F, row.names = F, sep = "\t", quote = F)


#prey.dat
b <- proteinGroups |>
  filter(Peptides > 1,
         Unique.peptides > 0,
         !grepl("CON", Protein.IDs),
         !grepl("REV", Protein.IDs))|>
  #left_join(length_info, by = "ID") |>
  mutate(Gene = vapply(strsplit(Gene.names, ";"), "[", 1, FUN.VALUE = character(1)),
         Protein = vapply(strsplit(Majority.protein.IDs, ";"), "[", 1, FUN.VALUE = character(1))) |>
  dplyr::select(Protein, Sequence.length, Gene)

b |> write.table(file = "prey.txt", 
                 col.names = F, row.names = F, sep = "\t", quote = F)


#bait.dat
c <- data.frame(
  "IP" = unique(a$name)
) |>
  mutate(Bait = vapply(strsplit(IP, "_"), "[", 1, FUN.VALUE = character(1)),
         Control = if_else(Bait == "WT", "C", "T"))

c |> write.table(file = "bait.txt", 
                 col.names = F, row.names = F, sep = "\t", quote = F)
