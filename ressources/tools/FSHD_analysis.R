library(dplyr)
library(tidyr)


# R-script for analyzing blastn-out of FSHD-regions (blastn with -outfmt 6) 


titles <- c("read-id", "region", "percent-identity", "region-length", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "e-value", "bitscore")

args <- commandArgs(trailingOnly = TRUE)

if(length(args)==0){
  print("Please include arguments!")
  stop("Requires command line argument.")
}

mat <- read.table(args[1], header = FALSE, col.names = titles, fill=TRUE,stringsAsFactor=FALSE)
pas <- read.table(args[2], header = TRUE, fill=TRUE,stringsAsFactor=FALSE)
setwd(args[3])

##

region <- unique(mat$region)

mat <- mat %>%
  mutate(mat, region = recode(region,
                              "4qB_pLAM." = "B_pLAM_4qB",
                              "10q26_D4Z4." = "c10_D4Z4",
                              "4q35_D4Z4." = "c4_D4Z4",
                              "4qA_pLAM." = "A_pLAM_4qA",
                              "10q26_pLAM." = "c10_pLAM",
                              "4qB_D4F104S1." = "B_D4F104S1_4qB",
                              "4qA_D4F104S1." = "A_D4F104S1_4qA",
                              "10qA_D4F104S1." = "c10_D4F104S1_10qA",
                              "CLUHP4-201_exon1." = "CLUHP4_201_exon1",
                              "DUX4_end" = "DUX4_end",
                              "chr10_ctrl" = "c10_ctrl")) 



mat <- mat[order(mat$read.id, mat$qstart), ]

# relabel B-4qB lower 90% to pLAM_4qB_low_pid
mat <- mat %>%
  mutate(region = if_else(region == "B_pLAM_4qB" & percent.identity < 90,
                          "pLAM_4qB_low_pid", region))




matsub_fin <- mat %>%
  group_by(read.id) %>%
  # identify start positions
  mutate(
    control = ifelse(any(region == "CLUHP4_201_exon1"), qstart[region == "CLUHP4_201_exon1"][1], NA),
    A_start = ifelse(any(region == "A_D4F104S1_4qA" | region == "B_D4F104S1_4qB"), qstart[region == "A_D4F104S1_4qA" | region == "B_D4F104S1_4qB"][1], NA)
  ) %>%
  # delete c10_D4Z4 before CLUHP4 but only if the order CLUHP4 - D4F104S1 is there
  filter(!(region %in% "c10_D4Z4" & 
             !is.na(control) & !is.na(A_start) & control < A_start & qstart < control)) %>%
  # same for reverse read (Reverse-Read)
  filter(!(region %in% "c10_D4Z4" & 
             !is.na(control) & !is.na(A_start) & control > A_start & qstart > control)) %>%
  # label all c4_D4Z4 to c4_ctrl
  mutate(
    region = case_when(
      region == "c4_D4Z4" & 
        !is.na(control) & !is.na(A_start) & control > A_start & qstart > control ~ "c4_ctrl",  
      TRUE ~ region  # all other regions are unchanged
    )
  ) %>%
  mutate(
    region = case_when(
      region == "c4_D4Z4" & 
        !is.na(control) & !is.na(A_start) & control < A_start & qstart < control ~ "c4_ctrl",  
      TRUE ~ region  # all other regions are unchanged
    )
  ) %>%
  ungroup() %>%
  select(-control, -A_start)

expected_offsets <- c(-41983)
tolerance <- 1000

matsub_updated <- matsub_fin %>%
  group_by(read.id) %>%
  mutate(
    D4F104S1_pos = qstart[region %in% c("A_D4F104S1_4qA", "B_D4F104S1_4qB")][1],
    CLUHP4_pos = qstart[region == "CLUHP4_201_exon1"][1],
    read_orientation = ifelse(!is.na(D4F104S1_pos) & !is.na(CLUHP4_pos) &
                                CLUHP4_pos < D4F104S1_pos, "forward",
                              ifelse(CLUHP4_pos > D4F104S1_pos, "reverse", NA_character_)),
  ) %>%
  ungroup() %>%
  group_by(read.id) %>%
  mutate(
    has_ctrl = any(region == "c4_ctrl")
  ) %>%
  rowwise() %>%
  mutate(
    expected_positions = list(
      if (!is.na(read_orientation) && read_orientation == "forward") {
        D4F104S1_pos + expected_offsets
      } else {
        NA_real_
      }
    ),
    should_insert = read_orientation == "forward" &&
      !has_ctrl &&
      !any(is.na(expected_positions)) &&
      all(expected_positions >= (D4F104S1_pos + expected_offsets - tolerance) & expected_positions <= (D4F104S1_pos + expected_offsets) & expected_positions > 0)
  ) %>%
  ungroup()

no_hit_reads <- matsub_updated %>%
  filter(read_orientation == "forward", should_insert) %>%
  distinct(read.id, D4F104S1_pos) %>%
  rowwise() %>%
  mutate(
    qstart = list(D4F104S1_pos + expected_offsets),
    region = list(rep("no_chr4_hit", length(expected_offsets))),
    percent.identity = list(rep(1, length(expected_offsets)))  
  ) %>%
  unnest(c(qstart, region, percent.identity)) %>%
  ungroup()  

matsub_final <- matsub_updated %>%
  bind_rows(no_hit_reads) %>%
  arrange(read.id, qstart) %>%
  group_by(read.id, region) %>%
  filter(
    (region != "c10_ctrl") | (length == max(length))
  ) %>%
  ungroup()


### Count RU

chr4.d4z4count <- matsub_final %>%
  filter(region == "c4_D4Z4" & length > 1500) %>%
  count(region, read.id)
  
# Keep original values for other regions
original_values <- matsub_final %>%
  filter(region != "c10_D4Z4", region != "c4_D4Z4") %>%
  group_by(read.id, region) %>%
  summarise(percent.identity = max(percent.identity)) 

# Add values for region, whereas the count is used
combined <- original_values %>%
  bind_rows(chr4.d4z4count %>% rename(percent.identity = n))

# Fill in table

matneu <- combined %>%
  pivot_wider(names_from = region, values_from = percent.identity) %>%
  rename(RU_count = c4_D4Z4)

# filter reads

plam_check_cols <- intersect(c("pLAM_4qB_low_pid", "B_pLAM_4qB"), colnames(matneu))
missing_check_cols <- intersect(
  c("A_D4F104S1_4qA", "A_pLAM_4qA", "B_D4F104S1_4qB", "c10_D4F104S1_10qA",
    "RU_count", "DUX4_end", "CLUHP4_201_exon1", "c4_ctrl", "c10_ctrl"),
  colnames(matneu)
)

false_4qB <- matneu %>%
  filter(
    if_any(all_of(plam_check_cols), ~ !is.na(.)) &
      if_all(all_of(missing_check_cols), is.na)
  )

false_D4F104S1 <- matneu %>%
  filter(if_all(any_of(c("A_D4F104S1_4qA", "B_D4F104S1_4qB", "c10_D4F104S1_10qA"))) < 90 & if_all(any_of(c("B_pLAM_4qB", "pLAM_4qB_low_pid", "A_pLAM_4qA", "RU_count", "DUX4_end", "CLUHP4_201_exon1", "c4_ctrl")), is.na))

false_c10_ctrl <- matneu %>%
  filter(if_all(any_of(c("c10_ctrl"))) & if_all(any_of(c("B_pLAM_4qB", "A_D4F104S1_4qA", "A_pLAM_4qA", "B_D4F104S1_4qB", "c10_D4F104S1_10qA", "RU_count", "DUX4_end")), is.na))

matfilt <- matneu %>%
  filter(!(read.id %in% false_4qB$read.id)) %>%
  filter(!(read.id %in% false_D4F104S1$read.id))

D4Z4_only <- matneu %>%
  filter(!is.na(RU_count)) %>%
  filter(if_all(any_of(c("A_D4F104S1_4qA", "A_pLAM_4qA", "B_D4F104S1_4qB", "c10_D4F104S1_10qA", "B_pLAM_4qB", "pLAM_4qB_low_pid", "CLUHP4_201_exon1", "DUX4_end", "c10_ctrl")), is.na))

matfilt_final <- matfilt %>%
  filter(!(read.id %in% D4Z4_only$read.id))

##sort matfilt2
cols_to_keep <- c("read.id", "c10_ctrl", "pLAM_4qB_low_pid", "c4_ctrl", "CLUHP4_201_exon1",
                  "A_D4F104S1_4qA", "B_D4F104S1_4qB", "c10_D4F104S1_10qA",
                  "A_pLAM_4qA", "B_pLAM_4qB", "c10_pLAM", "DUX4_end",
                  "RU_count")

for (col in cols_to_keep) {
  if (!col %in% names(matfilt_final)) {
    matfilt_final[[col]] <- NA
  }
}
matfilt_final <- matfilt_final[, cols_to_keep]



# create repeat-order

repeat_order <- matsub_final %>%
  select(read.id, region, percent.identity, qstart, length) %>%
  mutate(type = case_when(
    grepl("c10_D4F104S1_10qA", region) ~ "10qA_D4F104S1",
    grepl("B_D4F104S1_4qB", region) ~ "4qB_D4F104S1",
    grepl("A_D4F104S1_4qA", region) ~ "4qA_D4F104S1",
    grepl("c10_D4Z4", region) ~ "c10",
    grepl("c4_D4Z4", region) ~ "c4",
    grepl("c10_pLAM", region) ~ "10qA_pLAM",
    grepl("A_pLAM_4qA", region) ~ "4qA_pLAM",
    grepl("pLAM_4qB_low_pid", region) ~ "4qB_pLAM_low-pid",
    grepl("B_pLAM_4qB", region) ~ "4qB_pLAM",
    grepl("CLUHP4_201_exon1", region) ~ "CLUHP4-201_exon1",
    grepl("DUX4_end", region) ~ "DUX4_end",
    grepl("c4_ctrl_L", region) ~ "c4_ctrl_L",
    grepl("c4_ctrl_S", region) ~ "c4_ctrl_S",
    grepl("c4_ctrl", region) ~ "c4_ctrl",
    grepl("no_chr4_hit", region) ~ "no_chr4_hit",
    grepl("c10_ctrl", region) ~ "c10_ctrl",
    TRUE ~ NA_character_
  )) %>%
  filter(!(region == "c4_D4Z4" & length <= 1500)) %>%
  filter(!(region == "c10_D4Z4" & length <= 1500)) %>%
  arrange(read.id, qstart) %>%
  group_by(read.id) %>%
  mutate(
    qstart = as.numeric(qstart),
    cluster = cumsum(
      (qstart - lag(qstart, default = -Inf)) > 100
    )
  ) %>%
  group_by(read.id, cluster) %>%
  filter(percent.identity == max(percent.identity, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-cluster) %>%
  arrange(read.id, qstart)

# test if gaps are between c4/c10 and insert gap_size
gaps <- repeat_order %>%
  arrange(read.id, qstart) %>%
  group_by(read.id) %>%
  mutate(
    this_is_D4Z4 = type %in% c("c4", "c10"),
    next_is_D4Z4 = lead(type) %in% c("c4", "c10"),
    next_qstart = lead(qstart),
    gap_size = next_qstart - qstart
  ) %>%
  # only when both neighbouring RU-types c4 and c10 and a gap >5kb
  filter(this_is_D4Z4 & next_is_D4Z4 & gap_size > 5000) %>%
  mutate(
    qstart = qstart + 1,  # gap starts directly after current RU
    region = paste0("gap_", gap_size, "bp"),
    percent.identity = NA,
    length = NA,
    type = region
  ) %>%
  select(read.id, region, percent.identity, qstart, length, type)

# insert gap-regions
repeat_order <- repeat_order %>%
  bind_rows(gaps) %>%
  arrange(read.id, qstart)

# label most distal RU before `_pLAM` with S or L
repeat_order <- repeat_order %>%
  group_by(read.id) %>%
  mutate(
    # Test if previous or the following element has a *_pLAM 
    last_repeat = lead(type) %in% c("4qA_pLAM") |
      lag(type) %in% c("4qA_pLAM"),
    type = if_else(last_repeat & type %in% c("c10", "c4"), 
                   paste0(type, if_else(length < 2500 & length > 1000, "S", "L")), 
                   type)
  ) %>%
  ungroup()

# create repeat-sequence
all_repeat_order <- repeat_order %>%
  filter(!is.na(type)) %>%
  group_by(read.id) %>%
  summarise(repeat_sequence = paste(type, collapse = " - ")) %>%
  ungroup()

# implement repeat-sequence

matfilt_final <- matfilt_final %>%
  left_join(all_repeat_order, by = c("read.id"))
matfilt_final <- matfilt_final %>%
  arrange(desc(RU_count))
matfilt_final <- matfilt_final %>%
  left_join(pas, by = "read.id")

### D4Z4-analysis

# count partial-repeats from reads with repeats only

D4Z4_onlysort <- D4Z4_only[order(D4Z4_only$RU_count, decreasing = TRUE),]
D4Z4_only_count <- table(D4Z4_onlysort$RU_count)

matD <- matsub_final %>%
  filter(read.id %in% D4Z4_onlysort$read.id)


D4Z4_RU <- matD %>%
  select(read.id, region, percent.identity, qstart, length) %>%
  mutate(type = case_when(
    grepl("c10_D4F104S1_10qA", region) ~ "10qA_D4F104S1",
    grepl("B_D4F104S1_4qB", region) ~ "4qB_D4F104S1",
    grepl("A_D4F104S1_4qA", region) ~ "4qA_D4F104S1",
    grepl("c4_D4Z4", region) ~ "c4",
    grepl("c10_D4Z4", region) ~ "c10",
    grepl("c10_pLAM", region) ~ "10qA_pLAM",
    grepl("A_pLAM_4qA", region) ~ "4qA_pLAM",
    grepl("pLAM_4qB_low_pid", region) ~ "4qB_pLAM_low-pid",
    grepl("CLUHP4_201_exon1", region) ~ "CLUHP4-201_exon1",
    grepl("B_pLAM_4qB", region) ~ "4qB_pLAM",
    grepl("DUX4_end", region) ~ "DUX4_end",
    grepl("c4_ctrl", region) ~ "c4_ctrl",
    grepl("no_chr4_hit", region) ~ "no_chr4_hit",
    grepl("c10_ctrl", region) ~ "c10_ctrl",
    TRUE ~ NA_character_
  )) %>%
  filter(!(region == "c4_D4Z4" & length <= 1500)) %>%
  filter(!(region == "c10_D4Z4" & length <= 1500)) %>%
  arrange(read.id, qstart) %>%
  group_by(read.id) %>%
  mutate(
    qstart = as.numeric(qstart),
    cluster = cumsum(
      (qstart - lag(qstart, default = -Inf)) > 100
    )
  ) %>%
  group_by(read.id, cluster) %>%
  filter(percent.identity == max(percent.identity, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-cluster) %>%
  arrange(read.id, qstart)

# create repeat-sequence
D4Z4_RU_order <- D4Z4_RU %>%
  filter(!is.na(type)) %>%
  group_by(read.id) %>%
  summarise(repeat_sequence = paste(type, collapse = " - ")) %>%
  ungroup()

# implement repeat-sequence
D4Z4_onlysort <- as.data.frame(D4Z4_onlysort) %>%
  left_join(D4Z4_RU_order, by = c("read.id"))

# count c4/c10 per read  
D4Z4_onlysort <- D4Z4_onlysort %>%
  mutate(
    c4_count = lengths(regmatches(repeat_sequence, gregexpr("\\bc4\\b", repeat_sequence))),
    c10_count = lengths(regmatches(repeat_sequence, gregexpr("\\bc10\\b", repeat_sequence))),
    chr_assignment = case_when(
      c4_count > c10_count ~ "chr4",
      c10_count > c4_count ~ "chr10",
      c4_count > 0 & c4_count == c10_count ~ "mixed",
      TRUE ~ "none"
    )
  )

# chr4 reads
D4Z4_chr4 <- D4Z4_onlysort %>%
  filter(chr_assignment == "chr4")
D4Z4_chr4_only_count <- table(D4Z4_chr4$RU_count)
D4Z4_chr4_max <- max(D4Z4_chr4$RU_count, na.rm = TRUE)
D4Z4_chr4_max <- sprintf("%s RU", D4Z4_chr4_max)

# chr10 reads
D4Z4_chr10 <- D4Z4_onlysort %>%
  filter(chr_assignment == "chr10")
D4Z4_chr10_only_count <- table(D4Z4_chr10$RU_count)
D4Z4_chr10_max <- max(D4Z4_chr10$RU_count, na.rm = TRUE)
D4Z4_chr10_max <- sprintf("%s RU", D4Z4_chr10_max)

# mixed reads
D4Z4_mixed <- D4Z4_onlysort %>%
  filter(chr_assignment == "mixed")
mixed_max <- max(pmax(D4Z4_mixed$c4_count, D4Z4_mixed$c10_count, na.rm = TRUE))
D4Z4_mixed_max <- sprintf("mixed (max %s RU)", mixed_max)


# create repeat-order

matfilt2 <- data.frame(matfilt_final)
matsub <- matsub_final


## chr10

chr10 <- matfilt2 %>%
  filter(c10_D4F104S1_10qA != "NULL" | (if_all(any_of(c("c10_pLAM")))) != "NULL") %>%
  filter(c10_D4F104S1_10qA > A_D4F104S1_4qA | (if_all(any_of(c("c10_pLAM")))) > (if_all(any_of(c("A_pLAM_4qA")))))  

susp_chr4 <- matfilt2 %>%
  filter(!(if_all(any_of(c("pLAM_4qB_low_pid", "c4_ctrl")), is.na)))
         
chr10 <- chr10 %>%
  filter(!(read.id %in% susp_chr4$read.id)) %>%
  arrange(desc(RU_count))

chr10 <- chr10[order(chr10$RU_count, decreasing = TRUE),]

chr10_blast <- matsub_final %>%
  filter(read.id %in% chr10$read.id) 

chr10_complete <- chr10 %>%
  filter(if_any(c("c10_pLAM", "DUX4_end"), ~ !is.na(.)) &
           if_any(c("c10_D4F104S1_10qA", "CLUHP4_201_exon1"), ~ !is.na(.)))

chr10 <- chr10 %>%
  mutate(status = if_else(read.id %in% chr10_complete$read.id, "complete", "partial")) %>%
  relocate(status, .before = repeat_sequence) 

chr10_complete <- chr10_complete %>%
  mutate(status = "complete") %>%
  relocate(status, .before = repeat_sequence)

## chr4

chr4 <- matfilt2 %>%
  filter(A_D4F104S1_4qA != "NULL" | (if_all(any_of(c("A_pLAM_4qA")))) != "NULL" | if_all(any_of(c("B_pLAM_4qB"))) != "NULL" | if_all(any_of(c("pLAM_4qB_low_pid"))) != "NULL") %>%
  filter(A_D4F104S1_4qA > c10_D4F104S1_10qA | (if_all(any_of(c("A_pLAM_4qA")))) > (if_all(any_of(c("c10_pLAM")))) | if_all(any_of(c("B_pLAM_4qB"))) > 89.9)

chr4 <- bind_rows(chr4, susp_chr4)
chr4 <- chr4 %>%
  distinct(read.id, .keep_all = TRUE)
chr4 <- chr4[order(chr4$RU_count, decreasing = TRUE),]
#chr4$c10_D4Z4 <- NULL


#chr4: Filter the 4qA Haplotypes

A_4qA <- chr4 %>%
  filter((if_all(any_of(c("A_pLAM_4qA")))) > (if_all(any_of(c("c10_pLAM")))) | A_D4F104S1_4qA > pmax(B_D4F104S1_4qB, c10_D4F104S1_10qA, na.rm = TRUE)) %>%
  arrange(desc(RU_count)) %>%
  mutate(
    warning = if_else(
      grepl("c10_ctrl", repeat_sequence) | grepl("no_chr4_hit", repeat_sequence),
      "Chr10 markers detected, possible chr10 read (f.ex. 10A176T/10A180T)!",
      NA_character_
    )
  )

A_4qA <- A_4qA %>%
  mutate(
    S_or_L = case_when(
      grepl("c4S|c10S", repeat_sequence) ~ "S",
      grepl("c4L|c10L", repeat_sequence) ~ "L",
      TRUE ~ NA_character_
    )
  )
A_4qA <- A_4qA %>%
  mutate(S_or_L = ifelse(is.na(S_or_L), "", S_or_L))

A_4qA_blast <- matsub %>%
  filter(read.id %in% A_4qA$read.id)

A_4qA_complete <- A_4qA %>%
  filter(if_any(c("A_pLAM_4qA", "DUX4_end"), ~ !is.na(.)) &
  if_any(c("A_D4F104S1_4qA", "CLUHP4_201_exon1", "pLAM_4qB_low_pid"), ~ !is.na(.)))

A_4qA <- A_4qA %>%
  mutate(status = if_else(read.id %in% A_4qA_complete$read.id, "complete", "partial")) %>%
  relocate(S_or_L, .before = repeat_sequence) %>%
  relocate(status, .after = S_or_L)

A_4qA_complete <- A_4qA_complete %>%
  mutate(status = "complete") %>%
  relocate(S_or_L, .before = repeat_sequence) %>%
  relocate(status, .after = S_or_L)  
    
#chr4: Filter the 4qB Haplotypes


B_4qB <- chr4 %>%
  filter((if_any(c("B_pLAM_4qB"), ~ !is.na(.))) | (if_all(any_of(c("B_pLAM_4qB")))) > (if_all(any_of(c("A_pLAM_4qA")))) | B_D4F104S1_4qB > pmax(A_D4F104S1_4qA, c10_D4F104S1_10qA, na.rm = TRUE) | (if_all(any_of(c("pLAM_4qB_low_pid")))) > (if_all(any_of(c("A_pLAM_4qA"))))) %>%
  arrange(desc(RU_count)) %>%
  mutate(
    warning = if_else(
      grepl("c10_ctrl", repeat_sequence) | grepl("no_chr4_hit", repeat_sequence),
      "Chr 10 markers detected, possible chr 10 read (f.ex. 10B161T)!",
      NA_character_
    )
  )

B_4qB_blast <- matsub %>%
  filter(read.id %in% B_4qB$read.id)

B_4qB_complete <- B_4qB %>%
  filter(if_any(c("B_pLAM_4qB", "DUX4_end"), ~ !is.na(.)) &
  if_any(c("B_D4F104S1_4qB", "CLUHP4_201_exon1", "pLAM_4qB_low_pid"), ~ !is.na(.)))

B_4qB <- B_4qB %>%
  mutate(status = if_else(read.id %in% B_4qB_complete$read.id, "complete", "partial")) %>%
  relocate(status, .before = repeat_sequence) 

B_4qB_complete <- B_4qB_complete %>%
  mutate(status = "complete") %>%
  relocate(status, .before = repeat_sequence)  

chr4_undefined <- chr4 %>%
  filter(!(read.id %in% A_4qA$read.id)) %>%
  filter(!(read.id %in% B_4qB$read.id))

chr4_undefined_complete <- chr4_undefined %>%
  filter(if_any(c("B_pLAM_4qB", "A_pLAM_4qA", "c10_pLAM", "DUX4_end"), ~ !is.na(.)) &
           if_any(c("B_D4F104S1_4qB", "c10_D4F104S1_10qA", "A_D4F104S1_4qA", "CLUHP4_201_exon1", "pLAM_4qB_low_pid"), ~ !is.na(.)))

chr4_undefined <- chr4_undefined %>%
  mutate(status = if_else(read.id %in% chr4_undefined_complete$read.id, "complete", "partial")) %>%
  relocate(status, .before = repeat_sequence) 

chr4_undefined_complete <- chr4_undefined_complete %>%
  mutate(status = "complete") %>%
  relocate(status, .before = repeat_sequence) 



## check for CHIMERIC reads:

chimeric_reads <- bind_rows(A_4qA_complete, chr10_complete, B_4qB_complete, chr4_undefined_complete)

chimeric_reads <- chimeric_reads[
  (!is.na(chimeric_reads$A_D4F104S1_4qA) & !is.na(chimeric_reads$c10_D4F104S1_10qA) & 
     chimeric_reads$A_D4F104S1_4qA != chimeric_reads$c10_D4F104S1_10qA) |
    (!is.na(chimeric_reads$B_D4F104S1_4qB) & !is.na(chimeric_reads$c10_D4F104S1_10qA) & 
       chimeric_reads$B_D4F104S1_4qB != chimeric_reads$c10_D4F104S1_10qA) |
    (!is.na(chimeric_reads$B_D4F104S1_4qB) & !is.na(chimeric_reads$A_D4F104S1_4qA) & 
       chimeric_reads$B_D4F104S1_4qB != chimeric_reads$A_D4F104S1_4qA) 
  , ]

chimeric_reads <- chimeric_reads %>% distinct(read.id, .keep_all = TRUE)

front_reads <- chimeric_reads %>%
  select(read.id, contains("D4F104S1")) %>%
  pivot_longer(cols = c("A_D4F104S1_4qA", "B_D4F104S1_4qB", "c10_D4F104S1_10qA"), 
               names_to = "read.id_type", values_to = "percent.identity") %>%
  mutate(type = case_when(
    grepl("^A", read.id_type) ~ "A",
    grepl("^c10", read.id_type) ~ "c10",
    grepl("^B", read.id_type) ~ "B",
    TRUE ~ NA_character_)) %>%
  group_by(read.id) %>%
  mutate(highest_identity = ifelse(all(is.na(percent.identity)), NA_real_, max(percent.identity, na.rm = TRUE))) %>%
  filter(percent.identity == highest_identity) %>%
  slice_max(order_by = percent.identity, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(read.id, front_type = type)

back_reads <- chimeric_reads %>%
  select(read.id, contains("pLAM")) %>%
  pivot_longer(cols = c("A_pLAM_4qA", "B_pLAM_4qB", "c10_pLAM"), names_to = "read.id_type", values_to = "percent.identity") %>%
  mutate(type = case_when(
    grepl("^A", read.id_type) ~ "A",
    grepl("^c10", read.id_type) ~ "c10",
    grepl("^B", read.id_type) ~ "B",
    TRUE ~ NA_character_
  )) %>%
  group_by(read.id) %>%
  mutate(highest_identity = ifelse(all(is.na(percent.identity)), NA_real_, max(percent.identity, na.rm = TRUE))) %>%
  filter(percent.identity == highest_identity) %>%
  slice_max(order_by = percent.identity, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(read.id, back_type = type)  # Nur relevante Spalten behalten

# Keep only reads where front and back types differ
filtered_reads <- front_reads %>%
  left_join(back_reads, by = "read.id") %>%
  filter(front_type != back_type)

#  Extract the filtered read.id values and keep only those rows from the original chimeric_reads_t table
final_filtered_table <- chimeric_reads %>%
  filter(read.id %in% filtered_reads$read.id)

# Add the front_type and back_type columns to the final filtered table
final_filtered_table <- final_filtered_table %>%
  left_join(filtered_reads %>% select(read.id, front_type, back_type), by = "read.id")

final_filtered_table <- final_filtered_table %>%
  rename(
    D4F104S1_type = front_type,
    pLAM_type = back_type
  )

chimeric_reads <- final_filtered_table %>% distinct(read.id, .keep_all = TRUE)

chimeric_reads <- chimeric_reads %>% 
  mutate(
    D4F104S1_type = case_when(
      D4F104S1_type == "c10" ~ "10qA",
      D4F104S1_type == "A" ~ "4qA",
      D4F104S1_type == "B" ~ "4qB",
      TRUE ~ D4F104S1_type  # Keep other values unchanged
    ),
    pLAM_type = case_when(
      pLAM_type == "c10" ~ "10qA",
      pLAM_type == "A" ~ "4qA",
      pLAM_type == "B" ~ "4qB",
      TRUE ~ pLAM_type  # Keep other values unchanged
    )
  )

chimeric_reads <- chimeric_reads %>%
  mutate(
    count_c4 = lengths(regmatches(repeat_sequence, gregexpr("\\bc4(S|L)?\\b", repeat_sequence))),
    count_c10 = lengths(regmatches(repeat_sequence, gregexpr("\\bc10\\b", repeat_sequence))),
    warning = case_when(
      # Original rule for chr4-Read
      D4F104S1_type == "10qA" & 
        (!is.na(c4_ctrl) | !is.na(pLAM_4qB_low_pid)) ~ 
        "very likely chr4-read and not chimeric as chr4_identifier (pLAM_4qB_low_pid and/or c4_ctrl) were detected.",
      
      D4F104S1_type %in% c("4qA", "4qB") & 
        pLAM_type == "10qA" & 
        grepl("no_chr4_hit", repeat_sequence) & grepl("c10_ctrl", repeat_sequence) ~ 
        "very likely chr10-read and not chimeric as chr10-identifier (c10_ctrl) and lack of chr4 hit (no_chr4_hit) were both detected.",
      
      # Original rule for chr10-Read
      D4F104S1_type %in% c("4qA", "4qB") & 
        pLAM_type == "10qA" & 
        grepl("c10_ctrl", repeat_sequence) ~ 
        "very likely chr10-read and not chimeric as a chr10-identifier was detected.",
      
      # Original rule for chr10-Read
      D4F104S1_type %in% c("4qA", "4qB") & 
        pLAM_type == "10qA" & 
        grepl("no_chr4_hit", repeat_sequence) ~ 
        "possibly chr10-read and not chimeric as no chr4-identifier was detected on the position (with 1000bp tolerance), marked by no_chr4_hit",
      
      # Wenn D4F104S1_type == 10qA und nur c4-Repeats vorkommen
      D4F104S1_type == "10qA" &
        count_c4 > 0 & count_c10 == 0 ~ 
        "possibly chr4-read and not chimeric as only chr4_D4Z4 were detected. Check for possible 10A/BT haplotype (non-standard-chromosome.",
      
      # Wenn D4F104S1_type == 4qA oder 4qB und nur c10-Repeats vorkommen
      D4F104S1_type %in% c("4qA", "4qB") &
        count_c10 > 0 & count_c4 == 0 ~ 
        "possibly chr10-read and not chimeric as only chr10_D4Z4 were detected.",
      
      TRUE ~ NA_character_
    )
  )

chimeric_reads_type <- chimeric_reads %>%
  select(read.id, D4F104S1_type, pLAM_type) %>%
  distinct()

chimeric_blast <- matsub %>%
  filter(read.id %in% chimeric_reads_type$read.id)

chimeric_repeat_order <- chimeric_reads %>%
  select("read.id", "D4F104S1_type", "pLAM_type", "RU_count", "repeat_sequence", "warning")


# filter out chimeric-reads from 4qA and chr10-reads 

A_4qA_complete <- A_4qA_complete %>%
  anti_join(chimeric_reads, by = "read.id")

chr10_complete <- chr10_complete %>%
  anti_join(chimeric_reads, by = "read.id")

A_4qA <- A_4qA %>%
  anti_join(chimeric_reads, by = "read.id")

chr10 <- chr10 %>%
  anti_join(chimeric_reads, by = "read.id")

B_4qB <- B_4qB %>%
  anti_join(chimeric_reads, by = "read.id")

B_4qB_complete <- B_4qB_complete %>%
  anti_join(chimeric_reads, by = "read.id")

chr4_undefined <- chr4_undefined %>%
  anti_join(chimeric_reads, by = "read.id")

# extract read-ids from 4qB to a .txt

B_4qB_ID <- data.frame(B_4qB$read.id)
write.table(B_4qB_ID, file = "4qB_all-reads-ID.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

B_4qB_ID_complete <- data.frame(B_4qB_complete$read.id)
write.table(B_4qB_ID_complete, file = "4qB_complete-reads-ID.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# extract read-ids from 4qA to a .txt

A_4qA_ID <- data.frame(A_4qA$read.id)
write.table(A_4qA_ID, file = "4qA_all-reads-ID.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

A_4qA_ID_complete <- data.frame(A_4qA_complete$read.id)
write.table(A_4qA_ID_complete, file = "4qA_complete-reads-ID.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# extract read-ids from chimeric_reads

chimeric_reads_ID <- data.frame(chimeric_reads$read.id)
write.table(chimeric_reads_ID, file = "chimeric-reads-ID.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# write D4Z4-Infos to table
D4Z4_chr4_ID <- data.frame(D4Z4_chr4$read.id)
write.table(D4Z4_chr4_ID, file = "D4Z4-only_chr4-reads-ID.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

D4Z4_chr10_ID <- data.frame(D4Z4_chr10$read.id)
write.table(D4Z4_chr10_ID, file = "D4Z4-only_chr10-reads-ID.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

chr4_undefined_ID <- data.frame(chr4_undefined$read.id)
write.table(chr4_undefined_ID, file = "chr4-undefined_all-reads-ID.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

chr10_ID <- data.frame(chr10$read.id)
write.table(chr10_ID, file = "chr10_all-reads-ID.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

chr10_ID_complete <- data.frame(chr10_complete$read.id)
write.table(chr10_ID_complete, file = "chr10_complete-reads-ID.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# create statistic-table 

B_4qB_partial <- B_4qB %>%
  filter(!(read.id %in% B_4qB_complete$read.id))

A_4qA_partial <- A_4qA %>%
  filter(!(read.id %in% A_4qA_complete$read.id))

chr10_total_reads <- as.numeric(count(chr10) + count(D4Z4_chr10))
chr10_reads_HP_marker <- as.numeric(count(chr10))
chr10_complete_reads <- as.numeric(count(chr10_complete))
q10A_repeats <- chr10_complete$RU_count
q10A_repeats <- sort(q10A_repeats)

chr4_total_reads <- as.numeric(count(A_4qA) + count(B_4qB) + count(chr4_undefined) + count(D4Z4_chr4))
chr4_reads_HP_marker <- as.numeric(count(A_4qA) + count(B_4qB))
chr4_undefined_reads <- as.numeric(count(chr4_undefined))

q4A_reads <- as.numeric(count(A_4qA))
q4A_partial_repeats <- A_4qA_partial$RU_count
q4A_partial_repeats[is.na(q4A_partial_repeats)] <- 0
q4A_partial_repeats <- as.numeric(sort(q4A_partial_repeats))
q4A_complete_reads <- as.numeric(count(A_4qA_complete))
q4A_repeats <- A_4qA_complete$RU_count
q4A_repeats <- sort(q4A_repeats)

q4B_reads <- as.numeric(count(B_4qB))
q4B_partial_repeats <- B_4qB_partial$RU_count
q4B_partial_repeats[is.na(q4B_partial_repeats)] <- 0
q4B_partial_repeats <- as.numeric(sort(q4B_partial_repeats))
q4B_complete_reads <- as.numeric(count(B_4qB_complete))
q4B_repeats <- B_4qB_complete$RU_count
q4B_repeats <- sort(q4B_repeats)
chr4_partial_repeats <- D4Z4_only_count

chimeric_total_reads <- as.numeric(count(chimeric_reads))
chimeric_reads_repeats <- chimeric_reads$RU_count
chimeric_reads_repeats <- sort(chimeric_reads_repeats)

D4Z4_chr4_read_count <- as.numeric(count(D4Z4_chr4))
D4Z4_chr10_read_count <- as.numeric(count(D4Z4_chr10))

# group infos
# read-counts_chr10
read_types_chr10 <- c("chr10_total_reads", "chr10_reads_HP_marker", "chr10_complete_reads", "10qA-complete-RU", "chr10_reads_with_D4Z4_only", "chr10_reads_D4Z4_only longest RU")
count_chr10 <- c(chr10_total_reads, chr10_reads_HP_marker, chr10_complete_reads, paste(q10A_repeats, collapse = ", "), D4Z4_chr10_read_count, D4Z4_chr10_max)

#chr4
read_types_chr4 <- c("chr4_total_reads", "chr4_reads_HP_marker", "chr4_undefined_reads", "chr4_reads_with_D4Z4_only", "chr4_reads_D4Z4_only longest RU")
count_chr4 <- c(chr4_total_reads, chr4_reads_HP_marker, chr4_undefined_reads, D4Z4_chr4_read_count, D4Z4_chr4_max)

# 4qA-repeats
reads_4qA <- c("4qA-all-reads", "4qA-partial-RU", "4qA-complete-reads", "4qA-complete-RU")
count_4qA <- c(q4A_reads,  paste(q4A_partial_repeats, collapse = ", "), q4A_complete_reads, 
           paste(q4A_repeats, collapse = ", "))

# 4qB-repeats
reads_4qB <- c("4qB-all-reads", "4qB-partial-RU", "4qB-complete-reads", "4qB-complete-RU")
count_4qB <- c(q4B_reads, paste(q4B_partial_repeats, collapse = ", "), q4B_complete_reads, 
           paste(q4B_repeats, collapse = ", "))

# chimeric-repeats
reads_chim <- c("chimeric-reads", "chimeric-RU")
count_chim <- c(chimeric_total_reads, paste(chimeric_reads_repeats, collapse = ", "))

# write statistics to csv.files

stats_10 <- data.frame(read_types_chr10, count_chr10, stringsAsFactors = FALSE)
stats_4 <- data.frame(read_types_chr4, count_chr4, stringsAsFactors = FALSE)
stats_A <- data.frame(reads_4qA, count_4qA, stringsAsFactors = FALSE)
stats_B <- data.frame(reads_4qB, count_4qB, stringsAsFactors = FALSE)
stats_chim <- data.frame(reads_chim, count_chim, stringsAsFactors = FALSE)

if (NROW(A_4qA) > 0) {write.csv2(A_4qA, file = "A_4qA_all-reads.csv", row.names = FALSE)}
if (NROW(B_4qB) > 0) {write.csv2(B_4qB, file = "B_4qB_all-reads.csv", row.names = FALSE)}
if (NROW(B_4qB_complete) > 0) {write.csv2(B_4qB_complete, file = "B_4qB_complete-reads.csv", row.names = FALSE)}
if (NROW(A_4qA_complete) > 0) {write.csv2(A_4qA_complete, file = "A_4qA_complete-reads.csv", row.names = FALSE)}
if (NROW(chimeric_reads) > 0) {write.csv2(chimeric_reads, file = "chimeric_reads.csv", row.names = FALSE)}
if (NROW(chr10_complete) > 0) {write.csv2(chr10_complete, file = "chr10_complete-reads.csv", row.names = FALSE)}
if (NROW(chr10) > 0) {write.csv2(chr10, file = "chr10_all-reads.csv", row.names = FALSE)}
if (NROW(D4Z4_onlysort) > 0) {write.csv2(D4Z4_onlysort, file = "D4Z4-only_all-reads.csv", row.names = FALSE)}
if (NROW(chr4_undefined) > 0) {write.csv2(chr4_undefined, file = "chr4_undefined_all-reads.csv", row.names = FALSE)}

D4Z4_chr4_only_count <- as.data.frame(D4Z4_chr4_only_count)
if (ncol(D4Z4_chr4_only_count) == 2) {
  colnames(D4Z4_chr4_only_count) <- c("Count of partial RU.", "Number of reads.")
} else {
  D4Z4_chr4_only_count <- data.frame(
    `Count of partial RU.` = character(0),
    `Number of reads.` = integer(0)
  )
}

D4Z4_chr10_only_count <- as.data.frame(D4Z4_chr10_only_count)
if (ncol(D4Z4_chr10_only_count) == 2) {
  colnames(D4Z4_chr10_only_count) <- c("Count of partial RU.", "Number of reads.")
} else {
  D4Z4_chr10_only_count <- data.frame(
    `Count of partial RU.` = character(0),
    `Number of reads.` = integer(0)
  )
}


sink("FSHD_overview-statistics.csv")

cat('FSHD-Analysis - overall statistics')
cat('\n')
cat('\n')
cat('This statistics gives an overview of the read-types and counts of reads and repeat-units (RU).')
cat('\n')
cat('The D4Z4-repeat-array is distinguishable for the chromosome and haplotye (4qA/4qB/10qA) at its proximal (D4F104S1 in DBET) and distal end (pLAM in DUX4 gene-body).')
cat('\n')
cat('Reads with HP marker: reads which have the regions incorporated for distinguishing the chromosome and haplotype.')
cat('\n')
cat('The most distal RU before pLAM is checked if it belongs to the short (S) or long (L) haplotype for 4qA.')
cat('\n')
cat('The workflow incorporates markers for chromosome 4, which further helps to distinguish the reads.')
cat('\n')
cat('\n')
cat('chr4-markers: pLAM_4qB_low_pid followed by two c4_ctrl (DUX4L9-201) hits upstream of the proximal start D4F104S1 of the D4Z4-array are hits only to be found on chromosome 4.')
cat('\n')
cat('chr10-markers: c10_ctrl hit upstream of the proximal start D4F104S1 of the D4Z4-array is only to be found on chromosome 10.')
cat('\n')
cat('The workflow also tests if the chr4-markers at the proximal end are present regarding the length of the read and gives a "no-chr4-hit" if not.')
cat('\n')
cat('complete reads: reads which include a complete D4Z4-repeat-array with proximal and distal end. The reads completely overlap the whole array.')
cat('\n')
cat('partial reads: reads which include a partial repeat-array with either the proximal or distal end. Therefore the reads only reach into the array either from the proximal or distal side.')
cat('\n')
cat('\n')
cat('\n')
cat('____________________________')
cat('\n')
cat('Statistics - Chromosome 4:')
cat('\n')
write.csv2(stats_4, row.names = FALSE)
cat('\n')
cat('Statistics - 4qA-Haplotype:')
cat('\n')
write.csv2(stats_A, row.names = FALSE)
cat('\n')
cat('Statistics - 4qB-Haplotype:')
cat('\n')
write.csv2(stats_B, row.names = FALSE)
cat('\n')
cat('____________________________')
cat('\n')
cat('Statistics - Chimeric reads (possible sub-haplotypes or (in rare events) translocations):')
cat('\n')
write.csv2(stats_chim, row.names = FALSE)
cat('\n')
cat('____________________________')
cat('\n')
cat('Statistics - Chromosome 10:')
cat('\n')
write.csv2(stats_10, row.names = FALSE)
cat('\n')
cat('____________________________')

cat('\n')
cat('\n')
cat('Count of reads from chromosome 4 which only contain D4Z4-repeats.')
cat('\n')
write.csv2(D4Z4_chr4_only_count, row.names = FALSE)
cat('\n')
cat('____________________________')

cat('\n')
cat('\n')
cat('Count of reads from chromosome 10 which only contain D4Z4-repeats.')
cat('\n')
write.csv2(D4Z4_chr10_only_count, row.names = FALSE)
cat('\n')
cat('____________________________')

cat('\n')
cat('\n')
cat('Chimeric reads:')
cat('\n')
cat('\n')
cat('Shows haplotype for proximal (D4F104S1) and distal (pLAM) end of the D4Z4-array of the chimeric reads.')
cat('\n')
cat('Attention: There are sub-haplotypes for 4qA/4qB which are nearly not distinguishable in their proximal-front (D4F104S1) region from chr10.')
cat('\n')
cat('The tool recognizes chr4-markers applied into the blast-workflow and gives warnings if chimeric reads might belong to a rare haplotype.')
cat('\n')
cat('\n')
write.csv2(chimeric_repeat_order, row.names = FALSE)
cat('\n')
cat('\n')
cat('ATTENTION: Please further investigate ambiguous reads to confirm the correct chromosome and allele (f.ex. SSLP-polymorphism)!')
cat('\n')
cat('4qA-sub-haplotypes not distinguishable in their proximal-front (D4F104S1) region from chr10: 4A166 (not permissive), 4A166H, 4A168, but distinguishable from the distal 4qA-pLAM.')
cat('\n')
cat('4qB-sub-haplotypes not distinguishable in their proximal-front (D4F104S1) region from chr10: 4B162, 4B168, 4B170H. But well distinguishable from the distal 4qB-pLAM.')
cat('\n')
cat('\n')
cat('\n')
cat('____________________________')
cat('\n')
cat('\n')
cat('DUCKS4 - Version 2.1.0')

sink()

# give blast-output

output_dir <- "blast_results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

chimeric_blast <- mat %>% filter(mat$read.id %in% chimeric_reads_ID$chimeric_reads.read.id)
write.csv2(chimeric_blast, file = file.path(output_dir, "chimeric_reads_blast.csv"), row.names = FALSE)

Aq4_blast <- mat %>% filter(mat$read.id %in% A_4qA_ID$A_4qA.read.id)
write.csv2(Aq4_blast, file = file.path(output_dir, "4qA_reads_blast.csv"), row.names = FALSE)

Bq4_blast <- mat %>% filter(mat$read.id %in% B_4qB_ID$B_4qB.read.id)
write.csv2(Bq4_blast, file = file.path(output_dir, "4qB_reads_blast.csv"), row.names = FALSE)

chr10_blast <- mat %>% filter(mat$read.id %in% chr10_ID$chr10.read.id)
write.csv2(chr10_blast, file = file.path(output_dir, "chr10_reads_blast.csv"), row.names = FALSE)

chr4_undefined_blast <- mat %>% filter(mat$read.id %in% chr4_undefined_ID$chr4_undefined.read.id)
write.csv2(chr4_undefined_blast, file = file.path(output_dir, "chr4_undefined_reads_blast.csv"), row.names = FALSE)


blast <- basename(args[1])
target <- file.path(output_dir, basename(blast))
file.rename(blast, target)

pas <- basename(args[2])
target <- file.path(output_dir, basename(pas))
file.rename(pas, target)


#################################################################################################################################
####                 FINISH :)             #####
#################################################################################################################################



