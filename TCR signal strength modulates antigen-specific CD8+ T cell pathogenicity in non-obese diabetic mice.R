library(magrittr)
library("readr")
library(readxl)
library(writexl)
library("rTCRBCRr")
library(tidyr)
library(dplyr)
library(ggplot2)

CD5hi_1_TRAV <- read_excel("CD5hi_CD5lo_TCR_Repertoires_TRAV_TRBV_CDR3.xlsx", sheet = "CD5hi_1_TRAV")
CD5hi_1_TRBV <- read_excel("CD5hi_CD5lo_TCR_Repertoires_TRAV_TRBV_CDR3.xlsx", sheet = "CD5hi_1_TRBV")

CD5hi_2_TRAV <- read_excel("CD5hi_CD5lo_TCR_Repertoires_TRAV_TRBV_CDR3.xlsx", sheet = "CD5hi_2_TRAV")
CD5hi_2_TRBV <- read_excel("CD5hi_CD5lo_TCR_Repertoires_TRAV_TRBV_CDR3.xlsx", sheet = "CD5hi_2_TRBV")

CD5lo_1_TRAV <- read_excel("CD5hi_CD5lo_TCR_Repertoires_TRAV_TRBV_CDR3.xlsx", sheet = "CD5lo_1_TRAV")
CD5lo_1_TRBV <- read_excel("CD5hi_CD5lo_TCR_Repertoires_TRAV_TRBV_CDR3.xlsx", sheet = "CD5lo_1_TRBV")

CD5lo_2_TRAV <- read_excel("CD5hi_CD5lo_TCR_Repertoires_TRAV_TRBV_CDR3.xlsx", sheet = "CD5lo_2_TRAV")
CD5lo_2_TRBV <- read_excel("CD5hi_CD5lo_TCR_Repertoires_TRAV_TRBV_CDR3.xlsx", sheet = "CD5lo_2_TRBV")

### Fig.7B-E ###
## TRA ## 
# Sum up the count for each gene in the V column
CD5hi_1_TRAV_counts <- CD5hi_1_TRAV %>%
  group_by(V) %>%
  summarise(total_count = sum(Count))
CD5hi_2_TRAV_counts <- CD5hi_2_TRAV %>%
  group_by(V) %>%
  summarise(total_count = sum(Count))
CD5lo_1_TRAV_counts <- CD5lo_1_TRAV %>%
  group_by(V) %>%
  summarise(total_count = sum(Count))
CD5lo_2_TRAV_counts <- CD5lo_2_TRAV %>%
  group_by(V) %>%
  summarise(total_count = sum(Count))
# Sum up the count for each gene in the J column
CD5hi_1_TRAJ_counts <- CD5hi_1_TRAV %>%
  group_by(J) %>%
  summarise(total_count = sum(Count))
CD5hi_2_TRAJ_counts <- CD5hi_2_TRAV %>%
  group_by(J) %>%
  summarise(total_count = sum(Count))
CD5lo_1_TRAJ_counts <- CD5lo_1_TRAV %>%
  group_by(J) %>%
  summarise(total_count = sum(Count))
CD5lo_2_TRAJ_counts <- CD5lo_2_TRAV %>%
  group_by(J) %>%
  summarise(total_count = sum(Count))

## TRB ## 
# Sum up the count for each gene in the V column
CD5hi_1_TRBV_counts <- CD5hi_1_TRBV %>%
  group_by(V) %>%
  summarise(total_count = sum(Count))
CD5hi_2_TRBV_counts <- CD5hi_2_TRBV %>%
  group_by(V) %>%
  summarise(total_count = sum(Count))
CD5lo_1_TRBV_counts <- CD5lo_1_TRBV %>%
  group_by(V) %>%
  summarise(total_count = sum(Count))
CD5lo_2_TRBV_counts <- CD5lo_2_TRBV %>%
  group_by(V) %>%
  summarise(total_count = sum(Count))
# Sum up the count for each gene in the J column
CD5hi_1_TRBJ_counts <- CD5hi_1_TRBV %>%
  group_by(J) %>%
  summarise(total_count = sum(Count))
CD5hi_2_TRBJ_counts <- CD5hi_2_TRBV %>%
  group_by(J) %>%
  summarise(total_count = sum(Count))
CD5lo_1_TRBJ_counts <- CD5lo_1_TRBV %>%
  group_by(J) %>%
  summarise(total_count = sum(Count))
CD5lo_2_TRBJ_counts <- CD5lo_2_TRBV %>%
  group_by(J) %>%
  summarise(total_count = sum(Count))

# Get the intersect of gene sets
CD5hi_TRAV_common_genes <- intersect(CD5hi_1_TRAV_counts$V, CD5hi_2_TRAV_counts$V)
CD5lo_TRAV_common_genes <- intersect(CD5lo_1_TRAV_counts$V, CD5lo_2_TRAV_counts$V)
TRAV_common_genes <- intersect(CD5hi_TRAV_common_genes, CD5lo_TRAV_common_genes)
CD5hi_TRAJ_common_genes <- intersect(CD5hi_1_TRAJ_counts$J, CD5hi_2_TRAJ_counts$J)
CD5lo_TRAJ_common_genes <- intersect(CD5lo_1_TRAJ_counts$J, CD5lo_2_TRAJ_counts$J)
TRAJ_common_genes <- intersect(CD5hi_TRAJ_common_genes, CD5lo_TRAJ_common_genes)

CD5hi_TRBV_common_genes <- intersect(CD5hi_1_TRBV_counts$V, CD5hi_2_TRBV_counts$V)
CD5lo_TRBV_common_genes <- intersect(CD5lo_1_TRBV_counts$V, CD5lo_2_TRBV_counts$V)
TRBV_common_genes <- intersect(CD5hi_TRBV_common_genes, CD5lo_TRBV_common_genes)
CD5hi_TRBJ_common_genes <- intersect(CD5hi_1_TRBJ_counts$J, CD5hi_2_TRBJ_counts$J)
CD5lo_TRBJ_common_genes <- intersect(CD5lo_1_TRBJ_counts$J, CD5lo_2_TRBJ_counts$J)
TRBJ_common_genes <- intersect(CD5hi_TRBJ_common_genes, CD5lo_TRBJ_common_genes)

# Reshape the data into long format and assign correct sample names
TRAV_all_counts <- bind_rows(
  CD5hi_1_TRAV_counts %>% mutate(Sample = "CD5hi_1_TRAV"),
  CD5hi_2_TRAV_counts %>% mutate(Sample = "CD5hi_2_TRAV"),
  CD5lo_1_TRAV_counts %>% mutate(Sample = "CD5lo_1_TRAV"),
  CD5lo_2_TRAV_counts %>% mutate(Sample = "CD5lo_2_TRAV")
)
TRBV_all_counts <- bind_rows(
  CD5hi_1_TRBV_counts %>% mutate(Sample = "CD5hi_1_TRBV"),
  CD5hi_2_TRBV_counts %>% mutate(Sample = "CD5hi_2_TRBV"),
  CD5lo_1_TRBV_counts %>% mutate(Sample = "CD5lo_1_TRBV"),
  CD5lo_2_TRBV_counts %>% mutate(Sample = "CD5lo_2_TRBV")
)
TRAJ_all_counts <- bind_rows(
  CD5hi_1_TRAJ_counts %>% mutate(Sample = "CD5hi_1_TRAJ"),
  CD5hi_2_TRAJ_counts %>% mutate(Sample = "CD5hi_2_TRAJ"),
  CD5lo_1_TRAJ_counts %>% mutate(Sample = "CD5lo_1_TRAJ"),
  CD5lo_2_TRAJ_counts %>% mutate(Sample = "CD5lo_2_TRAJ")
)
TRBJ_all_counts <- bind_rows(
  CD5hi_1_TRBJ_counts %>% mutate(Sample = "CD5hi_1_TRBJ"),
  CD5hi_2_TRBJ_counts %>% mutate(Sample = "CD5hi_2_TRBJ"),
  CD5lo_1_TRBJ_counts %>% mutate(Sample = "CD5lo_1_TRBJ"),
  CD5lo_2_TRBJ_counts %>% mutate(Sample = "CD5lo_2_TRBJ")
)

# Pivot the data frame to wide format
TRAV_all_counts <-
  pivot_wider(
    data = TRAV_all_counts,
    names_from = Sample,
    values_from = total_count
  )
TRBV_all_counts <-
  pivot_wider(
    data = TRBV_all_counts,
    names_from = Sample,
    values_from = total_count
  )
TRAJ_all_counts <-
  pivot_wider(
    data = TRAJ_all_counts,
    names_from = Sample,
    values_from = total_count
  )
TRBJ_all_counts <-
  pivot_wider(
    data = TRBJ_all_counts,
    names_from = Sample,
    values_from = total_count
  )

# Filter TRAV_all_counts to keep only common genes
TRAV_data_filtered <- TRAV_all_counts %>%
  filter(V %in% TRAV_common_genes)
write_xlsx(TRAV_data_filtered, "TRAV_all.xlsx", format_headers = TRUE)
# Plot Fig.7B in Prism

# Filter TRAJ_all_counts to keep only common genes
TRAJ_data_filtered <- TRAJ_all_counts %>%
  filter(J %in% TRAJ_common_genes)
write_xlsx(TRAJ_data_filtered, "TRAJ_all.xlsx", format_headers = TRUE)
# Plot Fig.7C in Prism

# Filter TRBV_all_counts to keep only common genes
TRBV_data_filtered <- TRBV_all_counts %>%
  filter(V %in% TRBV_common_genes)
write_xlsx(TRBV_data_filtered, "TRBV_all.xlsx", format_headers = TRUE)
# Plot Fig.7D in Prism

# Filter TRBJ_all_counts to keep only common genes
TRBJ_data_filtered <- TRBJ_all_counts %>%
  filter(J %in% TRBJ_common_genes)
write_xlsx(TRBJ_data_filtered, "TRBJ_all.xlsx", format_headers = TRUE)
# Plot Fig.7E in Prism

### Fig.7F ### 
# library("rTCRBCRr")
# library("readr")
# Define the file paths
input_files <- c(
  "CD5hi_1_trimmed.clonotypes.TRA.txt",
  "CD5hi_1_trimmed.clonotypes.TRB.txt",
  "CD5hi_2_trimmed.clonotypes.TRA.txt",
  "CD5hi_2_trimmed.clonotypes.TRB.txt",
  "CD5lo_1_trimmed.clonotypes.TRA.txt",
  "CD5lo_1_trimmed.clonotypes.TRB.txt",
  "CD5lo_2_trimmed.clonotypes.TRA.txt",
  "CD5lo_2_trimmed.clonotypes.TRB.txt"
)
# Read the MiXCR files into a list of data frames
mixcr_data <- lapply(input_files, read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
present_tool <- c("trust", "mixcr")[2]
CD5hi_1_TRA <- mixcr_data[[1]] %>%
  format_clonotype_to_immunarch_style(., clonotyping_tool = present_tool) %>%
  remove_nonproductive_CDR3aa %>%
  annotate_chain_name_and_isotype_name %>%
  merge_convergent_clonotype
CD5hi_1_TRB <- mixcr_data[[2]] %>%
  format_clonotype_to_immunarch_style(., clonotyping_tool = present_tool) %>%
  remove_nonproductive_CDR3aa %>%
  annotate_chain_name_and_isotype_name %>%
  merge_convergent_clonotype
CD5hi_2_TRA <- mixcr_data[[3]] %>%
  format_clonotype_to_immunarch_style(., clonotyping_tool = present_tool) %>%
  remove_nonproductive_CDR3aa %>%
  annotate_chain_name_and_isotype_name %>%
  merge_convergent_clonotype
CD5hi_2_TRB <- mixcr_data[[4]] %>%
  format_clonotype_to_immunarch_style(., clonotyping_tool = present_tool) %>%
  remove_nonproductive_CDR3aa %>%
  annotate_chain_name_and_isotype_name %>%
  merge_convergent_clonotype
CD5lo_1_TRA <- mixcr_data[[5]] %>%
  format_clonotype_to_immunarch_style(., clonotyping_tool = present_tool) %>%
  remove_nonproductive_CDR3aa %>%
  annotate_chain_name_and_isotype_name %>%
  merge_convergent_clonotype
CD5lo_1_TRB <- mixcr_data[[6]] %>%
  format_clonotype_to_immunarch_style(., clonotyping_tool = present_tool) %>%
  remove_nonproductive_CDR3aa %>%
  annotate_chain_name_and_isotype_name %>%
  merge_convergent_clonotype
CD5lo_2_TRA <- mixcr_data[[7]] %>%
  format_clonotype_to_immunarch_style(., clonotyping_tool = present_tool) %>%
  remove_nonproductive_CDR3aa %>%
  annotate_chain_name_and_isotype_name %>%
  merge_convergent_clonotype
CD5lo_2_TRB <- mixcr_data[[8]] %>%
  format_clonotype_to_immunarch_style(., clonotyping_tool = present_tool) %>%
  remove_nonproductive_CDR3aa %>%
  annotate_chain_name_and_isotype_name %>%
  merge_convergent_clonotype

TRA_all_CDR3 <- bind_rows(
  CD5hi_1_TRA %>% mutate(Sample = "CD5hi_TRA"),
  CD5hi_2_TRA %>% mutate(Sample = "CD5hi_TRA"),
  CD5lo_1_TRA %>% mutate(Sample = "CD5lo_TRA"),
  CD5lo_2_TRA %>% mutate(Sample = "CD5lo_TRA")
)
TRB_all_CDR3 <- bind_rows(
  CD5hi_1_TRB %>% mutate(Sample = "CD5hi_TRB"),
  CD5hi_2_TRB %>% mutate(Sample = "CD5hi_TRB"),
  CD5lo_1_TRB %>% mutate(Sample = "CD5lo_TRB"),
  CD5lo_2_TRB %>% mutate(Sample = "CD5lo_TRB")
)

## Plot the TRA_all_CDR3 ## 
ggplot(data = TRA_all_CDR3, aes(x = as.factor(nchar(CDR3.aa)), y = Proportion, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "CDR3 length (AA)", y = "Relative frequency (%)") + 
  scale_y_continuous(limits = c(0, max(TRA_all_CDR3$Proportion) + 0.0015), expand = c(0, 0)) +
  scale_x_discrete(limits = 3:30) +
  scale_fill_manual(values = c("CD5hi_TRA" = "orange", "CD5lo_TRA" = "navy")) +
  theme(legend.position = c(0.87, 0.83))
# > TRA_t_test_result
# Welch Two Sample t-test
# data:  TRA_cdr3_length_cd5hi and TRA_cdr3_length_cd5lo
# t = -0.91952, df = 3153, p-value = 0.3579
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.15286114  0.05525888
# sample estimates:
# mean of x mean of y 
# 13.58864  13.63744

# Plot the TRB_all_CDR3 ##
# geom_point(stat = "identity", position = position_dodge(width = 0.9), color = "white")
ggplot(data = TRB_all_CDR3, aes(x = as.factor(nchar(CDR3.aa)), y = Proportion, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "CDR3 length (AA)", y = "Relative frequency (%)") + 
  scale_y_continuous(limits = c(0, max(TRB_all_CDR3$Proportion) + 0.0015), expand = c(0, 0)) +
  scale_x_discrete(limits = 4:30) +
  scale_fill_manual(values = c("CD5hi_TRB" = "orange", "CD5lo_TRB" = "navy")) +
  theme(legend.position = c(0.87, 0.83))
# > TRB_t_test_result
# Welch Two Sample t-test
# data:  TRB_cdr3_length_cd5hi and TRB_cdr3_length_cd5lo
# t = -2.2129, df = 5835.2, p-value = 0.02694
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.158334018 -0.009581303
# sample estimates:
# mean of x mean of y 
# 13.74450  13.82845 

### Fig.7G ###
## TRA ##
## Extract unique CDR3 sequences from each dataset along with their clone counts ## 
TRA_unique_CDR3_hi_1 <- unique(CD5hi_1_TRA$CDR3.aa)
TRA_clone_count_CDR3_hi_1 <- CD5hi_1_TRA$Clones[match(TRA_unique_CDR3_hi_1, CD5hi_1_TRA$CDR3.aa)]
TRA_unique_CDR3_hi_2 <- unique(CD5hi_2_TRA$CDR3.aa)
TRA_clone_count_CDR3_hi_2 <- CD5hi_2_TRA$Clones[match(TRA_unique_CDR3_hi_2, CD5hi_2_TRA$CDR3.aa)]
TRA_unique_CDR3_lo_1 <- unique(CD5lo_1_TRA$CDR3.aa)
TRA_clone_count_CDR3_lo_1 <- CD5lo_1_TRA$Clones[match(TRA_unique_CDR3_lo_1, CD5lo_1_TRA$CDR3.aa)]
TRA_unique_CDR3_lo_2 <- unique(CD5lo_2_TRA$CDR3.aa)
TRA_clone_count_CDR3_lo_2 <- CD5lo_2_TRA$Clones[match(TRA_unique_CDR3_lo_2, CD5lo_2_TRA$CDR3.aa)]
# Function to calculate unique overlap count
calculate_overlap_count <- function(unique1, count1, unique2, count2) {
  intersected <- intersect(unique1, unique2)
  count_overlap <- sum(pmin(count1[match(intersected, unique1)], count2[match(intersected, unique2)]))
  return(count_overlap)
}
## Calculate unique overlap counts ##
TRA_overlap_count_hi_hi <- calculate_overlap_count(TRA_unique_CDR3_hi_1, TRA_clone_count_CDR3_hi_1, TRA_unique_CDR3_hi_2, TRA_clone_count_CDR3_hi_2)
TRA_overlap_count_hi_lo1 <- calculate_overlap_count(TRA_unique_CDR3_hi_1, TRA_clone_count_CDR3_hi_1, TRA_unique_CDR3_lo_1, TRA_clone_count_CDR3_lo_1)
TRA_overlap_count_hi_lo2 <- calculate_overlap_count(TRA_unique_CDR3_hi_1, TRA_clone_count_CDR3_hi_1, TRA_unique_CDR3_lo_2, TRA_clone_count_CDR3_lo_2)
TRA_overlap_count_hi2_lo1 <- calculate_overlap_count(TRA_unique_CDR3_hi_2, TRA_clone_count_CDR3_hi_2, TRA_unique_CDR3_lo_1, TRA_clone_count_CDR3_lo_1)
TRA_overlap_count_hi2_lo2 <- calculate_overlap_count(TRA_unique_CDR3_hi_2, TRA_clone_count_CDR3_hi_2, TRA_unique_CDR3_lo_2, TRA_clone_count_CDR3_lo_2)
TRA_overlap_count_lo_lo <- calculate_overlap_count(TRA_unique_CDR3_lo_1, TRA_clone_count_CDR3_lo_1, TRA_unique_CDR3_lo_2, TRA_clone_count_CDR3_lo_2)
# Create a 4x4 matrix to store the overlap counts
TRA_overlap_count_matrix <- matrix(0, nrow = 4, ncol = 4)
rownames(TRA_overlap_count_matrix) <- colnames(TRA_overlap_count_matrix) <- c("CD5hi_1_TRA", "CD5hi_2_TRA", "CD5lo_1_TRA", "CD5lo_2_TRA")
# Fill in the overlap counts in the matrix
TRA_overlap_count_matrix[1, 2] <- TRA_overlap_count_matrix[2, 1] <- TRA_overlap_count_hi_hi
TRA_overlap_count_matrix[1, 3] <- TRA_overlap_count_matrix[3, 1] <- TRA_overlap_count_hi_lo1
TRA_overlap_count_matrix[1, 4] <- TRA_overlap_count_matrix[4, 1] <- TRA_overlap_count_hi_lo2
TRA_overlap_count_matrix[2, 3] <- TRA_overlap_count_matrix[3, 2] <- TRA_overlap_count_hi2_lo1
TRA_overlap_count_matrix[2, 4] <- TRA_overlap_count_matrix[4, 2] <- TRA_overlap_count_hi2_lo2
TRA_overlap_count_matrix[3, 4] <- TRA_overlap_count_matrix[4, 3] <- TRA_overlap_count_lo_lo
# Print the overlap count matrix
print(TRA_overlap_count_matrix)
# This code calculates the unique overlap counts between datasets based on the Clones column. 
# The resulting matrix overlap_count_matrix contains the counts of shared clones between each pair of datasets.
## Plot in Prism
# > print(TRA_overlap_count_matrix)
#              CD5hi_1_TRA CD5hi_2_TRA CD5lo_1_TRA CD5lo_2_TRA
# CD5hi_1_TRA           0          50          37          37
# CD5hi_2_TRA          50           0          39          34
# CD5lo_1_TRA          37          39           0          27
# CD5lo_2_TRA          37          34          27           0

## TRB ##
## Extract unique CDR3 sequences from each dataset along with their clone counts ##
TRB_unique_CDR3_hi_1 <- unique(CD5hi_1_TRB$CDR3.aa)
TRB_clone_count_CDR3_hi_1 <- CD5hi_1_TRB$Clones[match(TRB_unique_CDR3_hi_1, CD5hi_1_TRB$CDR3.aa)]
TRB_unique_CDR3_hi_2 <- unique(CD5hi_2_TRB$CDR3.aa)
TRB_clone_count_CDR3_hi_2 <- CD5hi_2_TRB$Clones[match(TRB_unique_CDR3_hi_2, CD5hi_2_TRB$CDR3.aa)]
TRB_unique_CDR3_lo_1 <- unique(CD5lo_1_TRB$CDR3.aa)
TRB_clone_count_CDR3_lo_1 <- CD5lo_1_TRB$Clones[match(TRB_unique_CDR3_lo_1, CD5lo_1_TRB$CDR3.aa)]
TRB_unique_CDR3_lo_2 <- unique(CD5lo_2_TRB$CDR3.aa)
TRB_clone_count_CDR3_lo_2 <- CD5lo_2_TRB$Clones[match(TRB_unique_CDR3_lo_2, CD5lo_2_TRB$CDR3.aa)]
# Calculate unique overlap counts
TRB_overlap_count_hi_hi <- calculate_overlap_count(TRB_unique_CDR3_hi_1, TRB_clone_count_CDR3_hi_1, TRB_unique_CDR3_hi_2, TRB_clone_count_CDR3_hi_2)
TRB_overlap_count_hi_lo1 <- calculate_overlap_count(TRB_unique_CDR3_hi_1, TRB_clone_count_CDR3_hi_1, TRB_unique_CDR3_lo_1, TRB_clone_count_CDR3_lo_1)
TRB_overlap_count_hi_lo2 <- calculate_overlap_count(TRB_unique_CDR3_hi_1, TRB_clone_count_CDR3_hi_1, TRB_unique_CDR3_lo_2, TRB_clone_count_CDR3_lo_2)
TRB_overlap_count_hi2_lo1 <- calculate_overlap_count(TRB_unique_CDR3_hi_2, TRB_clone_count_CDR3_hi_2, TRB_unique_CDR3_lo_1, TRB_clone_count_CDR3_lo_1)
TRB_overlap_count_hi2_lo2 <- calculate_overlap_count(TRB_unique_CDR3_hi_2, TRB_clone_count_CDR3_hi_2, TRB_unique_CDR3_lo_2, TRB_clone_count_CDR3_lo_2)
TRB_overlap_count_lo_lo <- calculate_overlap_count(TRB_unique_CDR3_lo_1, TRB_clone_count_CDR3_lo_1, TRB_unique_CDR3_lo_2, TRB_clone_count_CDR3_lo_2)
# Create a 4x4 matrix to store the overlap counts
TRB_overlap_count_matrix <- matrix(0, nrow = 4, ncol = 4)
rownames(TRB_overlap_count_matrix) <- colnames(TRB_overlap_count_matrix) <- c("CD5hi_1_TRB", "CD5hi_2_TRB", "CD5lo_1_TRB", "CD5lo_2_TRB")
# Fill in the overlap counts in the matrix
TRB_overlap_count_matrix[1, 2] <- TRB_overlap_count_matrix[2, 1] <- TRB_overlap_count_hi_hi
TRB_overlap_count_matrix[1, 3] <- TRB_overlap_count_matrix[3, 1] <- TRB_overlap_count_hi_lo1
TRB_overlap_count_matrix[1, 4] <- TRB_overlap_count_matrix[4, 1] <- TRB_overlap_count_hi_lo2
TRB_overlap_count_matrix[2, 3] <- TRB_overlap_count_matrix[3, 2] <- TRB_overlap_count_hi2_lo1
TRB_overlap_count_matrix[2, 4] <- TRB_overlap_count_matrix[4, 2] <- TRB_overlap_count_hi2_lo2
TRB_overlap_count_matrix[3, 4] <- TRB_overlap_count_matrix[4, 3] <- TRB_overlap_count_lo_lo
# Print the overlap count matrix
print(TRB_overlap_count_matrix)
# This code calculates the unique overlap counts between datasets based on the Clones column. 
# The resulting matrix overlap_count_matrix contains the counts of shared clones between each pair of datasets.
## Plot in Prism
# > print(TRB_overlap_count_matrix)
#             CD5hi_1_TRB CD5hi_2_TRB CD5lo_1_TRB CD5lo_2_TRB
# CD5hi_1_TRB           0          55          29          46
# CD5hi_2_TRB          55           0          26          34
# CD5lo_1_TRB          29          26           0          27
# CD5lo_2_TRB          46          34          27           0

### Fig.7I ###
# library(dplyr)
# Define hydropathy scale dictionary
hydropathy_scale <- list(I = 4.5, V = 4.2, L = 3.8, F = 2.8, C = 2.5, 
                         M = 1.9, A = 1.8, G = -0.4, T = -0.7, S = -0.8,
                         W = -0.9, Y = -1.3, P = -1.6, H = -3.2, E = -3.5,
                         Q = -3.5, D = -3.5, N = -3.5, K = -3.9, R = -4.5)

# Define the calculateHydrophobicity function
calculateHydrophobicity <- function(cdr3_seq, hydropathy_scale) {
  cdr3_ls <- unlist(strsplit(cdr3_seq, ''))
  sum_hdp <- sum(sapply(cdr3_ls, function(aa) hydropathy_scale[[aa]]))
  result <- sum_hdp / length(cdr3_ls)
  return(result)
}

# Apply the calculateHydrophobicity function to the CDR3.aa column
CD5hi_1_TRA$hydrophobicity <- sapply(CD5hi_1_TRA$CDR3.aa, calculateHydrophobicity, hydropathy_scale)
CD5hi_2_TRA$hydrophobicity <- sapply(CD5hi_2_TRA$CDR3.aa, calculateHydrophobicity, hydropathy_scale)
CD5lo_1_TRA$hydrophobicity <- sapply(CD5lo_1_TRA$CDR3.aa, calculateHydrophobicity, hydropathy_scale)
CD5lo_2_TRA$hydrophobicity <- sapply(CD5lo_2_TRA$CDR3.aa, calculateHydrophobicity, hydropathy_scale)
CD5hi_1_TRB$hydrophobicity <- sapply(CD5hi_1_TRB$CDR3.aa, calculateHydrophobicity, hydropathy_scale)
CD5hi_2_TRB$hydrophobicity <- sapply(CD5hi_2_TRB$CDR3.aa, calculateHydrophobicity, hydropathy_scale)
CD5lo_1_TRB$hydrophobicity <- sapply(CD5lo_1_TRB$CDR3.aa, calculateHydrophobicity, hydropathy_scale)
CD5lo_2_TRB$hydrophobicity <- sapply(CD5lo_2_TRB$CDR3.aa, calculateHydrophobicity, hydropathy_scale)

# List of sample names
samples <- c("CD5hi_1_TRA", "CD5hi_2_TRA", "CD5lo_1_TRA", "CD5lo_2_TRA", "CD5hi_1_TRB", "CD5hi_2_TRB", "CD5lo_1_TRB",  "CD5lo_2_TRB")
# Function to calculate average hydrophobicity score
calculate_avg_hydrophobicity <- function(sample_name) {
  # CD5hi_1_TRA, CD5hi_1_TRB, etc. are data frames with columns hydrophobicity and Clones
  sum_hydrophobicity <- sum(get(sample_name)$hydrophobicity)
  sum_clones <- sum(get(sample_name)$Clones)
  avg_hydrophobicity <- sum_hydrophobicity / sum_clones
  return(avg_hydrophobicity)
}
# Apply the function to each sample
avg_hydrophobicity_scores <- lapply(samples, calculate_avg_hydrophobicity)
# Print the results
names(avg_hydrophobicity_scores) <- samples
print(avg_hydrophobicity_scores)
# Plot in Prism
# $CD5hi_1_TRA
# [1] 0.1565368
# $CD5hi_2_TRA
# [1] 0.1646369
# $CD5lo_1_TRA
# [1] 0.1410456
# $CD5lo_2_TRA
# [1] 0.1362827
# $CD5hi_1_TRB
# [1] -0.1712868
# $CD5hi_2_TRB
# [1] -0.1582817
# $CD5lo_1_TRB
# [1] -0.1583921
# $CD5lo_2_TRB
# [1] -0.1748877

### Fig.7J ###
CD5hi_1_TRA_repertoire_metrics <- CD5hi_1_TRA %>% 
  compute_repertoire_metrics_by_chain_name
CD5hi_1_TRB_repertoire_metrics <- CD5hi_1_TRB %>% 
  compute_repertoire_metrics_by_chain_name
CD5hi_2_TRA_repertoire_metrics <- CD5hi_2_TRA %>% 
  compute_repertoire_metrics_by_chain_name
CD5hi_2_TRB_repertoire_metrics <- CD5hi_2_TRB %>% 
  compute_repertoire_metrics_by_chain_name
CD5lo_1_TRA_repertoire_metrics <- CD5lo_1_TRA %>% 
  compute_repertoire_metrics_by_chain_name
CD5lo_1_TRB_repertoire_metrics <- CD5lo_1_TRB %>% 
  compute_repertoire_metrics_by_chain_name
CD5lo_2_TRA_repertoire_metrics <- CD5lo_2_TRA %>% 
  compute_repertoire_metrics_by_chain_name
CD5lo_2_TRB_repertoire_metrics <- CD5lo_2_TRB %>% 
  compute_repertoire_metrics_by_chain_name
## Plot in Prism ## 
# > CD5hi_1_TRA_repertoire_metrics
#      diversity  clonality richness  evenness      median
# IGH        NA         NA       NA        NA          NA
# IGK        NA         NA       NA        NA          NA
# IGL        NA         NA       NA        NA          NA
# TRA  6.751191 0.01489907      947 0.9851009 0.000748503
# TRB        NA         NA       NA        NA          NA
# TRD        NA         NA       NA        NA          NA
# TRG        NA         NA       NA        NA          NA

# > CD5hi_2_TRA_repertoire_metrics
#      diversity  clonality richness  evenness       median
# IGH        NA         NA       NA        NA           NA
# IGK        NA         NA       NA        NA           NA
# IGL        NA         NA       NA        NA           NA
# TRA 6.7423344 0.01420314      934 0.9857969 0.0007830854
# TRB        NA         NA       NA        NA           NA
# TRD 0.6931472 0.00000000        2 1.0000000 0.5000000000
# TRG        NA         NA       NA        NA           NA

# > CD5lo_1_TRA_repertoire_metrics
#      diversity  clonality richness  evenness       median
# IGH        NA         NA       NA        NA           NA
# IGK        NA         NA       NA        NA           NA
# IGL        NA         NA       NA        NA           NA
# TRA   6.54199 0.01454793      764 0.9854521 0.0009496676
# TRB        NA         NA       NA        NA           NA
# TRD   1.56071 0.03027610        5 0.9697239 0.1666666667
# TRG        NA         NA       NA        NA           NA

# > CD5lo_2_TRA_repertoire_metrics
#      diversity  clonality richness  evenness       median
# IGH        NA         NA       NA        NA           NA
# IGK        NA         NA       NA        NA           NA
# IGL        NA         NA       NA        NA           NA
# TRA 6.5125700 0.01544367      746 0.9845563 0.0009460738
# TRB        NA         NA       NA        NA           NA
# TRD 0.6931472 0.00000000        2 1.0000000 0.5000000000
# TRG        NA         NA       NA        NA           NA

# > CD5hi_2_TRB_repertoire_metrics
#      diversity clonality richness  evenness       median
# IGH        NA        NA       NA        NA           NA
# IGK        NA        NA       NA        NA           NA
# IGL        NA        NA       NA        NA           NA
# TRA        NA        NA       NA        NA           NA
# TRB  7.234592 0.0130706     1526 0.9869294 0.0004878049
# TRD        NA        NA       NA        NA           NA
# TRG        NA        NA       NA        NA           NA

# > CD5hi_2_TRB_repertoire_metrics
#      diversity clonality richness  evenness       median
# IGH        NA        NA       NA        NA           NA
# IGK        NA        NA       NA        NA           NA
# IGL        NA        NA       NA        NA           NA
# TRA        NA        NA       NA        NA           NA
# TRB  7.234592 0.0130706     1526 0.9869294 0.0004878049
# TRD        NA        NA       NA        NA           NA
# TRG        NA        NA       NA        NA           NA

# > CD5lo_1_TRB_repertoire_metrics
#      diversity  clonality richness  evenness       median
# IGH        NA         NA       NA        NA           NA
# IGK        NA         NA       NA        NA           NA
# IGL        NA         NA       NA        NA           NA
# TRA        NA         NA       NA        NA           NA
# TRB  7.060919 0.01438148     1292 0.9856185 0.0005602241
# TRD        NA         NA       NA        NA           NA
# TRG        NA         NA       NA        NA           NA

# > CD5lo_2_TRB_repertoire_metrics
#      diversity  clonality richness  evenness       median
# IGH        NA         NA       NA        NA           NA
# IGK        NA         NA       NA        NA           NA
# IGL        NA         NA       NA        NA           NA
# TRA        NA         NA       NA        NA           NA
# TRB  7.155223 0.01296586     1407 0.9870341 0.0005208333
# TRD        NA         NA       NA        NA           NA
# TRG        NA         NA       NA        NA           NA

### Supplemental material ###
# Extract V, D, and J columns
# Count occurrences of each combination
VJdata_CD5hi_1_TRAV <- CD5hi_1_TRAV %>%
  select(V, J) %>%
  group_by(V, J) %>%
  summarise(Count = n())
VJdata_CD5hi_2_TRAV <- CD5hi_2_TRAV %>%
  select(V, J) %>%
  group_by(V, J) %>%
  summarise(Count = n())
TRAVJdata_CD5hi <- bind_rows(VJdata_CD5hi_1_TRAV, VJdata_CD5hi_2_TRAV)
TRAVJdata_CD5hi <- 
  TRAVJdata_CD5hi %>%
  group_by(V, J) %>%
  summarise(Count = n())

VJdata_CD5lo_1_TRAV <- CD5lo_1_TRAV %>%
  select(V, J) %>%
  group_by(V, J) %>%
  summarise(Count = n())
VJdata_CD5lo_2_TRAV <- CD5lo_2_TRAV %>%
  select(V, J) %>%
  group_by(V, J) %>%
  summarise(Count = n())
TRAVJdata_CD5lo <- bind_rows(VJdata_CD5lo_1_TRAV, VJdata_CD5lo_2_TRAV)
TRAVJdata_CD5lo <- 
  TRAVJdata_CD5lo %>%
  group_by(V, J) %>%
  summarise(Count = n())

VJdata_CD5hi_1_TRBV <- CD5hi_1_TRBV %>%
  select(V, J) %>%
  group_by(V, J) %>%
  summarise(Count = n())
VJdata_CD5hi_2_TRBV <- CD5hi_2_TRBV %>%
  select(V, J) %>%
  group_by(V, J) %>%
  summarise(Count = n())
TRBVJdata_CD5hi <- bind_rows(VJdata_CD5hi_1_TRBV, VJdata_CD5hi_2_TRBV)
TRBVJdata_CD5hi <- 
  TRBVJdata_CD5hi %>%
  group_by(V, J) %>%
  summarise(Count = n())

VJdata_CD5lo_1_TRBV <- CD5lo_1_TRBV %>%
  select(V, J) %>%
  group_by(V, J) %>%
  summarise(Count = n())
VJdata_CD5lo_2_TRBV <- CD5lo_2_TRBV %>%
  select(V, J) %>%
  group_by(V, J) %>%
  summarise(Count = n())
TRBVJdata_CD5lo <- bind_rows(VJdata_CD5lo_1_TRBV, VJdata_CD5lo_2_TRBV)
TRBVJdata_CD5lo <- 
  TRBVJdata_CD5lo %>%
  group_by(V, J) %>%
  summarise(Count = n())

# Add a label column to each dataset indicating CD5hi or CD5lo
TRAVJdata_CD5hi$Label <- "CD5hi"
TRAVJdata_CD5lo$Label <- "CD5lo"
TRBVJdata_CD5hi$Label <- "CD5hi"
TRBVJdata_CD5lo$Label <- "CD5lo"

# Combine the datasets
TRAV_combined_data <- rbind(TRAVJdata_CD5hi, TRAVJdata_CD5lo)
TRBV_combined_data <- rbind(TRBVJdata_CD5hi, TRBVJdata_CD5lo)

# facet_wrap(~ Label, ncol = 1) +  # Side-by-side position
ggplot(TRAV_combined_data, aes(x = V, y = J, size = Count, color = Label)) +
  geom_point(alpha = ifelse(TRAV_combined_data$Label == "CD5lo", 0.5, 0.5)) +  # Set alpha to 1 for CD5lo, 0.5 for CD5hi
  scale_size_continuous(range = c(1, 4)) +
  labs(x = "TRAV genes", y = "TRAJ genes", title = "TRA repertoire") +
  scale_color_manual(values = c("CD5hi" = "red", "CD5lo" = "blue")) + 
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggplot(TRBV_combined_data, aes(x = V, y = J, size = Count, color = Label)) +
  geom_point(alpha = ifelse(TRBV_combined_data$Label == "CD5lo", 0.5, 0.5)) +  # Set alpha to 1 for CD5lo, 0.5 for CD5hi
  scale_size_continuous(range = c(1, 4)) +
  labs(x = "TRBV genes", y = "TRBJ genes", title = "TRB repertoire") +
  scale_color_manual(values = c("CD5hi" = "red", "CD5lo" = "blue")) + 
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
