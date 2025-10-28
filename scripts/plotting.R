library(DescTools)
library(tidyverse)
library(data.table)
library(taxonomizr)
library(dtplyr)
library(rlang)
library(ggforce)
library(openxlsx)
library(rioja)
library(patchwork)
library(wesanderson)
library(optparse)
library(yaml)
library(readxl)
library(networkD3)
library(RColorBrewer)
library(ggplot2)

source("/Users/bfj994/Documents/Thule/functions/mergebamfilterandmetaDMGA.R")
read_with_id <- function(file) {
  
  # Extract only the CGG number and remove everything after it
  filename <- basename(file)
  print(paste("Filename:", filename))  # Check the exact filename being passed
  
  # Use a more general regex to match the CGG ID part before any dot or extra text
  sample_id <- sub("\\..*$", "", filename)
  print(paste("Sample ID:", sample_id))  # Check the extracted sample ID
  
  # Read the file and skip the first two lines if they match the undesired pattern
  df <- tryCatch(
    {
      temp_df <- read_tsv(file, col_names = FALSE, col_types = cols(), quote = "", comment = "", progress = FALSE, skip = 2)
      
      # Check if the file has at least 5 columns
      if (ncol(temp_df) < 6) {
        message("Skipping file: ", file, " (not enough columns)")
        return(NULL)
      }
      
      temp_df |> 
        select(1:6) |>  # Keep only the first 5 columns
        rename(Read_ID = X1, Sequence = X2, Length = X3, Count = X4, Score = X5, LCA_Assignment = X6) |> 
        mutate(CGG_ID = sample_id)  # Append cleaned CGG_ID
    },
    error = function(e) {
      message("Skipping file: ", file, " (error reading)")
      return(NULL)
    }
  )
  
  return(df)
}
setwd("/Users/bfj994/Documents/Iceland/Rosaceae/")

meta <- read.csv("../metaData.txt")

origP <- list.files(path = "Pyrus_original", 
                       pattern = "lca$", 
                       full.names = TRUE)
origLCAP <- map_dfr(origP, read_with_id)

newP <- list.files(path = "Pyrus_new",
                  pattern = "lca.gz$",
                  full.names = TRUE)
newLCAP <- map_dfr(newP, read_with_id)
origM <- list.files(path = "Malus_original", 
                    pattern = "lca$", 
                    full.names = TRUE)
origLCAM <- map_dfr(origM, read_with_id)

newM <- list.files(path = "Malus_new",
                   pattern = "lca.gz$",
                   full.names = TRUE)
newLCAM <- map_dfr(newM, read_with_id)

origLCA <- rbind(origLCAM, origLCAP)
newLCA <- rbind(newLCAM, newLCAP)


# Step 1: Merge the dataframes
merged_df <- origLCA %>%
  inner_join(newLCA, by = c("Read_ID", "CGG_ID", "Sequence"), suffix = c("_orig", "_new")) %>%
  inner_join(meta, by = "CGG_ID")

# Count the frequency of each combination of old and new LCA assignments
heatmap_data_count <- merged_df %>%
  mutate(
    LCA_Assignment_orig = str_extract(LCA_Assignment_orig, '(?<=:").*?(?=":)'),
    LCA_Assignment_new  = str_extract(LCA_Assignment_new,  '(?<=:").*?(?=":)')
  ) %>%
  count(LCA_Assignment_orig, LCA_Assignment_new, CGG_ID, Date) %>%
  filter(n >= 10)

ggplot(heatmap_data_count, aes(x = LCA_Assignment_orig, y = LCA_Assignment_new, fill = n)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10", breaks = scales::breaks_log(10)) +
  labs(
    x = "LCA assignment for the original DB",
    y = "LCA assignment for a broader Rosaceae DB",
    fill = "Count"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 16),     # Increased x-axis label size
    axis.text.y = element_text(size = 16),                            # Increased y-axis label size
    axis.title.x = element_text(size = 16),                           # Increased x-axis title size
    axis.title.y = element_text(size = 16),                           # Increased y-axis title size
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.key.height = unit(1.5, "cm"),
    strip.text = element_text(size = 16)                              # Increased facet label size
  ) +
  facet_wrap(~Date, labeller = labeller(Date = function(x) paste("Year CE", x)), ncol = 5)

ggsave("heatmap_Rosaceae.png", width = 13, height = 17, dpi = 350)

setwd("/Users/bfj994/Documents/Iceland/Poaceae/")
origAe <- list.files(path = "Aegilops_original", 
                    pattern = "lca$", 
                    full.names = TRUE)
origLCAAe <- map_dfr(origAe, read_with_id)
newAe <- list.files(path = "Aegilops_new",
                   pattern = "lca.gz$",
                   full.names = TRUE)
newLCAAe <- map_dfr(newAe, read_with_id)

origAv <- list.files(path = "Avena_original", 
                    pattern = "lca$", 
                    full.names = TRUE)
origLCAAv <- map_dfr(origAv, read_with_id)
newAv <- list.files(path = "Avena_new",
                   pattern = "lca.gz$",
                   full.names = TRUE)
newLCAAv <- map_dfr(newAv, read_with_id)

origB <- list.files(path = "Brachypodium_original", 
                     pattern = "lca$", 
                     full.names = TRUE)
origLCAB <- map_dfr(origB, read_with_id)
newB <- list.files(path = "Brachypodium_new",
                    pattern = "lca.gz$",
                    full.names = TRUE)
newLCAB <- map_dfr(newB, read_with_id)

origHe <- list.files(path = "Helictochloa_original", 
                     pattern = "lca$", 
                     full.names = TRUE)
origLCAHe <- map_dfr(origHe, read_with_id)
newHe <- list.files(path = "Helictochloa_new",
                    pattern = "lca.gz$",
                    full.names = TRUE)
newLCAHe <- map_dfr(newHe, read_with_id)

origHo <- list.files(path = "Hordelymus_original", 
                     pattern = "lca$", 
                     full.names = TRUE)
origLCAHo <- map_dfr(origHo, read_with_id)
newHo <- list.files(path = "Hordelymus_new",
                    pattern = "lca.gz$",
                    full.names = TRUE)
newLCAHo <- map_dfr(newHo, read_with_id)

origK <- list.files(path = "Koeleria_original", 
                     pattern = "lca$", 
                     full.names = TRUE)
origLCAK <- map_dfr(origK, read_with_id)
newK <- list.files(path = "Koeleria_new",
                    pattern = "lca.gz$",
                    full.names = TRUE)
newLCAK <- map_dfr(newK, read_with_id)

origPh <- list.files(path = "Phalaris_original", 
                     pattern = "lca$", 
                     full.names = TRUE)
origLCAPh <- map_dfr(origPh, read_with_id)
newPh <- list.files(path = "Phalaris_new",
                    pattern = "lca.gz$",
                    full.names = TRUE)
newLCAPh <- map_dfr(newPh, read_with_id)

origT <- list.files(path = "Triticum_original", 
                     pattern = "lca$", 
                     full.names = TRUE)
origLCAT <- map_dfr(origT, read_with_id)
newT <- list.files(path = "Triticum_new",
                    pattern = "lca.gz$",
                    full.names = TRUE)
newLCAT <- map_dfr(newT, read_with_id)


origLCA <- rbind(origLCAAe, origLCAAv, origLCAB, origLCAHe, origLCAHo, origLCAK, origLCAPh, origLCAT)
newLCA <- rbind(newLCAAe, newLCAAv, newLCAB, newLCAHe, newLCAK, newLCAHo, newLCAPh, newLCAT)
merged_df <- origLCA %>%
  inner_join(newLCA, by = c("Read_ID", "CGG_ID", "Sequence"), suffix = c("_orig", "_new")) %>%
  inner_join(meta, by = "CGG_ID")


heatmap_data_count <- merged_df %>%
  mutate(
    LCA_Assignment_orig = str_extract(LCA_Assignment_orig, '(?<=:").*?(?=":)'),
    LCA_Assignment_new  = str_extract(LCA_Assignment_new,  '(?<=:").*?(?=":)')
  ) %>%
  count(LCA_Assignment_orig, LCA_Assignment_new, CGG_ID, Date) %>%
  filter(n >= 100)

# Create the heatmap using ggplot2
ggplot(heatmap_data_count, aes(x = LCA_Assignment_orig, y = LCA_Assignment_new, fill = n)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10", breaks = scales::breaks_log(10)) +  # Log scale for color mapping, real values in legend
  labs(x = "LCA assignment for the original DB",
       y = "LCA assignment for a broader Poaceae DB",
       fill = "Count") +  # Legend still shows real counts
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 16),     # Increased x-axis label size
    axis.text.y = element_text(size = 16),                            # Increased y-axis label size
    axis.title.x = element_text(size = 16),                           # Increased x-axis title size
    axis.title.y = element_text(size = 16),                           # Increased y-axis title size
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.key.height = unit(1.5, "cm"),
    strip.text = element_text(size = 16)                              # Increased facet label size
  ) +
  facet_wrap(~Date, labeller = labeller(Date = function(x) paste("Year CE", x)), ncol = 7)  # Add "Year" before each Date
ggsave("heatmap_Poaceae.png", width = 18, height = 20, dpi = 350)


