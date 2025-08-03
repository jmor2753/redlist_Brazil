
## Jaccard distance
# Load libraries
library(tidyverse)
library(vegan)      # for Jaccard distances
library(reshape2)   # for pivoting if needed

# Read the data
df <- read.csv("spps_dataset.csv", sep = ";")

# View structure

df %>%
  count(scientificName, list) %>%
  filter(n > 1)
## removing duplicates
df_clean <- df %>%
  distinct(scientificName, list, .keep_all = TRUE)

df_clean %>%
  count(scientificName, list) %>%
  filter(n > 1)
# Create presence/absence matrix: species x list
presence_matrix <- df %>%
  mutate(present = 1) %>%
  distinct(scientificName, list, .keep_all = TRUE) %>%
  pivot_wider(names_from = list, values_from = present, values_fill = 0)

# Remove scientificName column
# Step 1: Build the transposed matrix with no NAs
mat <- presence_matrix %>%
  select(scientificName,`MMA 2022`:`IUCN 2024`) %>%
  t() %>% # Replace any NA with 0
  janitor::row_to_names(row_number = 1) %>%
  as.matrix() 


# Make sure it's numeric
storage.mode(mat) <- "numeric"

# Compute Jaccard distance
jaccard_dist <- vegdist(mat, method = "jaccard")
as.matrix(jaccard_dist)

## similarity
similarity_matrix <- 1 - as.matrix(jaccard_dist)




# Full pairwise bootstrap function
library(tidyverse)
library(vegan)

bootstrap_jaccard <- function(mat, reps = 1000) {
  list_names <- rownames(mat)
  results <- expand.grid(list1 = list_names, list2 = list_names, stringsAsFactors = FALSE)
  results <- results %>% filter(list1 < list2)  # only unique pairs
  
  boot_data <- results %>%
    rowwise() %>%
    mutate(
      ci = list({
        v1 <- mat[list1, ]
        v2 <- mat[list2, ]
        reps_dist <- replicate(reps, {
          idx <- sample(seq_along(v1), replace = TRUE)
          as.numeric(vegdist(rbind(v1[idx], v2[idx]), method = "jaccard"))
        })
        list(mean = mean(reps_dist), lower = quantile(reps_dist, 0.025), upper = quantile(reps_dist, 0.975))
      }),
      mean = ci$mean,
      lower = ci$lower,
      upper = ci$upper
    ) %>%
    select(-ci) %>%
    ungroup()
  
  return(boot_data)
}

# Run it
set.seed(23123)
ci_matrix <- bootstrap_jaccard(mat, reps = 1000)
ci_matrix

## similarity
ci_sim_matrix <- ci_matrix %>%
  mutate(
    mean_sim  = 1 - mean,
    lower_sim = 1 - upper,  # invert upper for lower
    upper_sim = 1 - lower   # invert lower for upper
  )
#write.csv(ci_sim_matrix, "Similarity_CI.csv")
## plotting

similarity_plot <- ggplot(ci_matrix, aes(x = list1, y = list2, fill = 1 - mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f\n[%.2f–%.2f]", 1- mean, 1-upper, 1-lower)), size = 5) +
  scale_fill_gradient(
    high = "darkolivegreen2",
    low = "pink1",
    name = "Similarity index\n[95% CI]"
  ) +
  theme_minimal() +
  labs(
    title = "",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid = element_blank(),
    aspect.ratio = 1.2
  )

similarity_plot

### Graphs ###
df_scatter <- read.delim("scatter_plot_data.csv", header = TRUE, encoding = "Latin1", sep=";")

df_scatter_prop <- df_scatter %>%
  group_by(class, order, list) %>%
  summarise(sum_exp = sum(exp),
            sum_obs = sum(obt)) %>%
  ungroup()


df_scatter_total <- df_scatter_prop %>%
  group_by(class, list) %>%
  summarise(total_exp = sum(sum_exp),
            total_obs = sum(sum_obs)) %>%
  left_join(., df_scatter_prop) %>%
  mutate(prop_exp = (sum_exp/total_exp),
         prop_obs = (sum_obs/total_obs)) %>%
  ungroup() %>%
  mutate(order_2 = ifelse(prop_exp < 0.005 | prop_obs < 0.005,
                          "Other (<0.5%)",
                          order))




# Assume your dataframe is called df
df_long <- df_scatter_total %>%
  pivot_longer(cols = c(prop_exp, prop_obs),
               names_to = "type",
               values_to = "proportion") %>%
  mutate(type = recode(type,
                       prop_exp = "Expected",
                       prop_obs = "Observed"))


df_long <- df_long %>%
  unite("class_type", class, type, sep = "_", remove = FALSE) %>%
  mutate(list_type = paste(list, type, sep = "_")) %>%
  mutate(x_label = paste(type, list, sep = "\n"))

## adjusting colours
orders_without_other <- unique(df_long$order_2)
orders_without_other <- orders_without_other[orders_without_other != "Other (<0.5%)"]

library(RColorBrewer)
library(colorspace)

# Step 1: Get qualitative palettes
qual_pals <- c("Set1", "Set2", "Set3", "Dark2", "Accent", "Paired")
col_vector <- unique(unlist(lapply(qual_pals, function(pal) {
  brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
})))

# Step 2: Remove gray-like colors
rgb_matrix <- col2rgb(col_vector)
is_gray <- apply(rgb_matrix, 2, function(x) max(x) - min(x) < 15)
col_vector <- col_vector[!is_gray]

# Step 3: Convert to LAB and cluster
lab_colors <- as(hex2RGB(col_vector), "LAB")@coords
set.seed(42)
k <- length(orders_without_other)
km <- kmeans(lab_colors, centers = k)

# Step 4: Choose representative color from each cluster
final_colors <- sapply(1:k, function(i) {
  col_vector[which(km$cluster == i)[1]]
})

# Step 5: Assign names and add "Other"
color_palette <- final_colors
names(color_palette) <- orders_without_other
color_palette <- c(color_palette, "Other (<0.5%)" = "black")

## reordering
# Reorder order_2 by mean proportion across all bars
order_levels <- df_long %>%
  group_by(order_2) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE)) %>%
  arrange(desc(mean_prop)) %>%
  pull(order_2)

# Apply this order to the factor
df_long$order_2 <- factor(df_long$order_2, levels = order_levels)

prop_plot <- ggplot(df_long, aes(x = type, y = proportion, fill = order_2)) +
  geom_bar(stat = "identity") +
  facet_grid(class ~ list) +
  scale_fill_manual(values = color_palette) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = "Order",
    title = ""
  ) +
  theme_minimal(base_size = 14) +  # ← increase text size
  theme(
    panel.grid.major.x = element_blank(),        # Remove grid
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 13),
    strip.text = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black", size = 14))





#### scatterplot


df_scatter_prop_perlist <- df_scatter %>%
  group_by(class, order, list, iucn_category) %>%
  summarise(sum_exp = sum(exp),
            sum_obs = sum(obt)) %>%
  ungroup()



df_scatter_total_perlist <- df_scatter_prop_perlist %>%
  group_by(class, order, list) %>%
  summarise(total_exp = sum(sum_exp),
            total_obs = sum(sum_obs)) %>%
  left_join(., df_scatter_prop_perlist) %>%
  mutate(prop_exp = (sum_exp/total_exp),
         prop_obs = (sum_obs/total_obs)) %>%
  ungroup() 




# Most critical to least worrying (including DD last)
iucn_order <- c("CR", "EN", "VU", "NT", "LC", "DD")

# reorder
df_scatter_total_perlist <- df_scatter_total_perlist %>%
  mutate(iucn_category = factor(iucn_category, levels = iucn_order))


library(ggh4x) ## for nested facets

df_scatter_total_perlist$prop_obs <- ifelse(is.nan(df_scatter_total_perlist$prop_obs), 0, df_scatter_total_perlist$prop_obs)

# Add labels to selected points
df_labeled <- df_scatter_total_perlist %>%
  group_by(list, iucn_category ) %>%
  mutate(diff = abs(sum_exp - sum_obs),
         label_text = ifelse(diff == 1.000 & sum_obs > 10, order, NA)) %>%
  mutate(label_text = ifelse(class == "Insecta", label_text, NA)) # initial labeling where expected ≈ observed

# Select 2 points where expected > observed (largest gaps)
set.seed(2310)
top_over_exp <- df_scatter_total_perlist %>%
  mutate(gap = sum_exp - sum_obs) %>%
  filter(gap > 50, iucn_category %in% c("CR", "EN", "VU", "NT")) %>%
  group_by(list, iucn_category) %>%
  slice_max(order_by = gap, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(order, iucn_category) %>%
  slice_max(order_by = gap, n = 1, with_ties = FALSE) %>%  # keep top 1 per order
  ungroup() %>%
  slice_sample(n = 6)

# Combine both sets
df_labels_combined <- bind_rows(
  df_labeled %>% filter(!is.na(label_text)),
  top_over_exp %>% mutate(label_text = order)
) %>%
  distinct()

# Combine both sets with color flag
df_labels_combined <- bind_rows(
  df_labeled %>%
    filter(!is.na(label_text)) %>%
    mutate(label_color = "black"),
  top_over_exp %>%
    mutate(label_text = order,
           label_color = "red")
) %>%
  distinct()


## final plot

diagplot_labeled <- ggplot(df_scatter_total_perlist, aes(x = sqrt(sum_exp), y = sqrt(sum_obs),
                                                         fill = class, shape = list)) + 
  geom_point(size = 2.5, stroke = 0.5, color = "black", alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.5) +
  geom_text_repel(data = df_labels_combined,
                  aes(x = sqrt(sum_exp), y = sqrt(sum_obs),
                      label = label_text, color = label_color),
                  size = 3,
                  max.overlaps = 15,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = "gray30",
                  inherit.aes = FALSE,
                  show.legend = FALSE)  + # hide label color from legend 
  labs(
    x = expression(sqrt(Expected~count)),
    y = expression(sqrt(Observed~count))) + 
  ggh4x::facet_nested(iucn_category ~ class, scales = "free") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 15, color = "black"),
    axis.text = element_text(size = 13, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(),
    strip.text = element_text(color = "black"),
    strip.text.x.top = element_text(size = 17, color = "black"),
    strip.background = element_rect(fill = NA, color = "black"),
    strip.background.x = element_rect(fill = "grey95"),
    strip.background.y = element_rect(fill = "grey92"),
    strip.text.x.bottom = element_text(size = 13, color = "black"),
    strip.text.y.right = element_text(size = 15, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(color = "black", size = 15),
    legend.title = element_text(color = "black", size = 15)
  ) + 
  scale_fill_discrete(name = "Class") +
  scale_shape_manual(
    name = "List",
    values = c("ICMBIO" = 21, "MMA" = 24, "IUCN" = 22)
  ) + guides(fill = "none") +
  scale_color_manual(
    values = c("red" = "firebrick3", "black" =" black")
  ) 

diagplot_labeled



#write.csv(df_labels_combined[1:10], "TableS1_Lopesetal.csv")


#####

library(dplyr)
library(DescTools)



ccc_results <- df_scatter_total_perlist %>%
  group_by(class, list) %>%
  filter(n() >= 2) %>%  # require at least 2 points to compute CCC
  reframe(
    ccc = tryCatch(
      {
        if (sd(sum_exp) > 0 & sd(sum_obs) > 0) {
          CCC(sum_exp, sum_obs)$rho.c
        } else {
          NA
        }
      },
      error = function(e) NA
    ),
    n = n()
  )

ccc_results_df <- ccc_results %>%
  mutate(ccc_est = ccc_results$ccc$est,
         ccc_lwr = ccc_results$ccc$lwr.ci,
         ccc_upr = ccc_results$ccc$upr.ci) %>%
  select(-ccc)


## final table 
library(tidyr)
library(dplyr)

ccc_wide <- ccc_results_df %>%
  pivot_wider(
    id_cols = class,
    names_from = list,
    values_from = c(ccc_est, ccc_lwr, ccc_upr),
    names_glue = "{list}_{.value}"
  )

#write.csv(ccc_wide, "CCC_CI.csv")

# STEP 2: Plot CCC results
cccplot <- ccc_results_df %>%
  filter(!is.na(ccc_est)) %>%
  ggplot(aes(x = reorder(class, ccc_est), y = ccc_est, fill = class)) +
  geom_errorbar(aes(ymin = ccc_lwr, ymax = ccc_upr),
                width = 0.1,
                color = "black",
                position = position_dodge(width = 0.8)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, col = "black", alpha = 0.5) + 
  #geom_point(size = 1.5, stroke = 0.5, color = "black", position = position_dodge(width = 0.8), pch = 21, alpha = 0.5) +
  ggh4x::facet_nested(~ list, scales = "free_y") +
  scale_y_continuous(limits = c(-0.1, 1.1), expand = c(0, 0)) +
  labs(
    x = "Class",
    y = "Concordance correlation coefficient (CCC)\n[95% CI]",
    fill = "Class",
    shape = "Class"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.2),
    axis.title = element_text(color = "black", size = 13),
    strip.text = element_text(color = "black", size = 13, hjust = 0.5),
    strip.background = element_rect(fill = "grey95", color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    aspect.ratio = 3
  )

cccplot
similarity_plot
prop_plot 
diagplot 



### supplementary Table S2



ccc_results_percat <- df_scatter_total_perlist %>%
  group_by(class, iucn_category) %>%
  filter(n() >= 2) %>%  # require at least 2 points to compute CCC
  reframe(
    ccc = tryCatch(
      {
        if (sd(sum_exp) > 0 & sd(sum_obs) > 0) {
          CCC(sum_exp, sum_obs)$rho.c
        } else {
          NA
        }
      },
      error = function(e) NA
    ),
    n = n()
  )

ccc_results_df_percat <- ccc_results_percat %>%
  mutate(ccc_est = ccc_results_percat$ccc$est,
         ccc_lwr = ccc_results_percat$ccc$lwr.ci,
         ccc_upr = ccc_results_percat$ccc$upr.ci) %>%
  select(-ccc)

#write.csv(ccc_results_df_percat, "TableS2_Lopeetal.csv")
