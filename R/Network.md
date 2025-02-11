# Networks :microbe:

**Goal:** Traditional microbial ecology typically focuses on microbial shifts of particular bacteria. Networks are an opportunity to look at how the microbiome is organized and the interactions or co-occurence of bacteria.


**Determine the effect of replicate size on network metrics:**
```{r}
# Load required libraries
library(phyloseq)
library(NetCoMi)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)  # For data reshaping

# Set the directory where phyloseq objects are stored
save_dir <- "/Users/nicholas.macknight/Desktop/Ofav SCTLD/R/NetCoMi/SubsettedNetworks"
phyloseq_objects <- ComparingSusceptibilityList  # Placeholder for actual loading of phyloseq objects

# Define the range of sample sizes
sample_sizes <- seq(10, 200, by = 10)

# Initialize an empty list to store metrics for all objects
all_metrics <- list()

# Loop through each phyloseq object
for (physeq_obj_name in names(phyloseq_objects)) {
  physeq_obj <- get(physeq_obj_name)  # This assumes that 'get()' fetches the actual Phyloseq object
  
  # Skip if the number of samples is less than the smallest sample size
  if (nsamples(physeq_obj) < min(sample_sizes)) {
    next
  }
  
  # Initialize a data frame to store the metrics for this object
  metrics_df <- data.frame(
    Sample_Size = integer(),
    Edge_Density = numeric(),
    Average_Degree = numeric(),
    Modularity = numeric(),
    Average_Geodesic = numeric(),
    Network_Diameter = numeric(),
    Average_Closeness = numeric(),
    Average_Betweenness = numeric(),
    Vertices = integer(),
    stringsAsFactors = FALSE
  )
  
  # Check and process each sample size
  for (n in sample_sizes) {
    if (nsamples(physeq_obj) >= n) {
      # Sample and calculate metrics
      sampled_physeq_obj <- prune_samples(sample(sample_names(physeq_obj), n, replace = FALSE), physeq_obj)
      physeq_order <- tax_glom(sampled_physeq_obj, taxrank = "Order")
      physeq_order_renamed <- renameTaxa(physeq_order, numDupli = "Order")

      net_order <- netConstruct(physeq_order_renamed, taxRank = "Order", measure = "pearson",
                                zeroMethod = "multRepl", normMethod = "clr", sparsMethod = "threshold", thresh = 0.5)
      graph_net <- igraph::graph_from_adjacency_matrix(net_order$adjaMat1, weighted = TRUE, mode = "undirected")

      # Calculate metrics
      metrics <- data.frame(
        Sample_Size = n,
        Edge_Density = igraph::edge_density(graph_net),
        Average_Degree = mean(degree(graph_net)),
        Modularity = max(cluster_fast_greedy(graph_net)$modularity),
        Network_Diameter = igraph::diameter(graph_net, directed = FALSE, weights = NA),
        Average_Closeness = mean(closeness(graph_net, mode = "all"), na.rm = TRUE),
        Average_Betweenness = mean(betweenness(graph_net, directed = FALSE, weights = NA), na.rm = TRUE),
        Vertices = vcount(graph_net)
      )

      # Append the results to the data frame
      metrics_df <- rbind(metrics_df, metrics)
    }
  }
  
  # Store the compiled metrics in the list
  all_metrics[[physeq_obj_name]] <- metrics_df
}

# Combine all the metrics into a single dataframe
combined_metrics <- bind_rows(all_metrics, .id = "Phyloseq_Object")

# Optionally, write combined metrics to a CSV file
write.csv(combined_metrics, file = paste0(save_dir, "/iterated_network_metrics_10-200.csv"), row.names = FALSE)


# Reshape for plotting
long_metrics <- pivot_longer(combined_metrics, cols = c("Average_Degree","Edge_Density", "Modularity","Network_Diameter","Average_Closeness","Average_Betweenness"), 
                             names_to = "Metric", values_to = "Value")

# Plotting
ggplot(long_metrics, aes(x = Sample_Size, y = Value, color = Phyloseq_Object, group = interaction(Phyloseq_Object, Metric))) +
  geom_line() +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Network Metrics Across Different Sample Sizes", x = "Sample Size", y = "Metric Value") +
  theme(legend.position = "right")


# Find the last sample size for label placement
last_points <- long_metrics %>%
  group_by(Phyloseq_Object, Metric) %>%
  slice(which.max(Sample_Size))

# Plotting
ggplot(long_metrics, aes(x = Sample_Size, y = Value, color = Phyloseq_Object, group = interaction(Phyloseq_Object, Metric))) +
  geom_line() +
  geom_text(data = last_points, aes(label = Phyloseq_Object), nudge_x = 1, hjust = 0, vjust = 0, size = 3, show.legend = FALSE) +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Network Metrics Across Different Sample Sizes", x = "Sample Size", y = "Metric Value") +
  theme(legend.position = "right")
# Save the plot
#ggsave("network_metrics_plot.png", width = 10, height = 8, units = "in")


# Just Modularity
# Reshape for plotting
long_metrics <- pivot_longer(combined_metrics, cols = "Modularity", 
                             names_to = "Metric", values_to = "Value")

# Find the last sample size for label placement
last_points <- long_metrics %>%
  group_by(Phyloseq_Object, Metric) %>%
  slice(which.max(Sample_Size))

# Plotting
ggplot(long_metrics, aes(x = Sample_Size, y = Value, color = Phyloseq_Object, group = Phyloseq_Object)) +
  geom_line() +
  geom_text(data = last_points, aes(label = Phyloseq_Object), nudge_x = 1, hjust = 0, vjust = 0, size = 3, show.legend = FALSE) +
  theme_minimal() +
  labs(title = "Modularity Across Different Sample Sizes", x = "Sample Size", y = "Modularity") +
  theme(legend.position = "right")

```
<div align="center">
<img width="708" alt="image" src="https://github.com/user-attachments/assets/a6da00ee-e4b2-4fc5-8c3e-40020cd8ba20" />
</div>

> Conclusion: There is definitely a strong effect of sample size and the metrics do begin to plateua, certainly by 100 samples, arguable by 50 we see saturation dependent on the metric. So because some of these phyloseq objects have 8-14 samples, we will need to curate larger phyloseq objects and revisit the questions we can ask.



# Comparing Susceptibility
> ps_C_Resistant has 39 samples as our minimum sample size.
> So if we were to sample 39 samples from each of the four datasets, here is what would happen:
```{r}
# Load required libraries
library(phyloseq)
library(NetCoMi)
library(igraph)
library(dplyr)  # For data manipulation
library(ggplot2)  # For plotting

# Set the directory where plots will be saved
save_dir <- "/Users/nicholas.macknight/Desktop/Ofav SCTLD/R/NetCoMi/SubsettedNetworks"

# Ensure the directory exists and create it if it does not
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# Original character vector
ComparingSusceptibility <- c("ps_C_Resistant", "ps_C_Intermediate", "ps_C_Susceptible", "ps_C_HSusceptible")

# Convert to a named list with the same values for simplicity
ComparingSusceptibilityList <- setNames(as.list(ComparingSusceptibility), ComparingSusceptibility)

# Print the new list structure
print(ComparingSusceptibilityList)

# Initialize an empty data frame to store the metrics for all phyloseq objects
network_metrics_ComparingSusceptibility <- data.frame(
  Phyloseq_Object = character(),
  Edge_Density = numeric(),
  Average_Degree = numeric(),
  Modularity = numeric(),
  Connectance = numeric(),
  Average_Geodesic = numeric(),
  Network_Diameter = numeric(),
  Average_Closeness = numeric(),
  Average_Betweenness = numeric(),
  Vertices = integer(),
  stringsAsFactors = FALSE
)

# Loop through each phyloseq object specified in ComparingSusceptibility
for (physeq_obj_name in ComparingSusceptibility) {
  # Get the phyloseq object
  physeq_obj <- get(physeq_obj_name)
  
  # Randomly sample 39 samples
  sampled_physeq_obj <- prune_samples(sample(sample_names(physeq_obj), 39, replace = FALSE), physeq_obj)
  
  # Agglomerate to Order level and rename taxa to ensure unique identifiers
  physeq_order <- tax_glom(sampled_physeq_obj, taxrank = "Order")
  physeq_order_renamed <- renameTaxa(physeq_order, 
                                     pat = "<name>", 
                                     substPat = "<name>_<subst_name>(<subst_R>)",
                                     numDupli = "Order")
  
  # Perform Network Construction
  net_order <- netConstruct(physeq_order_renamed,
                            taxRank = "Order",
                            measure = "pearson",
                            zeroMethod = "multRepl",
                            normMethod = "clr",
                            sparsMethod = "threshold",
                            thresh = 0.4,  # Using fixed threshold
                            verbose = 3)
  
  # Convert to igraph object
  graph_net <- igraph::graph_from_adjacency_matrix(net_order$adjaMat1, weighted = TRUE, mode = "undirected")

  # Calculate various network metrics
  edge_density <- igraph::graph.density(graph_net)
  network_diameter <- igraph::diameter(graph_net, directed = FALSE, weights = NA)
  avg_geodesic <- mean(igraph::distances(graph_net, mode = "all"), na.rm = TRUE)
  avg_degree <- mean(degree(graph_net))
  avg_closeness <- mean(closeness(graph_net, mode = "all"), na.rm = TRUE)
  avg_betweenness <- mean(betweenness(graph_net, directed = FALSE, weights = NA), na.rm = TRUE)
  connectance <- igraph::graph.density(graph_net)
  modularity <- max(cluster_fast_greedy(graph_net)$modularity)
  vertices <- vcount(graph_net)

  # Append the results to the data frame
  network_metrics_ComparingSusceptibility <- network_metrics_ComparingSusceptibility %>%
    add_row(
      Phyloseq_Object = physeq_obj_name,
      Edge_Density = edge_density,
      Average_Degree = avg_degree,
      Modularity = modularity,
      Connectance = connectance,
      Average_Geodesic = avg_geodesic,
      Network_Diameter = network_diameter,
      Average_Closeness = avg_closeness,
      Average_Betweenness = avg_betweenness,
      Vertices = vertices
    )
  
  # Print completion message
  print(paste("Completed network analysis for:", physeq_obj_name))
}

# Print completion message for all metrics collection
print("All network metrics collected.")

# Optionally, write the combined metrics to a CSV file
write.csv(network_metrics_ComparingSusceptibility, file = paste0(save_dir, "/network_metrics_ComparingSusceptibility.csv"), row.names = FALSE)

# Print the path to the saved CSV
print(paste("Combined network metrics saved to:", paste0(save_dir, "/network_metrics_ComparingSusceptibility.csv")))

# Return or print the final data frame for verification in the console
print(network_metrics_ComparingSusceptibility)

# Plotting each metric in a separate plot
metrics_to_plot <- c("Edge_Density", "Average_Degree", "Average_Betweenness", "Average_Closeness", "Network_Diameter", "Vertices", "Modularity", "Average_Geodesic", "Connectance")

network_metrics_ComparingSusceptibility$Phyloseq_Object <- factor(network_metrics_ComparingSusceptibility$Phyloseq_Object,
                                                                  levels = ComparingSusceptibility)

for (metric in metrics_to_plot) {
  plot <- ggplot(network_metrics_ComparingSusceptibility, aes(x = Phyloseq_Object, y = .data[[metric]], group = Phyloseq_Object, color = Phyloseq_Object)) +
    geom_line() +
    geom_point() +
    labs(title = paste("Comparison of", metric, "across Susceptiblity Categories"),
         x = "Stage of Susceptibility", y = metric) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")  # No need for legend if colors match the objects
  
  # Optionally, print the plot in the console
  print(plot)
}

```
Results of sampling 39 samples per resilience category:
<div align="center">
<img width="300" alt="image" src="https://github.com/user-attachments/assets/ee55878a-1440-4e1d-826c-f5fcbf4b3fda" /> <img width="300" alt="image" src="https://github.com/user-attachments/assets/b0e7c866-d504-4253-81ae-259defc45c11" /> <img width="300" alt="image" src="https://github.com/user-attachments/assets/d4f47b1f-68da-4fe6-98ce-513de6f3b2a0" /> <img width="300" alt="image" src="https://github.com/user-attachments/assets/8c6af032-e4ea-4dd1-9baa-064e4a809bad" /> <img width="300" alt="image" src="https://github.com/user-attachments/assets/76b0eff2-fda3-4d7e-ac25-044ef3a74772" /> <img width="300" alt="image" src="https://github.com/user-attachments/assets/1c8a2a75-cca1-4292-a9d4-f02848251424" /> <img width="300" alt="image" src="https://github.com/user-attachments/assets/3399a2a1-53fe-473e-b291-8c5d2de963db" /> <img width="300" alt="image" src="https://github.com/user-attachments/assets/e7ae46b1-13be-40f6-90bb-6f4ae7fe1127" />
</div>

> The single dots are not aesthetic, but networks are not iterative and do not produce variation when normally applied.

# Comparing Susceptibility - 10 Iterations
> 10 iterations each of subsampled phyloseq objects at 39 samples
```{r}
# Load required libraries
library(phyloseq)
library(NetCoMi)
library(igraph)
library(dplyr)
library(ggplot2)

# Set directories
save_dir <- "/Users/nicholas.macknight/Desktop/Ofav SCTLD/R/NetCoMi/SubsettedNetworks"
full_metrics_path <- "/Users/nicholas.macknight/Desktop/Ofav SCTLD/R/NetCoMi/network_metrics_newcalc.csv"

# Ensure directory exists
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# Load full dataset metrics for comparison
full_metrics <- read.csv(full_metrics_path)

# List of phyloseq object names
phyloseq_objects <- ComparingSusceptibility

# Initialize a list to store metrics data frames for all objects
all_metrics <- list()

# Loop through each phyloseq object
for (physeq_obj_name in phyloseq_objects) {
  # Retrieve the phyloseq object
  physeq_obj <- get(physeq_obj_name)
  
  # Initialize a data frame to store metrics for current object
  metrics_df <- data.frame(
  Phyloseq_Object = character(),
  Edge_Density = numeric(),
  Average_Degree = numeric(),
  Modularity = numeric(),
  Average_Geodesic = numeric(),
  Network_Diameter = numeric(),
  Average_Closeness = numeric(),
  Average_Betweenness = numeric(),
  Vertices = integer(),
  stringsAsFactors = FALSE
)
  
  # Perform subsampling and calculate metrics 10 times
  for (i in 1:10) {
    sampled_physeq_obj <- prune_samples(sample(sample_names(physeq_obj), 39, replace = FALSE), physeq_obj)
    physeq_order <- tax_glom(sampled_physeq_obj, taxrank = "Order")
    physeq_order_renamed <- renameTaxa(physeq_order, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)", numDupli = "Order")
    
    net_order <- netConstruct(physeq_order_renamed, taxRank = "Order", measure = "pearson", zeroMethod = "multRepl", normMethod = "clr", sparsMethod = "threshold", thresh = 0.5)
    graph_net <- igraph::graph_from_adjacency_matrix(net_order$adjaMat1, weighted = TRUE, mode = "undirected")
    
    # Calculate network metrics
    metrics_df <- rbind(metrics_df, data.frame(
      Iteration = i,
      Edge_Density = igraph::edge_density(graph_net),
      Average_Degree = mean(degree(graph_net)),
      Modularity = max(cluster_fast_greedy(graph_net)$modularity),
      Average_Geodesic = mean(distances(graph_net, mode = "all"), na.rm = TRUE),
      Network_Diameter <- igraph::diameter(graph_net, directed = FALSE, weights = NA),
      Average_Closeness <- mean(closeness(graph_net, mode = "all"), na.rm = TRUE),
      Average_Betweenness <- mean(betweenness(graph_net, directed = FALSE, weights = NA), na.rm = TRUE),
      Vertices <- vcount(graph_net)
    ))
  }
  
  
  # Store metrics dataframe in list
  all_metrics[[physeq_obj_name]] <- metrics_df
}

# Combine all metrics into a single dataframe for plotting
combined_metrics <- bind_rows(all_metrics, .id = "Phyloseq_Object")

# Clean up the column names
names(combined_metrics) <- gsub("\\..*$", "", names(combined_metrics))

# Now your column names should be cleaned up
print(names(combined_metrics))

library(dplyr)
library(ggplot2)
library(ggpubr)


# Set the order of the Phyloseq Objects as desired
combined_metrics$Phyloseq_Object <- factor(combined_metrics$Phyloseq_Object, levels = c("ps_C_Resistant", "ps_C_Intermediate", "ps_C_Susceptible", "ps_C_HSusceptible"))

# Plotting with boxplots grouped by Phyloseq Object
ggplot(combined_metrics, aes(x = Phyloseq_Object, y = Modularity, fill = Phyloseq_Object)) +
  geom_boxplot() +
  theme(legend.position = "none")
  

# > Conclusion: sample size seems to have an effect on network metrics. Next we will see when there is a sufficient sample size to reach a saturation point in network metrics. 



# Load required libraries
library(ggplot2)
library(tidyr)  # For pivoting data

# Assuming 'combined_metrics' contains other metrics as well
# Reshape the data to a long format where each row is a unique metric for a given Phyloseq_Object
long_format_metrics <- pivot_longer(
  combined_metrics,
  cols = -Phyloseq_Object,  # Exclude the Phyloseq_Object column from the transformation
  names_to = "Metric",
  values_to = "Value"
)

# Plotting with boxplots grouped by Phyloseq Object for each metric
p <- ggplot(long_format_metrics, aes(x = Phyloseq_Object, y = Value, fill = Phyloseq_Object)) +
  geom_boxplot() +
  facet_wrap(~ Metric, scales = "free_y") +  # Create a separate plot for each metric
  theme_minimal() +
  labs(title = "Variation in Network Metrics Across Phyloseq Objects", x = "Stage of Susceptibility", y = "Metric Value") +
  theme(legend.position = "none")  # Remove the legend

# Print the plot
print(p)

```

### Calculating Statistical signficance between susceptibilty categories for each network metric. Manually
```{r}
# Load necessary library
library(dplyr)

# Names of metrics to analyze
metrics <- colnames(combined_metrics)[-c(1, 2,6)]  # Exclude 'Phyloseq_Object' and 'Iteration'

# Initialize a list to store results
results_list <- list()

# Loop through each metric
for (metric in metrics) {
  # ANOVA to test differences across Phyloseq objects
  anova_results <- aov(reformulate("Phyloseq_Object", response = metric), data = combined_metrics)
  summary_anova <- summary(anova_results)
  
  # Extract the p-value correctly from the ANOVA summary
  p_value <- summary_anova[[1]]["Phyloseq_Object", "Pr(>F)"]

  # Check if the ANOVA is significant
  if (p_value < 0.05) {
    # If significant, proceed with Tukey HSD test
    tukey_results <- TukeyHSD(anova_results)
    
    # Convert to dataframe and filter for significant comparisons
    significant_comparisons <- as.data.frame(tukey_results[["Phyloseq_Object"]])
    significant_comparisons <- significant_comparisons %>% 
      filter(`p adj` < 0.05) %>%
      mutate(Metric = metric)  # Add the metric name for clarity

    # Store results
    results_list[[metric]] <- significant_comparisons
  } else {
    results_list[[metric]] <- data.frame(Metric = metric, Comparison = NA, Difference = NA, `Lower Bound` = NA, `Upper Bound` = NA, `p adj` = NA)
  }
}

results_list
```
### Manually adding signficance annotations to plots
```{r}
# Set directory for saving plots
plot_save_dir <- "/Users/nicholas.macknight/Desktop/Ofav SCTLD/R/NetCoMi/SubsettedNetworks/"

# Ensure the directory exists, create if it does not
if (!dir.exists(plot_save_dir)) {
  dir.create(plot_save_dir, recursive = TRUE)
}

max_y <- max(combined_metrics$Average_Degree, na.rm = TRUE)
increment <- (max_y - min(combined_metrics$Average_Degree, na.rm = TRUE)) * 0.1  # Adjust increment as needed
Average_Degree_plot <- ggplot(combined_metrics, aes(x = Phyloseq_Object, y = Average_Degree, fill = Phyloseq_Object)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_signif(comparisons = list(c("ps_D_pre", "ps_D_ini"),c("ps_D_pre", "ps_D_fin"),c("ps_D_ini", "ps_D_ear"),c("ps_D_fin", "ps_D_ear")),y_position = c(max_y + increment, max_y + 2 * increment, max_y + 3 * increment, max_y + 4 * increment),tip_length = 0.05) +
  theme_minimal() +
  labs(title = "Average Degree Throughout Disease Exposure", subtitle="The number of connections that each bacteria in the network has with other bacteria", x = "Exposure Timeline", y = "Average Degree") +
  theme(legend.position = "none")
Average_Degree_plot
ggsave(filename = paste0(plot_save_dir, "/ComparingSusceptibility_Average_Degree_Plot.png"), plot = Average_Degree_plot, width = 4, height = 4, units = "in")
ggsave(filename = paste0(plot_save_dir, "/ComparingSusceptibility_Average_Degree_Plot.pdf"), plot = Average_Degree_plot, width = 4, height = 4, units = "in")


max_y <- max(combined_metrics$Edge_Density, na.rm = TRUE)
increment <- (max_y - min(combined_metrics$Edge_Density, na.rm = TRUE)) * 0.1  # Adjust increment as needed
Edge_Density_plot <- ggplot(combined_metrics, aes(x = Phyloseq_Object, y = Edge_Density, fill = Phyloseq_Object)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_signif(comparisons = list(c("ps_D_pre", "ps_D_ini"),c("ps_D_pre", "ps_D_fin"),c("ps_D_ini", "ps_D_ear"),c("ps_D_fin", "ps_D_ear")),y_position = c(max_y + increment, max_y + 2 * increment, max_y + 3 * increment, max_y + 4 * increment),tip_length = 0.05) +
  theme_minimal() +
  labs(title = "Edge Density Throughout Disease Exposure", subtitle="Actual connections relative to possible connections.", x = "Exposure Timeline", y = "Edge Density") +
  theme(legend.position = "none")
Edge_Density_plot
ggsave(filename = paste0(plot_save_dir, "/ComparingSusceptibility_Edge_Density_Plot.png"), plot = Edge_Density_plot, width = 4, height = 4, units = "in")
ggsave(filename = paste0(plot_save_dir, "/ComparingSusceptibility_Edge_Density_Plot.pdf"), plot = Edge_Density_plot, width = 4, height = 4, units = "in")


max_y <- max(combined_metrics$Modularity, na.rm = TRUE)
increment <- (max_y - min(combined_metrics$Modularity, na.rm = TRUE)) * 0.1  # Adjust increment as needed
Modularity_plot <- ggplot(combined_metrics, aes(x = Phyloseq_Object, y = Modularity, fill = Phyloseq_Object)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_signif(comparisons = list(c("ps_D_pre", "ps_D_ear"),c("ps_D_pre", "ps_D_ini"),c("ps_D_pre", "ps_D_fin"),c("ps_D_ini", "ps_D_ear"),c("ps_D_fin", "ps_D_ear")),y_position = c(max_y + increment, max_y + 2 * increment, max_y + 3 * increment, max_y + 4 * increment),tip_length = 0.05) +
  theme_minimal() +
  labs(title = "Modularity Throughout Disease Exposure", subtitle="Higher modularity means higher compartmentalization, or less globally connected",x = "Exposure Timeline", y = "Modularity") +
  theme(legend.position = "none")
Modularity_plot
ggsave(filename = paste0(plot_save_dir, "/ComparingSusceptibility_Modularity_Plot.png"), plot = Modularity_plot, width = 4, height = 4, units = "in")
ggsave(filename = paste0(plot_save_dir, "/ComparingSusceptibility_Modularity_Plot.pdf"), plot = Modularity_plot, width = 4, height = 4, units = "in")

max_y <- max(combined_metrics$Network_Diameter, na.rm = TRUE)
increment <- (max_y - min(combined_metrics$Network_Diameter, na.rm = TRUE)) * 0.1  # Adjust increment as needed
Network_Diameter_plot <- ggplot(combined_metrics, aes(x = Phyloseq_Object, y = Network_Diameter, fill = Phyloseq_Object)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  #geom_signif(comparisons = list(c("ps_D_pre", "ps_D_ear"),c("ps_D_pre", "ps_D_ini"),c("ps_D_pre", "ps_D_fin"),c("ps_D_ini", "ps_D_ear"),c("ps_D_fin", "ps_D_ear")),y_position = c(max_y + increment, max_y + 2 * increment, max_y + 3 * increment, max_y + 4 * increment),tip_length = 0.05) +
  theme_minimal() +
  labs(title = "Network Diameter Throughout Disease Exposure", subtitle="The shortest path between the two farthest points in a network",x = "Exposure Timeline", y = "Network Diameter") +
  theme(legend.position = "none")
Network_Diameter_plot
ggsave(filename = paste0(plot_save_dir, "/ComparingSusceptibility_Network_Diameter_Plot.png"), plot = Network_Diameter_plot, width = 4, height = 4, units = "in")
ggsave(filename = paste0(plot_save_dir, "/ComparingSusceptibility_Network_Diameter_Plot.pdf"), plot = Network_Diameter_plot, width = 4, height = 4, units = "in")

max_y <- max(combined_metrics$Average_Closeness, na.rm = TRUE)
increment <- (max_y - min(combined_metrics$Average_Closeness, na.rm = TRUE)) * 0.1  # Adjust increment as needed
Average_Closeness_plot <- ggplot(combined_metrics, aes(x = Phyloseq_Object, y = Average_Closeness, fill = Phyloseq_Object)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_signif(comparisons = list(c("ps_D_pre", "ps_D_ini"),c("ps_D_pre", "ps_D_fin"),c("ps_D_ini", "ps_D_ear"),c("ps_D_fin", "ps_D_ear")),y_position = c(max_y + increment, max_y + 2 * increment, max_y + 3 * increment, max_y + 4 * increment),tip_length = 0.05) +
  theme_minimal() +
  labs(title = "Average Closeness Throughout Disease Exposure", subtitle= "How close a node is to all other nodes in the network", x = "Exposure Timeline", y = "Average Closeness") +
  theme(legend.position = "none")
Average_Closeness_plot
ggsave(filename = paste0(plot_save_dir, "/ComparingSusceptibility_Average_Closeness_Plot.png"), plot = Average_Closeness_plot, width = 4, height = 4, units = "in")
ggsave(filename = paste0(plot_save_dir, "/ComparingSusceptibility_Average_Closeness_Plot.pdf"), plot = Average_Closeness_plot, width = 4, height = 4, units = "in")

max_y <- max(combined_metrics$Average_Betweenness, na.rm = TRUE)
increment <- (max_y - min(combined_metrics$Average_Betweenness, na.rm = TRUE)) * 0.1  # Adjust increment as needed
Average_Betweenness_plot <- ggplot(combined_metrics, aes(x = Phyloseq_Object, y = Average_Betweenness, fill = Phyloseq_Object)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_signif(comparisons = list(c("ps_D_pre", "ps_D_ini"),c("ps_D_pre", "ps_D_fin"),c("ps_D_ini", "ps_D_ear"),c("ps_D_fin", "ps_D_ear")),y_position = c(max_y + increment, max_y + 2 * increment, max_y + 3 * increment, max_y + 4 * increment),tip_length = 0.05) +
  theme_minimal() +
  labs(title = "Average Betweenness Throughout Disease Exposure", subtitle= "When a node lies on paths between other nodes and acts as bridge points", x = "Exposure Timeline", y = "Average Betweenness") +
  theme(legend.position = "none")
Average_Betweenness_plot
ggsave(filename = paste0(plot_save_dir, "/ComparingSusceptibility_Average_Betweenness_Plot.png"), plot = Average_Betweenness_plot, width = 4, height = 4, units = "in")
ggsave(filename = paste0(plot_save_dir, "/ComparingSusceptibility_Average_Betweenness_Plot.pdf"), plot = Average_Betweenness_plot, width = 4, height = 4, units = "in")


```
Results of 10 iterations of sampling 39 replicates:

> The problrm here is that the Resilient dataset is being completely sampled each time and has no variation which is a bias when performing...an analysis of variation! So to solve this, less replicates than 39 will need to be sampled, but which number should be chosen that isnt arbitrary? 

<div align="center">
<img width="400" alt="image" src="https://github.com/user-attachments/assets/f189ce45-a560-4a14-b9b0-7468904d2f7d" /> <img width="400" alt="image" src="https://github.com/user-attachments/assets/7d41b242-7321-4b14-885c-e6d1fa6591be" /> <img width="400" alt="image" src="https://github.com/user-attachments/assets/a1caba71-b0a3-4cc4-bfdd-56f8afb826b3" /> <img width="400" alt="image" src="https://github.com/user-attachments/assets/8ed3a353-a9c8-4a2c-9895-70a268a3b343" /> <img width="400" alt="image" src="https://github.com/user-attachments/assets/b7dbc19f-2fa8-408c-97da-53e6333e3943" /> <img width="400" alt="image" src="https://github.com/user-attachments/assets/4dc764f5-d80f-468d-aa33-7824062f5c21" />
</div>



> Concluding thoughts. in each network method appraoch, all samples from each dataset, using the limiting datasets number of samples (39 replicates from Resistant), and using the iterations, the resistant networks are always distinct. This gives me hope that there is a true biological distinction here and we need to figure out a way to confidently remove technical bias in how the networks were constructed. time to dig into the literature. 






