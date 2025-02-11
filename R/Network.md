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


