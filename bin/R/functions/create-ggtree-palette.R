create_tree_palette <- function(input,lin){

  df <- input |> 
    count(d5.10.15) |> # count unique superclusters
    arrange(desc(n)) |> # sort high to low
    select(supercluster=d5.10.15,n) |> # rename to supercluster
    filter(supercluster != "NA") |> # remove NAs (ungrouped)
    mutate(d5.10.15=supercluster) |> # duplicate supercluster ID
    separate_wider_delim(cols = d5.10.15, delim = ".",
                        names = c("d5", "d10", "d15"),
                        too_many = "merge") |>
    separate_wider_delim(cols = d15, delim = "-",
                        names = c("d15"),
                        too_many = "drop")

# Generate a unique color for each supercluster using randomcoloR
unique_clusters <- unique(df$supercluster)
n_clusters <- length(unique_clusters)
palette <- distinctColorPalette(n_clusters)  # Generates distinct random colors

# Create a color palette for each supercluster
supercluster_palette <- tibble(
  supercluster = unique_clusters,
  color = palette
)

# Reshape the original `df` to long format and join with `supercluster_palette`
component_mapping_df <- df %>% 
  pivot_longer(cols = starts_with("d"), names_to = "distance", values_to = "component") %>%
  left_join(supercluster_palette, by = "supercluster")
  
# Prioritize the color of the largest supercluster for each component-distance combination
component_palette <- component_mapping_df %>%
  group_by(component, distance) %>%
  slice_max(n, with_ties = FALSE) %>% # Keep the row with the largest supercluster
  ungroup() %>%
  select(component, distance, component_color = color) %>%
  distinct() # Ensure unique rows
  
# Return the lineage ID to the component cluster ID
component_palette <- component_palette %>%
  mutate(component = paste0(component, "-", lin))  # Append "-L{lineage}"
  
# Create a combined palette dataframe for reference, if needed
combined_palette <- bind_rows(
  supercluster_palette %>% select(name=supercluster, palette_color=color),
  component_palette %>% select(name=component, palette_color=component_color)
)
  
# Convert combined_palette to a named vector
color_palette <- setNames(combined_palette$palette_color, combined_palette$name)

# Return the color palette
return(color_palette)
  
}