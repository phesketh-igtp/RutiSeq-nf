create_tree_palette <- function(input, lin){

  df <- input |>
    dplyr::count(merged_clusterID) |> # count unique superclusters
    arrange(desc(n)) |> # sort high to low# rename to supercluster
    filter(merged_clusterID != "NA")

  # Generate a unique color for each supercluster using randomcoloR
  unique_clusters <- unique(df$merged_clusterID)
  n_clusters <- length(unique_clusters)
  palette <- randomcoloR::distinctColorPalette(n_clusters)

  #Create a color palette for each supercluster
  supercluster_palette <- tibble(
    merged_clusterID = unique_clusters,
    color = palette
  )

  # Reshape the original `df` to long format and join with
  ## `supercluster_palette`
  component_mapping_df <- input |>
    dplyr::select(!Tip_lable) |>
    tidyr::pivot_longer(cols = starts_with("t"),
                 names_to = "distance",
                 values_to = "component") |>
    left_join(supercluster_palette, by = "merged_clusterID")

  # Prioritize the color of the largest supercluster for each
  ## component-distance combination
  component_palette <- component_mapping_df |>
    dplyr::select(component, distance, component_color = color) |>
    dplyr::distinct() # Ensure unique rows

  # Create a combined palette dataframe for reference, if needed
  combined_palette <- dplyr::bind_rows(
    supercluster_palette |> 
      dplyr::select(name = merged_clusterID, palette_color = color),
    component_palette |> 
      dplyr::select(name = component, palette_color = component_color)
  )

  # Convert combined_palette to a named vector
  color_palette <- setNames(combined_palette$palette_color,
                            combined_palette$name)

  # Return the color palette
  return(color_palette)

}