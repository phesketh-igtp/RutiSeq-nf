dictionary_rename <- function(df,dict_path) { 

# import dictionary
dict <- read.csv(dict_path)

# create vector name
dict_names <- dict |> 
    select(new.name,old.name) |>
    deframe()

# Create vector for cols to keep (using 'new.name' since the cols will be renamed)
cols.to.keep <- dict |> 
    filter(final == "Y") |>
    select(new.name) |>
    deframe()

# Rename cols
df <- df |> 
    rename(all_of(dict_names)) |> # Rename cols using dict_names
    select(all_of(cols.to.keep)) # Keep only cols.to.keep
    
} 
