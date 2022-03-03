# Takes the filtered data, removes everything but the markers, keeping chip/side

library(tidyverse)

# Read in file with cells as rows and info columns for batches
data = read_csv("cleaned_filtered_data.csv")

# Filters out the QCHeLa and NA values in the Chip column
first_filter = filter(data, Chip != "_QCHeLa_")
second_filter = filter(first_filter, Chip != "NA")

# Transposes it to proteins as rows, then sets rownames to columns
all_t = t(as.matrix(second_filter))
all_t = data.frame(all_t)
all_t = rownames_to_column(all_t, var = "Accession")

# Saves the Chip/Side and other info columns and transposes it
info_rows = all_t[618:622,]
rownames(info_rows) = NULL
info_rows = column_to_rownames(info_rows, var = "Accession")
info_t = data.frame(t(as.matrix(info_rows)))

# Saves known markers for the phases
markers = read_csv("phase_markers.csv")

# Only saves the marker proteins, then makes the Accession into rownames
joined = inner_join(all_t, markers)
ordered = arrange(joined, Phase)
ordered = ordered[,1:159]
ordered = column_to_rownames(ordered, var = "Accession")

# Transposes the data back to cells as rows
final = data.frame(t(as.matrix(ordered)))

# Makes the values numeric 
final = lapply(final, as.numeric)

# Adds back on the info columns 
final_final = data.frame(final, info_t)

# Writes it to a *.csv
write_csv(final_final, "Markers_only_chip_and_side.csv")
