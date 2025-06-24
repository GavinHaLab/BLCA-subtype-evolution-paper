
##########
# Author: Allie S. Kreitman
# Usage: Rscript clone_tree_generation_for_machina_BLCA_TAN.R 
# Date: 2024-07-03, 2024-07-05
# Goal: This script is to format the lychee data from BCLA Rapid autopsy project to generate the clone tree input for running tool MACHINA
##########

########## DOCUMENT SETUP ----------------------------
# load packages
library(tidyverse)
library(stringr)



# USE THIS FOR SINGLE FILES
# Read the .dot file into R as a text file
filepath <- "../../data/19-004_patient_sample_lichee_tree.dot"


#for (filepath in dot_files_list){
message("Processing file: ", filepath) 
dot_file <- readLines(filepath, warn = FALSE)
patient_id <- str_extract(filepath, "\\d+-\\d+") %>% str_replace("-", "_")


########## READ DOT FILE INTO R AND CONVERT TO TABLE OF 'child' AND 'parent' FOR EACH CONNECTION
# Initialize an empty list to store edges
edges <- list()

# Start processing from line 4
for (i in 4:length(dot_file)) {
  line <- dot_file[i]
  
  # Extract edges using regex
  if (grepl("->", line)) {
    # Split the line by "->"
    nodes <- str_split(line, "->")[[1]]
    
    # Trim leading and trailing spaces from nodes
    node1 <- trimws(nodes[1])
    node2 <- trimws(gsub("\\[.*?\\];?$", "", nodes[2]))  # Remove any annotations in square brackets & trailing semicolon
    
    # Add to the edges list
    edges <- c(edges, list(c(node1, node2)))
  }
}

# Convert the list of edges into a data frame
edges_df <- do.call(rbind, edges)
edges_df <- data.frame(edges_df)
colnames(edges_df) <- c("parent", "child")

########## READ DOT FILE AND PULL THE SAMPLE ID AND ITS LEAF NUMBER 

# Step 1: Extract relevant lines containing numbers and labels
extracted_lines <- grep("\\[image=", readLines(textConnection(dot_file)), value = TRUE)

# Step 2: Extract number and label from each line
extracted_data <- lapply(extracted_lines, function(line) {
  # Extract number
  number <- as.numeric(sub("^([0-9]+)\\s.*", "\\1", line))
  
  # Extract label
  label <- sub('.*label="([^"]+)".*', '\\1', line)
  
  list(number = number, label = label)
})

# Step 3: Convert list to data frame
leaf_labels_fromDot <- do.call(rbind.data.frame, extracted_data)


########## BUILD OUT ALL CONNECTIONS IN THE TREE ---------------
# remove 0 node from edges
edges_df <- edges_df %>% filter(parent != 0)

# Identify the ends (leafs) by using values in the "child" col that are not in the "parent" col 
leafs_list <- setdiff(edges_df$child, edges_df$parent)
# leafs_list <- rev(leafs_list)
# Identify the top node by by using values in the "parent" col that are not in the "child" col
top_node <- setdiff(edges_df$parent, edges_df$child)

# list of unique nodes/ numbers in the sample 
nodes_list <- unique(c(edges_df$child, edges_df$parent))

# make a df to input the growing connections data
# headers will be "leafs" and all the nodes that are not leafs
headers <- c("leaf", paste("node_", sep = "", setdiff(nodes_list, leafs_list)))
# make a dataframe
tree_connections_df <- setNames(data.frame(matrix(ncol = length(headers), nrow = 0)), headers)
tree_connections_df$leaf <- as.numeric(tree_connections_df$leaf)
# tree_connections_df <- setNames(data.frame(matrix(ncol = length(headers), nrow = length(leafs_list))), headers)
# tree_connections_df$leafs <- leafs_list

for (i in 1:length(leafs_list)){ # loop through all the leafs in the tree
  # reset the tree location to the middle
  loop_location = "middle_of_tree"
  
  # make a temp leaf value to cycle through
  current_leaf <- as.numeric(leafs_list[i])
  
  # make a temp df to fill in for that leaf
  tree_connections_df_temp <- setNames(data.frame(matrix(ncol = length(headers), nrow = 0)), headers) %>% 
    add_row(leaf = current_leaf) 
  
  while(loop_location != "origin"){ # make a loop within each leaf for all the nodes until we reach the top of the tree
    
    # For the leafs, find the "parent" values that led there
    next_connected_nodes <- edges_df %>% 
      filter(child %in% current_leaf) %>% 
      pull(parent) %>% 
      unique()
    
    #format the nodes as they are written in the df headers
    next_connected_nodes_headers <- paste("node_", sep = "", next_connected_nodes)
    
    # input the parent values into a df of connections
    #tree_connections_df_temp$leaf <- as.numeric(leafs_list[i])
    tree_connections_df_temp[,next_connected_nodes_headers] <- TRUE
    
    # save the "parent" values as new temporary leafs
    current_leaf <- next_connected_nodes
    
    #check whether we are at the top of the tree by seeing if the only new leafs are the top node
    loop_location <- ifelse(top_node %in% current_leaf & length(current_leaf) == 1, "origin", loop_location)
    
    # repeat the above process until you get to the top_node  
    
    rm(next_connected_nodes_headers)
  }
  
  # add the dataframe for the last leaf into the larger dataframe
  tree_connections_df <- rbind(tree_connections_df, tree_connections_df_temp)
  
  # clean up from the loop
  rm(next_connected_nodes, loop_location, current_leaf, tree_connections_df_temp)
  
  # repeat the process for all the samples (leafs) 
}


# FORMAT THE OUTPUT FOR MACHINA ----------------
machina_format <- tree_connections_df %>% 
  pivot_longer( # pivot the table so it is just leafs and nodes and logical node output
    cols = starts_with("node"),
    names_to = "node",
    values_to = "node_logical"
  ) %>% 
  filter(!is.na(node_logical)) # remove entries where the logical is false (so that node does not connect to the leaf)

# Create a helper table to map unique values to letters
letter_map <- data.frame(
  node = unique(machina_format$node),
  node_letter = letters[seq_along(unique(machina_format$node))]
) %>% 
  mutate(node_number = gsub("node_", "", node))


# Join the helper table to df and select relevant columns
machina_format <- machina_format %>%
  left_join(letter_map, by = "node") %>%
  group_by(node_letter) %>% 
  mutate(numeric_leaf_value = row_number()) %>% 
  ungroup() %>% 
  #mutate(numeric_leaf_value = dense_rank(leaf)) %>%  #add leaf number from 1:lenght(leafs) 
  mutate(node_leaf_combo = paste(node_letter, sep = "", numeric_leaf_value))

# edges to letters
edges_letter <- left_join(edges_df, letter_map, by = join_by(parent == node_number)) %>% # join parent number with letter key
  rename(parent_letter = node_letter) %>%  #rename key col to parent_letter
  filter(!child %in% leafs_list) %>% #remove the leafs because they have been separately delt with above
  select(-node) %>% #remove node col
  left_join(letter_map, by = join_by(child == node_number)) %>% # join child numbers with letter key
  rename(child_letter = node_letter) %>%  #rename node letter with child letter to distinguish from parent letter
  select(parent_letter, child_letter)
  
# make a machina input file
machina_input_tree <- machina_format %>% 
  select(node_letter, node_leaf_combo) %>% 
  rename(parent_letter = node_letter) %>% 
  rename(child_letter = node_leaf_combo) %>% 
  rbind(edges_letter %>% 
          select(parent_letter, child_letter))

# machina leaf labels
machina_labels <- machina_format %>% 
  left_join(leaf_labels_fromDot, by = join_by(leaf == number)) %>% 
  select(node_leaf_combo, label)

# color labeling file
color_file <- leaf_labels_fromDot %>% 
  select(label) %>% 
  mutate(color = row_number())
  
# save machina outputs
write_tsv(machina_labels, paste0("../../results/patient_sample_lichee_tree", patient_id,"-", "labeling.tsv"), col_names = FALSE)
write_tsv(machina_input_tree, paste0("../../results/patient_sample_lichee_tree", patient_id,"-", "tree.tsv"), col_names = FALSE)
write_tsv(color_file, paste0("../../results/patient_sample_lichee_tree", patient_id,"-", "colorfile.tsv"), col_names = FALSE)

#}


