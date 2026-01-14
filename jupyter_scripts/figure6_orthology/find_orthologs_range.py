import networkx as nx
import pickle
from collections import Counter

def find_one_to_one_orthologs(graph_file="ortholog_graph.pkl", 
                              output_file="one_to_one_orthologs.txt", 
                              num_species=7,
                              min_species=5):
    """
    Identifies 1:1 ortholog groups from the homology graph.
    
    Parameters:
    - graph_file: Path to the input graph pickle file
    - output_file: Path to save the results
    - num_species: Total number of species in the dataset
    - min_species: Minimum number of species required in a 1:1 ortholog group
    """
    print(f"Loading graph from {graph_file}...")
    try:
        with open(graph_file, 'rb') as f_in:
            G = pickle.load(f_in)
    except FileNotFoundError:
        print(f"Error: Graph file {graph_file} not found.")
        return
    
    print(f"Finding connected components and filtering for 1:1 orthologs with at least {min_species} species...")
    one_to_one_groups = []
    components = list(nx.connected_components(G))
    print(f"Found {len(components)} connected components.")
    
    # Track statistics for reporting
    groups_by_species_count = {i: 0 for i in range(min_species, num_species+1)}
    
    for i, component in enumerate(components):
        # Count species in this component
        species_in_component = [node.split("::")[0] for node in component]  # Assumes node name format "Species::Gene::Exon"
        species_counts = Counter(species_in_component)
        
        # Check if there are at least min_species distinct species, each appearing exactly once
        if (len(species_counts) >= min_species and 
            len(species_counts) <= num_species and 
            all(count == 1 for count in species_counts.values())):
            
            # Add this component to our ortholog groups
            one_to_one_groups.append(sorted(list(component)))  # Sort for consistent output
            
            # Track statistics
            groups_by_species_count[len(species_counts)] += 1
        
        if (i + 1) % 1000 == 0:
            print(f"Processed {i+1}/{len(components)} components...")
    
    # Report statistics
    print(f"Identified {len(one_to_one_groups)} potential 1:1 ortholog groups:")
    for species_count, group_count in groups_by_species_count.items():
        if group_count > 0:
            print(f"  - {group_count} groups with {species_count} species")
    
    # Save the groups to a file, one group per line, exons separated by tabs
    with open(output_file, 'w') as f_out:
        # First write a header line with the count of species in each group
        f_out.write("species_count\tortholog_group\n")
        
        # Then write each group with its species count
        for group in one_to_one_groups:
            species_count = len(set(node.split("::")[0] for node in group))
            # 修复这一行，避免在f-string中使用反斜杠
            group_str = ",".join(group)
            f_out.write(f"{species_count}\t{group_str}\n")
    
    print(f"1:1 ortholog groups saved to {output_file}")
    return one_to_one_groups

# --- Main Execution ---
find_one_to_one_orthologs(num_species=9, min_species=1)  # Set num_species to total species count, min_species to 5
