import networkx as nx
import pickle
from collections import Counter

def find_one_to_one_orthologs(graph_file="ortholog_graph.pkl", output_file="one_to_one_orthologs.txt", num_species=7):
    """Identifies 1:1 ortholog groups from the homology graph."""
    print(f"Loading graph from {graph_file}...")
    try:
        with open(graph_file, 'rb') as f_in:
            G = pickle.load(f_in)
    except FileNotFoundError:
        print(f"Error: Graph file {graph_file} not found.")
        return

    print("Finding connected components and filtering for 1:1 orthologs...")
    one_to_one_groups = []
    components = list(nx.connected_components(G))
    print(f"Found {len(components)} connected components.")

    for i, component in enumerate(components):
        if len(component) == num_species:
            # Check if all nodes (exons) are from distinct species
            species_in_component = [node.split("::")[0] for node in component] # Assumes node name format "Species::Gene::Exon"
            species_counts = Counter(species_in_component)

            # Check if there are exactly num_species distinct species, each appearing once
            if len(species_counts) == num_species and all(count == 1 for count in species_counts.values()):
                one_to_one_groups.append(sorted(list(component))) # Sort for consistent output

        if (i + 1) % 1000 == 0:
             print(f"Processed {i+1}/{len(components)} components...")


    print(f"Identified {len(one_to_one_groups)} potential 1:1 ortholog groups.")

    # Save the groups to a file, one group per line, exons separated by tabs
    with open(output_file, 'w') as f_out:
        for group in one_to_one_groups:
            f_out.write("\t".join(group) + "\n")

    print(f"1:1 ortholog groups saved to {output_file}")

# --- Main Execution ---
find_one_to_one_orthologs(num_species=8) # Set num_species to 7 for your case
