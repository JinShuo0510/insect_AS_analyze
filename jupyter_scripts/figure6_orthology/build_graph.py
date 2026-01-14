import pandas as pd
import networkx as nx
import os
import glob
import pickle # To save the graph object

def build_homology_graph(species_list, pairs_dir="reciprocal_results", graph_output_file="ortholog_graph.pkl"):
    """Builds a graph where nodes are exons and edges represent reciprocal homology."""
    G = nx.Graph()
    print("Building homology graph...")

    # Find all reciprocal pairs files
    pair_files = glob.glob(os.path.join(pairs_dir, "*_*_reciprocal_pairs.tsv"))

    if not pair_files:
        print("Error: No reciprocal pair files found in", pairs_dir)
        return None

    processed_pairs = 0
    for f in pair_files:
        try:
            df = pd.read_csv(f, sep='\t')
            for index, row in df.iterrows():
                exon1 = row['name_s1_orig']
                exon2 = row['name_s2_orig']
                # Add nodes implicitly by adding edge
                G.add_edge(exon1, exon2)
                processed_pairs += 1
        except pd.errors.EmptyDataError:
            print(f"Skipping empty file: {f}")
        except Exception as e:
            print(f"Error reading file {f}: {e}")

    print(f"Graph built with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges from {processed_pairs} pairs.")

    # Save the graph object for later use
    with open(graph_output_file, 'wb') as f_out:
        pickle.dump(G, f_out)
    print(f"Graph saved to {graph_output_file}")

    return G

# --- Main Execution ---
#species_list = ['speciesA', 'speciesB', 'speciesC', 'speciesD', 'speciesE', 'speciesF', 'speciesG']
species_list = [
        "Acyrthosiphon_pisum","Aedes_aegypti", "Apis_mellifera", "Blattella_germanica",
        "Bombyx_mori", "Drosophila_mojavensis", "Gryllus_bimaculatus",
        "Helicoverpa_armigera", "Tribolium_castaneum"
    ]
build_homology_graph(species_list, pairs_dir="reciprocal_results")
