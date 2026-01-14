import pandas as pd
import sys
import os
import glob
import numpy as np
from datetime import datetime
import warnings

# Suppress specific pandas warnings if desired (optional)
warnings.filterwarnings('ignore', category=FutureWarning)

def load_exon_lengths(species, exon_dir="exons"):
    """
    Loads the original exon annotation BED file for a species, calculates
    the length (end - start) for each exon, and returns a Series
    mapping exon name to its length.
    Assumes BED file has at least 6 columns: chr, start, end, name, score, strand.
    """
    exon_file = os.path.join(exon_dir, f"{species}_exons.bed")
    if not os.path.exists(exon_file):
        print(f"Error: Original exon file not found: {exon_file}")
        return None
    
    try:
        # Read only necessary columns: start (col 1), end (col 2), name (col 3) (0-based indexing)
        # Use low_memory=False for potentially large files
        df_exons = pd.read_csv(exon_file, sep='\t', header=None, 
                               usecols=[1, 2, 3], names=['start', 'end', 'name'],
                               low_memory=False)
        
        # Check for empty file
        if df_exons.empty:
            print(f"Warning: Original exon file is empty: {exon_file}")
            return pd.Series(dtype=int) # Return empty Series

        # Calculate length
        df_exons['length'] = df_exons['end'] - df_exons['start']
        
        # Handle potential non-positive lengths (optional, based on data quality)
        df_exons = df_exons[df_exons['length'] > 0]

        # Set name as index and keep only the length Series
        exon_lengths = df_exons.set_index('name')['length']
        
        # Check for duplicate exon names - keeping the first occurrence
        if exon_lengths.index.has_duplicates:
            print(f"Warning: Duplicate exon names found in {exon_file}. Keeping first occurrence.")
            exon_lengths = exon_lengths[~exon_lengths.index.duplicated(keep='first')]
            
        return exon_lengths
        
    except pd.errors.EmptyDataError:
        print(f"Warning: Original exon file is empty: {exon_file}")
        return pd.Series(dtype=int)
    except Exception as e:
        print(f"Error reading or processing original exon file {exon_file}: {e}")
        return None

def calculate_reciprocal_pairs(species1, species2, overlap_dir="overlaps", exon_dir="exons", threshold=0.6):
    """
    Identifies reciprocal best overlapping exons between two species
    satisfying the intersection/(original union) ratio threshold.
    Uses original exon lengths for calculation.
    """
    s1_to_s2_file = os.path.join(overlap_dir, f"{species1}_to_{species2}_overlaps.txt")
    s2_to_s1_file = os.path.join(overlap_dir, f"{species2}_to_{species1}_overlaps.txt")
    output_file = os.path.join("reciprocal_results", f"{species1}_{species2}_reciprocal_pairs.tsv") # Ensure output dir exists
    
    if not os.path.exists(s1_to_s2_file) or not os.path.exists(s2_to_s1_file):
        print(f"Warning: Overlap files missing for pair {species1}-{species2}. Skipping.")
        return None

    print(f"Processing reciprocal overlaps for {species1} and {species2}...")

    # --- Load Original Exon Lengths ---
    s1_lengths = load_exon_lengths(species1, exon_dir)
    s2_lengths = load_exon_lengths(species2, exon_dir)

    if s1_lengths is None or s2_lengths is None:
        print(f"Error loading original exon lengths for {species1} or {species2}. Skipping pair.")
        return None
    if s1_lengths.empty or s2_lengths.empty:
        print(f"Warning: Original exon lengths missing for {species1} or {species2}. Skipping pair.")
        return None

    # Define column indices from bedtools intersect -wo output (0-based)
    # We need: name_s1_orig (col 3), name_s2_orig (col 10), overlap_bp (col 14)
    cols_indices_needed = [3, 10, 14]
    cols_names_needed = ['name_s1_orig', 'name_s2_orig', 'overlap_bp']

    try:
        # --- Read A->B overlaps ---
        df_a_to_b = pd.read_csv(s1_to_s2_file, sep='\t', header=None, 
                                usecols=cols_indices_needed, names=cols_names_needed,
                                low_memory=False)
        
        if df_a_to_b.empty:
            print(f"Warning: Empty overlap file {s1_to_s2_file}.")
            # Fall through to handle empty reciprocal_df later
        else:
             # Convert overlap to int
            df_a_to_b['overlap_bp'] = pd.to_numeric(df_a_to_b['overlap_bp'], errors='coerce').fillna(0).astype(int)
            # Keep only positive overlaps
            df_a_to_b = df_a_to_b[df_a_to_b['overlap_bp'] > 0]
            # Merge with original lengths
            df_a_to_b = pd.merge(df_a_to_b, s1_lengths.rename('len_s1_orig'), left_on='name_s1_orig', right_index=True, how='inner')
            df_a_to_b = pd.merge(df_a_to_b, s2_lengths.rename('len_s2_orig'), left_on='name_s2_orig', right_index=True, how='inner')
            df_a_to_b = df_a_to_b[['name_s1_orig', 'name_s2_orig', 'len_s1_orig', 'len_s2_orig', 'overlap_bp']].drop_duplicates()
            # Ensure lengths are integers
            df_a_to_b['len_s1_orig'] = df_a_to_b['len_s1_orig'].astype(int)
            df_a_to_b['len_s2_orig'] = df_a_to_b['len_s2_orig'].astype(int)


        # --- Read B->A overlaps ---
        df_b_to_a = pd.read_csv(s2_to_s1_file, sep='\t', header=None, 
                                usecols=cols_indices_needed, names=cols_names_needed,
                                low_memory=False)
        
        if df_b_to_a.empty:
             print(f"Warning: Empty overlap file {s2_to_s1_file}.")
             # Create an empty DataFrame with expected columns for merge check if df_a_to_b wasn't empty
             df_b_to_a_check = pd.DataFrame(columns=['name_s1_orig', 'name_s2_orig'])
        else:
            # Convert overlap to int
            df_b_to_a['overlap_bp'] = pd.to_numeric(df_b_to_a['overlap_bp'], errors='coerce').fillna(0).astype(int)
            # Keep only positive overlaps
            df_b_to_a = df_b_to_a[df_b_to_a['overlap_bp'] > 0]
            
            # IMPORTANT: In B->A file, col 3 is S2 name, col 10 is S1 name
            # Rename immediately after reading for clarity
            df_b_to_a = df_b_to_a.rename(columns={'name_s1_orig': 'name_s2_orig', # Original S2 name
                                               'name_s2_orig': 'name_s1_orig'}) # Original S1 name
            
            # Keep only columns needed for the reciprocal check
            df_b_to_a_check = df_b_to_a[['name_s1_orig', 'name_s2_orig']].drop_duplicates()


        # --- Reciprocal Check ---
        # Handle case where one or both input dataframes might be empty after filtering/loading
        if df_a_to_b.empty or df_b_to_a_check.empty:
            print(f"No valid overlaps found in one or both directions for {species1}-{species2}.")
            reciprocal_df = pd.DataFrame(columns=['name_s1_orig', 'name_s2_orig', 'len_s1_orig', 'len_s2_orig', 'overlap_bp']) # Empty DF
        else:
            reciprocal_df = pd.merge(
                df_a_to_b, # Contains overlaps A->B with correct original lengths
                df_b_to_a_check, # Contains pairs found in B->A direction
                on=['name_s1_orig', 'name_s2_orig'],
                how='inner'
            )

        total_reciprocal_found = len(reciprocal_df) # Count before thresholding

        # --- Calculate Ratio using CORRECT lengths ---
        if not reciprocal_df.empty:
            # Calculate union length: len(A) + len(B) - overlap
            reciprocal_df['union_len'] = reciprocal_df['len_s1_orig'] + reciprocal_df['len_s2_orig'] - reciprocal_df['overlap_bp']

            # Define safe ratio calculation function
            def calculate_ratio(row):
                # Ensure union is positive and overlap is valid relative to original lengths
                if row['union_len'] > 0 and \
                   row['overlap_bp'] >= 0 and \
                   row['overlap_bp'] <= row['len_s1_orig'] and \
                   row['overlap_bp'] <= row['len_s2_orig']:
                    return row['overlap_bp'] / row['union_len']
                return 0.0 # Return 0 for invalid cases

            # Apply ratio calculation
            reciprocal_df['ratio'] = reciprocal_df.apply(calculate_ratio, axis=1)
            
            # Filter invalid ratios (e.g., due to potential data inconsistencies)
            reciprocal_df = reciprocal_df[reciprocal_df['ratio'] >= 0] 

            # --- Filter by Threshold ---
            final_pairs = reciprocal_df[reciprocal_df['ratio'] >= threshold].copy()
        else:
            # If reciprocal_df was empty initially
            final_pairs = reciprocal_df.copy() # Still empty, but has columns

        # --- Prepare Output and Stats ---
        if final_pairs.empty:
            print(f"No pairs passed the threshold ({threshold}) for {species1} to {species2}")
            passed_threshold_count = 0
            avg_ratio, max_ratio, min_ratio, median_ratio = 0, 0, 0, 0
            # Save empty file with header
            header_cols = ['species1', 'species2', 'name_s1_orig', 'name_s2_orig', 
                           'len_s1_orig', 'len_s2_orig', 'overlap_bp', 'union_len', 'ratio']
            pd.DataFrame(columns=header_cols).to_csv(output_file, sep='\t', index=False, header=True)

        else:
            passed_threshold_count = len(final_pairs)
            final_pairs['species1'] = species1
            final_pairs['species2'] = species2
            
            # Select and order final columns
            final_columns = ['species1', 'species2', 'name_s1_orig', 'name_s2_orig', 
                             'len_s1_orig', 'len_s2_orig', 'overlap_bp', 'union_len', 'ratio']
            final_pairs = final_pairs[final_columns]
            
            final_pairs.to_csv(output_file, sep='\t', index=False, header=True)
            print(f"Saved {len(final_pairs)} reciprocal pairs (ratio >= {threshold}) to {output_file}")

            # Calculate stats on the pairs *that passed* the threshold
            avg_ratio = final_pairs['ratio'].mean()
            max_ratio = final_pairs['ratio'].max()
            min_ratio = final_pairs['ratio'].min()
            median_ratio = final_pairs['ratio'].median()
            
        # Return statistics
        return {
            'species1': species1,
            'species2': species2,
            'total_reciprocal': total_reciprocal_found, # Total pairs found before threshold
            'passed_threshold': passed_threshold_count, # Pairs after threshold
            'avg_ratio': avg_ratio if passed_threshold_count > 0 else 0,
            'max_ratio': max_ratio if passed_threshold_count > 0 else 0,
            'min_ratio': min_ratio if passed_threshold_count > 0 else 0,
            'median_ratio': median_ratio if passed_threshold_count > 0 else 0
        }

    except pd.errors.EmptyDataError:
        print(f"Warning: One or both overlap files for {species1}-{species2} resulted in empty data after initial read. Skipping.")
        return None
    except KeyError as e:
         print(f"Error processing pair {species1}-{species2}: Missing expected column or index key: {e}. Check input file formats and column definitions.")
         return None
    except Exception as e:
        print(f"Error processing pair {species1}-{species2}: {e}")
        # Consider adding more specific error logging here if needed
        import traceback
        traceback.print_exc() # Print traceback for debugging
        return None

# (Keep the get_species_from_files and generate_detailed_report functions as they were)
def get_species_from_files(overlap_dir="overlaps"):
    """
    Extracts the list of species from existing overlap files
    """
    species_set = set()
    
    # Pattern: species1_to_species2_overlaps.txt
    overlap_files = glob.glob(os.path.join(overlap_dir, "*_to_*_overlaps.txt"))
    
    for file_path in overlap_files:
        filename = os.path.basename(file_path)
        parts = filename.split('_to_')
        if len(parts) == 2:
            species1 = parts[0]
            species2 = parts[1].replace('_overlaps.txt', '')
            species_set.add(species1)
            species_set.add(species2)
    
    return sorted(list(species_set))

def generate_detailed_report(stats, threshold, output_dir="reciprocal_results"):
    """
    Generate a detailed Markdown report.
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Convert stats to DataFrame
    if not stats:
        print("No data available for report generation.")
        return
        
    stats_df = pd.DataFrame(stats)
    
    # Save detailed stats to TSV
    stats_tsv_file = os.path.join(output_dir, "reciprocal_summary.tsv")
    stats_df.round(4).to_csv(stats_tsv_file, sep='\t', index=False)
    print(f"Summary statistics saved to {stats_tsv_file}")
    
    # Create markdown report
    report_file = os.path.join(output_dir, "reciprocal_analysis_report.md")
    
    with open(report_file, 'w') as f:
        # Title and introduction
        f.write(f"# Reciprocal Exon Overlap Analysis Report\n\n")
        f.write(f"*Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*\n\n")
        f.write(f"This report summarizes the results of reciprocal exon overlap analysis between pairs of species, ")
        f.write(f"using an **Intersection / Union ratio threshold of {threshold:.2f}**. ")
        f.write(f"The union length is calculated based on the **original exon lengths**.\n\n")
        
        # Summary statistics
        f.write("## Overall Summary\n\n")
        # Calculate overall stats carefully, handling potential NaNs if some pairs had no passing results
        valid_stats = stats_df[stats_df['passed_threshold'] > 0]
        num_species = len(set(stats_df['species1']).union(set(stats_df['species2'])))
        avg_reciprocal = stats_df['total_reciprocal'].mean()
        avg_passed = stats_df['passed_threshold'].mean()
        avg_ratio_overall = valid_stats['avg_ratio'].mean() if not valid_stats.empty else 0

        f.write(f"- Total species analyzed: {num_species}\n")
        f.write(f"- Total species pairs processed: {len(stats_df)}\n")
        f.write(f"- Average reciprocal pairs found (before threshold): {avg_reciprocal:.2f}\n")
        f.write(f"- Average pairs passing threshold (ratio >= {threshold:.2f}): {avg_passed:.2f}\n")
        f.write(f"- Average similarity ratio (for pairs passing threshold): {avg_ratio_overall:.4f}\n\n")
        
        # Detailed results table
        f.write(f"## Detailed Results (Threshold: {threshold:.2f})\n\n")
        f.write("| Species 1 | Species 2 | Total Reciprocal Found | Passing Threshold | Avg Ratio (Passed) | Min Ratio (Passed) | Max Ratio (Passed) | Median Ratio (Passed) |\n")
        f.write("|-----------|-----------|------------------------|-------------------|--------------------|--------------------|--------------------|-----------------------|\n")
        
        # Format output for better readability
        stats_df_sorted = stats_df.sort_values(by=['species1', 'species2'])
        for _, row in stats_df_sorted.iterrows():
             # Format ratios only if threshold passed > 0
            avg_r = f"{row['avg_ratio']:.4f}" if row['passed_threshold'] > 0 else "N/A"
            min_r = f"{row['min_ratio']:.4f}" if row['passed_threshold'] > 0 else "N/A"
            max_r = f"{row['max_ratio']:.4f}" if row['passed_threshold'] > 0 else "N/A"
            med_r = f"{row['median_ratio']:.4f}" if row['passed_threshold'] > 0 else "N/A"

            f.write(f"| {row['species1']} | {row['species2']} | {row['total_reciprocal']} | ")
            f.write(f"{row['passed_threshold']} | {avg_r} | {min_r} | ")
            f.write(f"{max_r} | {med_r} |\n")
        
        f.write("\n\n")
        
        # Recommendations and next steps
        f.write("## Recommendations for Further Analysis\n\n")
        f.write("1. **Phylogenetic Context**: Correlate the number/ratio of reciprocal pairs with the phylogenetic distance between species pairs.\n")
        f.write("2. **Functional Enrichment**: Analyze the GO terms or pathways associated with genes containing reciprocal exons, especially those with high ratios.\n")
        f.write("3. **Sequence Level Conservation**: Perform sequence alignments for the identified reciprocal exon pairs to assess nucleotide/amino acid identity.\n")
        f.write("4. **Exon Gain/Loss**: Investigate species pairs with very low numbers of reciprocal exons to potentially identify lineage-specific exon evolution.\n")
        f.write("5. **Parameter Sensitivity**: Evaluate how changing the ratio threshold affects the number and characteristics of identified pairs.\n\n")
        
        f.write("## Methods\n\n")
        f.write("Exon coordinates from each source species were projected onto the target species' genome using `halLiftover`. ")
        f.write("Overlaps between lifted exons and the target species' annotated exons were identified using `bedtools intersect -wo`. ")
        f.write("Reciprocal overlaps were found by requiring an exon pair (A, B) to appear in both the A-to-B and B-to-A overlap results. ")
        f.write(f"A similarity ratio was calculated for each reciprocal pair as `Overlap Length / (Original Length A + Original Length B - Overlap Length)`. ")
        f.write(f"Pairs with a ratio ≥ {threshold:.2f} were retained as significant reciprocal best overlapping exons.\n")
    
    print(f"Detailed report generated at {report_file}")


# --- Main Execution ---
def main():
    # Set threshold value (can be made configurable via command line)
    threshold = 0.6 # Example: Lower threshold
    
    # Define directories
    exon_dir = "exons"
    overlap_dir = "overlaps"
    output_dir = "reciprocal_results" # Directory for final pairs and reports
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Define species list (or get dynamically if preferred)
    # species_list = get_species_from_files(overlap_dir) # Option to get from files
    species_list = [
        "Acyrthosiphon_pisum","Aedes_aegypti", "Apis_mellifera", "Blattella_germanica", 
        "Bombyx_mori", "Drosophila_mojavensis", "Gryllus_bimaculatus", 
        "Helicoverpa_armigera", "Tribolium_castaneum"
    ]

    if not species_list:
         print("Error: No species found or defined. Exiting.")
         sys.exit(1)

    print(f"Processing {len(species_list)} species: {', '.join(species_list)}")
    print(f"Using overlap ratio threshold: {threshold}")
    
    # Store statistics for summary report
    stats = []
    
    # Process each unique species pair
    processed_pairs = set()
    for i in range(len(species_list)):
        for j in range(len(species_list)): # Check all pairs, including self if needed? No, skip self.
            if i == j:
                continue
                
            s1 = species_list[i]
            s2 = species_list[j]
            
            # Ensure pair is processed only once (e.g., A-B, not B-A again)
            pair = tuple(sorted((s1, s2)))
            if pair in processed_pairs:
                continue
            processed_pairs.add(pair)

            # Calculate reciprocal pairs for the canonical pair (e.g., alphabetically sorted)
            result = calculate_reciprocal_pairs(pair[0], pair[1], 
                                                overlap_dir=overlap_dir, 
                                                exon_dir=exon_dir, 
                                                threshold=threshold)
            if result:
                stats.append(result)
            else:
                 # Log that the pair was skipped or failed
                 print(f"Calculation failed or skipped for pair: {pair[0]}-{pair[1]}")

    
    # Generate detailed report if any stats were collected
    if stats:
        generate_detailed_report(stats, threshold, output_dir)
    else:
        print("No statistics collected, skipping report generation.")

    print("\nFinished processing all pairs.")

if __name__ == "__main__":
    main()
