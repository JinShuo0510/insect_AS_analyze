import pandas as pd
import sys
import os
import glob
import numpy as np
from datetime import datetime

def calculate_reciprocal_pairs(species1, species2, overlap_dir="overlaps", exon_dir="exons", threshold=0.6):
    """
    Identifies reciprocal best overlapping exons between two species
    satisfying the intersection/union ratio threshold.
    """
    s1_to_s2_file = os.path.join(overlap_dir, f"{species1}_to_{species2}_overlaps.txt")
    s2_to_s1_file = os.path.join(overlap_dir, f"{species2}_to_{species1}_overlaps.txt")
    output_file = os.path.join(overlap_dir, f"{species1}_{species2}_reciprocal_pairs.tsv")
    
    if not os.path.exists(s1_to_s2_file) or not os.path.exists(s2_to_s1_file):
        print(f"Warning: Overlap files missing for pair {species1}-{species2}. Skipping.")
        return None
    
    print(f"Processing reciprocal overlaps for {species1} and {species2}...")
    
    # Define column names based on BED6+1 input and bedtools intersect -wo output
    cols_intersect = [
        'chr_A_lifted', 'start_A_lifted', 'end_A_lifted', 'name_s1_orig', 'score_A', 'strand_A', 'len_s1_orig',
        'chr_B_target', 'start_B_target', 'end_B_target', 'name_s2_orig', 'score_B', 'strand_B', 'len_s2_orig',
        'overlap_bp'
    ]
    
    try:
        # Read A->B overlaps
        df_a_to_b = pd.read_csv(s1_to_s2_file, sep='\t', header=None, names=cols_intersect,
                               converters={'len_s1_orig': int, 'len_s2_orig': int, 'overlap_bp': int})
        
        # Check for empty file
        if df_a_to_b.empty:
            print(f"Warning: Empty overlap file {s1_to_s2_file}. Skipping.")
            return None
        
        df_a_to_b = df_a_to_b[['name_s1_orig', 'name_s2_orig', 'len_s1_orig', 'len_s2_orig', 'overlap_bp']].drop_duplicates()
        
        # Read B->A overlaps
        df_b_to_a = pd.read_csv(s2_to_s1_file, sep='\t', header=None, names=cols_intersect,
                               converters={'len_s1_orig': int, 'len_s2_orig': int, 'overlap_bp': int})
        
        # Check for empty file
        if df_b_to_a.empty:
            print(f"Warning: Empty overlap file {s2_to_s1_file}. Skipping.")
            return None
        
        df_b_to_a = df_b_to_a[['name_s1_orig', 'name_s2_orig', 'len_s1_orig', 'len_s2_orig', 'overlap_bp']].drop_duplicates()
        
        # Rename columns for consistency before merge
        df_b_to_a = df_b_to_a.rename(columns={'name_s1_orig': 'name_s2_orig', 'name_s2_orig': 'name_s1_orig',
                                           'len_s1_orig': 'len_s2_orig', 'len_s2_orig': 'len_s1_orig'})
        
        # --- Reciprocal Check ---
        reciprocal_df = pd.merge(
            df_a_to_b,
            df_b_to_a[['name_s1_orig', 'name_s2_orig']], 
            on=['name_s1_orig', 'name_s2_orig'],
            how='inner'
        )
        
        # --- Calculate Ratio ---
        reciprocal_df['union_len'] = reciprocal_df['len_s1_orig'] + reciprocal_df['len_s2_orig'] - reciprocal_df['overlap_bp']
        
        # Create safe ratio calculation function
        def calculate_ratio(row):
            if row['union_len'] > 0:
                return row['overlap_bp'] / row['union_len']
            return 0
            
        # Apply ratio calculation
        reciprocal_df['ratio'] = reciprocal_df.apply(calculate_ratio, axis=1)
        
        # --- Filter by Threshold ---
        final_pairs = reciprocal_df[reciprocal_df['ratio'] >= threshold].copy()
        
        # 关键修复：检查是否有满足条件的对
        if final_pairs.empty:
            print(f"No pairs passed the threshold for {species1} to {species2}")
            # 创建一个空的DataFrame，但预先定义好所有列
            columns = ['species1', 'species2', 'name_s1_orig', 'name_s2_orig', 
                      'len_s1_orig', 'len_s2_orig', 'overlap_bp', 'union_len', 'ratio']
            final_pairs = pd.DataFrame(columns=columns)
            # 写入空文件
            final_pairs.to_csv(output_file, sep='\t', index=False, header=True)
            print(f"Saved 0 reciprocal pairs to {output_file}")
            
            # 返回统计信息
            return {
                'species1': species1,
                'species2': species2,
                'total_reciprocal': len(reciprocal_df),
                'passed_threshold': 0,
                'avg_ratio': 0,
                'max_ratio': 0,
                'min_ratio': 0,
                'median_ratio': 0
            }
        else:
            # 有满足条件的对，继续处理
            final_pairs['species1'] = species1
            final_pairs['species2'] = species2
            
            # 选择和保存最终结果
            final_columns = ['species1', 'species2', 'name_s1_orig', 'name_s2_orig', 
                           'len_s1_orig', 'len_s2_orig', 'overlap_bp', 'union_len', 'ratio']
            final_pairs = final_pairs[final_columns]
            
            final_pairs.to_csv(output_file, sep='\t', index=False, header=True)
            print(f"Saved {len(final_pairs)} reciprocal pairs to {output_file}")
            
            # 返回统计信息
            return {
                'species1': species1,
                'species2': species2,
                'total_reciprocal': len(reciprocal_df),
                'passed_threshold': len(final_pairs),
                'avg_ratio': final_pairs['ratio'].mean(),
                'max_ratio': final_pairs['ratio'].max(),
                'min_ratio': final_pairs['ratio'].min(),
                'median_ratio': final_pairs['ratio'].median()
            }
        
    except pd.errors.EmptyDataError:
        print(f"Warning: One or both overlap files for {species1}-{species2} are empty. Skipping.")
        return None
    except Exception as e:
        print(f"Error processing pair {species1}-{species2}: {e}")
        return None

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
    Generate a detailed HTML report with interactive tables and plots
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Convert stats to DataFrame
    if not stats:
        print("No data available for report generation.")
        return
        
    stats_df = pd.DataFrame(stats)
    
    # Save detailed stats to CSV
    stats_df.to_csv(f"{output_dir}/reciprocal_summary.tsv", sep='\t', index=False)
    
    # Create markdown report
    report_file = f"{output_dir}/reciprocal_analysis_report.md"
    
    with open(report_file, 'w') as f:
        # Title and introduction
        f.write(f"# Reciprocal Exon Overlap Analysis Report\n\n")
        f.write(f"*Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*\n\n")
        f.write(f"This report summarizes the results of reciprocal exon overlap analysis between pairs of species, ")
        f.write(f"using a similarity threshold of {threshold}.\n\n")
        
        # Summary statistics
        f.write("## Overall Summary\n\n")
        f.write(f"- Total species analyzed: {len(set(stats_df['species1']).union(set(stats_df['species2'])))}\n")
        f.write(f"- Total species pairs: {len(stats_df)}\n")
        f.write(f"- Average reciprocal pairs found: {stats_df['total_reciprocal'].mean():.2f}\n")
        f.write(f"- Average pairs passing threshold: {stats_df['passed_threshold'].mean():.2f}\n")
        f.write(f"- Average similarity ratio: {stats_df['avg_ratio'].mean():.4f}\n\n")
        
        # Detailed results table
        f.write("## Detailed Results\n\n")
        f.write("| Species 1 | Species 2 | Total Reciprocal | Passing Threshold | Avg Ratio | Min Ratio | Max Ratio | Median Ratio |\n")
        f.write("|-----------|-----------|------------------|-------------------|-----------|-----------|-----------|-------------|\n")
        
        for _, row in stats_df.iterrows():
            f.write(f"| {row['species1']} | {row['species2']} | {row['total_reciprocal']} | ")
            f.write(f"{row['passed_threshold']} | {row['avg_ratio']:.4f} | {row['min_ratio']:.4f} | ")
            f.write(f"{row['max_ratio']:.4f} | {row['median_ratio']:.4f} |\n")
        
        f.write("\n\n")
        
        # Recommendations and next steps
        f.write("## Recommendations for Further Analysis\n\n")
        f.write("1. **Phylogenetic Analysis**: Combine these results with phylogenetic distances to explore evolutionary patterns.\n")
        f.write("2. **Functional Enrichment**: Analyze the functional categories of genes containing highly conserved exons.\n")
        f.write("3. **Sequence Conservation**: Examine sequence-level conservation within these reciprocal exon pairs.\n")
        f.write("4. **Species-Specific Features**: Identify exons unique to specific species or taxonomic groups.\n\n")
        
        f.write("## Methods\n\n")
        f.write("Reciprocal best overlaps were identified by comparing exon coordinates projected between species pairs in both directions. ")
        f.write(f"A similarity ratio was calculated as (overlap length) / (union length) for each exon pair, ")
        f.write(f"and pairs with a ratio ≥ {threshold} were considered significant matches.\n")
    
    print(f"Detailed report generated at {report_file}")

# --- Main Execution ---
def main():
    # Set threshold value (can be made configurable via command line)
    threshold = 0.6
    
    # Create output directory
    output_dir = "reciprocal_results"
    os.makedirs(output_dir, exist_ok=True)
    
    # Define species list
    species_list = [
        "Aedes_aegypti", "Apis_mellifera", "Blattella_germanica", 
        "Bombyx_mori", "Drosophila_mojavensis", "Gryllus_bimaculatus", 
        "Helicoverpa_armigera", "Tribolium_castaneum"
    ]
    
    print(f"Found {len(species_list)} species: {', '.join(species_list)}")
    
    # Store statistics for summary report
    stats = []
    
    # Process each species pair
    for i in range(len(species_list)):
        for j in range(i + 1, len(species_list)):
            s1 = species_list[i]
            s2 = species_list[j]
            result = calculate_reciprocal_pairs(s1, s2, overlap_dir="overlaps", threshold=threshold)
            if result:
                stats.append(result)
    
    # Generate detailed report
    generate_detailed_report(stats, threshold, output_dir)
    
    print("Finished processing all pairs.")

if __name__ == "__main__":
    main()
