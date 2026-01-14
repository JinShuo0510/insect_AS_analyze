import pandas as pd
import math
import sys
from io import StringIO

# --- Configuration ---
PWM_FILE = 'PWM.txt' # Path to your PWM file
RESULTS_FILE = 'all_tissues_enriched_in_downstream.csv' # Path to your results CSV file
OUTPUT_FILE = 'all_tissues_enriched_in_downstream_annotated_results.csv' # Name for the output file
RELATIVE_SCORE_THRESHOLD = 0.5 # Threshold from the paper (p >= 0.5 * P_max)
PSEUDOCOUNT = 1e-9 # Small value to avoid log(0) or division by zero

# --- Helper Functions ---

def parse_pwm_file(filepath):
    """Parses the custom PWM file format."""
    motifs = {}
    current_motif_data = None
    pwm_lines = []
    metadata = {}

    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: # Skip empty lines
                    continue

                if line.startswith('RBP'):
                    # Start of a new RBP block, save the previous one if exists
                    if current_motif_data and pwm_lines:
                        try:
                            # Convert PWM lines to matrix (list of dicts)
                            header = pwm_lines[0].split()
                            if header[0].upper() != 'POS':
                                print(f"Warning: Unexpected PWM header format for {metadata.get('Motif', 'Unknown')}: {pwm_lines[0]}. Skipping motif.", file=sys.stderr)
                            else:
                                matrix = []
                                for pwm_line in pwm_lines[1:]:
                                    parts = pwm_line.split()
                                    pos_data = {}
                                    # Handle potential extra columns if any, focus on A, C, G, U
                                    base_index = 1 # Start index for ACGU values
                                    base_map = { 'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3 } # Map U/T to same index
                                    
                                    # Ensure we have enough parts
                                    if len(parts) < 5:
                                         print(f"Warning: Malformed PWM line for {metadata.get('Motif', 'Unknown')}: {pwm_line}. Skipping line.", file=sys.stderr)
                                         continue

                                    vals = {}
                                    try:
                                        vals['A'] = float(parts[base_index + base_map['A']])
                                        vals['C'] = float(parts[base_index + base_map['C']])
                                        vals['G'] = float(parts[base_index + base_map['G']])
                                        # Use U value, treat T the same as U
                                        vals['U'] = float(parts[base_index + base_map['U']])
                                        vals['T'] = vals['U'] # Allow matching T in hexamer
                                    except (ValueError, IndexError) as e:
                                        print(f"Error parsing PWM values for {metadata.get('Motif', 'Unknown')}: {pwm_line} - {e}. Skipping line.", file=sys.stderr)
                                        continue

                                    # Normalize if frequencies don't sum to 1 (add pseudocount first)
                                    total = sum(vals[b] for b in ['A', 'C', 'G', 'U']) + 4 * PSEUDOCOUNT
                                    matrix.append({b: (vals[b] + PSEUDOCOUNT) / total for b in ['A', 'C', 'G', 'T', 'U']}) # Store T explicitly

                                if matrix: # Only add if PWM parsing was successful
                                    current_motif_data['pwm'] = matrix
                                    motifs[metadata['Motif']] = current_motif_data
                                else:
                                     print(f"Warning: Could not parse PWM matrix for {metadata.get('Motif', 'Unknown')}. Skipping motif.", file=sys.stderr)


                        except Exception as e:
                            print(f"Error processing motif {metadata.get('Motif', 'Unknown')}: {e}", file=sys.stderr)

                    # Reset for the new motif
                    current_motif_data = {}
                    metadata = {}
                    pwm_lines = []
                    parts = line.split(maxsplit=1) # Split only on the first space
                    if len(parts) == 2:
                         metadata['RBP_ID'] = parts[1]
                    else:
                         metadata['RBP_ID'] = "Unknown_RBP_ID"
                         print(f"Warning: Could not parse RBP ID line: {line}", file=sys.stderr)
                    current_motif_data['rbp_id'] = metadata['RBP_ID']


                elif line.startswith('RBP Name'):
                     metadata['RBP_Name'] = line.split(maxsplit=2)[-1] if len(line.split()) > 2 else 'Unknown'
                     current_motif_data['rbp_name'] = metadata['RBP_Name']
                elif line.startswith('Gene'):
                     metadata['Gene'] = line.split(maxsplit=1)[-1] if len(line.split()) > 1 else 'Unknown'
                     current_motif_data['gene'] = metadata['Gene']
                elif line.startswith('Motif'):
                     metadata['Motif'] = line.split(maxsplit=1)[-1] if len(line.split()) > 1 else 'Unknown_Motif_ID'
                     current_motif_data['motif_id'] = metadata['Motif'] # Use the specific motif ID as key later
                elif line.startswith('Family'):
                     current_motif_data['family'] = line.split(maxsplit=1)[-1] if len(line.split()) > 1 else 'Unknown'
                elif line.startswith('Species'):
                     current_motif_data['species'] = line.split(maxsplit=1)[-1] if len(line.split()) > 1 else 'Unknown'
                elif line.startswith('Pos') or line[0].isdigit():
                    pwm_lines.append(line)
                # else: Ignore other lines for now

            # Save the last motif after the loop ends
            if current_motif_data and pwm_lines and 'motif_id' in current_motif_data:
                try:
                    header = pwm_lines[0].split()
                    if header[0].upper() == 'POS':
                        matrix = []
                        for pwm_line in pwm_lines[1:]:
                           parts = pwm_line.split()
                           if len(parts) < 5: continue # Skip malformed
                           vals = {}
                           base_index = 1
                           base_map = { 'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3 }
                           try:
                                vals['A'] = float(parts[base_index + base_map['A']])
                                vals['C'] = float(parts[base_index + base_map['C']])
                                vals['G'] = float(parts[base_index + base_map['G']])
                                vals['U'] = float(parts[base_index + base_map['U']])
                                vals['T'] = vals['U']
                           except (ValueError, IndexError) as e:
                                print(f"Error parsing PWM values for {metadata.get('Motif', 'Unknown')}: {pwm_line} - {e}. Skipping line.", file=sys.stderr)
                                continue

                           total = sum(vals[b] for b in ['A', 'C', 'G', 'U']) + 4 * PSEUDOCOUNT
                           matrix.append({b: (vals[b] + PSEUDOCOUNT) / total for b in ['A', 'C', 'G', 'T', 'U']})

                        if matrix:
                            current_motif_data['pwm'] = matrix
                            motifs[current_motif_data['motif_id']] = current_motif_data
                        else:
                            print(f"Warning: Could not parse PWM matrix for {metadata.get('Motif', 'Unknown')}. Skipping motif.", file=sys.stderr)
                except Exception as e:
                    print(f"Error processing last motif {metadata.get('Motif', 'Unknown')}: {e}", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: PWM file not found at {filepath}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while reading PWM file: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Parsed {len(motifs)} motifs from {filepath}")
    return motifs


def calculate_max_pwm_prob(pwm_matrix, length=6):
    """Calculates the maximum probability sequence score for a PWM window."""
    if not pwm_matrix or len(pwm_matrix) < length:
        return 0.0 # Cannot calculate for invalid or too short PWM

    max_prob = 0.0
    # Slide a window of 'length' over the PWM
    for start_pos in range(len(pwm_matrix) - length + 1):
        current_window_max_prob = 1.0
        for i in range(length):
            pos_matrix = pwm_matrix[start_pos + i]
            # Find max probability at this position (handle T/U)
            max_p = max(pos_matrix.get('A', 0), pos_matrix.get('C', 0),
                        pos_matrix.get('G', 0), pos_matrix.get('U', 0)) # PWM uses U
            current_window_max_prob *= max(max_p, PSEUDOCOUNT) # Avoid multiplying by zero

        max_prob = max(max_prob, current_window_max_prob)

    return max_prob

def calculate_pwm_hexamer_prob(pwm_matrix, hexamer):
    """Calculates the probability of a PWM window generating a hexamer."""
    hexamer = hexamer.upper().replace('T', 'U') # Standardize hexamer to use U
    length = len(hexamer)
    if not pwm_matrix or len(pwm_matrix) < length or length != 6:
        return 0.0

    max_prob = 0.0
    # Slide a window of 'length' over the PWM
    for start_pos in range(len(pwm_matrix) - length + 1):
        current_window_prob = 1.0
        valid_hex = True
        for i in range(length):
            nucleotide = hexamer[i]
            pos_matrix = pwm_matrix[start_pos + i]
            if nucleotide in pos_matrix:
                current_window_prob *= max(pos_matrix[nucleotide], PSEUDOCOUNT) # Use value from PWM
            else:
                 #print(f"Warning: Invalid nucleotide '{nucleotide}' in hexamer '{hexamer}'.", file=sys.stderr)
                 # If hexamer has non-ACGTU, probability is 0 for this window
                 current_window_prob = 0.0
                 valid_hex = False
                 break # Exit inner loop
        if valid_hex:
            max_prob = max(max_prob, current_window_prob)

    return max_prob


# --- Main Execution ---

# 1. Parse PWMs
motifs_data = parse_pwm_file(PWM_FILE)
if not motifs_data:
    print("Error: No motifs were parsed successfully. Exiting.", file=sys.stderr)
    sys.exit(1)

# 2. Pre-calculate Max Probabilities for each motif
print("Pre-calculating maximum motif scores...")
motif_max_probs = {}
for motif_id, data in motifs_data.items():
    if 'pwm' in data:
        motif_max_probs[motif_id] = calculate_max_pwm_prob(data['pwm'], length=6)
    else:
        motif_max_probs[motif_id] = 0.0
        print(f"Warning: Motif {motif_id} missing PWM matrix for max prob calculation.", file=sys.stderr)
print("Done pre-calculating.")


# 3. Read Results File
try:
    results_df = pd.read_csv(RESULTS_FILE)
    print(f"Read {len(results_df)} rows from {RESULTS_FILE}")
    if 'hexamer' not in results_df.columns:
         raise ValueError("Results file must contain a 'hexamer' column.")
except FileNotFoundError:
    print(f"Error: Results file not found at {RESULTS_FILE}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"Error reading results file {RESULTS_FILE}: {e}", file=sys.stderr)
    sys.exit(1)


# 4. Annotate Hexamers
print("Annotating hexamers...")
annotations = []
total_rows = len(results_df)

for index, row in results_df.iterrows():
    hexamer = str(row['hexamer']).strip() # Ensure it's a string and remove whitespace

    # Progress indicator
    if (index + 1) % 100 == 0 or index == total_rows - 1:
        print(f"Processing row {index + 1}/{total_rows}...", end='\r')

    best_match = {
        'rbp_name': None,
        'motif_id': None,
        'match_score': 0.0,
        'max_score': 0.0,
        'relative_score': 0.0,
        'is_match': False
    }

    if len(hexamer) != 6 or not all(c in 'ACGTU' for c in hexamer.upper()):
         # print(f"Warning: Skipping invalid hexamer '{hexamer}' at index {index}.", file=sys.stderr)
         annotations.append(best_match)
         continue


    current_best_relative_score = -1.0

    for motif_id, data in motifs_data.items():
        if 'pwm' not in data:
            continue # Skip if PWM is missing

        pwm = data['pwm']
        max_score = motif_max_probs.get(motif_id, 0.0)

        if max_score <= PSEUDOCOUNT * 10: # Avoid division by zero or near-zero
            continue

        match_score = calculate_pwm_hexamer_prob(pwm, hexamer)
        relative_score = match_score / max_score

        if relative_score >= RELATIVE_SCORE_THRESHOLD:
             # Check if this match is better than the current best for this hexamer
             if relative_score > current_best_relative_score:
                  current_best_relative_score = relative_score
                  best_match['rbp_name'] = data.get('rbp_name', 'Unknown')
                  best_match['motif_id'] = motif_id
                  best_match['match_score'] = match_score
                  best_match['max_score'] = max_score
                  best_match['relative_score'] = relative_score
                  best_match['is_match'] = True
             # Optional: Handle ties or collect all matches above threshold here if needed

    annotations.append(best_match)

print("\nAnnotation complete.")

# 5. Add annotations to DataFrame
print("Adding annotation columns to results...")
anno_df = pd.DataFrame(annotations)
anno_df.columns = ['Matching_RBP_Name', 'Matching_Motif_ID', 'Matching_Score',
                   'Max_Motif_Score', 'Relative_Score', 'Is_Match']

# Ensure indices align if results_df was filtered or modified
results_df.reset_index(drop=True, inplace=True)
anno_df.reset_index(drop=True, inplace=True)

annotated_df = pd.concat([results_df, anno_df], axis=1)

# 6. Write Output File
print(f"Writing annotated results to {OUTPUT_FILE}...")
try:
    annotated_df.to_csv(OUTPUT_FILE, index=False)
    print("Done.")
except Exception as e:
    print(f"Error writing output file {OUTPUT_FILE}: {e}", file=sys.stderr)
