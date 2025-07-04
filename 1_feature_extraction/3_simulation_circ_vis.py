#!/usr/bin/env python3
import os
import sys
import pandas as pd


def convert_tsv_to_csv(input_tsv_path, output_csv_path):
    """
    Convert TSV to CSV and modify column names
    """
    # Check if input file exists
    if not os.path.exists(input_tsv_path):
        raise FileNotFoundError(f"Input file {input_tsv_path} does not exist.")
    
    try:
        # Read TSV file
        df = pd.read_csv(input_tsv_path, sep='\t')
        
        # Rename columns
        column_mapping = {
            'type': 'circrna_type',
            'exon_count': 'block_count'
        }
        df = df.rename(columns=column_mapping)
        
        # Save to CSV
        df.to_csv(output_csv_path, index=False)
        print(f"Converted {input_tsv_path} to {output_csv_path}")
        return df
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

def process_circRNA_metadata(df):
    """
    Process circRNA metadata for visualization
    """
    # Rename 'splice_type' to 'splice_site_type'
    df = df.rename(columns={'splice_type': 'splice_site_type'})
    
    # Define canonical splice site pairs
    canonical_pairs = {
        'GT-AG': True,
        'GC-AG': True,
        'AT-AC': True
    }
    
    # Add 'is_canonical' column
    df['is_canonical'] = df['splice_site_type'].apply(lambda x: canonical_pairs.get(x, False))
    
    # Modify non-standard splice sites
    df.loc[~df['splice_site_type'].isin(canonical_pairs.keys()), 'splice_site_type'] = 'Non-canonical'
    
    # Define donor and acceptor site mapping
    def get_donor_acceptor(splice_type):
        donor_acceptor_map = {
            'GT-AG': ('GT', 'AG'),
            'GC-AG': ('GC', 'AG'),
            'AT-AC': ('AT', 'AC'),
            'Non-canonical': ('Unknown', 'Unknown')
        }
        return donor_acceptor_map.get(splice_type, ('Unknown', 'Unknown'))
    
    # Add donor and acceptor site columns
    df[['donor_site', 'acceptor_site']] = df['splice_site_type'].apply(get_donor_acceptor).apply(pd.Series)
    
    # Add genomic span column
    df['genomic_span'] = df['end'] - df['start']
    
    # Add block_sizes and block_starts columns (placeholder)
    df['block_sizes'] = df['mature_length'].apply(lambda x: f"[{x}]")
    df['block_starts'] = '[0]'
    
    # Ensure proper block_count calculation and grouping
    # Group exon counts to match the visualization script's expectations
    def group_exon_count(count):
        if count == 1:
            return 1
        elif count == 2:
            return 2
        elif count == 3:
            return 3
        elif count == 4:
            return 4
        else:
            return 5  # 5+ category
    
    df['block_count_grouped'] = df['block_count'].apply(group_exon_count)
    return df

def run_visualization_script(processed_csv_path, output_dir):
    """
    Run the visualization script with processed CSV
    """
    # Use direct import of the module by importing the script directly
    vis_script_path = '/scratch/aerusakovich/~/benchmark/Scripts/1_feature_extraction/2_database_vis.py'
    
    # Use importlib to import the module from the file path
    import importlib.util
    spec = importlib.util.spec_from_file_location("database_vis", vis_script_path)
    database_vis = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(database_vis)
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate visualizations
    database_vis.generate_visualizations(processed_csv_path, output_dir)

def main():
    # Input and output paths
    input_tsv_path = '/scratch/aerusakovich/sim_ciri_long_jobim/all/circRNA_metadata.tsv'
    processed_csv_path = '/scratch/aerusakovich/sim_ciri_long_jobim/all/circRNA_processed.csv'
    output_dir = '/scratch/aerusakovich/sim_ciri_long_jobim/all/visualization_output'
    
    # Convert TSV to CSV
    df = convert_tsv_to_csv(input_tsv_path, processed_csv_path)
    
    # Process metadata
    processed_df = process_circRNA_metadata(df)
    processed_df.to_csv(processed_csv_path, index=False)
    print(f"Processed metadata saved to {processed_csv_path}")
    
    # Run visualization
    run_visualization_script(processed_csv_path, output_dir)

if __name__ == "__main__":
    main()