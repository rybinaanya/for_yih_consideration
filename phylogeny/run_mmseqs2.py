#!/usr/bin/env python
import pandas as pd
import os
from Bio import SeqIO
import shutil
from pathlib import Path

def run_mmseqs2(query_fasta, target_fasta, output_dir, db_id, num_threads=4):
    """
    Runs MMseqs2 to align query protein sequences against target protein sequences.
    
    Parameters:
        query_fasta (str): Path to the query FASTA file.
        target_fasta (str): Path to the target FASTA file.
        output_dir (str): Directory to store MMseqs2 results.
        num_threads (int): Number of threads for MMseqs2.
    
    Returns:
        str: Path to the MMseqs2 result file.
    """
    os.makedirs(output_dir, exist_ok=True)
    db_name = os.path.join(output_dir, db_id)
    result_file = os.path.join(output_dir, f"{db_id}_results.m8")
    tmp_dir=os.path.join(output_dir, 'tmp')
    
    # Create databases
    os.system(f"mmseqs createdb {query_fasta} {db_name}_query --dbtype 2")
    os.system(f"mmseqs createdb {target_fasta} {db_name}_target --dbtype 2")
    
    # Run alignment
    os.system(f"mmseqs search {db_name}_query {db_name}_target {db_name}_result {tmp_dir} --threads {num_threads} --search-type 2")
    
    # Convert results to BLAST-like format
    os.system(f"mmseqs convertalis {db_name}_query {db_name}_target {db_name}_result {result_file} --search-type 2")

    os.system(f"rm -rf {tmp_dir}")
    
    return result_file

def parse_mmseqs_results(result_file):
    """
    Parse MMseqs2 results and extract relevant information.
    
    Parameters:
        result_file (str): Path to MMseqs2 result file.
    
    Returns:
        pd.DataFrame: DataFrame containing parsed results.
    """
    columns = [
        'query', 'target', 'identity', 'alignment_length', 'mismatches', 
        'gap_openings', 'query_start', 'query_end', 'target_start', 
        'target_end', 'evalue', 'bit_score'
    ]
    
    df = pd.read_csv(result_file, sep='\t', header=None, names=columns)
    return df

def download_github_folder(repo, folder_path, local_dir):
    """
    Download a folder from a GitHub repo
    """
    api_url = f"https://api.github.com/repos/{repo}/contents/{folder_path}"
    response = requests.get(api_url)
    response.raise_for_status()
    
    os.makedirs(local_dir, exist_ok=True)
    
    for item in response.json():
        if item['type'] == 'file':
            file_response = requests.get(item['download_url'])
            file_path = os.path.join(local_dir, item['name'])
            
            with open(file_path, 'wb') as f:
                f.write(file_response.content)
            print(f"Downloaded: {item['name']}")
        elif item['type'] == 'dir':
            subfolder_path = os.path.join(local_dir, item['name'])
            download_github_folder(repo, f"{folder_path}/{item['name']}", subfolder_path)

# Download Genomes_fna folder fro github
repo = "zsfrbkv/ShigellaProject"
folder_path = "1Tree/Data/Genomes_fna"
local_dir = "Genomes_fna"
download_github_folder(repo, folder_path, local_dir)


# Get all genome files
genomes_dir = local_dir 
genome_files = [
    os.path.join(genomes_dir, f) 
    for f in os.listdir(genomes_dir) 
    if f.endswith('.fna')
]

query_fasta='yih_genes_k12.fasta'
output_dir='yih_homologous_genes'
mmseqs2_results = []
for f in genome_files:
    db_id='_'.join(os.path.basename(f).split('_')[:2])
    mmseqs_result_file = run_mmseqs2(Path(query_fasta), Path(f), output_dir, db_id)
    mmseqs_result = parse_mmseqs_results(mmseqs_result_file)
    mmseqs_result['genome_id'] = db_id
    mmseqs2_results.append(mmseqs_result)


pd.concat(mmseqs2_results).to_csv('yih_genes_mmseq.tsv', sep='\t', index=False)