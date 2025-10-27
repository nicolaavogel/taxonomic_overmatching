import os
import time
import random
from Bio import Entrez

def download_refseq_sequences(family, email, output_dir, retries=5, sleep_time=30, max_downloads_before_sleep=10):
    Entrez.email = email
    search_statement_refseq = f'"{family}"[organism] AND (refseq[filter] AND ("100000"[SLEN] : "50000000000000"[SLEN]))'
    
    print(f"Searching NCBI for: {search_statement_refseq}")
    search_handle = Entrez.esearch(db="nucleotide", term=search_statement_refseq, retmax=100000)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    ids = search_results["IdList"]
    if not ids:
        print("No sequences found.")
        return
    
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Found {len(ids)} sequences. Downloading...")
    
    download_count = 0
    
    for seq_id in ids:
        attempt = 0
        while attempt < retries:
            try:
                fetch_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
                fasta_data = fetch_handle.read()
                fetch_handle.close()
                
                # Write data to file
                with open(os.path.join(output_dir, f"{seq_id}.fasta"), "w") as f:
                    f.write(fasta_data)
                
                print(f"Downloaded: {seq_id}")
                download_count += 1
                
                # Check if we should sleep after a certain number of downloads
                if download_count % max_downloads_before_sleep == 0:
                    sleep_duration = random.randint(10, 60)  # Random sleep time between 10 and 60 seconds
                    print(f"Sleeping for {sleep_duration} seconds to avoid overload...")
                    time.sleep(sleep_duration)
                
                break  # Break the retry loop if download is successful
            except Exception as e:
                attempt += 1
                print(f"Error downloading {seq_id}, attempt {attempt} of {retries}: {e}")
                
                if attempt < retries:
                    sleep_duration = random.randint(5, 20)  # Random sleep before retrying
                    print(f"Retrying in {sleep_duration} seconds...")
                    time.sleep(sleep_duration)
                else:
                    print(f"Failed to download {seq_id} after {retries} attempts.")
    
    print("Download complete.")

# Usage example
download_refseq_sequences("Rosaceae", "n.alexandra.vogel@gmail.com", "/projects/caeg/data/pp_analysis/Iceland/Rosaceae_db")

#download_refseq_sequences("Poaceae", "n.alexandra.vogel@gmail.com", "/projects/caeg/data/pp_analysis/Iceland/Poaceae_db")
