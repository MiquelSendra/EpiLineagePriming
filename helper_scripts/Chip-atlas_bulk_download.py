import os
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor

# 1. Define paths
main_folder = r"S:\LAB_MT\RESULTADOS\MiQuel\Paper Endocardio\Sendra_2023\new manuscript\DifferentiationPaper\ES_ATAC"
queries_folder = os.path.join(main_folder, "queries")
data_folder = os.path.join(main_folder, "data")

# 2. Create 'data' folder if not exists
os.makedirs(data_folder, exist_ok=True)

# 3. Find and merge all TSV files
all_tsv_files = [f for f in os.listdir(queries_folder) if f.endswith(".tsv")]
dfs = []

for f in all_tsv_files:
    df = pd.read_csv(os.path.join(queries_folder, f), sep="\t", dtype=str)
    dfs.append(df)

# 4. Concatenate and drop duplicates by the first column ("SRX ID")
merged_df = pd.concat(dfs, ignore_index=True)
merged_df = merged_df.drop_duplicates(subset=merged_df.columns[0])

# 5. Save merged file
merged_queries_path = os.path.join(main_folder, "merged_queries.tsv")
merged_df.to_csv(merged_queries_path, sep="\t", index=False)

# 6. Download .bw files in parallel
def download_bw(srx_id):
    url = f"https://chip-atlas.dbcls.jp/data/mm10/eachData/bw/{srx_id}.bw"
    output_path = os.path.join(data_folder, f"{srx_id}.bw")
    if not os.path.exists(output_path):  # Skip if already downloaded
        try:
            response = requests.get(url, stream=True, timeout=10)
            response.raise_for_status()
            with open(output_path, "wb") as f_out:
                for chunk in response.iter_content(chunk_size=8192):
                    f_out.write(chunk)
            print(f"Downloaded: {srx_id}")
        except Exception as e:
            print(f"Error downloading {srx_id}: {e}")

# 7. Download in parallel using ThreadPoolExecutor
srx_ids = merged_df.iloc[:, 0].tolist()
with ThreadPoolExecutor(max_workers=10) as executor:
    executor.map(download_bw, srx_ids)

print("All downloads completed!")
