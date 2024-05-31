import pandas as pd


csv_files = [
    "/home/olivia/Downloads/tns_search.csv",
    "/home/olivia/Downloads/tns_search(1).csv",
    "/home/olivia/Downloads/tns_search(2).csv",
    "/home/olivia/Downloads/tns_search(3).csv",
    "/home/olivia/Downloads/tns_search(4).csv",
]

dfs = []

for file in csv_files:
    df = pd.read_csv(file)
    dfs.append(df)

merged_df = pd.concat(dfs, ignore_index=True)

output_file_path = "/home/olivia/Downloads/tns_search_merged.csv"
merged_df.to_csv(output_file_path, index=False)
