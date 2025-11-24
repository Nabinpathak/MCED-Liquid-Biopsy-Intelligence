import pandas as pd
import os

input_folder = r"D:\MCED_DP\mced_dataset\healthy_CpG"
output_folder = r"D:\MCED_DP\mced_dataset\healthy_bedgraph"

# Create output folder if not exists
os.makedirs(output_folder, exist_ok=True)

for file in os.listdir(input_folder):
    if file.lower().endswith(".bed"):
        input_path = os.path.join(input_folder, file)
        
        # Read BED file
        df = pd.read_csv(input_path, sep="\t", header=None)
        
        # Select first 4 columns
        df_bg = df[[0, 1, 2, 3]]
        
        # Extract filename without extension
        basename = os.path.splitext(file)[0]
        
        # Output file path with same name
        output_path = os.path.join(output_folder, basename + ".bedgraph")
        
        # Save file
        df_bg.to_csv(output_path, sep="\t", header=False, index=False)
        
        print(f"Converted: {file}  -->  {basename}.bedgraph")

print("DONE! All .bed files converted to .bedgraph.")
