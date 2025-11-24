import os

# Input and output folders
input_folder = r"D:\MCED_DP\mced_dataset\tumor_CpG"
output_folder = r"D:\MCED_DP\mced_dataset\tumor_processed"

os.makedirs(output_folder, exist_ok=True)

for filename in os.listdir(input_folder):
    if not filename.endswith(".bedGraph") and not filename.endswith(".bedGraph.gz"):
        continue
    
    infile = os.path.join(input_folder, filename)
    outfile = os.path.join(output_folder, filename.replace(".bedGraph", "_norm.bedGraph"))

    with open(infile, "r", errors="ignore") as fin, open(outfile, "w") as fout:
        for line in fin:
            if line.startswith("track") or line.startswith("#"):
                continue
            
            parts = line.strip().split()
            
            # Tumor data must have at least 6 columns
            if len(parts) < 6:
                continue
            
            chrom = parts[0]
            start = parts[1]
            end = parts[2]
            
            coverage = float(parts[4])
            methyl_reads = float(parts[5])
            
            # avoid division by zero
            if coverage == 0:
                methyl_fraction = 0
            else:
                methyl_fraction = methyl_reads / coverage
            
            fout.write(f"{chrom}\t{start}\t{end}\t{methyl_fraction:.6f}\n")

    print(f"Converted: {filename} → {os.path.basename(outfile)}")

print("\nDONE. All tumor files converted!")
