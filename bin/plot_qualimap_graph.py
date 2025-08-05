import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv("qualimap_summary.csv")

df["number_of_mapped_reads(%)"] = pd.to_numeric(df["number_of_mapped_reads(%)"], errors='coerce')
df["number_of_unmapped_reads(%)"] = pd.to_numeric(df["number_of_unmapped_reads(%)"], errors='coerce')


plt.figure(figsize=(10, 6))

plt.bar(df["filename"], df["number_of_mapped_reads(%)"], label="Mapped_reads(%)", alpha=0.7, color='red')
plt.bar(df["filename"], df["number_of_unmapped_reads(%)"], label="Unmapped_reads(%)", alpha=0.7, color='blue')

plt.xlabel("File Name")
plt.ylabel("Number of reads (%)")
plt.title("Comparison of Mapped reads and Unmapped reads")
plt.xticks(rotation=90)
plt.legend()
plt.grid(axis="y", linestyle="--", alpha=0.7)

plt.savefig("qualimap_comparison.png", dpi=300, bbox_inches="tight")

