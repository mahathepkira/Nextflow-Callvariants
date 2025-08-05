import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv("summary.csv")

df["Total_Sequences_BeforeQC"] = pd.to_numeric(df["Total_Sequences_BeforeQC"], errors='coerce')
df["Total_Sequences_AfterQC"] = pd.to_numeric(df["Total_Sequences_AfterQC"], errors='coerce')


plt.figure(figsize=(10, 6))

plt.bar(df["file_name"], df["Total_Sequences_BeforeQC"], label="Before QC", alpha=0.7, color='red')
plt.bar(df["file_name"], df["Total_Sequences_AfterQC"], label="After QC", alpha=0.7, color='blue')

plt.xlabel("File Name")
plt.ylabel("Total Sequences")
plt.title("Comparison of Total Sequences Before and After QC")
plt.xticks(rotation=90)
plt.legend()
plt.grid(axis="y", linestyle="--", alpha=0.7)

plt.savefig("qc_comparison.png", dpi=300, bbox_inches="tight")

