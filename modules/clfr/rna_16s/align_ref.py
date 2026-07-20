import subprocess
import pandas as pd
import re
from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for server
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--bam", type=str, required=True, help="Input BAM file")
parser.add_argument("--ref_fasta", type=str, required=True, help="Reference FASTA file")
parser.add_argument("--outdir", type=str, required=True, help="Output directory")
parser.add_argument("--samtools", type=str, default="samtools", help="Path to samtools")
args = parser.parse_args()

# ---------------- 配置部分 ----------------
bam_file = args.bam
ref_fasta = args.ref_fasta
outdir = args.outdir
samtools = args.samtools
Path(outdir).mkdir(parents=True, exist_ok=True)

output_csv = f"{outdir}/abundance_align_ref.csv"
output_png = f"{outdir}/abundance_align_ref.png"

# Zymo 官方理论丰度 (16S Only 列，细菌部分，百分比)
theoretical = {
    "Pseudomonas aeruginosa": 4.2,
    "Escherichia coli": 10.1,
    "Salmonella enterica": 10.4,
    "Lactobacillus fermentum": 18.4,
    "Enterococcus faecalis": 9.9,
    "Staphylococcus aureus": 15.5,
    "Listeria monocytogenes": 14.1,
    "Bacillus subtilis": 17.4
    # Saccharomyces 和 Cryptococcus 在 16S 中理论上接近 0，不列入
}

# ---------------- 函数：从 FASTA header 提取种名 ----------------
def extract_species_from_fasta(fasta_path):
    species_map = {}
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()[1:]  # 去掉 >
                # 假设 header 如 ">Pseudomonas_aeruginosa_16S" 或类似
                # 根据你的 FASTA 调整正则
                match = re.search(r'^([A-Za-z]+_[a-z]+)', header)  # 取前两个词，如 Pseudomonas_aeruginosa
                if match:
                    species = match.group(1).replace("_", " ")
                    species_map[header] = species
    return species_map

# ---------------- 主流程 ----------------
print("Step 1: Running samtools idxstats...")
idxstats_output = subprocess.check_output([samtools, "idxstats", bam_file]).decode().strip()

# 解析 idxstats 输出
# 格式: ref_name   ref_length   mapped_reads   unmapped_reads
data = []
for line in idxstats_output.splitlines():
    if not line.strip():
        continue
    parts = line.split("\t")
    if len(parts) == 4:
        ref_name, length, mapped, unmapped = parts
        data.append({
            "ref_header": ref_name,
            "length": int(length),
            "mapped_reads": int(mapped),
            "unmapped_reads": int(unmapped)
        })

df = pd.DataFrame(data)

# Step 2: 映射到种名
species_map = extract_species_from_fasta(ref_fasta)
df["species"] = df["ref_header"].map(species_map)

# 过滤掉未映射到的（通常是 18S 或无关）
df = df[df["species"].notna()]

# Step 3: 按种聚合（如果一个种有多条 16S 序列，如多拷贝）
df_agg = df.groupby("species").agg({
    "mapped_reads": "sum",
    "length": "sum"  # 可选：总长度，用于加权
}).reset_index()

# 计算总 mapped reads（只算细菌部分）
total_mapped = df_agg["mapped_reads"].sum()

# Step 4: 计算观测相对丰度 (%)
df_agg["observed_percent"] = (df_agg["mapped_reads"] / total_mapped * 100).round(2)

# Step 5: 加入理论值 & 计算偏差
df_agg["theoretical_percent"] = df_agg["species"].map(theoretical).round(2)
df_agg["difference"] = (df_agg["observed_percent"] - df_agg["theoretical_percent"]).round(2)
df_agg["fold_change"] = (df_agg["observed_percent"] / df_agg["theoretical_percent"]).round(2)  # >1 表示过高

# 排序（按理论丰度降序）
df_agg = df_agg.sort_values("theoretical_percent", ascending=False)

# 输出
print("\nZymo 16S abundance：")
print(df_agg.to_string(index=False))

# 保存到 CSV
df_agg.to_csv(output_csv, index=False)
# print(f"\nsave to: {output_csv}")

# 可视化
df_agg.plot(x="species", y=["theoretical_percent", "observed_percent"], kind="bar")
plt.title("Observed vs Theoretical Abundance (16S)")
plt.ylabel("Relative Abundance (%)")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(output_png)
# print(f"save fig: {output_png}")
