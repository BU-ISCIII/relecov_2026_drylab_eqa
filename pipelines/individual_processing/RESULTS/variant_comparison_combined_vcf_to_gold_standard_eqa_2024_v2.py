import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import re
import json
import os.path

# --- Argumentos ---
parser = argparse.ArgumentParser(description="Localizar variants_long_table CSV y Excel, etiquetarlos y guardarlos.")
parser.add_argument("-b", "--base-dir", type=str, default=".", help="Carpeta base donde buscar")
parser.add_argument("-o", "--output-dir", type=str, default=".", help="Carpeta para guardar los CSVs")
args = parser.parse_args()

# --- Función para añadir columna EQA ---
def add_eqa_column(df, samples_list, sample_col="SAMPLE"):
    """
    Añade la columna EQA a la df con la identificación oportuna
    Ordena por muestra y posición
    """
    df = df.copy()
    df["EQA"] = pd.Series(pd.NA, index=df.index, dtype="object")

    for eqaid in samples_list:
        df.loc[df[sample_col].astype(str).str.contains(eqaid, na=False), "EQA"] = eqaid

    if "POS" in df.columns:
        df = df.sort_values(by=["EQA", "POS"], ascending=[True, True])

    return df

# --- Función flexible para comparar ALT ---
def comparar_alt(row, col_ref='REF_gold', col_gold='ALT_gold', col_comp='ALT_comp'):
    """
    Clasifica diferencias entre gold standard y muestra reportada.

    Reglas:
    - wrong_nt: misma POS y misma REF, ALT distinto, ambos de longitud 1
    - insertions: inserción no presente en gold
    - deletions: deleción no presente en gold
    - missing: variante presente en gold y ausente en muestra
    - denovo: SNP de novo de un solo nucleótido en muestra, ausente en gold
    """
    ref = str(row[col_ref]) if pd.notna(row[col_ref]) else ''
    alt_gold = str(row[col_gold]) if pd.notna(row[col_gold]) else ''
    alt_comp = str(row[col_comp]) if pd.notna(row[col_comp]) else ''

    has_gold = bool(alt_gold)
    has_comp = bool(alt_comp)

    # Caso 1: está en gold pero no en muestra
    if has_gold and not has_comp:
        return "Missing variant"

    # Caso 2: está en muestra pero no en gold
    if has_comp and not has_gold:
        # Inserción de novo
        if len(alt_comp) > len(ref):
            return "Insertion relative to gold standard"
        # Deleción de novo
        elif len(ref) > len(alt_comp):
            return "Deletion relative to gold standard"
        # SNP de novo
        elif len(ref) == 1 and len(alt_comp) == 1:
            return "De novo reported variant"
        else:
            return "De novo reported variant"

    # Caso 3: están ambos, comparar
    if has_gold and has_comp:
        # iguales, no debería entrar aquí normalmente
        if alt_gold == alt_comp:
            return "Match"

        # mismo REF, distinto ALT
        # SNP incorrecto
        if len(ref) == 1 and len(alt_gold) == 1 and len(alt_comp) == 1:
            return "Wrong nucleotide"

        # inserción relativa a gold
        if len(alt_comp) > len(alt_gold):
            return "Insertion relative to gold standard"

        # deleción relativa a gold
        if len(alt_comp) < len(alt_gold):
            return "Deletion relative to gold standard"

        # fallback
        return "Wrong nucleotide"

    return "Unknown"


def comparar_gold_vs(df_comp, col_suffix, extra_cols_gold=None, extra_cols_comp=None):
    """
    df_comp: DataFrame de muestras a comparar con var_gold
    col_suffix: sufijo para las columnas de df_comp ('enviados' o 'isciii')
    extra_cols_gold: lista de columnas del gold standard a incluir
    extra_cols_comp: lista de columnas del comparado a incluir
    """
    if extra_cols_gold is None:
        extra_cols_gold = []
    if extra_cols_comp is None:
        extra_cols_comp = []

    if df_comp.empty:
        merged = var_gold.copy()
        merged['SAMPLE_ID'] = pd.NA
        for col in [f'ALT_{col_suffix}', f'CHROM_{col_suffix}'] + extra_cols_comp:
            merged[col] = pd.NA
    else:
        eqa_validos = df_comp['EQA'].dropna().unique()

        gold_filtered = var_gold[var_gold['EQA'].isin(eqa_validos)].copy()
        df_comp = df_comp.copy()

        gold_filtered = gold_filtered.rename(columns={
            'SAMPLE': 'SAMPLE_gold',
            'CHROM': 'CHROM_gold',
            'REF': 'REF_gold',
            'ALT': 'ALT_gold'
        })

        df_comp = df_comp.rename(columns={
            'SAMPLE': f'SAMPLE_{col_suffix}',
            'CHROM': f'CHROM_{col_suffix}',
            'REF': f'REF_{col_suffix}',
            'ALT': f'ALT_{col_suffix}'
        })

        merged = gold_filtered.merge(
            df_comp,
            left_on=['EQA', 'POS', 'REF_gold'],
            right_on=['EQA', 'POS', f'REF_{col_suffix}'],
            how='outer',
            indicator=True
        )

        merged['REF_gold'] = merged['REF_gold'].combine_first(merged[f'REF_{col_suffix}'])
        merged['SAMPLE_ID'] = merged.get(f'SAMPLE_{col_suffix}').combine_first(merged.get('SAMPLE_gold'))

    diff = merged[
        (merged['_merge'] != 'both') |
        (merged['ALT_gold'] != merged.get(f'ALT_{col_suffix}', pd.NA))
    ].copy()

    for col in [f'ALT_{col_suffix}', f'CHROM_{col_suffix}']:
        if col not in diff.columns:
            diff[col] = pd.NA

    for col in extra_cols_gold:
        col_gold = f"{col}_gold"
        if col_gold not in diff.columns:
            if col in diff.columns:
                diff[col_gold] = diff[col]
            else:
                diff[col_gold] = pd.NA

    for col in extra_cols_comp:
        col_comp = f"{col}_{col_suffix}"
        if col_comp not in diff.columns:
            if col in diff.columns:
                diff[col_comp] = diff[col]
            else:
                diff[col_comp] = pd.NA

    diff['POS'] = pd.to_numeric(diff['POS'], errors='coerce').astype('Int64')

    diff_final = diff[
        ['SAMPLE_ID', 'EQA', 'POS', 'REF_gold', 'ALT_gold', f'ALT_{col_suffix}', f'CHROM_{col_suffix}'] +
        [f"{col}_gold" for col in extra_cols_gold] +
        [f"{col}_{col_suffix}" for col in extra_cols_comp]
    ].copy()

    diff_final[f'RESULTADOS_{col_suffix}'] = diff_final.apply(
        comparar_alt,
        axis=1,
        col_ref='REF_gold',
        col_gold='ALT_gold',
        col_comp=f'ALT_{col_suffix}'
    )

    return diff_final.sort_values(by=['EQA', 'POS'])


def calculate_values_eqa(merged_all: pd.DataFrame, vlt_lab, vlt_gold):
    variants_dict = {}
    resultados_map = {
              "Wrong nucleotide": "wrong_nt",
              "Insertion relative to gold standard": "insertions",
              "Deletion relative to gold standard": "deletions",
              "Missing variant": "missing",
              "De novo reported variant": "denovo"
            }
    
    # Set up the initial values
    df = merged_all.copy()
    vlt_df = vlt_lab.copy()
    vlt_gold_df = vlt_gold.copy()
    af_threshold = 0.75

    # Find number of successful hits
    for sample, vlt_df_sample in vlt_df.groupby("SAMPLE"):
        set1 = set(map(tuple, vlt_df_sample[["POS", "REF", "ALT"]].values))
        set2 = set(map(tuple, vlt_gold_df[vlt_gold_df["SAMPLE"] == sample][["POS", "REF", "ALT"]].values))

        num_matches = len(set1 & set2)
        if sample not in variants_dict:
            variants_dict[sample] = {}
        variants_dict[sample]["successful_hits"] = num_matches
    


    for sample, group in vlt_df.groupby("SAMPLE"):
        if sample not in variants_dict:
            variants_dict[sample] = {}

        # Check presence of low- and high-frequency alleles
        af_values = pd.to_numeric(group["AF"], errors="coerce")

        has_low_freq = (af_values < af_threshold).any()
        has_high_freq = (af_values > af_threshold).any()

        variants_dict[sample]["high_and_low_freq"] = bool(has_high_freq and has_low_freq)
        variants_dict[sample]["high_freq_only"] = bool(has_high_freq and not has_low_freq)
        variants_dict[sample]["low_freq_only"] = bool(has_low_freq and not has_high_freq)

        # Keep only variants with AF > 0.75
        group_high_af = group[af_values > af_threshold].copy()

        # Number of variants in consensus considering only AF > 0.75
        n_variants_consensus = len(group_high_af)
        variants_dict[sample]["number_of_variants_in_consensus_vcf"] = n_variants_consensus
        
        # Number of variants with effect considering only AF > 0.75
        df_variants_effect = group_high_af[group_high_af["EFFECT"] == "missense_variant"]
        n_variants_effect = len(df_variants_effect)
        variants_dict[sample]["number_of_variants_with_effect_vcf"] = n_variants_effect

        if df.empty:
            for key, json_equivalent in resultados_map.items():
                variants_dict[sample][json_equivalent] = 0

    for sample, reported_df in df.groupby("EQA"):
        resultados_counts = reported_df["RESULTADOS_enviados"].value_counts(dropna=False).to_dict()
        for key, json_equivalent in resultados_map.items():
            variants_dict[sample][json_equivalent] = int(resultados_counts.get(key,0))

    if os.path.isfile("variants_report.json"):
        try:
            with open("variants_report.json", "r") as f:
                existing = json.load(f)
            if existing:
                existing.update(variants_dict)
                variants_dict = existing
        except json.JSONDecodeError:
            pass
    with open("variants_report.json", "w") as f:
        json.dump(variants_dict, f)

def get_pattern_vlt(gold_standard_filename):
    if "FLU" in gold_standard_filename:
        organism = "Influenza"
    else:
        organism = "coronavirus_2"
    if "1" in gold_standard_filename:
        tech = "Illumina"
    else:
        tech = "Nanopore"
    return f"*{organism}*{tech}*csv"

BASE_DIR = Path(args.base_dir).resolve()
BASE_NAME = BASE_DIR.name  # Nombre para los archivos de salida

COD_DIR = Path(f"{args.base_dir}/../").resolve().name

VCF_DIR = BASE_DIR / "variants_long_table"
GOLD_DIR = BASE_DIR / "gold_standard/variants_long_table"

# Carpeta fija para resultados
OUTPUT_DIR = BASE_DIR / "variant_analysis_result"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --- Archivos de entrada ---
for gold_standard in GOLD_DIR.rglob("*.csv"):
    if "FLU" in gold_standard.name:
        continue
    var_gold = pd.read_csv(gold_standard)
    # Renombrar DP y AF antes de cualquier merge
    if 'DP' in var_gold.columns:
        var_gold.rename(columns={'DP': 'DP_gold'}, inplace=True)
    if 'AF' in var_gold.columns:
        var_gold.rename(columns={'AF': 'AF_gold'}, inplace=True)

    var_enviados = pd.DataFrame()

    # --- Lista de muestras ---
    samples = list(set(var_gold.SAMPLE.values))

    # --- Buscar archivos CSV en la carpeta base ---
    # Primero, buscar el pattern del long table basandonos en el gold standard
    variant_long_table_pattern = get_pattern_vlt(gold_standard.name)

    for file in VCF_DIR.rglob(variant_long_table_pattern):
        if "full" in file.name:
            continue
        else:
            var_enviados = pd.read_csv(file)
    if var_enviados.empty:
        continue

    # --- Añadir columna EQA ---
    var_gold = add_eqa_column(var_gold, samples)
    var_enviados = add_eqa_column(var_enviados, samples)

    # --- Asegurar tipos SI NO ESTÄN VACIOS ---
    for df in (var_gold, var_enviados):
        if not df.empty:
            for col in ['EQA','POS','ALT','REF','CHROM']:
                if col in df.columns:
                    df[col] = df[col].astype(str)


    # --- Comparar y guardar ---

    diff_enviados = comparar_gold_vs(var_enviados, 'enviados',
                                    extra_cols_gold=['DP','AF'],
                                    extra_cols_comp=['DP','AF'])

                                    
    ##########################################

    # --- 🧮 2. Resolver duplicados DP_gold y AF_gold ---
    def combine_cols(col_x, col_y):
        """Combina columnas duplicadas prefiriendo valores no nulos o idénticos."""
        return np.where(
            pd.notna(col_x) & pd.notna(col_y) & (col_x == col_y),
            col_x,
            np.where(pd.notna(col_x), col_x, col_y)
        )

    merged_all = diff_enviados

    # --- 🧾 4. Añadir COD_LAB desde BASE_NAME ---
    match = re.search(r'(COD-\d+)', COD_DIR)
    sample_tag = match.group(1) if match else BASE_NAME
    #FIXME THIS WILL WORK FOR LATER PURPOSES - JUST NOT NOW, SO HARCODING
    merged_all['COD_LAB'] = sample_tag
    # --- 🗃️ 5. Ordenar columnas ---
    final_cols = [
        'COD_LAB', #'SAMPLE_ID',
        'EQA', 'POS',
        'REF_gold', 'ALT_gold',
        'ALT_enviados',
        'DP_gold','DP_enviados',
        'AF_gold','AF_enviados',
        'RESULTADOS_enviados',
    ]

    for col in final_cols:
        if col not in merged_all.columns:
            merged_all[col] = pd.NA

    merged_all = merged_all[final_cols]

    # --- 🔤 6. Renombrar columnas para presentación ---
    merged_all.rename(columns={
        'REF_gold': 'REF (NC_045512.2)',
        'ALT_gold': 'Gold_Standard'
    }, inplace=True)

    # --- 📊 7. Ordenar por EQA y posición ---
    merged_all = merged_all.sort_values(by=['EQA','POS'])

    # --- 7.5 Calcular stats y guardar ---
    calculate_values_eqa(merged_all, var_enviados, var_gold)

    # --- 💾 8. Guardar CSV final ---
    path_combined = OUTPUT_DIR / f"{sample_tag}_diferencias_{gold_standard.name}_COMBINADO_v2.csv"
    merged_all.to_csv(path_combined, index=False)
    print(f"✅ Archivo combinado guardado: {path_combined}")

