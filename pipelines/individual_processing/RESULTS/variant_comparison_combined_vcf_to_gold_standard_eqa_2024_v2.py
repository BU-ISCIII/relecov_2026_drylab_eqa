import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import re

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
    df["EQA"] = np.nan
    for eqaid in samples_list:
        print(df.columns)
        df.loc[df[sample_col].str.contains(eqaid, na=False), "EQA"] = eqaid
    if "POS" in df.columns:
        df = df.sort_values(by=["EQA", "POS"], ascending=[True, True])
    return df

# --- Función flexible para comparar ALT ---
def comparar_alt(row, col_gold='ALT_gold', col_comp='ALT_comp'):
    """
    Comparación de los valores en columnas ALT donde se recogen 
    los cambios con respecto a la refrencia
    """
    alt_gold = str(row[col_gold]) if pd.notna(row[col_gold]) else ''
    alt_comp = str(row[col_comp]) if pd.notna(row[col_comp]) else ''

    if alt_gold and not alt_comp:
        return f"Cambio en gold standard no encontrado en archivo VCF (revisar Ns, wildtype)"
    if alt_comp and not alt_gold:
        return f"Cambio en archivo VCF no encontado en gold standard (revisar Ns, variantes)"
    if alt_gold and alt_comp:
        if len(alt_gold) == len(alt_comp):
            return "Wrong nucleotide"
        elif len(alt_gold) > len(alt_comp):
            return "Deletion relative to gold standard"
        elif len(alt_gold) < len(alt_comp):
            return "Insertion relative to gold standard"
    return "Unknown"

def comparar_gold_vs(df_comp, col_suffix, extra_cols_gold=None, extra_cols_comp=None):
    """
    df_comp: DataFrame de muestras a comparar con var_gold
    col_suffix: sufijo para las columnas de df_comp ('enviados' o 'isciii')
    extra_cols_gold: lista de columnas del gold standard a incluir (ej: ['DP', 'AF'])
    extra_cols_comp: lista de columnas del comparado a incluir (ej: ['DP', 'AF'])
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
        eqa_validos = df_comp['EQA'].unique()

        gold_filtered = var_gold[var_gold['EQA'].isin(eqa_validos)].copy()
        gold_filtered['REF_gold'] = gold_filtered['REF']
        gold_filtered['ALT_gold'] = gold_filtered['ALT']

        df_comp = df_comp.copy()
        df_comp[f'REF_{col_suffix}'] = df_comp['REF']
        df_comp[f'ALT_{col_suffix}'] = df_comp['ALT']

        merged = gold_filtered.merge(
            df_comp,
            on=['EQA','POS', 'REF', 'ALT'],
            how='outer',
            suffixes=('_gold', f'_{col_suffix}'),
            indicator=True
        )

        if f'REF_{col_suffix}' in merged.columns:
            merged['REF_gold'] = merged['REF_gold'].combine_first(merged[f'REF_{col_suffix}'])
        merged['SAMPLE_ID'] = merged.get(f'SAMPLE_{col_suffix}').combine_first(merged.get('SAMPLE_gold'))

    diff = merged[(merged['_merge'] != 'both') | (merged['ALT_gold'] != merged.get(f'ALT_{col_suffix}', None))]

    # Columnas obligatorias
    for col in [f'ALT_{col_suffix}', f'CHROM_{col_suffix}']:
        if col not in diff.columns:
            diff[col] = pd.NA

    # Columnas extra del gold
    for col in extra_cols_gold:
        col_gold = f"{col}_gold"
        # Solo crear si NO existe ya
        if col_gold not in diff.columns:
            if col in diff.columns:
                diff[col_gold] = diff[col]
            else:
                diff[col_gold] = pd.NA

    # Columnas extra del comparado
    for col in extra_cols_comp:
        if col in diff.columns:
            diff[f"{col}_{col_suffix}"] = diff[col]
        else:
            diff[f"{col}_{col_suffix}"] = pd.NA

    diff['POS'] = diff['POS'].astype(float).astype('Int64')

    diff_final = diff[['SAMPLE_ID','EQA','POS','REF_gold','ALT_gold',f'ALT_{col_suffix}',f'CHROM_{col_suffix}'] +
                      [f"{col}_gold" for col in extra_cols_gold] +
                      [f"{col}_{col_suffix}" for col in extra_cols_comp]].copy()

    diff_final[f'RESULTADOS_{col_suffix}'] = diff_final.apply(
        comparar_alt,
        axis=1,
        col_gold='ALT_gold',
        col_comp=f'ALT_{col_suffix}'
    )

    return diff_final.sort_values(by=['EQA','POS'])


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

VCF_DIR = BASE_DIR / "variants_long_table"
GOLD_DIR = BASE_DIR / "gold_standard/variants_long_table"

# Carpeta fija para resultados
OUTPUT_DIR = BASE_DIR / "variant_analysis_result"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --- Archivos de entrada ---
for gold_standard in GOLD_DIR.rglob("*.csv"):
    if "FLU" in gold_standard:
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
            print(file.name)
            var_enviados = pd.read_csv(file)

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
    match = re.search(r'(COD-\d+)', BASE_NAME)
    sample_tag = match.group(1) if match else BASE_NAME
    #FIXME THIS WILL WORK FOR LATER PURPOSES - JUST NOT NOW, SO HARCODING
    merged_all['COD_LAB'] = "COD-2400"

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

    # --- 💾 8. Guardar CSV final ---
    path_combined = OUTPUT_DIR / f"{BASE_NAME}_diferencias_{gold_standard.name}_COMBINADO_v2.csv"
    merged_all.to_csv(path_combined, index=False)
    print(f"✅ Archivo combinado guardado: {path_combined}")
