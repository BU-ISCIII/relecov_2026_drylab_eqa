{# =========================
  JINJA2 TEMPLATE: RELECOV EQA REPORT
  Context expected:
    - general: dict (from general.json)
    - labdata: dict (from lab_<LAB_COD>.json)
========================= #}

{# ---------- Helpers / macros ---------- #}
{% macro pct(x, decimals=1) -%}
  {{ "%.{}f".format(decimals)|format(x) }}%
{%- endmacro %}
{% macro render_figure(path, caption=None, figure_class=None, has_panels=False) -%}
{% if path -%}
{% set _figure_style = figure_cfg.style if figure_cfg is defined and figure_cfg.style else ('max-width: 100%;' if has_panels else 'max-width: 80%;') %}
{% set _auto_width = '100%' if has_panels else '80%' %}
{% set _disable_max_width = false %}
{% if 'max-width:' in _figure_style and 'width:' not in _figure_style|replace('max-width', '') %}
{% set _auto_width = _figure_style|replace('max-width:', '')|replace(';', '')|trim %}
{% if '%' in _auto_width and (_auto_width|replace('%', '')|float) > 100 %}
{% set _disable_max_width = true %}
{% endif %}
{% elif 'width:' in _figure_style %}
{% set _auto_width = _figure_style|replace('width:', '')|replace(';', '')|trim %}
{% endif %}
<figure{% if figure_class %} class="{{ figure_class }}"{% endif %}>
<img src="{{ path }}" alt="{{ caption|default('Figure') }}" style="width: {{ _auto_width }}; {% if _disable_max_width %}max-width: none;{% else %}{{ _figure_style }}{% endif %}"/>
</figure>
{%- endif %}
{%- endmacro %}
{% macro software_label(name, version=None, db_version=None) -%}
{%- if name -%}{{ name }}{{ " (" ~ version ~ ")" if version else "" }}{{ "; DB " ~ db_version if db_version else "" }}{%- else -%}NA{%- endif -%}
{%- endmacro %}
{% macro discrepancy_label(key) -%}
{%- set labels = {
  "wrong_nt": "Wrong nucleotide",
  "ambiguity2nt": "Ambiguity instead of nucleotide",
  "nt2ambiguity": "Nucleotide instead of ambiguity",
  "ns2nt": "Stretch of Ns instead of nucleotide stretch",
  "nt2ns": "Nucleotide stretch instead of stretch of Ns",
  "insertions": "Insertion relative to gold standard",
  "deletions": "Deletion relative to gold standard",
  "missing": "Missing variant",
  "denovo": "De novo variant"
} -%}
{{ labels.get(key, key if key is not none else "NA") }}
{%- endmacro %}

{% set fig_counter = namespace(value=0) %}
{% set table_counter = namespace(value=0) %}
{% set appendix_fig_counter = namespace(value=0) %}
{% set appendix_table_counter = namespace(value=0) %}
{% set figure_cfg = namespace(style=None) %}
{% set consensus_appendix_entries = namespace(value=[]) %}
{% set variant_sars_appendix_entries = namespace(value=[]) %}
{% set variant_flu_appendix_entries = namespace(value=[]) %}
{% set classification_appendix_entries = namespace(value=[]) %}
{% set qc_appendix_entries = namespace(value=[]) %}
{% set benchmark_appendix_entries = namespace(value=[]) %}
{% set lab_consensus_appendix_entries = namespace(value=[]) %}
{% set lab_variant_figure_appendix_entries = namespace(value=[]) %}
{% set lab_classification_figure_appendix_entries = namespace(value=[]) %}
{% set lab_workflow_figure_appendix_entries = namespace(value=[]) %}
{% set lab_qc_figure_appendix_entries = namespace(value=[]) %}
{% set lab_metadata_metrics_appendix_entries = namespace(value=[]) %}
{% set metadata_metric_labels = {
  "per_genome_greater_10x": "% Genome > 10x",
  "depth_of_coverage_value": "Depth of coverage mean value",
  "per_Ns": "% Ns",
  "per_reads_virus": "% Reads virus",
  "per_reads_host": "% Reads host"
} %}

# RELECOV 2.0 - Consolidation of WGS and RT-PCR activities for SARS-CoV-2 in Spain towards sustainable use and integration of enhanced infrastructure and capacities in the RELECOV network

##### Sarai Varona, Enrique Sapena, Pablo Mata, Alejandro Bernabéu, Pau Pascual, Magdalena Matito, Juan Ledesma, Sara Monzón, Isabel Cuesta

## Table of Contents

- [Executive Summary](#executive-summary)
- [1. Introduction](#1-introduction)
- [2. Scope of the EQA](#2-scope-of-the-eqa)
- [3. Dataset Design and Sample Selection Criteria](#3-dataset-design-and-sample-selection-criteria)
    - [3.1. Rationale for Dataset Selection](#31-rationale-for-dataset-selection)
    - [3.2. SARS-CoV-2 Dataset Selection](#32-sars-cov-2-dataset-selection)
    - [3.3. Influenza Dataset Selection](#33-influenza-dataset-selection)
- [4. Methodology of Evaluation](#4-methodology-of-evaluation)
    - [4.1. Submission Completeness](#41-submission-completeness)
    - [4.2. Evaluation of Consensus Genome Reconstruction Performance](#42-evaluation-of-consensus-genome-reconstruction-performance)
    - [4.3. Evaluation of Variant Detection Accuracy](#43-evaluation-of-variant-detection-accuracy)
    - [4.4. Evaluation of Lineage, Subtype and Clade Assignment](#44-evaluation-of-lineage-subtype-and-clade-assignment)
    - [4.5. Evaluation of Metadata Completeness and Compliance](#45-evaluation-of-metadata-completeness-and-compliance)
    - [4.6. Pipeline Benchmarking and Comparative Performance](#46-pipeline-benchmarking-and-comparative-performance)
- [5. General Results](#5-general-results)
    - [5.1. Submission Completeness](#51-submission-completeness)
    - [5.2. Consensus Genome Reconstruction Performance](#52-consensus-genome-reconstruction-performance)
    - [5.3. Variant Detection Accuracy](#53-variant-detection-accuracy)
    - [5.4. Lineage, Subtype and Clade Assignment](#54-lineage-subtype-and-clade-assignment)
    - [5.5. Metadata Completeness and Compliance](#55-metadata-completeness-and-compliance)
    - [5.6. Pipeline Benchmarking and Comparative Performance](#56-pipeline-benchmarking-and-comparative-performance)
- [6. Component-Specific Results](#6-component-specific-results)
    - [6.1. SARS1 (SARS-CoV-2, Illumina)](#61-sars1-sars-cov-2-illumina)
    - [6.2. SARS2 (SARS-CoV-2, Oxford Nanopore Technologies)](#62-sars2-sars-cov-2-oxford-nanopore-technologies)
    - [6.3. FLU1 (Influenza virus, Illumina)](#63-flu1-influenza-virus-illumina)
    - [6.4. FLU2 (Influenza virus, Oxford Nanopore Technologies)](#64-flu2-influenza-virus-oxford-nanopore-technologies)
- [7. Discussion](#7-discussion)
    - [7.1. Consensus Genome Reconstruction](#71-consensus-genome-reconstruction)
    - [7.2. Variant Detection and Reporting](#72-variant-detection-and-reporting)
    - [7.3. Classification and QC Interpretation](#73-classification-and-qc-interpretation)
    - [7.4. Workflow Diversity and Reporting Constraints](#74-workflow-diversity-and-reporting-constraints)
    - [7.5. Metadata Reporting and Schema Compliance](#75-metadata-reporting-and-schema-compliance)
    - [7.6. Implications for RELECOV 2.0](#76-implications-for-relecov-20)
- [8. Conclusions](#8-conclusions)
{% if labdata %}
- [9. Individual Laboratory Technical Report](#9-individual-laboratory-technical-report)
{% endif %}
- [Appendix](#appendix)

## Executive Summary

The 2026 RELECOV Dry-Lab EQA provides the first network-wide dry-lab assessment focused specifically on bioinformatic analytical performance across respiratory virus surveillance workflows. The exercise evaluated 20 distributed datasets grouped into four analytical components, comprising Illumina and Nanopore for SARS-CoV-2 and influenza. Participating laboratories were assessed on consensus genome reconstruction, variant reporting, lineage/type and clade assignment, metadata completeness, and the reproducibility of their declared analytical workflows relative to curated gold standards and network-wide distributions.

Nineteen laboratories participated, corresponding to {{ pct(general.total_participants_pct, 2) }} of invited laboratories, with high submission rates for expected analytical outputs: {{ pct(general.submission_rates_pct.fasta, 2) }} for consensus genome files and {{ pct(general.submission_rates_pct.vcf, 2) }} for VCF files.

Across the network, consensus genome reconstruction performed best in the Illumina-based components, with a combined median genome identity of {{ pct(general.general_results.consensus.median_identity_illumina_pct, 2) }}, compared with {{ pct(general.general_results.consensus.median_identity_nanopore_pct, 2) }} in the Nanopore-based components. However, broad identity ranges in SARS2 and FLU2 indicate that outlier submissions remained present, particularly in contexts where masking, coverage thresholds, and consensus-generation choices differed across laboratories.

Variant reporting showed clear methodological heterogeneity. For SARS-CoV-2, the median number of discrepancies relative to the curated reference variant sets was {{ general.general_results.sars_variants.median_discrepancy_illumina }} in the Illumina component and {{ general.general_results.sars_variants.median_discrepancy_nanopore }} in the Nanopore component. Influenza submissions were more heterogeneous structurally: the median number of high-frequency variants reported in metadata was {{ general.general_results.influenza_variants.median_variants_in_consensus }}, whereas the corresponding median obtained from submitted VCF files after filtering to AF >=75% was {{ general.general_results.influenza_variants.median_variants_in_consensus_vcf }}. The median discrepancy between these two representations was {{ general.general_results.influenza_variants.median_discrepancies_in_reported_variants }}, indicating that metadata declarations and filtered VCF content were not always directly concordant.

Classification performance was consistently higher for lineage/type assignment than for clade assignment. SARS-CoV-2 lineage concordance reached {{ pct(general.general_results.classification.sars_lineage_concordance_pct) }}, compared with {{ pct(general.general_results.classification.sars_clade_concordance_pct) }} for clade assignment, while influenza type/subtype concordance reached {{ pct(general.general_results.classification.influenza_type_concordance_pct) }}, compared with {{ pct(general.general_results.classification.influenza_clade_concordance_pct) }} for clade assignment. Review of submitted files further indicated that part of the excess discordance in clade assignment reflected field completion and nomenclature issues, including missing clade entries and lineage/type-like values entered in the clade field.

Metadata completeness and reporting remain major priorities for harmonisation. The median metadata completeness rate across participating laboratories was {{ pct(general.metadata_completeness.median_pct) }}, with values ranging from {{ pct(general.metadata_completeness.min_pct) }} to {{ pct(general.metadata_completeness.max_pct) }}. Although software names were reported for {{ pct(general.metadata_completeness.software_names_pct) }} of expected fields, only {{ pct(general.metadata_completeness.software_version_pct) }} of software-version fields, {{ pct(general.metadata_completeness.coverage_threshold_pct) }} of coverage thresholds, {{ pct(general.metadata_completeness.variant_calling_params_pct) }} of variant-calling parameter fields, and {{ pct(general.metadata_completeness.reference_genome_pct) }} of reference genome identifiers were completed. A total of {{ general.metadata_completeness.total_workflows }} distinct workflows were identified, together with substantial diversity in consensus, variant calling, and classification software.

Overall, the EQA shows that RELECOV laboratories already have substantial bioinformatic capacity, but also that inter-laboratory comparability is limited by heterogeneous thresholds, parameter reporting, reference selection, and uneven completion of metadata and QC fields. These findings support RELECOV 2.0 priorities centred on minimum performance standards, stronger metadata requirements, clearer reporting rules for consensus and variants, and component-aware benchmarking rather than a single cross-context pipeline ranking.

## 1. Introduction

The RELECOV Network aims to strengthen genomic surveillance of respiratory viruses by developing and harmonising analytical capacities across the participating laboratories. In this context, it was essential to **assess the consistency, reproducibility and maturity of the bioinformatics workflows implemented within the network**.

To this end, an **external quality assessment (EQA) exercise in dry lab format** was conducted, inspired by the ECDC’s 2024 dry-lab EQA. The exercise focused on the bioinformatic characterisation of respiratory viruses, covering key analytical tasks including viral genome reconstruction, variant identification, and lineage and clade assignment.

A central component of this initiative was to evaluate the range of analytical pipelines used across the RELECOV Network, identify their relative performance, and determine which approaches were best suited to support a genomic surveillance standard within the network. This evaluation contributes directly to **Objective 2.1** of RELECOV 2.0, which focuses on _improving deep knowledge of the capacities and methodologies of the laboratories belonging to the network, as well as identifying a common methodology adapted to them and to the needs of the platform_. Furthermore, the EQA provides the practical evidence base required for **Task T6.1**, which aims to  _identify the most suitable bioinformatic analysis method for each sequencing platform, through an intercomparison exercise with simulated data for bioinformaticians_, in order to define the workflow that should be integrated into the RELECOV analytical platform.

The exercise was also aligned with **Milestone M6.3**, which pertains to _define sequencing and analysis protocols for each of the sequencing platforms_. In addition, the exercise provided operational insights relevant to **Task T6.5**, which addresses _the adaptation and improvement of the analysis pipeline for the different sequencing platforms used by the laboratories of the network_. It also contributed to **Task T6.4**, related to _sequence metadata annotation with ontologies, schema generation, parsing and validation_, by highlighting practical issues affecting metadata completeness, controlled-vocabulary use, and the consistency of reported analytical parameters.

The overall objective of the exercise was to **assess the bioinformatic performance of the participating laboratories, identify areas for improvement, and promote the adoption of consistent and comparable analytical practices across the network**. The outcomes presented in this report strengthen RELECOV’s preparedness and response capacity in routine surveillance and public health emergencies, support the quality and robustness of genomic analyses performed throughout the network, and contribute directly to the fulfilment of key project objectives, milestones and deliverables.

## 2. Scope of the EQA

The 2026 RELECOV Dry-Lab External Quality Assessment (EQA) was designed to evaluate the bioinformatic analytical performance of laboratories participating in the RELECOV Network in the context of respiratory virus genomic surveillance.

Participating laboratories were provided with raw sequencing datasets corresponding to four independent analytical components:

- **SARS1**: Five SARS-CoV-2 samples sequenced using paired-end Illumina technology from the 2024 ECDC ESIB EQA
- **SARS2**: Five SARS-CoV-2 samples sequenced using Oxford Nanopore Technologies from the 2024 ECDC ESIB EQA
- **FLU1**: Five influenza virus samples sequenced using paired-end Illumina technology, 3 generated in-silico and 2 from the 2024 ECDC ESIB EQA.
- **FLU2**: Five influenza virus samples sequenced using Oxford Nanopore Technologies, 3 generated in-silico and 2 from the 2024 ECDC ESIB EQA.

Datasets were distributed as raw sequencing reads (.fastq files), and laboratories were free to analyse any subset of components according to their technical capacity and routine workflow.
Laboratories were requested to submit:

- For each analysed sample:
    - One consensus genome sequence in .fasta format, containing exclusively the target viral genome reconstructed from the provided reads.
    - One or more variant call files in .vcf format, listing detected nucleotide variants relative to the reference genome selected by the laboratory.
- A completed harmonised metadata template, documenting analytical tools, software versions, reference genomes used, parameter settings, coverage thresholds, Lineage, Subtype or clade assignment tools, and file paths to submitted outputs.

The primary objective of the exercise was to assess the consistency, reproducibility, and comparability of bioinformatic workflows currently implemented across the network. The evaluation focused on core analytical tasks that are essential for routine genomic surveillance and public health response, including:

- **Viral genome reconstruction**: Generation of high-quality consensus genome sequences from raw sequencing reads produced using Illumina and Oxford Nanopore Technologies platforms.
- **Variant identification and reporting**: Detection and annotation of nucleotide variants relative to a chosen reference genome, including evaluation of filtering criteria, allele frequency thresholds, and variant file standardisation.
- **Lineage, Subtype and clade assignment**: Accurate classification of reconstructed genomes using established nomenclature systems and version-controlled databases.
- **Metadata reporting and interoperability**: Completion of a harmonised metadata template capturing software versions, analytical parameters, reference genome selection, and file traceability, ensuring compatibility with automated validation and integration into the RELECOV analytical platform.

## 3. Dataset Design and Sample Selection Criteria

### 3.1. Rationale for Dataset Selection

The 2026 RELECOV Dry-Lab EQA was specifically designed for laboratories operating in a clinical and hospital-based diagnostic context, where routine genomic surveillance primarily involves human respiratory samples.

Sample selection followed three guiding principles:

- Representation of realistic genomic surveillance scenarios.
- Inclusion of defined analytical challenges.
- Ensuring methodological benchmarking robustness.

Datasets were derived from two sources:

- Reused datasets from the 2024 ECDC ESIB Dry-Lab EQA.
- Newly generated in-silico datasets constructed to simulate seasonal human influenza circulation.

The integration of both sources allowed alignment with internationally validated materials while tailoring the exercise to the operational reality of RELECOV clinical laboratories.

### 3.2. SARS-CoV-2 Dataset Selection

SARS-CoV-2 datasets were selected from the 2024 ECDC ESIB EQA to ensure comparability with internationally benchmarked material. Both Illumina and Nanopore panels included samples representing:

- High-quality baseline genomes.
- Low read-depth scenarios.
- Samples with numerous mixed sites.
- Contamination with non-target viral reads.
- Lineages of epidemiological relevance (e.g., recombinant or XBB-related lineages).

Only samples generated using the same ARTIC primer scheme (v4.1) were selected to avoid introducing variability associated with enrichment panel differences. This ensured that observed performance differences reflect analytical workflow characteristics rather than primer design heterogeneity.

{% set table_counter.value = table_counter.value + 1 %}

_**Table {{ table_counter.value }}**. Overview of SARS-CoV-2 datasets used in the RELECOV 2026 Dry-Lab EQA.
The table details sample origin, sequencing technology (Illumina paired-end or Oxford Nanopore Technologies), amplicon primer scheme version, and specific analytical characteristics intentionally selected to assess workflow robustness under challenging conditions._

| Sample | Source             | Platform | Amplicon primers version | Ref sample | Key Feature                                       | FASTQ files | Read layout | Clade Assignment | Lineage Assignment | Quality check |
|--------|--------------------|----------|--------------------------|------------|---------------------------------------------------|-------------|-------------|------------------|--------------------|---------------|
| SARS1  | ECDC-ESIB EQA 2024 | Illumina | ARTIC v4.1               | SARS2.04   | Influenza virus sample with some SARS-CoV-2 reads | 2           | Paired-end  | -                | -                  | Bad           |
| SARS2  | ECDC-ESIB EQA 2024 | Illumina | ARTIC v4.1               | SARS2.01   | High-quality baseline sample                      | 2           | Paired-end  | 21K              | BA.1.13            | Ok            |
| SARS3  | ECDC-ESIB EQA 2024 | Illumina | ARTIC v4.1               | SARS2.16   | XBB sample / insertion challenge                  | 2           | Paired-end  | 23A              | XBB.1.5            | Ok            |
| SARS4  | ECDC-ESIB EQA 2024 | Illumina | ARTIC v4.1               | SARS2.20   | Very low read depth                               | 2           | Paired-end  | 22E              | BQ.1.1             | Bad           |
| SARS5  | ECDC-ESIB EQA 2024 | Illumina | ARTIC v4.1               | SARS2.13   | >10 mixed sites                                   | 2           | Paired-end  | 21K              | BA.1.1             | Bad           |
| SARS6  | ECDC-ESIB EQA 2024 | Nanopore | ARTIC v4.1               | SARS1.01   | High-quality baseline sample                      | 1           | Single-end  | recombinant      | XCH.1              | Ok            |
| SARS7  | ECDC-ESIB EQA 2024 | Nanopore | ARTIC v4.1               | SARS1.09   | XBB sample / ambiguity next to a deletion         | 1           | Single-end  | 23A              | XBB.1.5.24         | Ok            |
| SARS8  | ECDC-ESIB EQA 2024 | Nanopore | ARTIC v4.1               | SARS1.15   | >10 mixed sites                                   | 1           | Single-end  | 23D              | XBB.1.9.1          | Bad           |
| SARS9  | ECDC-ESIB EQA 2024 | Nanopore | ARTIC v4.1               | SARS1.12   | Influenza virus sample with some SARS-CoV-2 reads | 1           | Single-end  | -                | -                  | Bad           |
| SARS10 | ECDC-ESIB EQA 2024 | Nanopore | ARTIC v4.1               | SARS1.05   | Very low read depth                               | 1           | Single-end  | -                | -                  | Bad           |

### 3.3. Influenza Dataset Selection

The influenza datasets provided in the 2024 ECDC ESIB EQA predominantly correspond to zoonotic influenza strains of animal origin, including H5N1, H5N6, and reassortant genomes.

While these datasets are valuable for specialised surveillance contexts, they do not represent the routine analytical scenario encountered by most RELECOV laboratories, which primarily process:

- Seasonal human Influenza A/H1N1
- Seasonal human Influenza A/H3N2

Given that the objective of this EQA is to benchmark bioinformatic workflows in a clinical hospital environment, it was considered methodologically necessary to include representative seasonal human influenza strains.

Therefore, selected ECDC influenza samples were complemented with newly generated in-silico datasets designed to simulate:

- Seasonal H1N1 circulation.
- Seasonal H3N2 circulation.

This design ensures that evaluation reflects the analytical demands of the RELECOV network.

#### In-Silico Influenza Dataset Construction

The following seasonal clades were selected as reference backbones:

- H1N1 clade D.3.1.1
- H1N1 clade C.1.9.3
- H3N2 clade K
- H3N2 clade J.2.2

To simulate realistic clinical complexity, the in-silico design incorporated:

- Reconstruction of minority variant consensus sequences
- Controlled mixing of major and minor variants at defined proportions
- Simulation of human background reads
- Introduction of contamination (e.g., SARS-CoV-2 or rhinovirus reads in selected samples).
- Segment-specific coverage dropouts (e.g., HA or NA depletion).
- Platform-specific read simulation using [ART (Illumina) v2016.06.05](https://surveillance.cancer.gov/genetic-simulation-resources/packages/art/) and [Badread (Nanopore) v0.4.1](https://github.com/rrwick/Badread).

This approach allowed precise control over:

- Variant frequency structure
- Segment coverage distribution
- Contamination levels
- Platform-dependent error profiles

{% set table_counter.value = table_counter.value + 1 %}

_**Table {{ table_counter.value }}**. Viral, host and contaminant composition design of in-silico influenza datasets used for benchmarking._

| Sample | Influenza reads | Host reads | Additional Viral reads  | Total reads | Analytical Challenge                |
|--------|-----------------|------------|-------------------------|-------------|-------------------------------------|
| FLU2   | 1378764         | 462520     | 0                       | 1841284     | Baseline performance assessment     |
| FLU4   | 181626          | 300000     | 200000 SARS-CoV-2 reads | 681626      | Contamination with SARS-CoV-2       |
| FLU5   | 1088000         | 100000     | 0                       | 1188000     | NA segment dropout                  |
| FLU7   | 5677            | 100        | 255 Rhinovirus reads    | 6032        | Cross-virus contamination challenge |
| FLU8   | 5380            | 300        | 0                       | 5680        | Baseline performance assessment     |
| FLU9   | 19989           | 500        | 0                       | 20489       | HA segment dropout                  |

{% set table_counter.value = table_counter.value + 1 %}

_**Table {{ table_counter.value }}**. Influenza virus samples used in the RELECOV 2026 Dry-Lab EQA, including sequencing platform, enrichment strategy, primer scheme, and key analytical features._

| Sample | Source    | Platform | Enrichment Strategy | Primer Scheme                                   | Read Layout | Ref_sample        | Type   | Clade HA  | Legacy Clade       | Key Feature                             | Quality check |
|--------|-----------|----------|---------------------|-------------------------------------------------|-------------|-------------------|--------|-----------| ------------------ | ----------------------------------------|---------------|
| FLU1   | ESIB 2024 | Illumina | Amplicon            | CommonUni12/13 (Van den Hoecke 2015)            | Paired-end  | INFL2.07          | A/H5N1 | 2.3.4.4b  | -                  | High-quality baseline sample (zoonotic) | Ok            |
| FLU2   | In-silico | Illumina | Amplicon            | Zhou 2009 single-reaction genomic amplification | Paired-end  | In-silico Sample1 | A/H1N1 | D.3.1.1   | 6B.1A.5a.2a.1      | High-quality baseline sample (human)    | Ok            |
| FLU3   | ESIB 2024 | Illumina | No enrichment       | —                                               | Paired-end  | INFL2.04          | —      | —         | -                  | No influenza (Rhinovirus only)          | Bad           |
| FLU4   | In-silico | Illumina | Amplicon            | Zhou 2009 single-reaction genomic amplification | Paired-end  | In-silico Sample3 | A/H3N2 | K         | 3C.2a1b.2a.2a.3a.1 | Contamination with SARS-CoV-2           | Ok            |
| FLU5   | In-silico | Illumina | Amplicon            | Zhou 2009 single-reaction genomic amplification | Paired-end  | In-silico Sample4 | A/H3N2 (A/H3 or A/H3Nx) | J.2.2     | 3C.2a1b.2a.2a.3a.1 | NA segment dropout                      | Bad           |
| FLU6   | ESIB 2024 | Nanopore | No enrichment       | —                                               | Single-end  | INFL1.02          | A/H5N6 | 2.3.4.4h  | -                  | High-quality baseline sample (zoonotic) | Ok            |
| FLU7   | In-silico | Nanopore | Amplicon            | Zhou 2009 single-reaction genomic amplification | Single-end  | In-silico Sample2 | A/H1N1 | C.1.9.3   | 6B.1A.5a.2a        | Contamination with Rhinovirus           | Ok            |
| FLU8   | In-silico | Nanopore | Amplicon            | Zhou 2009 single-reaction genomic amplification | Single-end  | In-silico Sample3 | A/H3N2 | K         | 3C.2a1b.2a.2a.3a.1 | High-quality baseline sample (human)    | Ok            |
| FLU9   | In-silico | Nanopore | Amplicon            | Zhou 2009 single-reaction genomic amplification | Single-end  | In-silico Sample1 | A/H1N1 (A/N1 or A/HxN1) | D.3.1.1   | 6B.1A.5a.2a.1      | HA segment dropout                      | Bad           |
| FLU10  | ESIB 2024 | Nanopore | Amplicon            | CommonUni12/13 (Van den Hoecke 2015)            | Single-end  | INFL1.08          | A/H5N1 | 2.3.4.4b  | -                  | High-quality baseline sample (zoonotic) | Ok            |

## 4. Methodology of Evaluation

The evaluation framework was designed to ensure objective, reproducible, and comparable assessment of analytical performance across participating laboratories. Submitted outputs were benchmarked against curated gold standard datasets from ECDC ESIB or generated in-silico.

The evaluation was structured into five independent analytical domains:

- Submission Completeness
- Consensus genome reconstruction performance
- Variant detection accuracy
- Lineage, Subtype and Clade Assignment
- Metadata completeness and compliance
- Pipeline Benchmarking and Comparative Performance

Each domain was assessed using predefined quantitative metrics to allow cross-laboratory comparison and pipeline benchmarking. Participation metrics were calculated at both component and laboratory level.

All the scripts and templates used for evaluation and to generate reports and plots is publicly available [in github](https://github.com/BU-ISCIII/relecov_2026_drylab_eqa)

### 4.1. Submission Completeness

Submission completeness was evaluated to quantify the extent to which participating laboratories provided the expected analytical outputs for the components they chose to analyse.

This assessment focused exclusively on:

- Number of components analysed per laboratory
- Number of consensus genome files (.fasta) submitted
- Number of variant call files (.vcf) submitted

Laboratories were free to analyse any subset of the four available components (SARS1, SARS2, FLU1, FLU2).

Network-level participation was summarised using:

- Total number of participating laboratories
- Number of laboratories per component
- Median number of components analysed per laboratory

For each analysed component, laboratories were expected to submit:

- One consensus genome file (.fasta) per sample
- One variant call file (.vcf) per sample

Submission completeness was calculated as:

$$
\text{File submission rate} =
\frac{\text{Number of submitted files}}
{\text{Total number of expected files}}
$$

Missing files were recorded but not penalised beyond descriptive reporting, as laboratories were allowed to participate selectively according to local analytical capacity.

### 4.2. Evaluation of Consensus Genome Reconstruction Performance

For each sample, a curated gold standard consensus genome was provided by the ECDC or generated internally. These reference sequences served as the ground truth for comparative analysis.

All submitted consensus sequences (.fasta) were:

- Aligned against the corresponding gold standard sequence using [Mafft v7.475](https://mafft.cbrc.jp/alignment/software/).
- Compared position-by-position relative to the declared reference genome coordinate system.

Differences between submitted sequences and gold standard sequences were categorised into the following classes:

- **Wrong nucleotide**: A nucleotide different from the allowed reference or ambiguity code.
- **Ambiguity instead of nucleotide**: Ambiguity codes introduced where a defined nucleotide was expected.
- **Nucleotide instead of ambiguity**: Defined nucleotide provided where an ambiguity code was expected.
- **Stretch of Ns instead of nucleotide stretch**: Continuous region of Ns where defined bases were expected.
- **Nucleotide stretch instead of stretch of Ns**: Defined bases provided where Ns were expected.
- **Insertion relative to gold standard**
- **Deletion relative to gold standard**

Each insertion, deletion, or contiguous stretch of Ns was counted as a single event.

For each laboratory and sample, the following summary metrics were compiled:

- Total number of nucleotide discrepancies
- Percentage genome identity relative to the curated gold standard reference sequence

The proportional contribution of each discrepancy category was calculated relative to the total number of discrepancies observed per component.

### 4.3. Evaluation of Variant Detection Accuracy

For influenza virus datasets, direct position-by-position comparison of reported variants against the curated reference variant set was not feasible under the same framework applied to SARS-CoV-2.

Unlike SARS-CoV-2, where laboratories predominantly use a shared and globally standardised reference genomes (either MN908947.3 or NC_045512.2), influenza virus analyses exhibited substantial heterogeneity in reference genome selection. As a result:

- Variant coordinates were reported relative to different reference accessions.
- Segment boundaries and numbering schemes varied.
- Insertions and deletions were represented inconsistently across reference backbones.

This heterogeneity prevented robust coordinate harmonisation across submissions without introducing alignment-dependent artifacts and interpretation bias.

#### 4.3.1. SARS-CoV-2

A curated reference variant set was generated for each SARS-CoV-2 sample. Variant positions were standardized relative to a defined coordinate system referred to the references used by Nextclade.

Submitted .vcf files were:

- Converted to a standardised long table format for coordinate comparison
- Compared position-by-position with the reference variant set

Differences between submitted variants and reference variant set were categorised into the following classes:

- **Wrong nucleotide**: A nucleotide different from the allowed reference or ambiguity code.
- **Insertion relative to gold standard**
- **Deletion relative to gold standard**
- **Missing variant**: Variants present in the reference but missing in the sample.
- **De novo**: Variants present in the sample but missing in the reference set.

Each insertion, deletion, or contiguous stretch of Ns was counted as a single event.

For each laboratory and sample, the total number of nucleotide discrepancies was calculated.

The proportional contribution of each discrepancy category was calculated relative to the total number of discrepancies observed per component.

Metadata describing the following analytical settings were collected to support result interpretation:

- Allele frequency thresholds
- Minimum coverage thresholds
- Reference genome selection

#### 4.3.2. Descriptive and Structural Variant Reporting Metrics

In addition to nucleotide-level discrepancy analysis for SARS-CoV-2, both SARS-CoV-2 and influenza submissions were evaluated using descriptive and structural reporting metrics to characterise reporting behaviour and methodological heterogeneity across laboratories.

For both viruses, the following reporting practice metrics were collected:

- Number of laboratories reporting high-frequency variants only.
- Number of laboratories reporting both high- and low-frequency variants.
- Number of laboratories reporting exclusively low-frequency variants.
- Total number of distinct reference genomes employed for variant calling or mapping.

For influenza virus, additional structural summary metrics were calculated because direct coordinate-harmonised comparison of all submitted variants was not methodologically robust across segment-specific references:

- Number of variants with an allele frequency higher than 75%.
- Number of variants with an allele frequency higher than 75% derived from VCF files.
- Total number of variants present in the submitted VCF.
- Discrepancies between variants with an allele frequency higher than 75% reported in the metadata and in the VCF file.

These metrics provide insight into:

- Variant reporting practices across laboratories.
- Heterogeneity in allele frequency thresholds.
- Diversity of reference genome usage.
- Internal consistency between consensus outputs and submitted VCF files.
- Degree of methodological standardisation within the network.

This evaluation approach allows characterisation of variant reporting behaviour while acknowledging the need for harmonisation in inherently reference-dependent analyses, particularly for segmented influenza genomes.

### 4.4. Evaluation of Lineage, Subtype and Clade Assignment

Classification outputs were evaluated separately according to virus type.

#### SARS-CoV-2

For each SARS-CoV-2 sample:

- Lineage assignment was compared to the gold standard lineage designation from the ECDC in 2024.
- Clade assignment was compared to the gold standard clade classification from the ECDC in 2024.

#### Influenza virus

For influenza samples, evaluation included:

- Virus type and subtype identification (e.g., Influenza A/HxNy) compared to the gold standard subtype of the ECDC or the reference genome’s subtype used to generate in-silico reads.
- Clade assignment compared to the gold standard subtype of the ECDC or the reference genome’s subtype used to generate in-silico reads.

#### Both viruses

For SARS-CoV-2 and Influenza viruses, concordance was assessed as:

- **Match**, when lineage/subtype or clade was correct.
- **Discrepancy**, when one of the classifications was incorrect.

If a laboratory did not report a lineage/subtype or clade assignment, that missing classification was also counted as a **Discrepancy** for evaluation purposes. Although these metadata fields were not mandatory in the submission template, participating RELECOV laboratories are expected to be able to determine both classification dimensions for analysed samples.

Potential contributors considered during result interpretation included:

- Differences in database versioning
- Differences in software versioning
- Reporting practices and field completion
- The possible relationship between consensus discrepancies and lineage/type assignment performance

Failure to identify virus presence in positive samples, or misclassification of negative samples, was recorded separately.

### 4.5. Evaluation of Metadata Completeness and Compliance

Metadata assessment focused on analytical transparency and interoperability rather than biological correctness. Before the start of the exercise, a metadata template with controlled-vocabulary dropdowns was distributed among the laboratories to review the available options and incorporation of missing software tools into the schema.

For each submitted sample, metadata completeness was calculated as:

$$
\text{Sample metadata completeness} =
\frac{\text{Number of correctly populated fields}}
{\text{Total number of applicable fields}}
$$

Laboratory-level and component-level completeness summaries were then derived from these sample-level values. Fields were evaluated for:

- Completion: Each sample has a list of minimum **recommended** fields, based on the sample characteristics. For each component/sample/lab the total number of completed minimum **recommended** fields was evaluated. Both mandatory and optional analytical fields were included in the completeness assessment, while fields not applicable to a laboratory’s selected components were excluded from scoring.
- Compliance with controlled vocabularies. Metadata entries were considered non-compliant when:
    - Controlled vocabulary options were bypassed
    - Free-text substitutions replaced defined values
- Valid file path reporting

This evaluation allowed quantification of metadata standardisation and reproducibility readiness across the network.

#### Evaluation of Sample Quality Control Assessment

The evaluation of sample quality control (QC) assessment was designed to determine whether participating laboratories correctly interpreted overall analytical quality status for each sample. Participating laboratories were required to report their own QC evaluation for each analysed sample within the metadata template. For every sample included in the exercise, a gold standard quality control classification was predefined based on the original ECDC dataset evaluation or the in-silico design specifications. Each sample was categorised as:

- Pass
- Fail

For each laboratory and sample, the reported QC classification was compared to the predefined gold standard QC status. Results were categorised as:

- **Match**: Laboratory-reported QC status identical to the gold standard classification.
- **Discrepancy**: Laboratory-reported QC status different from the gold standard classification.

For each laboratory, component, and the overall network, the following metrics were calculated:

- Total number of QC evaluations performed
- Number of Matches
- QC concordance rate, where:

$$
\text{QC concordance rate} =
\frac{\text{Number of Matches}}
{\text{Total QC evaluations}}
$$

QC evaluations were calculated only for samples analysed by the laboratory.

The QC assessment evaluation was limited to concordance analysis. The exercise did not attempt to infer the internal QC criteria applied by laboratories, but rather assessed agreement with the predefined gold standard QC status to evaluate interpretative consistency across the network.

### 4.6. Pipeline Benchmarking and Comparative Performance

The pipeline benchmarking analysis was designed to evaluate analytical performance at the pipeline and software level, rather than solely at the individual laboratory level. The objective was to identify which analytical workflows most consistently generate results that closely match the curated gold standard datasets.

For each declared pipeline or analytical workflow (including software combinations and parameter configurations), performance was aggregated across all laboratories using that approach.

The primary benchmarking criterion was based on these performance indicators:

- Median consensus genome identity relative to the curated gold standard.
- Median number of discrepancies relative to the curated gold standard.
- Exact lineage/type and clade classification concordance.
- Median metadata completeness

These metrics were analysed to determine whether pipelines achieving high consensus similarity also demonstrated consistent downstream analytical accuracy.

Benchmarking results were interpreted to identify:

- Pipelines demonstrating consistently low divergence from gold standards
- Parameter configurations associated with systematic discrepancies
- The impact of software versioning and reference genome selection

The benchmarking framework therefore provides an empirical basis for:

- Identifying best-performing analytical workflows
- Defining minimum performance criteria for network harmonisation
- Informing recommendations for standardisation within the RELECOV analytical platform

## 5. General Results

A total of 52 laboratories within the RELECOV network were invited to participate. Of these, {{ general.total_participants }} laboratories {{ pct(general.total_participants_pct) }} submitted results for one or more components with the following distribution:

- SARS1 (SARS-CoV-2, Illumina): {{ general.participation_per_component.SARS1 }} laboratories.
- SARS2 (SARS-CoV-2, Oxford Nanopore Technologies): {{ general.participation_per_component.SARS2 }} laboratories.
- FLU1 (Influenza virus, Illumina): {{ general.participation_per_component.FLU1 }} laboratories.
- FLU2 (Influenza virus, Oxford Nanopore Technologies): {{ general.participation_per_component.FLU2 }} laboratories.

The median number of components analysed per participating laboratory was {{ "%.0f"|format(general.median_components_analysed_per_lab) }}.

The results presented below are interpreted according to the evaluation framework described in [Section 4](#4-methodology-of-evaluation).

### 5.1. Submission Completeness

Across all components:

- {{ pct(general.submission_rates_pct.fasta) }} of laboratories submitted consensus genome files (.fasta), where applicable.
- {{ pct(general.submission_rates_pct.vcf) }} submitted variant call files (.vcf), where applicable.

Component-level submission totals are presented in Section 6 and reflect both the number of participating laboratories and the expected output files for each dataset.

### 5.2. Consensus Genome Reconstruction Performance

Across the two Illumina-based components, the combined median genome identity was {{ pct(general.general_results.consensus.median_identity_illumina_pct, 2) }}, compared with {{ pct(general.general_results.consensus.median_identity_nanopore_pct, 2) }} across the two Nanopore-based components. When grouped by virus, the combined median genome identity was {{ pct(general.general_results.consensus.median_identity_sars_pct, 2) }} across the SARS-CoV-2 components and {{ pct(general.general_results.consensus.median_identity_influenza_pct, 2) }} across the influenza components. Notably, the influenza median was also slightly lower than the combined Nanopore-based median, indicating that virus-specific analytical complexity likely contributed in addition to platform-related differences. Nanopore-based datasets also showed broader overall identity ranges, where low-identity outliers were present.

Dominant discrepancy patterns differed by component:

- In SARS1, the most frequent discrepancy category was stretches of Ns in the submitted consensus where defined nucleotides were present in the gold standard.
- In SARS2, the most frequent discrepancy category was defined nucleotides in the submitted consensus where stretches of Ns were present in the gold standard.
- FLU1 and FLU2 were both dominated by deletions relative to the gold standard.


Across components, many discrepancy categories had medians of zero, indicating that errors tended to be concentrated in a smaller number of laboratories or samples rather than being uniformly distributed across the network.

{% set fig_counter.value = fig_counter.value + 1 %}

{% set figure_cfg.style = "max-width: 98%;" %}
{{ render_figure(general.figures.consensus_summary, "Network-level consensus reconstruction performance summary.", has_panels=True ) }}

**_Figure {{ fig_counter.value }}_. Consensus genome reconstruction performance across components**. Panel **A** shows the distribution of nucleotide discrepancies relative to the gold standard across components, and panel **B** shows the corresponding distribution of genome identity values. In both panels, the central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

### 5.3. Variant Detection Accuracy

#### 5.3.1. SARS-CoV-2

For SARS-CoV-2 components (SARS1 and SARS2), variant detection accuracy was assessed against curated reference variant sets. Overall, submitted VCFs showed a median number of discrepancies of {{ "%.0f"|format(general.general_results.sars_variants.median_discrepancy_illumina) }} for Illumina component and a median number of {{ "%.0f"|format(general.general_results.sars_variants.median_discrepancy_nanopore) }} for Nanopore component, discrepancies relative to the reference variant set.

{% set fig_counter.value = fig_counter.value + 1 %}
Variant detection performance differed across components (Figure {{ fig_counter.value }}). Contextual factors documented in the metadata that may contribute to these differences included:

- Allele frequency thresholds used for incorporation into vcf files
- Variant normalization practices (variant caller software and params)


{% set figure_cfg.style = "max-width: 70%;" %}
{{ render_figure(general.figures.variant_summary, "Network-level variant detection performance summary." ) }}
**_Figure {{ fig_counter.value }}_. SARS-CoV-2 network-level variant detection performance summary**. Boxplots represent the number of variant discrepancies per SARS-CoV-2 component across participating laboratories. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

Variant evaluation included structural reporting characteristics and methodological heterogeneity. At network level:

- {{ general.general_results.sars_variants.high_and_low_freq_pct }}% of laboratories reported both high- and low-frequency variants.
- {{ general.general_results.sars_variants.low_freq_only_pct }}% reported exclusively low-frequency variants.
- {{ general.general_results.sars_variants.high_freq_only_pct }}% reported only high-frequency variants.

Additionally, a total of {{ general.general_results.sars_variants.total_distinct_references }} distinct reference genomes were employed for variant calling across SARS-CoV-2 components ({{ general.general_results.sars_variants.distinct_references | join(", ") }}).

{% set fig_counter.value = fig_counter.value + 1 %}

{% set figure_cfg.style = "max-width: 70%;" %}
{{ render_figure(
general.figures.sars_variant_reporting_summary,
"SARS-CoV-2 variant reporting practices across the network."
) }}

**_Figure {{ fig_counter.value }}_. SARS-CoV-2 variant reporting characteristics across the network**. Summarise the proportion of laboratories reporting high- and/or low-frequency variants.

#### 5.3.2. Influenza virus

For influenza virus components (FLU1 and FLU2), variant evaluation focused on structural reporting characteristics and methodological heterogeneity.

At network level:

- {{ general.general_results.influenza_variants.high_and_low_freq_pct }}% of laboratories reported both high- and low-frequency variants.
- {{ general.general_results.influenza_variants.low_freq_only_pct }}% reported exclusively low-frequency variants.
- {{ general.general_results.influenza_variants.high_freq_only_pct }}% reported only high-frequency variants.

Additionally, an estimated total of {{ "%.0f"|format(general.general_results.influenza_variants.total_distinct_references) }} distinct reference genomes were employed for variant calling or mapping across influenza components. This value was rounded to the nearest whole genome by dividing the total number of distinct fragment references ({{ general.general_results.influenza_variants.total_distinct_fragments }}) by 8 influenza genome segments.

Structural summary metrics derived from submitted influenza consensus sequences and VCF files are presented in Table {{ table_counter.value + 1 }}. These metrics capture the overall magnitude of reported variants in the metadata file and the discrepancy between reported variants with an allele frequency >= 75% in the metadata file and the VCF file, rather than direct nucleotide-level accuracy against a unified reference coordinate system.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level structural summary of influenza variant reporting.**

| Metric | Network median | Min-max |
|---|---:|---:|
| Variants with AF>=75% | {{ "%.0f"|format(general.general_results.influenza_variants.median_variants_in_consensus) }} | {{ "%.0f"|format(general.general_results.influenza_variants.min_variants_in_consensus) }}–{{ "%.0f"|format(general.general_results.influenza_variants.max_variants_in_consensus) }} |
| Variants with AF>=75% in VCF | {{ "%.0f"|format(general.general_results.influenza_variants.median_variants_in_consensus_vcf) }} | {{ "%.0f"|format(general.general_results.influenza_variants.min_variants_in_consensus_vcf) }}–{{ "%.0f"|format(general.general_results.influenza_variants.max_variants_in_consensus_vcf) }} |
| Discrepancies in reported variants | {{ "%.0f"|format(general.general_results.influenza_variants.median_discrepancies_in_reported_variants) }} | {{ "%.0f"|format(general.general_results.influenza_variants.min_discrepancies_in_reported_variants) }}–{{ "%.0f"|format(general.general_results.influenza_variants.max_discrepancies_in_reported_variants) }} |
| Total variants in VCF | {{ "%.0f"|format(general.general_results.influenza_variants.median_variants_in_vcf) }} | {{ "%.0f"|format(general.general_results.influenza_variants.min_variants_in_vcf) }}–{{ "%.0f"|format(general.general_results.influenza_variants.max_variants_in_vcf) }} |

{% set fig_counter.value = fig_counter.value + 1 %}

{% set figure_cfg.style = "max-width: 96%;" %}
{{ render_figure(
general.figures.influenza_variant_reporting_summary,
"Influenza variant reporting practices across the network."
) }}

**_Figure {{ fig_counter.value }}_. Influenza variant reporting characteristics across the network**. Summarise the proportion of laboratories reporting high- and/or low-frequency variants.

Together, these results show marked heterogeneity in influenza variant reporting within the network (Table {{ table_counter.value }}, Figure {{ fig_counter.value }}).

### 5.4. Lineage, Subtype and Clade Assignment

Overall concordance rates were:

- SARS-CoV-2 lineage assignment: **{{ pct(general.general_results.classification.sars_lineage_concordance_pct) }}** concordance.
- Influenza type/subtype identification: **{{ pct(general.general_results.classification.influenza_type_concordance_pct) }}** concordance.
- SARS-CoV-2 clade assignment: **{{ pct(general.general_results.classification.sars_clade_concordance_pct) }}** concordance.
- Influenza clade assignment: **{{ pct(general.general_results.classification.influenza_clade_concordance_pct) }}** concordance.

Across components, lineage/type concordance was consistently higher than clade concordance. SARS-CoV-2 lineage assignment reached {{ pct(general.general_results.classification.sars_lineage_concordance_pct) }}, compared with {{ pct(general.general_results.classification.sars_clade_concordance_pct) }} for SARS-CoV-2 clade assignment, while influenza type/subtype identification reached {{ pct(general.general_results.classification.influenza_type_concordance_pct) }} compared with {{ pct(general.general_results.classification.influenza_clade_concordance_pct) }} for influenza clade assignment.

{% set fig_counter.value = fig_counter.value + 1 %}

{% set figure_cfg.style = "max-width: 98%;" %}
{{ render_figure(general.figures.classification_summary, "Distribution of classification outcomes across participating laboratories.", has_panels=True ) }}

**_Figure {{ fig_counter.value }}_. Distribution of classification outcomes across participating laboratories.** Panel **A** shows **lineage/type assignments**, and panel **B** shows **clade assignments**. Stacked bars represent the percentage of all possible sample-level classifications across participating laboratories for each component. Bars are partitioned into **Match** (correct assignments relative to the curated gold standard), **Discrepancy** (incorrect assignments), and **Not provided** (classification not reported).

### 5.5. Metadata completeness and compliance

The evaluation of metadata focused on analytical transparency, reproducibility, and interoperability within the RELECOV network, including controlled vocabulary adherence, logical consistency, and reporting of analytical parameters.

#### Overall Completeness

{% set fig_counter.value = fig_counter.value + 1 %}

Across all participating laboratories, the metadata template was completed at a median completeness rate of {{ pct(general.metadata_completeness.median_pct) }}, with values ranging from {{ pct(general.metadata_completeness.min_pct) }} to {{ pct(general.metadata_completeness.max_pct) }}. Component-level median completeness values were similar overall, but the observed ranges remained broad in all components (Figure {{ fig_counter.value }}). The leading incompleteness drivers were variant calling, pre-processing, and mapping fields, followed by QC metrics, de-hosting, and consensus analysis fields.

Most frequent incompleteness drivers across the network:
<ul class="compact-list">
{% for item in general.metadata_completeness.top_5_primary_incompleteness_driver_counts %}
<li>{{ item.driver }} (missing in {{ item.count }} laboratories)</li>
{% endfor %}
</ul>



{% set figure_cfg.style = "max-width: 80%;" %}
{{ render_figure(general.figures.metadata_completeness_distribution,
  "Distribution of metadata completeness across participating laboratories.") }}

**_Figure {{ fig_counter.value }}_. Distribution of metadata completeness across participating laboratories**. Boxplots represent the distribution of sample-level metadata completeness percentages across the different components. Completeness was calculated for each submitted sample as the proportion of filled metadata fields relative to the total number of maximum expected metadata fields. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

#### Reporting of Analytical Parameters

Although core pipeline tools were generally reported, variability was observed in the level of parameter detail provided.

- {{ pct(general.metadata_completeness.software_names_pct) }} of the maximum software-name fields were completed across submitted samples.
- {{ pct(general.metadata_completeness.software_version_pct) }} of the maximum software-version fields were completed across submitted samples.
- {{ pct(general.metadata_completeness.coverage_threshold_pct) }} specified minimum coverage thresholds.
- {{ pct(general.metadata_completeness.variant_calling_params_pct) }} reported variant calling parameters, containing the potential allele frequency thresholds.
- {{ pct(general.metadata_completeness.reference_genome_pct) }} reported the reference genome accession or identifier.

The incomplete reporting of parameters limited the ability to fully reconstruct or reproduce analytical workflows in {{ pct(general.metadata_completeness.incomplete_parameters_pct) }} of submissions.

#### Controlled Vocabulary Compliance

All (100%) laboratories required either clarification through e-mail contact or metadata correction during validation, as reflected by the high proportion of submissions with incomplete parameters or controlled-vocabulary corrections. Compliance with controlled vocabulary requirements was assessed to determine the degree of metadata standardisation achieved across participating laboratories. Only fields expected to contain predefined categorical values were considered in this analysis; free-text fields such as software versions, file names, and parameter descriptions were excluded.

- 26.32% of submissions (5 laboratories) were fully compliant with controlled vocabulary requirements.
- 73.68% (14 laboratories) required at least one manual correction due to the use of a non-standard value in a controlled field.

The most common compliance issues included:

- Use of free-text entries instead of predefined software names in dropdown-based metadata fields, which prevented direct validation against the harmonised template.
- Incorrect or inconsistent completion of lineage, clade, influenza type, or subtype fields.
- Missing mandatory fields requiring subsequent normalisation.

#### Sample Quality Control Assessment

Sample quality control (QC) classifications reported by laboratories (Pass/Fail) were compared against the predefined gold standard QC status for each sample (ECDC or in-silico). QC agreement was evaluated as a binary outcome:

- Match: laboratory QC classification equals the gold standard QC status
- Discrepancy: laboratory QC classification differs from the gold standard QC status

Overall, the network achieved {{ pct(general.qc.reported_match_rate_pct) }} QC concordance, corresponding to {{ general.qc.matches }} Matches and {{ general.qc.discrepancies }} Discrepancies across {{ general.qc.total_evaluations }} evaluated sample-level QC decisions.

{% set fig_counter.value = fig_counter.value + 1 %}
QC concordance differed across components, ranging from {{ pct(general.components.SARS1.qc.reported_match_rate_pct) }} in SARS1 to {{ pct(general.components.FLU2.qc.reported_match_rate_pct) }} in FLU2, based on reported QC information (Figure {{ fig_counter.value }}).


{% set figure_cfg.style = "max-width: 80%;" %}
{{ render_figure(
general.figures.qc_match_rate_by_component,
"QC concordance by component (Match, Discrepancy, and Not provided relative to the gold standard)."
) }}

**_Figure {{ fig_counter.value }}_. QC concordance by component relative to the gold standard.** Stacked bars represent the proportion of sample-level QC outcomes classified as Match, Discrepancy, or Not provided for each component across participating laboratories. Not provided values correspond to missing QC assessments and are shown separately from true discrepancies.

### 5.6. Pipeline Benchmarking and Comparative Performance

The benchmarking analysis was designed to assess whether differences in analytical software and parameterisation were associated with measurable variability in performance across participating laboratories.

The submitted metadata documented heterogeneity in:

- Choice of consensus reconstruction software
- Variant calling strategies
- Lineage and clade assignment tool's version and database versions.
- Reference genome selection
- Coverage and allele frequency thresholds

#### Diversity of Analytical Workflows

The metadata submissions allowed characterisation of the analytical landscape currently implemented across the RELECOV network.

A total of {{ general.metadata_completeness.total_workflows }} distinct analytical workflows were identified across participating laboratories, defined as unique combinations of software tools and versions declared in the metadata template.

Substantial diversity was observed in the selection of core analytical tools, based on distinct declared software identities in the submitted metadata (software name plus version where applicable):

- Consensus reconstruction software ( {{ general.metadata_completeness.total_consensus_softwares }} distinct declared software identities )
- Variant calling tools ( {{ general.metadata_completeness.total_variant_softwares }} distinct declared software identities )
- SARS-CoV-2 lineage assignment software ( {{ general.metadata_completeness.total_lineage_assignment_softwares }} distinct declared software identities )
- Clade assignment software ( {{ general.metadata_completeness.total_clade_assignment_softwares }} distinct declared software identities )
- Influenza type assignment software ( {{ general.metadata_completeness.total_type_assignment_softwares }} distinct declared software identities )
- Influenza subtype assignment software ( {{ general.metadata_completeness.total_subtype_assignment_softwares }} distinct declared software identities )

For lineage, clade, type, and subtype benchmarking in Section 6, these declared software identities are further stratified by database version when that information was reported, so the benchmarking categories may be more granular than the metadata diversity counts summarised here.

Comparative performance analyses stratified by component are presented in Section 6, where software-level differences are evaluated within homogeneous analytical contexts (SARS-CoV-2 Illumina, SARS-CoV-2 Nanopore, Influenza Illumina, Influenza Nanopore).

Because performance differed by component and by metric, software-level comparisons are presented in Section 6 within component-specific contexts rather than collapsed into a single cross-component ranking.

This diversity shows that multiple analytical configurations are currently in use across the RELECOV network. These findings highlight the importance of harmonising minimum analytical criteria while preserving methodological flexibility within the network.

## 6. Component-specific Results

This section presents the analytical results disaggregated by component, allowing a detailed assessment of performance within each dataset and sequencing technology. For each component, results are structured according to participation and submission metrics, consensus genome reconstruction performance, variant detection accuracy, and Lineage, Subtype or clade assignment concordance, as applicable.

Component-level analyses enable identification of platform-specific patterns, dataset-dependent challenges, and variability associated with particular sample characteristics. This approach facilitates a more granular interpretation of performance differences observed at the network level and supports targeted harmonisation recommendations.

All component-level results below are reported using the same evaluation framework described in [Section 4](#4-methodology-of-evaluation).

{% for comp_code, comp_net in general.components.items() %}

### 6.{{ loop.index }}. {{ comp_code }} ({{ comp_net.name }})

#### 6.{{ loop.index }}.1. Participation and Submissions

A total of {{ comp_net.total_labs }} laboratories submitted results for the {{ comp_code }} component:

- A total of {{ comp_net.total_fasta }} consensus genome sequences (.fasta) were submitted.
- A total of {{ comp_net.total_vcf }} variant call files (.vcf) were submitted.
- The metadata template completeness for {{ comp_code }} submissions had a median of {{ pct(comp_net.metadata_completeness_median) }}.
{% if comp_net.top_5_primary_incompleteness_driver_counts %}

Most frequent incompleteness drivers in {{ comp_code }}:
<ul class="compact-list">
{% for item in comp_net.top_5_primary_incompleteness_driver_counts %}
<li>{{ item.driver }} (missing in {{ item.count }} laboratories)</li>
{% endfor %}
</ul>
{% endif %}

#### 6.{{ loop.index }}.2. Consensus Genome Reconstruction Performance

Consensus sequences were evaluated against the corresponding curated gold standard references for each sample in the {{ comp_code }} component.

Overall, {{ comp_code }} showed a median genome identity of {{ pct(comp_net.consensus.median_identity, 2) }}, with a median of {{ comp_net.consensus.median_discrepancies }} nucleotide discrepancies per sample (range: {{ comp_net.consensus.min_discrepancies }}–{{ comp_net.consensus.max_discrepancies }}) (Figure {{ fig_counter.value + 1 }}).
{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_consensus_metrics_table_num = appendix_table_counter.value %}

Figure {{ fig_counter.value + 1 }} presents the distribution of nucleotide discrepancies per sample across participating laboratories for {{ comp_code }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 90%;" %}
{{ render_figure(
  comp_net.consensus.fig_discrepancies_boxplot_by_sample,
  "Consensus discrepancies per sample for " ~ comp_code ~ " relative to the curated gold standard.",
  has_panels=False
) }}


**Figure {{ fig_counter.value }}. Consensus reconstruction performance by sample for {{ comp_code }}.** Panel A shows the distribution of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample, and Panel B shows the corresponding distribution of genome identity values. In both panels, the central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. In Panel B, the y-axis is truncated to highlight differences among high-identity values.

Figure {{ fig_counter.value + 1 }} presents the distribution of nucleotide discrepancy types per sample across participating laboratories for {{ comp_code }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 80%;" %}
{{ render_figure(
  comp_net.consensus.fig_discrepancies_stacked_by_sample,
  "Consensus discrepancy types per sample for " ~ comp_code ~ " relative to the curated gold standard."
) }}


**Figure {{ fig_counter.value }}. Consensus discrepancy type composition per sample for {{ comp_code }}.** Stacked bars represent the number and type of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample.
{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_consensus_sample_table_num = appendix_table_counter.value %}
{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_consensus_type_table_num = appendix_table_counter.value %}
{% set appendix_fig_counter.value = appendix_fig_counter.value + 1 %}
{% set appendix_consensus_type_fig_num = appendix_fig_counter.value %}
{% set _ = consensus_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "metrics_table_num": appendix_consensus_metrics_table_num,
  "sample_table_num": appendix_consensus_sample_table_num,
  "type_table_num": appendix_consensus_type_table_num,
  "type_fig_num": appendix_consensus_type_fig_num
}) %}

The dominant discrepancy pattern observed in {{ comp_code }} was {{ discrepancy_label(comp_net.consensus.dominant_discrepancy_pattern) }} (Figure {{ fig_counter.value }}). Sample-level consensus reconstruction summary metrics are provided in Appendix Table {{ appendix_consensus_metrics_table_num }}. A full sample-level breakdown of discrepancy categories is provided in Appendix Table {{ appendix_consensus_sample_table_num }}, while the aggregated discrepancy composition by type and the corresponding category-wise boxplot can be found in Appendix Table {{ appendix_consensus_type_table_num }} and Appendix Figure {{ appendix_consensus_type_fig_num }}, respectively.

#### 6.{{ loop.index }}.3. Variant Detection Accuracy

{% if comp_net.variant.median_discrepancies is defined %}

Variant call files (.vcf) submitted for the {{ comp_code }} component were compared against the curated reference variant set corresponding to each sample in the {{ comp_code }} component.

Overall, {{ comp_code }} showed a median of {{ comp_net.variant.median_discrepancies }} variant discrepancies per sample (range: {{ comp_net.variant.min_discrepancies }}–{{ comp_net.variant.max_discrepancies }}), together with a median of {{ comp_net.variant.median_successful_hits if comp_net.variant.median_successful_hits is not none else "NA" }} successful hits per sample (Table {{ table_counter.value + 1 }}, Figure {{ fig_counter.value + 1 }}).

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 80%;" %}
{{ render_figure(
  comp_net.variant.fig_discrepancies_stacked_by_sample,
  "Variant discrepancies per sample for " ~ comp_code ~ " relative to the curated gold standard."
) }}

**Figure {{ fig_counter.value }}. Distribution of variant discrepancies per sample for {{ comp_code }}.** Stacked bars represent the number of nucleotide discrepancies and discrepancy types relative to the curated gold standard across participating laboratories for each sample.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level SARS-CoV-2 variant reporting metrics per sample for {{ comp_code }}.**

| Sample ID | Expected hits | Median successful hits (n={{ comp_net.variant_metadata_reporting.successful_hits_reported_n_labs }}) | Median variants >=75% AF in metadata (n={{ comp_net.variant_metadata_reporting.reported_n_labs }}) | Median variants >=75% AF in VCF (n={{ comp_net.variant_metadata_reporting.variants_in_consensus_vcf_reported_n_labs }}) | Median variants with effect in metadata (n={{ comp_net.variant_metadata_reporting.variants_with_effect_reported_n_labs }}) | Median variants with effect in VCF (n={{ comp_net.variant_metadata_reporting.variants_with_effect_vcf_reported_n_labs }}) | Median discrepancies metadata vs VCF | Median effect discrepancies metadata vs VCF |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
{% for s in comp_net.variant.samples %}
| {{ s.collecting_lab_sample_id }} | {{ s.expected_hits if s.expected_hits is not none else "NA" }} | {{ s.median_successful_hits if s.median_successful_hits is not none else "NA" }} | {{ s.variants_in_consensus.median if s.variants_in_consensus and s.variants_in_consensus.median is not none else "NA" }} | {{ s.variants_in_consensus_vcf.median if s.variants_in_consensus_vcf and s.variants_in_consensus_vcf.median is not none else "NA" }} | {{ s.variants_with_effect.median if s.variants_with_effect and s.variants_with_effect.median is not none else "NA" }} | {{ s.variants_with_effect_vcf.median if s.variants_with_effect_vcf and s.variants_with_effect_vcf.median is not none else "NA" }} | {{ s.discrepancies_in_reported_variants.median if s.discrepancies_in_reported_variants and s.discrepancies_in_reported_variants.median is not none else "NA" }} | {{ s.discrepancies_in_reported_variants_effect.median if s.discrepancies_in_reported_variants_effect and s.discrepancies_in_reported_variants_effect.median is not none else "NA" }} |
{% endfor %}

{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_variant_profile_table_num = appendix_table_counter.value %}
{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_variant_type_table_num = appendix_table_counter.value %}
{% set appendix_fig_counter.value = appendix_fig_counter.value + 1 %}
{% set appendix_variant_type_fig_num = appendix_fig_counter.value %}
{% set _ = variant_sars_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "profile_table_num": appendix_variant_profile_table_num,
  "type_table_num": appendix_variant_type_table_num,
  "type_fig_num": appendix_variant_type_fig_num
}) %}

At component level, {{ comp_net.variant_metadata_reporting.reported_n_labs }} laboratories reported the number of variants in the metadata, whereas {{ comp_net.variant_metadata_reporting.not_reported_n_labs }} did not report this field for any sample in {{ comp_code }}. The median number of variants with an allele frequency (AF) >=75% was {{ comp_net.variant.median_variants_in_consensus if comp_net.variant.median_variants_in_consensus is not none else "NA" }} in the metadata and {{ comp_net.variant.median_variants_in_consensus_vcf if comp_net.variant.median_variants_in_consensus_vcf is not none else "NA" }} in the submitted VCF files (Table {{ table_counter.value }}, Figure {{ fig_counter.value + 1 }}).

{% set fig_counter.value = fig_counter.value + 1 %}
Figure {{ fig_counter.value }} summarises the distribution of declared variant reporting modes across submitted sample outputs in {{ comp_code }}.

{% set figure_cfg.style = "max-width: 70%;" %}
{{ render_figure(
  comp_net.variant.fig_reporting_mode_by_component,
  "Variant reporting practices for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Variant reporting practices for {{ comp_code }}.** Bars represent the proportion of submitted sample outputs classified as high and low frequency reporting, high frequency only, or low frequency only, according to the metadata declarations associated with the variant outputs for this component.

The dominant discrepancy pattern observed in {{ comp_code }} was {{ discrepancy_label(comp_net.variant.dominant_discrepancy_pattern) }}. These patterns should be interpreted in the context of threshold choices, reference-genome selection, and declared pipeline configurations, all of which can shift the balance between successful hits, missing expected variants, and de novo calls even when laboratories analyse the same raw data. The full sample-level variant calling profile is provided in Appendix Table {{ appendix_variant_profile_table_num }}, while the aggregated discrepancy composition by type and the corresponding category-wise boxplot can be found in Appendix Table {{ appendix_variant_type_table_num }} and Appendix Figure {{ appendix_variant_type_fig_num }}, respectively.

{% else %}

For the {{ comp_code }} component, variant evaluation focused on the agreement between variants with allele frequency above 75% reported in the metadata template and those represented in the submitted VCF files, together with the overall number of variants present in the VCF output. At component level, {{ comp_net.variant_metadata_reporting.reported_n_labs }} laboratories reported the number of variants in the metadata, whereas {{ comp_net.variant_metadata_reporting.not_reported_n_labs }} did not report this field for any sample in {{ comp_code }}.

{% set fig_counter.value = fig_counter.value + 1 %}
Figure {{ fig_counter.value }} summarises the distribution of declared variant reporting modes across submitted sample outputs in {{ comp_code }}.

{% set figure_cfg.style = "max-width: 80%;" %}
{{ render_figure(
  comp_net.variant.fig_reporting_mode_by_component,
  "Variant reporting practices for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Variant reporting practices for {{ comp_code }}.** Bars represent the proportion of submitted sample outputs classified as high and low frequency reporting, high frequency only, or low frequency only, according to the metadata declarations associated with the variant outputs for this component.

Overall, {{ comp_code }} showed a median of {{ comp_net.variant.median_variants_in_consensus if comp_net.variant.median_variants_in_consensus is not none else "NA" }} variants with allele frequency above 75% reported in the metadata template, compared with {{ comp_net.variant.median_variants_in_consensus_vcf if comp_net.variant.median_variants_in_consensus_vcf is not none else "NA" }} corresponding variants represented in the consensus-derived VCF. The median number of discrepancies between both representations was {{ comp_net.variant.median_discrepancies_in_reported_variants if comp_net.variant.median_discrepancies_in_reported_variants is not none else "NA" }}, while the median total number of variants present in the submitted VCF files was {{ comp_net.variant.median_variants_in_vcf if comp_net.variant.median_variants_in_vcf is not none else "NA" }} (Table {{ table_counter.value + 1 }}, Figure {{ fig_counter.value + 1 }}).

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level influenza variant reporting metrics per sample for {{ comp_code }}.**

| Sample ID | Median variants >=75% AF in metadata (n={{ comp_net.variant_metadata_reporting.reported_n_labs }}) | Median variants >=75% AF in VCF (n={{ comp_net.variant_metadata_reporting.variants_in_consensus_vcf_reported_n_labs }}) | Median discrepancies between metadata and VCF | Median total variants in VCF (n={{ comp_net.variant_metadata_reporting.total_variants_in_vcf_reported_n_labs }}) |
|---|---:|---:|---:|---:|
{% for s in comp_net.variant.samples %}
| {{ s.collecting_lab_sample_id }} | {{ s.variants_in_consensus.median if s.variants_in_consensus and s.variants_in_consensus.median is not none else "NA" }} | {{ s.variants_in_consensus_vcf.median if s.variants_in_consensus_vcf and s.variants_in_consensus_vcf.median is not none else "NA" }} | {{ s.discrepancies_in_reported_variants.median if s.discrepancies_in_reported_variants and s.discrepancies_in_reported_variants.median is not none else "NA" }} | {{ s.variants_in_vcf.median if s.variants_in_vcf and s.variants_in_vcf.median is not none else "NA" }} |
{% endfor %}

{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_flu_aggregated_table_num = appendix_table_counter.value %}
{% set _ = variant_flu_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "aggregated_table_num": appendix_flu_aggregated_table_num
}) %}

These patterns indicate that influenza discrepancies reflect not only analytical differences in variant detection, but also differences in reporting conventions, allele-frequency thresholds, and reference selection. The aggregated structural summary for {{ comp_code }} is provided in Appendix Table {{ appendix_flu_aggregated_table_num }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 90%;" %}
{{ render_figure(
  comp_net.variant.fig_reporting_summary_by_sample,
  "Influenza variant reporting summary by sample for " ~ comp_code ~ "."
) }}


**Figure {{ fig_counter.value }}. Influenza variant reporting summary by sample for {{ comp_code }}.** Panel A shows, for each sample, the distribution across participating laboratories of the number of variants with allele frequency above 75% reported in the metadata template, the corresponding number represented in the consensus-derived VCF, and the discrepancies between both representations. Panel B shows the distribution across participating laboratories of the total number of variants present in the submitted VCF files for each sample. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range within the plotted scale, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

#### 6.{{ loop.index }}.4. Lineage, Subtype and Clade Assignment

Lineage, subtype and clade assignments submitted for the {{ comp_code }} component were evaluated for concordance with the curated gold standard classifications.

Across all participating laboratories, lineage/subtype concordance reached {{ pct(comp_net.typing.lineage_hit_pct) }}, whereas clade concordance reached {{ pct(comp_net.typing.clade_hit_pct) }}. The sample-level outcome distribution also shows that part of the observed discordance was associated with missing classifications or inconsistent completion of classification fields rather than with uniform analytical failure across all submissions.

{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_classification_table_num = appendix_table_counter.value %}
{% set _ = classification_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "table_num": appendix_classification_table_num
}) %}

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 98%;" %}
{{ render_figure(
  comp_net.typing.fig_stacked_bar_by_sample,
  "Classification outcome distribution per sample for " ~ comp_code ~ "."
) }}


**Figure {{ fig_counter.value }}. Classification outcome distribution per sample for {{ comp_code }}.** Panel A shows the proportion of lineage/subtype assignment Match, Discrepancy, and Not provided outcomes across participating laboratories for each sample. Panel B shows the corresponding proportions for clade assignments. Percentages are calculated over all participating laboratories in the component, so the Not provided segment captures samples for which lineage/subtype or clade information was not reported. Detailed sample-level concordance percentages are provided in Appendix Table {{ appendix_classification_table_num }}.

#### 6.{{ loop.index }}.5. Sample Quality Control Assessment

Sample-level QC in {{ comp_code }} was evaluated as concordance between the laboratory-reported Pass/Fail classification and the predefined gold standard status. QC concordance was heterogeneous across samples, and some laboratories did not report a formal QC assessment. Network-wide concordance for reported QC decisions was {{ pct(comp_net.qc.reported_match_rate_pct) }}.

{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_qc_table_num = appendix_table_counter.value %}
{% set _ = qc_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "table_num": appendix_qc_table_num
}) %}

{% set fig_counter.value = fig_counter.value + 1 %}

{% set figure_cfg.style = "max-width: 90%;" %}
{{ render_figure(
comp_net.qc.fig_qc_match_by_sample,
"Sample-level QC concordance for " ~ comp_code ~ " (Match, Discrepancy, and Not provided relative to the gold standard)."
) }}


**_Figure {{ fig_counter.value }}_. Sample-level QC concordance for {{ comp_code }} relative to the gold standard.** Bars represent the proportion of Match, Discrepancy, and Not provided outcomes per sample across participating laboratories. Higher discrepancy rates indicate samples for which laboratories more frequently diverged from the predefined QC status, whereas the Not provided segment captures missing QC assessments and is not interpreted as analytical disagreement. Detailed sample-level percentages and counts are provided in Appendix Table {{ appendix_qc_table_num }}.

#### 6.{{ loop.index }}.6. Pipeline Benchmarking and Comparative Performance

This section presents an exploratory comparative analysis of declared workflow configurations within {{ comp_code }}. Because laboratories differed in reference selection, software versions, parameterisation, reporting detail, and internal decision criteria, the results below should be interpreted as descriptive comparisons of observed performance patterns rather than as a controlled ranking of pipelines.

{% if comp_net.benchmarking.bioinformatics_protocol %}
##### Bioinformatics protocol

Based on metadata submissions, {{ comp_net.benchmarking.bioinformatics_protocol.total_number }} distinct bioinformatics protocols were reported for the {{ comp_code }} component. These summaries compare declared workflow configurations as they were used in practice across participating laboratories.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 80%;" %}
{{ render_figure(
  comp_net.benchmarking.bioinformatics_protocol.fig_discrepancy_boxplot,
  "Distribution of consensus discrepancies by pipeline configuration for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of consensus discrepancies by declared pipeline configuration for {{ comp_code }}.** This boxplot summarises sample-level consensus discrepancies stratified by bioinformatics protocol. The left y-axis shows discrepancy counts, while the right y-axis overlays lineage/type and clade classification accuracy for the same software configuration. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% set table_counter.value = table_counter.value + 1 %}

**Table {{ table_counter.value }}. Performance summary of declared bioinformatics protocols for {{ comp_code }}.**

| Bioinformatics protocol | Version | N labs | Median genome identity (%) | Median discrepancies | Median metadata completeness (%) | Clade concordance (%) | Lineage/type concordance (%) |
|---|---:|---:|---:|---:|---:|---:|---:|
{% for p in comp_net.benchmarking.bioinformatics_protocol.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ "%.2f"|format(p.median_identity_pct) }} | {{ p.median_discrepancies }} | {{ "%.1f"|format(p.median_metadata_completeness_pct) }} | {{ "%.1f"|format(p.clade_hit_pct) }} | {{ "%.1f"|format(p.lineage_hit_pct) }} |
{% endfor %}

The observed differences across configurations should be read in the context of heterogeneous laboratory practices, including differences in reference choice, parameterisation, and thresholding. The table and figures therefore help identify recurrent performance patterns within {{ comp_code }}, but they do not support a strong cross-laboratory ranking of bioinformatics protocols.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 96%;" %}
{{ render_figure(
  comp_net.benchmarking.bioinformatics_protocol.fig_metric_boxplots,
  "Distribution of performance metrics by pipeline configuration for " ~ comp_code ~ ".",
  "benchmark-figure landscape-benchmark-figure" if ((comp_net.benchmarking.bioinformatics_protocol.n_plot_groups | default(comp_net.benchmarking.bioinformatics_protocol.total_number)) >= 6 and (comp_net.benchmarking.bioinformatics_protocol.panel_count | default(99)) <= 2) else "benchmark-figure",
  has_panels=True
) }}

**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared pipeline configuration for {{ comp_code }}.** Multi-panel boxplots summarise sample-level performance stratified by bioinformatics protocols. Panel A displays genome identity (%), Panel B metadata completeness (%), and Panel C exact classification concordance (%). X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Where required, Panel A uses a truncated y-axis to highlight differences among high-identity values. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.dehosting %}
##### De-hosting software

{{ comp_net.benchmarking.dehosting.total_number }} distinct de-hosting software declarations were reported for the {{ comp_code }} component.
{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_dehosting_table_num = appendix_table_counter.value %}
{% set _ = benchmark_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "kind": "dehosting",
  "table_num": appendix_dehosting_table_num
}) %}

{% set fig_counter.value = fig_counter.value + 1 %}

The distribution below reflects only configurations with evaluable percentage of host reads values in the reported metadata, so the number of boxplots may be lower than the total number of declared de-hosting configurations. The full list of declared configurations and associated summary values is provided in Appendix Table {{ appendix_dehosting_table_num }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 96%;" %}
{{ render_figure(
  comp_net.benchmarking.dehosting.fig_metric_boxplots,
  "Distribution of percentage of host reads metrics by dehosting software version for " ~ comp_code ~ ".",
  "benchmark-figure landscape-benchmark-figure" if ((comp_net.benchmarking.dehosting.n_plot_groups | default(comp_net.benchmarking.dehosting.total_number)) >= 6 and (comp_net.benchmarking.dehosting.panel_count | default(99)) <= 2) else "benchmark-figure"
) }}

**Figure {{ fig_counter.value }}. Distribution of percentage of host reads by declared dehosting software version for {{ comp_code }}.** Boxplots summarise sample-level percentage of host reads stratified by dehosting software version. Only configurations with evaluable percentage of host reads values are displayed, so some declared software categories may be absent from the plot. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each version, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.preprocessing %}
##### Preprocessing software

{{ comp_net.benchmarking.preprocessing.total_number }} distinct pre-processing software configurations were reported for the {{ comp_code }} component.

{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_preprocessing_table_num = appendix_table_counter.value %}
{% set _ = benchmark_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "kind": "preprocessing",
  "table_num": appendix_preprocessing_table_num
}) %}

{% set fig_counter.value = fig_counter.value + 1 %}

Only pre-processing configurations with evaluable observations for the displayed metrics contribute to the figure, so some declared categories may not appear in the plot. The complete list of declared configurations and their summary values is provided in Appendix Table {{ appendix_preprocessing_table_num }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 96%;" %}
{{ render_figure(
  comp_net.benchmarking.preprocessing.fig_metric_boxplots,
  "Distribution of performance metrics by pre-processing software configuration for " ~ comp_code ~ ".",
  "benchmark-figure landscape-benchmark-figure" if ((comp_net.benchmarking.preprocessing.n_plot_groups | default(comp_net.benchmarking.preprocessing.total_number)) >= 6 and (comp_net.benchmarking.preprocessing.panel_count | default(99)) <= 2) else "benchmark-figure",
  has_panels=True
) }}

**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared pre-processing software configuration for {{ comp_code }}.** Multi-panel boxplots summarise sample-level performance stratified by pre-processing software. Panel A displays Number of reads sequenced and Panel B Reads passing filters. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.mapping %}
##### Mapping software

{{ comp_net.benchmarking.mapping.total_number }} distinct mapping software configurations were reported for the {{ comp_code }} component.

{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_mapping_table_num = appendix_table_counter.value %}
{% set _ = benchmark_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "kind": "mapping",
  "table_num": appendix_mapping_table_num
}) %}

{% set fig_counter.value = fig_counter.value + 1 %}

The mapping boxplots include only configurations for which the relevant performance metrics were available, which means that fewer categories may be plotted than were originally declared. Full configuration-level summaries are reported in Appendix Table {{ appendix_mapping_table_num }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 96%;" %}
{{ render_figure(
  comp_net.benchmarking.mapping.fig_metric_boxplots,
  "Distribution of performance metrics by mapping software configuration for " ~ comp_code ~ ".",
  "benchmark-figure landscape-benchmark-figure" if ((comp_net.benchmarking.mapping.n_plot_groups | default(comp_net.benchmarking.mapping.total_number)) >= 6 and (comp_net.benchmarking.mapping.panel_count | default(99)) <= 2) else "benchmark-figure"
) }}

**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared mapping software configuration for {{ comp_code }}.** Boxplots summarise sample-level performance stratified by mapping software. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.assembly %}
##### Assembly software

{{ comp_net.benchmarking.assembly.total_number }} distinct assembly software configurations were reported for the {{ comp_code }} component.
{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_assembly_table_num = appendix_table_counter.value %}
{% set _ = benchmark_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "kind": "assembly",
  "table_num": appendix_assembly_table_num
}) %}

{% set fig_counter.value = fig_counter.value + 1 %}

The assembly figures are restricted to configurations with evaluable values for the displayed metrics. As a result, some declared assembly categories may be absent from the plots; the full set of declared configurations and summary values is provided in Appendix Table {{ appendix_assembly_table_num }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 96%;" %}
{{ render_figure(
  comp_net.benchmarking.assembly.fig_metric_boxplots,
  "Distribution of performance metrics by assembly software configuration for " ~ comp_code ~ ".",
  "benchmark-figure landscape-benchmark-figure" if ((comp_net.benchmarking.assembly.n_plot_groups | default(comp_net.benchmarking.assembly.total_number)) >= 6 and (comp_net.benchmarking.assembly.panel_count | default(99)) <= 2) else "benchmark-figure",
  has_panels=True
) }}

**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared assembly software configuration for {{ comp_code }}.** Multi-panel boxplots summarise sample-level performance stratified by assembly software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers. Panel B uses a truncated y-axis to highlight differences among high-identity values.

{% endif %}

{% if comp_net.benchmarking.consensus_software %}
##### Consensus software

{{ comp_net.benchmarking.consensus_software.total_number }} distinct consensus software configurations were reported for the {{ comp_code }} component.
{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_consensus_software_table_num = appendix_table_counter.value %}
{% set _ = benchmark_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "kind": "consensus_software",
  "table_num": appendix_consensus_software_table_num
}) %}
{% set fig_counter.value = fig_counter.value + 1 %}

Only consensus software configurations with sufficient evaluable data are visualised in the figure below, so the plotted set may be smaller than the total set of declarations. All declared configurations and their associated summary values can be reviewed in Appendix Table {{ appendix_consensus_software_table_num }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 96%;" %}
{{ render_figure(
  comp_net.benchmarking.consensus_software.fig_metric_boxplots,
  "Distribution of performance metrics by consensus software configuration for " ~ comp_code ~ ".",
  "benchmark-figure landscape-benchmark-figure" if ((comp_net.benchmarking.consensus_software.n_plot_groups | default(comp_net.benchmarking.consensus_software.total_number)) >= 6 and (comp_net.benchmarking.consensus_software.panel_count | default(99)) <= 2) else "benchmark-figure",
  has_panels=True
) }}

**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared consensus software configuration for {{ comp_code }}.** Multi-panel boxplots summarise sample-level performance stratified by consensus software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers. For SARS1 and FLU1, Panel B uses a truncated y-axis to highlight differences among high-identity values.

{% endif %}

{% if comp_net.benchmarking.variant_calling %}
##### Variant calling software

{{ comp_net.benchmarking.variant_calling.total_number }} distinct variant calling software configurations were reported for the {{ comp_code }} component.
{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_variant_calling_table_num = appendix_table_counter.value %}
{% set _ = benchmark_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "kind": "variant_calling",
  "table_num": appendix_variant_calling_table_num
}) %}

{% set fig_counter.value = fig_counter.value + 1 %}

The plotted variant calling categories correspond only to configurations with evaluable observations for the displayed metrics. Consequently, the figure may show fewer configurations than were declared overall; the complete summaries are listed in Appendix Table {{ appendix_variant_calling_table_num }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 96%;" %}
{{ render_figure(
  comp_net.benchmarking.variant_calling.fig_metric_boxplots,
  "Distribution of performance metrics by variant calling software configuration for " ~ comp_code ~ ".",
  "benchmark-figure landscape-benchmark-figure" if ((comp_net.benchmarking.variant_calling.n_plot_groups | default(comp_net.benchmarking.variant_calling.total_number)) >= 6 and (comp_net.benchmarking.variant_calling.panel_count | default(99)) <= 2) else "benchmark-figure",
  has_panels=True
) }}

{% if comp_code[:3] == "FLU" %}
**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared variant calling software configuration for {{ comp_code }}.** Panel A is a stacked bar chart showing the number of evaluable samples assigned to each allele frequency reporting pattern for each software configuration. Boxplot Panel B displays the number of reported variants with AF >=75%, Panel C the number of variants with AF >=75% in the submitted VCF, Panel D the number of variants with effect, Panel E metadata-VCF discrepancies, and Panel F the total number of variants present in the submitted VCF files. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. In the boxplots, the central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.
{% else %}
**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared variant calling software configuration for {{ comp_code }}.** Panel A is a stacked bar chart showing the number of evaluable samples assigned to each allele frequency reporting pattern for each software configuration. Boxplot Panel B displays discrepancies in reported variants with AF >=75% in the submitted VCF, Panel C discrepancies in reported variants with effect, Panel D successful hits, and Panel E total discrepancies. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. In the boxplots, the central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.
{% endif %}

{% endif %}

{% if comp_net.benchmarking.clade_assignment %}
##### Clade Assignment Software

{{ comp_net.benchmarking.clade_assignment.total_number }} distinct clade assignment software configurations were reported for the {{ comp_code }} component. For this category, configurations were counted as unique combinations of software name, software version, and clade assignment database version when available.
{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_clade_assignment_table_num = appendix_table_counter.value %}
{% set _ = benchmark_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "kind": "clade_assignment",
  "table_num": appendix_clade_assignment_table_num
}) %}

{% set fig_counter.value = fig_counter.value + 1 %}

Because clade concordance could not be evaluated for every declared configuration, the boxplot includes only categories with usable observations. The full list of declared configurations and their summary values is available in Appendix Table {{ appendix_clade_assignment_table_num }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 96%;" %}
{{ render_figure(
  comp_net.benchmarking.clade_assignment.fig_metric_boxplots,
  "Distribution of performance metrics by clade assignment software configuration for " ~ comp_code ~ ".",
  "benchmark-figure landscape-benchmark-figure" if ((comp_net.benchmarking.clade_assignment.n_plot_groups | default(comp_net.benchmarking.clade_assignment.total_number)) >= 6 and (comp_net.benchmarking.clade_assignment.panel_count | default(99)) <= 2) else "benchmark-figure"
) }}

**Figure {{ fig_counter.value }}. Distribution of clade concordance by declared clade assignment software configuration for {{ comp_code }}.** This boxplot summarises sample-level clade concordance stratified by clade assignment software configuration, where each configuration corresponds to a unique combination of software name, software version, and clade assignment database version when available. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.lineage_assignment %}
##### Lineage Assignment Software Name

{{ comp_net.benchmarking.lineage_assignment.total_number }} distinct lineage assignment software configurations were reported for the {{ comp_code }} component. For this category, configurations were counted as unique combinations of software name, software version, and lineage assignment database version when available.
{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_lineage_assignment_table_num = appendix_table_counter.value %}
{% set _ = benchmark_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "kind": "lineage_assignment",
  "table_num": appendix_lineage_assignment_table_num
}) %}

{% set fig_counter.value = fig_counter.value + 1 %}

Lineage assignment configurations are shown only when concordance values were evaluable for the submitted observations, so the plotted categories may represent only a subset of the declarations. The complete configuration-level summary is provided in Appendix Table {{ appendix_lineage_assignment_table_num }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 96%;" %}
{{ render_figure(
  comp_net.benchmarking.lineage_assignment.fig_metric_boxplots,
  "Distribution of performance metrics by lineage assignment software configuration for " ~ comp_code ~ ".",
  "benchmark-figure landscape-benchmark-figure" if ((comp_net.benchmarking.lineage_assignment.n_plot_groups | default(comp_net.benchmarking.lineage_assignment.total_number)) >= 6 and (comp_net.benchmarking.lineage_assignment.panel_count | default(99)) <= 2) else "benchmark-figure"
) }}

**Figure {{ fig_counter.value }}. Distribution of lineage concordance by declared lineage assignment software configuration for {{ comp_code }}.** This boxplot summarises sample-level lineage concordance stratified by lineage assignment software configuration, where each configuration corresponds to a unique combination of software name, software version, and lineage assignment database version when available. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.type_assignment %}
##### Type Assignment Software Name

{{ comp_net.benchmarking.type_assignment.total_number }} distinct type assignment software configurations were reported for the {{ comp_code }} component. For this category, configurations were counted as unique combinations of software name, software version, and type assignment database version when available.
{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_type_assignment_table_num = appendix_table_counter.value %}
{% set _ = benchmark_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "kind": "type_assignment",
  "table_num": appendix_type_assignment_table_num
}) %}

{% set fig_counter.value = fig_counter.value + 1 %}

The type assignment plot is limited to configurations with evaluable concordance results, and therefore may contain fewer categories than the total number declared in metadata. The complete list of declared configurations and summary values is provided in Appendix Table {{ appendix_type_assignment_table_num }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 96%;" %}
{{ render_figure(
  comp_net.benchmarking.type_assignment.fig_metric_boxplots,
  "Distribution of performance metrics by type assignment software configuration for " ~ comp_code ~ ".",
  "benchmark-figure landscape-benchmark-figure" if ((comp_net.benchmarking.type_assignment.n_plot_groups | default(comp_net.benchmarking.type_assignment.total_number)) >= 6 and (comp_net.benchmarking.type_assignment.panel_count | default(99)) <= 2) else "benchmark-figure"
) }}

**Figure {{ fig_counter.value }}. Distribution of type concordance by declared type assignment software configuration for {{ comp_code }}.** This boxplot summarises sample-level type concordance stratified by type assignment software configuration, where each configuration corresponds to a unique combination of software name, software version, and type assignment database version when available. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.subtype_assignment %}

##### Subtype Assignment Software Name

{{ comp_net.benchmarking.subtype_assignment.total_number }} distinct subtype assignment software configurations were reported for the {{ comp_code }} component. For this category, configurations were counted as unique combinations of software name, software version, and subtype assignment database version when available.
{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set appendix_subtype_assignment_table_num = appendix_table_counter.value %}
{% set _ = benchmark_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp_net": comp_net,
  "kind": "subtype_assignment",
  "table_num": appendix_subtype_assignment_table_num
}) %}

{% set fig_counter.value = fig_counter.value + 1 %}

Subtype assignment configurations are plotted only when evaluable concordance data were available, so some declared categories may not be represented in the figure. Appendix Table {{ appendix_subtype_assignment_table_num }} contains the full list of declarations and their summary values.

{% set fig_counter.value = fig_counter.value + 1 %}
{% set figure_cfg.style = "max-width: 96%;" %}
{{ render_figure(
  comp_net.benchmarking.subtype_assignment.fig_metric_boxplots,
  "Distribution of performance metrics by subtype assignment software configuration for " ~ comp_code ~ ".",
  "benchmark-figure landscape-benchmark-figure" if ((comp_net.benchmarking.subtype_assignment.n_plot_groups | default(comp_net.benchmarking.subtype_assignment.total_number)) >= 6 and (comp_net.benchmarking.subtype_assignment.panel_count | default(99)) <= 2) else "benchmark-figure"
) }}

**Figure {{ fig_counter.value }}. Distribution of subtype concordance by declared subtype assignment software configuration for {{ comp_code }}.** This boxplot summarises sample-level subtype concordance stratified by subtype assignment software configuration, where each configuration corresponds to a unique combination of software name, software version, and subtype assignment database version when available. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% endfor %}

## 7. Discussion

The 2026 RELECOV Dry-Lab EQA provides the first network-wide dry-lab assessment focused specifically on bioinformatic performance across consensus reconstruction, variant reporting, classification, metadata reporting, and QC interpretation. By combining ECDC datasets with in-silico influenza material, the exercise captures both routine-use analytical behaviour and performance under heterogeneous reference and reporting conditions.

### 7.1. Consensus Genome Reconstruction

Consensus reconstruction results were strongest in the Illumina-based components overall, with a combined median genome identity of {{ pct(general.general_results.consensus.median_identity_illumina_pct, 2) }} compared with {{ pct(general.general_results.consensus.median_identity_nanopore_pct, 2) }} across the Nanopore-based components. When grouped by virus, the combined median genome identity was also higher for SARS-CoV-2 components ({{ pct(general.general_results.consensus.median_identity_sars_pct, 2) }}) than for influenza components ({{ pct(general.general_results.consensus.median_identity_influenza_pct, 2) }}). The fact that the influenza median remained slightly below the combined Nanopore-based median suggests that lower identity in influenza cannot be attributed to sequencing platform alone and was also influenced by virus-specific analytical complexity. At component level, SARS1 and SARS2 showed median identities of {{ pct(general.components.SARS1.consensus.median_identity_pct, 2) }} and {{ pct(general.components.SARS2.consensus.median_identity_pct, 2) }}, whereas FLU1 and FLU2 showed lower medians of {{ pct(general.components.FLU1.consensus.median_identity_pct, 2) }} and {{ pct(general.components.FLU2.consensus.median_identity_pct, 2) }}. These results suggest that influenza consensus reconstruction remains less standardised across the network than SARS-CoV-2 reconstruction and would benefit from clearer best-practice guidance on coverage thresholds, masking behaviour, indel handling, and segment-level quality criteria.

At the same time, the ranges observed across laboratories show that high medians did not eliminate outlier behaviour. In particular, the minimum identity values in SARS2 and FLU2 dropped to {{ pct(general.components.SARS2.consensus.identity_pct_min, 2) }} and {{ pct(general.components.FLU2.consensus.identity_pct_min, 2) }}, indicating that a subset of submissions diverged markedly from the curated gold standard.

The dominant discrepancy categories also differed by component. SARS1 was dominated by stretches of Ns in submitted consensuses where defined nucleotides were present in the gold standard (`ns2nt`), whereas SARS2 was dominated by defined nucleotides where stretches of Ns were present in the gold standard (`nt2ns`). These patterns are consistent with differences in masking behaviour and minimum coverage policies relative to the gold standard reconstruction criteria. This was especially clear in samples such as SARS4 in SARS1 and SARS8 in SARS2, where low depth or mixed-site complexity increased sensitivity to local masking decisions and likely exposed differences between ECDC reconstruction criteria and the thresholds applied within the network. In influenza, both FLU1 and FLU2 were dominated by deletions relative to the gold standard, suggesting that consensus generation parameters remain important contributors to inter-laboratory divergence. The most discrepancy-prone samples were in-silico datasets  with segment dropouts or contamination, which supports the interpretation that indel calling thresholds and segment-specific handling rules remain insufficiently harmonised in influenza workflows.

### 7.2. Variant Detection and Reporting

For SARS-CoV-2, variant detection performance did not follow a simple platform ranking. The median number of discrepancies relative to the curated variant set was {{ general.general_results.sars_variants.median_discrepancy_illumina }} in the Illumina component and {{ general.general_results.sars_variants.median_discrepancy_nanopore }} in the Nanopore component. This indicates that platform effects were present, but that they interacted with sample composition, reporting choices, and software configuration rather than determining performance on their own.

The submitted metadata documented substantial diversity in variant reporting behaviour. Across SARS-CoV-2 submissions, {{ general.general_results.sars_variants.high_and_low_freq_pct }} of laboratories reported both high- and low-frequency variants, whereas {{ general.general_results.sars_variants.high_freq_only_pct }} reported high-frequency variants only. Influenza reporting was more heterogeneous: {{ general.general_results.influenza_variants.high_and_low_freq_pct }} of laboratories reported both high- and low-frequency variants, {{ general.general_results.influenza_variants.low_freq_only_pct }} reported low-frequency variants only, and {{ general.general_results.influenza_variants.high_freq_only_pct }} reported high-frequency variants only.

Influenza results especially highlight the consequences of heterogeneous structural reporting. The network-level median number of variants with AF >=75% reported in metadata was {{ general.general_results.influenza_variants.median_variants_in_consensus }}, whereas the corresponding median derived from submitted VCF files was {{ general.general_results.influenza_variants.median_variants_in_consensus_vcf }}. The median discrepancy between these two representations was {{ general.general_results.influenza_variants.median_discrepancies_in_reported_variants }}, and the total number of variants present in submitted VCF files ranged from {{ general.general_results.influenza_variants.min_variants_in_vcf }} to {{ general.general_results.influenza_variants.max_variants_in_vcf }}. Together, these values indicate that influenza variant outputs were not directly comparable under a single harmonised coordinate framework and that reporting conventions differed markedly across laboratories. That heterogeneity is also visible in the mixed reporting modes, the use of multiple reference genomes, and the wide structural ranges in both total VCF content and metadata-VCF discrepancies.

In SARS-CoV-2, discrepancies between metadata-reported and VCF-derived values were generally limited, but the discrepancy profile still points to differences in filtering and interpretation rules. In SARS1, the most frequent discrepancy type was "De novo variants", and samples such as SARS3 and especially SARS4 showed the highest discrepancy burdens, consistent with the fact that the curated reference VCF retained only variants supported by a minimum total depth of 10x and a minimum allele frequency of 25%. Additional technical limitations also affected the calculation of VCF-derived variants with effect in this component: one submission was generated against an XBB reference genome rather than against the Wuhan-based reference framework used for the EQA gold standard, and another submission provided iVar TSV outputs instead of VCF files, which prevented standard annotation of VCF content for effect-based comparisons. This means that laboratories using different calling thresholds, or laboratories reporting only high-frequency variants, may inflate either "Missing expected variants" or "De novo variant" categories even when the underlying VCF is internally coherent. In SARS2, overall discrepancy levels also remained low, but "Missing expected variants" and "De novo variants" were again the most frequent discrepancy types, and the highest-burden samples were the low-quality materials SARS9 and SARS10. In that component, the curated reference VCF was more stringent, retaining only variants supported by a minimum total depth of 20x and a minimum allele frequency of 50%, which likely increased the sensitivity of concordance metrics to differences in filtering and reporting decisions. Technical limitations also affected the availability of some VCF-derived summary metrics in this component: one laboratory did not submit VCF files, and in another case the VCFs had been generated with Medaka, which does not provide allele-frequency information. As a result, it was not possible in those submissions to derive AF >=75% variant counts or to calculate annotated variants with effect from the VCF content.

For influenza, only a minority of laboratories reported metadata counts in a way that could be directly interpreted against the submitted VCFs, and some laboratories appear to have counted all VCF variants, including low-frequency calls, as if they belonged to the AF >=75% category. Others reported VCFs with very large total numbers of variants, indicating the inclusion of very low-frequency events supported by few reads. In practical terms, this means that the influenza discrepancies do not only reflect analytical differences in variant detection, but also differences in how laboratories interpreted software outputs and translated them into the metadata template. This reinforces the need for harmonised best practices defining which variants should be reported, under which AF thresholds, and how those thresholds should be represented in metadata.

### 7.3. Classification and QC Interpretation

Classification performance was acceptable overall but clearly stronger for lineage/type assignment than for clade assignment. SARS-CoV-2 lineage concordance reached {{ pct(general.general_results.classification.sars_lineage_concordance_pct) }}, compared with {{ pct(general.general_results.classification.sars_clade_concordance_pct) }} for SARS-CoV-2 clade assignment. Influenza type/subtype concordance reached {{ pct(general.general_results.classification.influenza_type_concordance_pct) }}, compared with {{ pct(general.general_results.classification.influenza_clade_concordance_pct) }} for influenza clade assignment. In the SARS-CoV-2 submissions, this difference is consistent with metadata-level reporting problems in the clade field itself: among the clade assignments reviewed in the submitted JSON files, some were left empty and others contained values that matched the lineage assignment or had lineage-like syntax rather than a clade designation. This suggests that part of the excess clade discordance reflects field completion and nomenclature issues in addition to true analytical misclassification.

QC interpretation showed additional between-component differences. Only 9 of the 19 participating laboratories reported at least one sample-level QC assessment in their submitted metadata. Among evaluable QC decisions, network-wide concordance was {{ pct(general.qc.reported_match_rate_pct) }}, and component-level concordance ranged from {{ pct(general.components.SARS1.qc.reported_match_rate_pct) }} in SARS1 to {{ pct(general.components.FLU2.qc.reported_match_rate_pct) }} in FLU2. This indicates that QC interpretation was not equally stable across all datasets, while also showing that many laboratories either did not apply or did not report a formal sample-level QC decision in the metadata template.

Because QC assignment depends on how each laboratory interprets coverage, ambiguity, contamination, and other internal acceptance criteria, these results are better understood as reflecting heterogeneity in QC interpretation and reporting rather than a simple measure of analytical correctness alone.

At sample level, the classification and QC discrepancies were concentrated in analytically difficult materials rather than being evenly distributed. In SARS1, the samples showing discordance against the ECDC gold standard were those flagged as low quality, such as SARS4 and SARS5, where low coverage or excess mixed sites likely shifted both QC status and downstream classification. In SARS2, SARS8 emerged as the most problematic lineage-assignment sample, with several laboratories reporting the parent lineage rather than the exact expected designation, again suggesting that small differences in masking and low-coverage treatment can propagate into classification discordance. At the same time, part of the observed clade discordance in SARS-CoV-2 was clearly metadata-driven: some laboratories populated the clade field with lineage values or left clade values empty, which inflated apparent clade disagreement independently of software performance or database version.

The influenza components showed a different pattern. A relatively high proportion of laboratories did not report clade values at all, although those that did generally reported concordant assignments; the residual heterogeneity mainly reflected differences between legacy clade labels and more specific subclade labels rather than outright analytical failure. Subtype discordance was more informative, especially in FLU1, where one A/H5N1 sample was reported as A/H5N6 by a laboratory using a BLAST-based assignment strategy, illustrating that subtype performance can still depend strongly on software choice and database interpretation. QC reporting for influenza was also sparse, particularly in FLU2, but among the small number of reported QC decisions the agreement with the expected status was generally high. This suggests that the main challenge in influenza QC was not incorrect interpretation as much as incomplete adoption of formal QC reporting.

### 7.4. Workflow Diversity and Reporting Constraints

The metadata confirms that RELECOV laboratories currently use a diverse analytical landscape. A total of {{ general.metadata_completeness.total_workflows }} distinct workflows were identified across participating laboratories, together with distinct tools or tool/version combinations for consensus genome generation ({{ general.metadata_completeness.total_consensus_softwares }}), variant calling ({{ general.metadata_completeness.total_variant_softwares }}), SARS-CoV-2 lineage assignment ({{ general.metadata_completeness.total_lineage_assignment_softwares }}), influeza type/subtype assignment ({{ general.metadata_completeness.total_subtype_assignment_softwares }}) and clade assignment ({{ general.metadata_completeness.total_clade_assignment_softwares }}).

This diversity is analytically valuable, but its interpretation is constrained by incomplete metadata reporting. Only {{ pct(general.metadata_completeness.software_version_pct) }} of software-version fields were completed, {{ pct(general.metadata_completeness.coverage_threshold_pct) }} of submitted samples specified a minimum coverage threshold, {{ pct(general.metadata_completeness.variant_calling_params_pct) }} reported variant calling parameters, and {{ pct(general.metadata_completeness.reference_genome_pct) }} reported a reference genome accession or identifier. For that reason, some plausible explanations for performance differences can only be discussed as contributing context rather than demonstrated causal effects.

The main incompleteness drivers were variant calling, pre-processing, and mapping fields, followed by QC metrics, de-hosting, consensus analysis, and classification-related metadata. This pattern suggests that laboratories were more consistent in declaring core tool identities than in documenting the exact thresholds and parameter sets that determine analytical behaviour.

The submitted metadata also supports the view that parameter heterogeneity contributed to consensus and variant calling variability. Across all submitted samples, laboratories reported at least 8 different conventions for minimum coverage thresholds, 13 distinct consensus parameter strings, 13 distinct mapping parameter strings, and 14 distinct variant calling parameter strings. These differences do not prove causality for any individual discrepancy, but they do show that laboratories were not applying a uniform set of masking, filtering, or coverage rules.

The benchmarking results also suggest that apparent top-performing configurations should be interpreted against the number of laboratories supporting them. In SARS1 and SARS2, several of the most favourable raw performance values were associated with single-laboratory configurations, whereas nf-core/viralrecon was supported by multiple laboratories and combined consistently high genome identity with low discrepancy burdens and strong metadata completeness. That pattern makes it a more informative indicator of reproducible network performance than a nominally better single-observation configuration. At the same time, the SARS results do not show a simple linear relationship between consensus discrepancy counts and classification concordance: most pipelines retained very high lineage or clade performance despite modest sequence-level variation, suggesting that SARS classification is robust to moderate reconstruction differences but still sensitive to reporting quality.

The influenza benchmarking profiles were more heterogeneous and showed stronger configuration-dependent effects. In FLU1, DRAGEN achieved the highest identity but only on a single observation, whereas custom pipelines and different IRMA versions showed broader performance ranges. INSaFLU performed substantially worse in both identity and discrepancy burden in that component, suggesting that workflow choice can have a larger practical impact in influenza than in SARS-CoV-2 under these datasets. In FLU2, IRMA v1.3.1 provided the most balanced profile across discrepancy burden and classification performance, whereas other versions, particularly IRMA v1.2.0, showed much poorer discrepancy behaviour despite acceptable classification fields. Taken together, these results indicate that influenza benchmarking is more sensitive to software versioning, parameterisation, and component-specific sample properties, and therefore requires more explicit best-practice recommendations rather than simple transfer of SARS-CoV-2 assumptions.


### 7.5. Metadata Reporting and Schema Compliance

Metadata quality affected not only interpretability but also interoperability. The exercise showed that laboratories were generally able to declare the core structure of their workflows, yet the level of detail required for reproducibility was often incomplete. This was especially evident for software-version fields, minimum coverage thresholds, variant calling parameters, and reference genome identifiers, all of which are necessary to reconstruct how a consensus or variant set was generated and to determine whether observed differences reflect analytical choice or true performance variation.

The validation process also showed that schema compliance remains an operational issue in its own right. As reported in the metadata results section, only 26.32% of submissions (5 laboratories) were fully compliant with controlled-vocabulary requirements, whereas 73.68% (14 laboratories) required at least one manual correction because non-standard values had been used in controlled fields. The most frequent problems involved free-text software names in dropdown-based fields, inconsistent completion of lineage, clade, type, or subtype assignments, and missing mandatory entries. These issues do not necessarily imply poor analytical practice, but they do reduce the value of the metadata for automated validation, cross-laboratory comparison, and downstream integration into the RELECOV platform.

A recurrent source of non-compliance was the use of free-text entries in fields for which predefined dropdown options were available, particularly for software names. This occurred even though the metadata template, including its controlled-vocabulary dropdowns, had been distributed two weeks before the start of the exercise to give laboratories time to review the available options and identify any missing software tools for possible inclusion in the schema, together with guidance on how mandatory fields should be completed when data were not available. In practice, this means that part of the harmonisation problem lies not only in analytical diversity itself, but in the difficulty of consistently mapping that diversity into a controlled metadata structure.

Another recurrent issue concerned database-version reporting for lineage and clade assignment tools. In some submissions, the value entered in the database-version field appears to correspond to the software version rather than the actual database version, particularly for tools such as Nextclade and Pangolin where both identifiers are distinct and analytically relevant. This suggests that the semantics of these metadata fields were not always clear to participants and that future versions of the template should include stricter validation rules and more explicit examples distinguishing software version from database release. More broadly, some laboratories also relied on relatively outdated software versions, such as older releases of iVar. This is not a minor technical detail, because software versioning can affect variant calling behaviour, allele-frequency estimation, and indel-depth interpretation, and therefore can alter both the analytical outputs themselves and the comparability of results across the network.

### 7.6. Implications for RELECOV 2.0

Taken together, the results support a harmonisation strategy centred on minimum performance and reporting standards rather than on enforcement of a single analytical pipeline. The data do not support a universal workflow ranking that would apply equally across all viruses, platforms, and tasks. Instead, they show that performance depends on the interaction between dataset characteristics, reporting conventions, software choice, and parameterisation.

## 8. Conclusions

The 2026 RELECOV Dry-Lab EQA shows that participating laboratories already have substantial bioinformatic capacity for respiratory virus genomic surveillance, but that performance and comparability still depend strongly on the analytical context in which each task is performed.

Consensus genome reconstruction was generally strongest in the Illumina-based components, while broader performance ranges in SARS2 and FLU2 indicate that a subset of submissions remained highly sensitive to masking behaviour, coverage thresholds, and consensus-generation choices. Variant analysis showed that direct SARS-CoV-2 comparison against curated reference sets is feasible, whereas influenza reporting remained much more heterogeneous because of mixed allele frequency reporting strategies, multiple reference backbones, and large discrepancies between metadata-reported and VCF-derived summaries.

Classification and QC interpretation further showed that harmonisation challenges are not limited to core sequence processing. Lineage/type assignment was more concordant than clade assignment, and part of the excess clade discordance in SARS-CoV-2 appears to reflect metadata completion and nomenclature problems in the clade field itself. QC interpretation was also unevenly reported, with only a subset of laboratories providing explicit sample-level QC assessments in the metadata template.

Overall, the results support RELECOV 2.0 priorities centred on:

- minimum performance standards for consensus reconstruction and variant reporting
- clearer rules for masking, coverage thresholds, and allele frequency reporting
- stronger metadata requirements for software versions, parameters, and reference genomes
- improved consistency in classification and QC field completion
- component-aware benchmarking rather than a single cross-context workflow ranking

Taken together, these findings provide a practical basis for harmonising analytical expectations across the network while preserving the methodological flexibility needed for different pathogens, sequencing platforms, and surveillance scenarios.

The EQA therefore provides a robust technical basis for harmonised, performance-driven genomic surveillance within RELECOV 2.0.

{% if labdata %}
{% set lab_code = labdata.lab.lab_cod | default(labdata.lab.submitting_institution_id) %}
<h2 id="9-individual-laboratory-technical-report" class="no-page-break">9. Individual Laboratory Technical Report</h2>

<h3 id="laboratory-{{ labdata.lab.lab_cod|lower }}" class="no-page-break">Laboratory: {{ labdata.lab.laboratory_name }} ({{ labdata.lab.lab_cod }})</h3>

This section provides a detailed technical assessment of the analytical results submitted by **{{ labdata.lab.lab_cod }}** within the 2026 RELECOV Dry-Lab EQA. Performance metrics are benchmarked against curated gold standards and contextualised relative to aggregated network-wide performance distributions. Network medians and interquartile ranges are provided for comparative interpretation, without disclosure of other laboratories’ identities.

The purpose of this section is to support technical optimisation, parameter harmonisation, and alignment with the analytical standards defined within RELECOV 2.0.

Only files, metadata fields, and derived analytical metrics actually provided by the laboratory are displayed in this individual report. If a file was not submitted, or a metadata field was not provided, the corresponding table entries, panels, or figures are omitted for that laboratory.

<h3 id="91-participation-overview" class="no-page-break">9.1. Participation Overview</h3>

The laboratory analysed **{{ labdata.components | length }}** out of 4 components. Network median components analysed per laboratory: **{{ general.median_components_analysed_per_lab }}**.

Analysed components:

{% for comp_code, comp_info in general.components.items() %}
- {{ comp_code }} ({{ comp_info.name }}): {{ "✔" if comp_code in labdata.components.keys() else "✖" }}
{% endfor %}

Regarding general metadata completeness:

- Metadata completeness for **{{ labdata.lab.lab_cod }}**: **{{ pct(labdata.metadata.completeness_pct) }}**
- Network median metadata completeness: **{{ pct(general.metadata_completeness.median_pct) }}**  
- Network range: **{{ pct(general.metadata_completeness.min_pct) }}–{{ pct(general.metadata_completeness.max_pct) }}**

{% if labdata.metadata.primary_incompleteness_drivers %}
Primary contributors to incompleteness for {{ labdata.lab.lab_cod }}:
<ul class="compact-list">
{% for d in labdata.metadata.primary_incompleteness_drivers %}
<li>{{ d }}</li>
{% endfor %}
</ul>
{% endif %}

{% for comp_code, comp in labdata.components.items() %}

### 9.{{ loop.index + 1 }}. {{ comp_code }} ({{ comp.display_name }})

The laboratory submitted results for the **{{ comp_code }}** component from {{ comp.sequencing_instrument_platform }} platform.

Number of ssubmitted outputs:

- `.fasta`: **{{ comp.metadata.fasta_submitted }} out of {{ comp.metadata.fasta_expected }} minimum expected**
- `.vcf`: **{{ comp.metadata.vcf_submitted }} out of {{ comp.metadata.vcf_expected }} minimum expected**

Sections, tables, and figures below are shown only when the corresponding files or metadata were provided for this component. Missing submissions or non-reported metadata fields are not displayed for **{{ labdata.lab.lab_cod }}**.

Regarding metadata completeness for {{ comp_code }}:

- Metadata completeness for **{{ comp.lab.lab_cod }}**: **{{ pct(comp.metadata.completeness_pct) }}**
- Network median metadata completeness: **{{ pct(general.components[comp_code].metadata_completeness_median) }}**  
- Network range: **{{ pct(general.components[comp_code].metadata_completeness_min_pct) }}–{{ pct(general.components[comp_code].metadata_completeness_max_pct) }}**

{% if comp.metadata.primary_incompleteness_drivers %}
Primary contributors to incompleteness for {{ comp_code }}:
<ul class="compact-list">
{% for d in comp.metadata.primary_incompleteness_drivers %}
<li>{{ d }}</li>
{% endfor %}
</ul>
{% endif %}

#### 9.{{ loop.index + 1 }}.1. Consensus Genome Reconstruction Performance

Consensus genome sequences (`.fasta`) submitted by **{{ labdata.lab.lab_cod }}** were compared against the curated gold standard reference for each sample included in the {{ comp_code }} component.

##### Per-sample summary metrics

{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set lab_consensus_metrics_table_num = appendix_table_counter.value %}
The detailed per-sample consensus reconstruction metrics for **{{ labdata.lab.lab_cod }}** are provided in Appendix Table {{ lab_consensus_metrics_table_num }}. The figure below summarises overall sequence similarity and discrepancy burden relative to the curated gold standard reference for {{ labdata.lab.lab_cod }} compared to the network.


{% set consensus_distribution_panel_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/consensus_distribution_panel.png" %}
{% if path_exists(consensus_distribution_panel_path) %}
{% set fig_counter.value = fig_counter.value + 1 %}

{% set figure_cfg.style = "max-width: 90%;" %}
{{ render_figure(
  consensus_distribution_panel_path,
  comp_code ~ ": distribution of consensus discrepancies and genome identity per sample across the network; black diamond indicates " ~ labdata.lab.lab_cod ~ "."
) }}

**Figure {{ fig_counter.value }}. Consensus reconstruction performance across participating laboratories ({{ comp_code }}).** Panel A shows the distribution of total consensus discrepancies per sample relative to the curated gold standard across the RELECOV network. Panel B shows the corresponding distribution of genome identity values per sample. In both panels, the central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. In Panel B, the y-axis is truncated to highlight differences among high-identity values. The black diamond corresponds to the results obtained by **{{ labdata.lab.lab_cod }}**.
{% endif %}

##### Discrepancy type breakdown per sample
{% set consensus_breakdown_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/consensus_discrepancy_breakdown_by_sample.png" %}
{% if path_exists(consensus_breakdown_path) %}
{% set fig_counter.value = fig_counter.value + 1 %}

{% set figure_cfg.style = "max-width: 80%;" %}
{{ render_figure(
  consensus_breakdown_path,
  comp_code ~ ": discrepancy type breakdown by sample for " ~ labdata.lab.lab_cod ~ "."
) }}

**Figure {{ fig_counter.value }}. Discrepancy type breakdown by sample for {{ labdata.lab.lab_cod }} ({{ comp_code }}).** Stacked bars show the contribution of each discrepancy category to the total consensus differences observed for each sample submitted by **{{ labdata.lab.lab_cod }}**.
{% endif %}

{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set lab_consensus_breakdown_table_num = appendix_table_counter.value %}
{% set _ = lab_consensus_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp": comp,
  "metrics_table_num": lab_consensus_metrics_table_num,
  "breakdown_table_num": lab_consensus_breakdown_table_num
}) %}
The full discrepancy-type breakdown per sample for **{{ labdata.lab.lab_cod }}** is provided in Appendix Table {{ lab_consensus_breakdown_table_num }}.

{% if comp.metadata.vcf_submitted >=1 %}

#### 9.{{ loop.index + 1 }}.2. Variant Detection Performance

{% if comp_code in ["SARS1", "SARS2"] %}
For SARS-CoV-2, variant call files (`.vcf`) submitted by **{{ labdata.lab.lab_cod }}** were compared against the curated reference variant set for each sample included in the {{ comp_code }} component.

The metrics presented in Table {{ table_counter.value }} summarise per-sample variant detection accuracy relative to the curated reference variant set and benchmark the laboratory’s results against the network median for the same sample. The laboratory-reported variant counts declared in the metadata were also compared against the values derived directly from the submitted VCF files for each sample.

{% set variant_detection_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/variant_metadata_vs_vcf_distribution.png" %}
{% if path_exists(variant_detection_path) %}
{% set fig_counter.value = fig_counter.value + 1 %}

{% set figure_cfg.style = "max-width: 96%;" %}
{{ render_figure(
  variant_detection_path,
  comp_code ~ ": distribution of variant detection metrics across the network; black diamond indicates " ~ labdata.lab.lab_cod ~ "."
) }}

**Figure {{ fig_counter.value }}. Variant detection performance across participating laboratories ({{ comp_code }}).** Panel A shows the distribution of total variant discrepancies per sample across the RELECOV network. Panel B shows the corresponding distribution of successful hits per sample. In both panels, the central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range across the network, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. The black diamond corresponds to the results obtained by **{{ labdata.lab.lab_cod }}**.
{% endif %}

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Per-sample variant detection performance metrics for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

| Sample ID | Reporting mode | Expected hits | {{ labdata.lab.lab_cod }} total discrepancies | Network median total discrepancies | {{ labdata.lab.lab_cod }} successful hits | Network median successful hits | Wrong variants | Insertions | Deletions | Missing expected variants | De novo variants |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{% for collecting_lab_sample_id, s in comp.samples.items() -%}
{% set ns = (general.components[comp_code].variant.samples | selectattr("collecting_lab_sample_id","equalto",collecting_lab_sample_id) | list | first) %}
| {{ collecting_lab_sample_id }} | {{ "High and low frequency" if s.variants.high_and_low_freq else ("High frequency only" if s.variants.high_freq_only else ("Low frequency only" if s.variants.low_freq_only else "NA")) }} | {{ ns.expected_hits if ns and ns.expected_hits is not none else "NA" }} | {{ s.variants.total_discrepancies if s.variants.total_discrepancies is not none else "NA" }} | {{ ns.median_discrepancies if ns else "NA" }} | {{ s.variants.successful_hits if s.variants.successful_hits is not none else "NA" }} | {{ ns.median_successful_hits if ns else "NA" }} | {{ s.variants.wrong_nt if s.variants.wrong_nt is not none else "NA" }} | {{ s.variants.insertions if s.variants.insertions is not none else "NA" }} | {{ s.variants.deletions if s.variants.deletions is not none else "NA" }} | {{ s.variants.missing if s.variants.missing is not none else "NA" }} | {{ s.variants.denovo if s.variants.denovo is not none else "NA" }} |
{% endfor %}

{% endif %}

{% if comp_code in ["FLU1", "FLU2"] %}
For influenza components, evaluation focused on structural reporting metrics and concordance between metadata-reported and VCF-derived variant counts for each sample.
{% else %}

{% endif %}

{% set table_counter.value = table_counter.value + 1 %}
{% if comp_code in ["SARS1", "SARS2"] %}
**Table {{ table_counter.value }}. Metadata-reported and VCF-derived variant metrics for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

| Sample ID | Reporting mode | Metadata: variants (AF >=75%) | VCF-derived variants (AF >=75%) | Metadata: variants with effect | VCF-derived variants with effect | Metadata-VCF discrepancies | Effect discrepancies |
|---|---|---:|---:|---:|---:|---:|---:|
{% for collecting_lab_sample_id, s in comp.samples.items() -%}
| {{ collecting_lab_sample_id }} | {{ "High and low frequency" if s.variants.high_and_low_freq else ("High frequency only" if s.variants.high_freq_only else ("Low frequency only" if s.variants.low_freq_only else "NA")) }} | {{ s.variants.number_of_variants_in_consensus if s.variants and s.variants.number_of_variants_in_consensus is not none else "NA" }} | {{ s.variants.number_of_variants_in_consensus_vcf if s.variants and s.variants.number_of_variants_in_consensus_vcf is not none else "NA" }} | {{ s.variants.number_of_variants_with_effect if s.variants and s.variants.number_of_variants_with_effect is not none else "NA" }} | {{ s.variants.number_of_variants_with_effect_vcf if s.variants and s.variants.number_of_variants_with_effect_vcf is not none else "NA" }} | {{ s.variants.discrepancies_in_reported_variants if s.variants and s.variants.discrepancies_in_reported_variants is not none else "NA" }} | {{ s.variants.discrepancies_in_reported_variants_effect if s.variants and s.variants.discrepancies_in_reported_variants_effect is not none else "NA" }} |
{% endfor %}
{% else %}
**Table {{ table_counter.value }}. Metadata-reported and VCF-derived variant metrics for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

| Sample ID | Reporting mode | Metadata: variants (AF >=75%) | VCF-derived variants (AF >=75%) | Metadata: variants with effect | Metadata-VCF discrepancies | Total variants in VCF |
|---|---|---:|---:|---:|---:|---:|
{% for collecting_lab_sample_id, s in comp.samples.items() -%}
| {{ collecting_lab_sample_id }} | {{ "High and low frequency" if s.variants.high_and_low_freq else ("High frequency only" if s.variants.high_freq_only else ("Low frequency only" if s.variants.low_freq_only else "NA")) }} | {{ s.variants.number_of_variants_in_consensus if s.variants and s.variants.number_of_variants_in_consensus is not none else "NA" }} | {{ s.variants.number_of_variants_in_consensus_vcf if s.variants and s.variants.number_of_variants_in_consensus_vcf is not none else "NA" }} | {{ s.variants.number_of_variants_with_effect if s.variants and s.variants.number_of_variants_with_effect is not none else "NA" }} | {{ s.variants.discrepancies_in_reported_variants if s.variants and s.variants.discrepancies_in_reported_variants is not none else "NA" }} | {{ s.variants.number_of_variants_in_vcf if s.variants and s.variants.number_of_variants_in_vcf is not none else "NA" }} |
{% endfor %}
{% endif %}
{% set variant_metrics_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/variant_metrics_distribution.png" %}
{% if path_exists(variant_metrics_path) %}
{% set appendix_fig_counter.value = appendix_fig_counter.value + 1 %}
{% set lab_variant_metrics_figure_num = appendix_fig_counter.value %}
{% set _ = lab_variant_figure_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp": comp,
  "figure_num": lab_variant_metrics_figure_num,
  "variant_metrics_path": variant_metrics_path
}) %}
The distribution of metadata-reported and VCF-derived variant metrics across participating laboratories for **{{ labdata.lab.lab_cod }}** is shown in Appendix Figure {{ lab_variant_metrics_figure_num }}.
{% endif %}
{% endif %}

#### 9.{{ loop.index + 1 }}.3. Lineage, Subtype and Clade Assignment

Lineage/type and clade assignments submitted by **{{ labdata.lab.lab_cod }}** were compared against the curated gold standard classifications for each sample included in the {{ comp_code }} component.
{% set classification_concordance_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/classification_dimension_concordance.png" %}
{% if path_exists(classification_concordance_path) %}
{% set appendix_fig_counter.value = appendix_fig_counter.value + 1 %}
{% set lab_classification_figure_num = appendix_fig_counter.value %}
{% set _ = lab_classification_figure_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp": comp,
  "figure_num": lab_classification_figure_num,
  "classification_concordance_path": classification_concordance_path
}) %}
The distribution of lineage/type and clade classification outcomes across participating laboratories for **{{ labdata.lab.lab_cod }}** is shown in Appendix Figure {{ lab_classification_figure_num }}.
{% endif %}

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Per-sample lineage/type and clade assignment results for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

| Sample ID | Expected lineage/type | Reported lineage/type | Expected clade | Reported clade | Number of matches | Number of discrepancies |
|---|---|---|---|---|---|---|
{% for collecting_lab_sample_id, s in comp.samples.items() -%}
| {{ collecting_lab_sample_id }} | {{ s.classification.expected_lineage }} | {{ s.classification.lineage_assignment }} | {{ s.classification.expected_clade }} | {{ s.classification.clade_assignment }} | {{ s.classification.number_matches }} | {{ s.classification.number_discrepancies }} |
{% endfor %}

#### 9.{{ loop.index + 1 }}.4. Pipeline Benchmarking and Comparative Performance

The analytical workflow declared by **{{ labdata.lab.lab_cod }}** was benchmarked against other workflows implemented across the RELECOV network for the {{ comp_code }} component.

Positioning was evaluated based on four primary performance indicators:

1. Total number of discrepancies
2. Median consensus genome identity relative to the curated gold standard.
3. Total number of lineage/type and clade classification matches.
4. Metadata completeness

{% set workflow_positioning_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/workflow_positioning_boxplots.png" %}
{% if path_exists(workflow_positioning_path) %}
{% set appendix_fig_counter.value = appendix_fig_counter.value + 1 %}
{% set lab_workflow_figure_num = appendix_fig_counter.value %}
{% set _ = lab_workflow_figure_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp": comp,
  "figure_num": lab_workflow_figure_num,
  "workflow_positioning_path": workflow_positioning_path
}) %}
The workflow positioning across the RELECOV network for **{{ labdata.lab.lab_cod }}** is shown in Appendix Figure {{ lab_workflow_figure_num }}.
{% endif %}

Table {{ table_counter.value }} summarises the software configuration declared by **{{ labdata.lab.lab_cod }}** for each analysed sample in {{ comp_code }}. Table {{ table_counter.value + 1 }} contextualises the performance of the declared workflow relative to aggregated network-level metrics. For all four indicators, the reported network median and min-max range correspond to laboratory-level summaries across participating laboratories for the same component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Declared workflow configuration for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

{% if comp_code in ["SARS1", "SARS2"] %}
| Sample ID | Bioinformatics protocol | Dehosting | Pre-processing | Mapping/Assembly | Variant calling | Consensus sequence | Lineage assignment | Clade assignment |
|---|---|---|---|---|---|---|---|---|
{% for collecting_lab_sample_id, s in comp.samples.items() -%}
{% set sb = s.software_benchmarking %}
| {{ collecting_lab_sample_id }} | {{ software_label(sb.bioinformatics_protocol_software_name, sb.bioinformatics_protocol_software_version) }} | {{ software_label(sb.dehosting_method_software_name, sb.dehosting_method_software_version) }} | {{ software_label(sb.preprocessing_software_name, sb.preprocessing_software_version) }} | {% if sb.mapping_software_name %}{{ software_label(sb.mapping_software_name, sb.mapping_software_version) }}{% elif sb.assembly %}{{ software_label(sb.assembly, sb.assembly_version) }}{% else %}NA{% endif %} | {{ software_label(sb.variant_calling_software_name, sb.variant_calling_software_version) }} | {{ software_label(sb.consensus_sequence_software_name, sb.consensus_sequence_software_version) }} | {{ software_label(sb.lineage_assignment_software_name, sb.lineage_assignment_software_version, sb.lineage_assignment_database_version) }} | {{ software_label(sb.clade_assignment_software_name, sb.clade_assignment_software_version, sb.clade_assignment_software_database_version) }} |
{% endfor %}
{% else %}
| Sample ID | Bioinformatics protocol | Dehosting | Pre-processing | Mapping/Assembly | Variant calling | Consensus sequence | Type assignment | Subtype assignment | Clade assignment |
|---|---|---|---|---|---|---|---|---|---|
{% for collecting_lab_sample_id, s in comp.samples.items() -%}
{% set sb = s.software_benchmarking %}
| {{ collecting_lab_sample_id }} | {{ software_label(sb.bioinformatics_protocol_software_name, sb.bioinformatics_protocol_software_version) }} | {{ software_label(sb.dehosting_method_software_name, sb.dehosting_method_software_version) }} | {{ software_label(sb.preprocessing_software_name, sb.preprocessing_software_version) }} | {% if sb.mapping_software_name %}{{ software_label(sb.mapping_software_name, sb.mapping_software_version) }}{% elif sb.assembly %}{{ software_label(sb.assembly, sb.assembly_version) }}{% else %}NA{% endif %} | {{ software_label(sb.variant_calling_software_name, sb.variant_calling_software_version) }} | {{ software_label(sb.consensus_sequence_software_name, sb.consensus_sequence_software_version) }} | {{ software_label(sb.type_assignment_software_name, sb.type_assignment_software_version, sb.type_assignment_software_database_version) }} | {{ software_label(sb.subtype_assignment_software_name, sb.subtype_assignment_software_version, sb.subtype_assignment_software_database_version) }} | {{ software_label(sb.clade_assignment_software_name, sb.clade_assignment_software_version, sb.clade_assignment_software_database_version) }} |
{% endfor %}
{% endif %}

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Workflow performance positioning for {{ labdata.lab.lab_cod }} within the network ({{ comp_code }}).**

| Metric | {{ labdata.lab.lab_cod }} workflow | Network median | Network min - max |
|---|---:|---:|---:|
| Total number of discrepancies | {{ comp.total_number_discrepancies_consensus if comp.total_number_discrepancies_consensus is not none else "NA" }} | {{ general.components[comp_code].workflow_total_discrepancies_median if general.components[comp_code].workflow_total_discrepancies_median is not none else "NA" }} | {{ general.components[comp_code].workflow_total_discrepancies_min if general.components[comp_code].workflow_total_discrepancies_min is not none else "NA" }} - {{ general.components[comp_code].workflow_total_discrepancies_max if general.components[comp_code].workflow_total_discrepancies_max is not none else "NA" }} |
| Median genome identity (%) | {{ pct(comp.median_genome_identity_pct, 4) if comp.median_genome_identity_pct is not none else "NA" }} | {{ pct(general.components[comp_code].workflow_median_identity_pct_median, 4) if general.components[comp_code].workflow_median_identity_pct_median is not none else "NA" }} | {{ pct(general.components[comp_code].workflow_median_identity_pct_min, 4) if general.components[comp_code].workflow_median_identity_pct_min is not none else "NA" }} - {{ pct(general.components[comp_code].workflow_median_identity_pct_max, 4) if general.components[comp_code].workflow_median_identity_pct_max is not none else "NA" }} |
| Total classification matches | {{ comp.total_classification_matches if comp.total_classification_matches is not none else "NA" }} | {{ general.components[comp_code].typing.total_classification_matches_median if general.components[comp_code].typing.total_classification_matches_median is not none else "NA" }} | {{ general.components[comp_code].typing.total_classification_matches_min if general.components[comp_code].typing.total_classification_matches_min is not none else "NA" }} - {{ general.components[comp_code].typing.total_classification_matches_max if general.components[comp_code].typing.total_classification_matches_max is not none else "NA" }}|
| Metadata completeness (%) | {{ pct(comp.metadata.completeness_pct, 2) if comp.metadata.completeness_pct is not none else "NA" }} | {{ pct(general.components[comp_code].metadata_completeness_median, 2) if general.components[comp_code].metadata_completeness_median is not none else "NA" }} | {{ general.components[comp_code].metadata_completeness_min_pct if general.components[comp_code].metadata_completeness_min_pct is not none else "NA" }} - {{ general.components[comp_code].metadata_completeness_max_pct if general.components[comp_code].metadata_completeness_max_pct is not none else "NA" }} |

#### 9.{{ loop.index + 1 }}.5. Metadata-Derived Analytical Metrics (per sample)

This section summarises selected quantitative analytical metrics declared in the metadata submission of **{{ labdata.lab.lab_cod }}**, disaggregated by sample within the {{ comp_code }} component.

Only metrics explicitly provided by the laboratory are included in the comparative assessment. Because laboratories may not complete all quantitative metadata fields for every sample, tables and panels below include only those metrics that were actually reported by **{{ labdata.lab.lab_cod }}**. Network-level medians and (min-max) ranges are shown for contextual interpretation.

##### Sample Quality Control Assessment

{{ labdata.lab.lab_cod }} QC evaluations (Pass/Fail) were compared against the predefined gold standard QC status for each sample within {{ comp_code }}. Samples without a laboratory-reported QC assessment are shown as `NA` in the table and are omitted from the comparative figure.
{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Sample-level QC assessment for {{ labdata.lab.lab_cod }} ({{ comp_code }}), benchmarked against network-level QC concordance.**

| Sample ID | Reported QC | Gold standard QC | Network % Match |
|---|---|---|---:|
{% for collecting_lab_sample_id, s in comp.samples.items() -%}
{% set ns = (general.components[comp_code].qc.samples | selectattr("collecting_lab_sample_id","equalto",collecting_lab_sample_id) | list | first) -%}
| {{ collecting_lab_sample_id }} | {{ s.qc_test if s.qc_test is not none else "NA" }} | {{ ns.gold_standard_qc if ns else "NA" }} | {{ pct(ns.reported_match_rate_pct) if ns and ns.reported_match_rate_pct is not none else "NA" }} |
{% endfor %}

{% set qc_tests_reported = (comp.samples.values() | selectattr("qc_test", "ne", none) | list | length) > 0 %}
{% set qc_match_rate_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/qc_match_rate.png" %}
{% if qc_tests_reported and path_exists(qc_match_rate_path) %}
{% set appendix_fig_counter.value = appendix_fig_counter.value + 1 %}
{% set lab_qc_figure_num = appendix_fig_counter.value %}
{% set _ = lab_qc_figure_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp": comp,
  "figure_num": lab_qc_figure_num,
  "qc_match_rate_path": qc_match_rate_path
}) %}

The sample-level QC concordance across the network for **{{ labdata.lab.lab_cod }}** is shown in Appendix Figure {{ lab_qc_figure_num }}.
{% else %}

No comparative QC concordance figure is shown for {{ comp_code }} because **{{ labdata.lab.lab_cod }}** did not report any sample-level QC assessment for this component.
{% endif %}

##### Other metrics

Additional metadata-derived analytical metrics were available for a subset of {{ comp_code }} samples, including genome coverage above 10x, mean depth of coverage, proportion of Ns, and the fraction of viral and host reads where reported. The comparative figure below summarises how the values reported by **{{ labdata.lab.lab_cod }}** relate to the network-wide distribution, while the full per-sample tables are provided in the appendix.
{% set metadata_metrics_reported = namespace(count=0) %}
{% set metadata_metrics_appendix_summary = namespace(first_table_num=None, last_table_num=None, first_sample_id=None, last_sample_id=None, sample_count=0) %}
{% for collecting_lab_sample_id, s in comp.samples.items() -%}
{% set m = s.metadata_metrics -%}
{% if m -%}
{% set ns = (general.components[comp_code].metadata_metrics.samples | selectattr("sample_id","equalto",collecting_lab_sample_id) | list | first) -%}
{% set sample_metrics = namespace(count=0) -%}
{% for metric_key in metadata_metric_labels.keys() -%}
{% if m.get(metric_key) is not none -%}
{% set sample_metrics.count = sample_metrics.count + 1 -%}
{% endif -%}
{% endfor -%}
{% if sample_metrics.count > 0 -%}
{% set metadata_metrics_reported.count = metadata_metrics_reported.count + sample_metrics.count -%}
{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set metadata_metrics_table_num = appendix_table_counter.value %}
{% set _ = lab_metadata_metrics_appendix_entries.value.append({
  "comp_code": comp_code,
  "comp": comp,
  "collecting_lab_sample_id": collecting_lab_sample_id,
  "sample": s,
  "network_sample_metrics": ns,
  "table_num": metadata_metrics_table_num
}) %}
{% if metadata_metrics_appendix_summary.first_table_num is none %}{% set metadata_metrics_appendix_summary.first_table_num = metadata_metrics_table_num %}{% endif %}
{% if metadata_metrics_appendix_summary.first_sample_id is none %}{% set metadata_metrics_appendix_summary.first_sample_id = collecting_lab_sample_id %}{% endif %}
{% set metadata_metrics_appendix_summary.last_table_num = metadata_metrics_table_num %}
{% set metadata_metrics_appendix_summary.last_sample_id = collecting_lab_sample_id %}
{% set metadata_metrics_appendix_summary.sample_count = metadata_metrics_appendix_summary.sample_count + 1 %}

{% endif %}
{% endif %}
{% endfor %}
{% if metadata_metrics_appendix_summary.sample_count > 0 %}
Appendix Table{% if metadata_metrics_appendix_summary.sample_count > 1 %}s{% endif %} {{ metadata_metrics_appendix_summary.first_table_num }}{% if metadata_metrics_appendix_summary.last_table_num != metadata_metrics_appendix_summary.first_table_num %}–{{ metadata_metrics_appendix_summary.last_table_num }}{% endif %} report the metadata-derived analytical metrics for sample{% if metadata_metrics_appendix_summary.sample_count > 1 %}s{% endif %} **{{ metadata_metrics_appendix_summary.first_sample_id }}{% if metadata_metrics_appendix_summary.last_sample_id != metadata_metrics_appendix_summary.first_sample_id %}–{{ metadata_metrics_appendix_summary.last_sample_id }}{% endif %}**.
{% endif %}

{% set metadata_metrics_panel_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/metadata_metrics_panel.png" %}
{% if metadata_metrics_reported.count > 0 and path_exists(metadata_metrics_panel_path) %}
{% set fig_counter.value = fig_counter.value + 1 %}

{% set figure_cfg.style = "max-width: 98%;" %}
{{ render_figure(
  metadata_metrics_panel_path,
  comp_code ~ ": distribution of metadata-derived analytical metrics across the network per sample; black diamond indicates " ~ labdata.lab.lab_cod ~ ".",
  has_panels=True
) }}

**Figure {{ fig_counter.value }}. Distribution of metadata-derived analytical metrics across participating laboratories ({{ comp_code }}).**
Panel A shows genome coverage above 10x, Panel B depth of coverage, Panel C proportion of Ns, Panel D viral reads, and Panel E host reads. Only metrics actually reported by **{{ labdata.lab.lab_cod }}** are shown, so only panels with evaluable data are displayed. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. The black diamond corresponds to the values reported by **{{ labdata.lab.lab_cod }}**.
{% else %}

No comparative metadata-derived analytical metrics figure is shown for {{ comp_code }} because **{{ labdata.lab.lab_cod }}** did not report any evaluable quantitative metadata metrics for this component.
{% endif %}

{% endfor %}

### Acknowledgement

We sincerely thank **{{ labdata.lab.lab_cod }}** for its participation in the 2026 RELECOV Dry-Lab EQA. The contribution of each laboratory is fundamental to maintaining analytical comparability, reproducibility, and interoperability across the network.

For any questions, technical clarifications, or follow-up discussions regarding this report, please contact the RELECOV WP.6 coordination team at [bioinformatica@isciii.es](mailto:bioinformatica@isciii.es).
{% endif %}

## Appendix

This appendix is reserved for supplementary material that may support interpretation of the report but is not essential to the main narrative. Additional figures, extended tables, sensitivity analyses, or other secondary outputs can be included here when relevant.

{# Use `appendix_fig_counter` and `appendix_table_counter` for supplementary material moved here.
   Refer to them from the main text as "Appendix Figure X" and "Appendix Table X". #}

{% set appendix_components = [
  ("SARS1", "SARS-CoV-2, Illumina"),
  ("SARS2", "SARS-CoV-2, Oxford Nanopore Technologies"),
  ("FLU1", "Influenza virus, Illumina"),
  ("FLU2", "Influenza virus, Oxford Nanopore Technologies")
] %}

{% for appendix_comp_code, appendix_comp_name in appendix_components %}
{% set appendix_ns = namespace(has_material=false) %}
{% for entry in consensus_appendix_entries.value %}{% if entry.comp_code == appendix_comp_code %}{% set appendix_ns.has_material = true %}{% endif %}{% endfor %}
{% for entry in variant_sars_appendix_entries.value %}{% if entry.comp_code == appendix_comp_code %}{% set appendix_ns.has_material = true %}{% endif %}{% endfor %}
{% for entry in variant_flu_appendix_entries.value %}{% if entry.comp_code == appendix_comp_code %}{% set appendix_ns.has_material = true %}{% endif %}{% endfor %}
{% for entry in classification_appendix_entries.value %}{% if entry.comp_code == appendix_comp_code %}{% set appendix_ns.has_material = true %}{% endif %}{% endfor %}
{% for entry in qc_appendix_entries.value %}{% if entry.comp_code == appendix_comp_code %}{% set appendix_ns.has_material = true %}{% endif %}{% endfor %}
{% for entry in benchmark_appendix_entries.value %}{% if entry.comp_code == appendix_comp_code %}{% set appendix_ns.has_material = true %}{% endif %}{% endfor %}

{% if appendix_ns.has_material %}
### {{ appendix_comp_code }} ({{ appendix_comp_name }})

{% for entry in consensus_appendix_entries.value %}
{% if entry.comp_code == appendix_comp_code %}
#### Consensus Genome Reconstruction Supplementary Material

**Appendix Table {{ entry.metrics_table_num }}. Network-level consensus reconstruction metrics per sample for {{ entry.comp_code }}.**

| Sample ID | Median genome identity (%) | Median discrepancies | Discrepancies min-max |
|---|---:|---:|---:|
{% for s in entry.comp_net.consensus.samples %}
| {{ s.collecting_lab_sample_id }} | {{ "%.2f"|format(s.median_identity_pct) }} | {{ s.median_discrepancies }} | {{ s.min }} – {{ s.max }} |
{% endfor %}

**Appendix Table {{ entry.sample_table_num }}. Network-level consensus discrepancy types per sample for {{ entry.comp_code }}.**

| Sample ID | Median of Wrong nucleotide | Median Ambiguity instead of nucleotide | Median Nucleotide instead of ambiguity | Median Stretch of Ns instead of nucleotide stretch | Median Nucleotide stretch instead of stretch of Ns | Median Insertion relative to gold standard | Median Deletion relative to gold standard |
|---|---:|---:|---:|---:|---:|---:|---:|
{% for s in entry.comp_net.consensus.samples %}
| {{ s.collecting_lab_sample_id }} | {{ s.wrong_nt }} | {{ s.ambiguity2nt }} | {{ s.nt2ambiguity }} | {{ s.ns2nt }} | {{ s.nt2ns }} | {{ s.insertions }} | {{ s.deletions }} |
{% endfor %}

**Appendix Table {{ entry.type_table_num }}. Network-level discrepancy composition by type for {{ entry.comp_code }}.**

| Discrepancy type | Network median per sample | Min-max occurrencies |
|---|---:|---:|
| Incorrect nucleotide | {{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.wrong_nt.median) }} | {{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.wrong_nt.min) }}–{{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.wrong_nt.max) }} |
| Ambiguity instead of nucleotide | {{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.ambiguity2nt.median) }} | {{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.ambiguity2nt.min) }}–{{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.ambiguity2nt.max) }} |
| Nucleotide instead of ambiguity | {{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.nt2ambiguity.median) }} | {{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.nt2ambiguity.min) }}–{{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.nt2ambiguity.max) }} |
| Stretch of Ns instead of nucleotide | {{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.ns2nt.median) }} | {{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.ns2nt.min) }}–{{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.ns2nt.max) }} |
| Nucleotide stretch instead of stretch of Ns | {{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.nt2ns.median) }} | {{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.nt2ns.min) }}–{{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.nt2ns.max) }} |
| Insertion relative to gold standard | {{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.insertions.median) }} | {{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.insertions.min) }}–{{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.insertions.max) }} |
| Deletion relative to gold standard | {{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.deletions.median) }} | {{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.deletions.min) }}–{{ "%.0f"|format(entry.comp_net.consensus.discrepancy_breakdown.deletions.max) }} |

Figure {{ entry.type_fig_num }} in the appendix summarises the contribution of each discrepancy category observed in {{ entry.comp_code }} relative to the curated gold standard.

{% set figure_cfg.style = "max-width: 96%;" %}
{{ render_figure(
  entry.comp_net.consensus.fig_discrepancy_type_boxplot,
  "Composition of consensus discrepancy types for " ~ entry.comp_code ~ " relative to the curated gold standard."
) }}

**Appendix Figure {{ entry.type_fig_num }}. Composition of consensus discrepancy types relative to the curated gold standard for {{ entry.comp_code }}.** Boxplots represent aggregated discrepancies across all submitted consensus sequences, stratified by discrepancy category. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.
{% endif %}
{% endfor %}

{% for entry in variant_sars_appendix_entries.value %}
{% if entry.comp_code == appendix_comp_code %}
<h4 class="appendix-landscape-heading">Variant Detection Accuracy Supplementary Material</h4>

**Appendix Table {{ entry.profile_table_num }}. Network-level SARS-CoV-2 variant calling profile per sample for {{ entry.comp_code }}.** The discrepancy-type columns correspond to the median count per sample across participating laboratories.

| Sample ID | Median successful hits | Median discrepancies | Discrepancies min-max | Median wrong nucleotide | Median insertions | Median deletions | Median missing | Median _de novo_ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
{% for s in entry.comp_net.variant.samples %}
| {{ s.collecting_lab_sample_id }} | {{ s.median_successful_hits if s.median_successful_hits is not none else "NA" }} | {{ s.median_discrepancies if s.median_discrepancies is not none else "NA" }} | {{ s.min if s.min is not none else "NA" }} – {{ s.max if s.max is not none else "NA" }} | {{ s.wrong_nt if s.wrong_nt is not none else "NA" }} | {{ s.insertions if s.insertions is not none else "NA" }} | {{ s.deletions if s.deletions is not none else "NA" }} | {{ s.missing if s.missing is not none else "NA" }} | {{ s.denovo if s.denovo is not none else "NA" }} |
{% endfor %}

**Appendix Table {{ entry.type_table_num }}. Network-level discrepancy composition by type for {{ entry.comp_code }}.** The discrepancy-type columns correspond to the median count per sample across participating laboratories.

| Discrepancy type | Network median per sample | Network min-max per sample |
|---|---:|---:|
| Incorrect nucleotide | {{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.wrong_nt.median) }} | {{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.wrong_nt.min) }}–{{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.wrong_nt.max) }} |
| Insertion relative to gold standard | {{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.insertions.median) }} | {{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.insertions.min) }}–{{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.insertions.max) }} |
| Deletions relative to gold standard | {{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.deletions.median) }} | {{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.deletions.min) }}–{{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.deletions.max) }} |
| Missing expected variants | {{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.missing.median) }} | {{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.missing.min) }}–{{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.missing.max) }} |
| De novo variants | {{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.denovo.median) }} | {{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.denovo.min) }}–{{ "%.0f"|format(entry.comp_net.variant.discrepancy_breakdown.denovo.max) }} |

Figure {{ entry.type_fig_num }} in the appendix summarises the contribution of each discrepancy category observed in {{ entry.comp_code }} relative to the curated gold standard.

{% set figure_cfg.style = "max-width: 90%;" %}
{{ render_figure(
  entry.comp_net.variant.fig_discrepancy_type_boxplot,
  "Composition of variant discrepancy types for " ~ entry.comp_code ~ " relative to the curated gold standard."
) }}

**Appendix Figure {{ entry.type_fig_num }}. Composition of variant discrepancy types relative to the curated gold standard for {{ entry.comp_code }}.** Boxplots represent aggregated discrepancies across all submitted variant calls, stratified by discrepancy category (incorrect nucleotide, excess ambiguous bases, and indels). Where required, a broken y-axis is used to preserve visual resolution in the lower discrepancy range while still displaying higher values above an empty interval. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.
{% endif %}
{% endfor %}

{% for entry in variant_flu_appendix_entries.value %}
{% if entry.comp_code == appendix_comp_code %}
<h4 class="appendix-landscape-heading">Variant Detection Accuracy Supplementary Material</h4>

**Appendix Table {{ entry.aggregated_table_num }}. Aggregated influenza variant reporting metrics for {{ entry.comp_code }}.**

| Metric | Network median | Network min-max |
|---|---:|---:|
| Variants >=75% AF in metadata | {{ entry.comp_net.variant.median_variants_in_consensus if entry.comp_net.variant.median_variants_in_consensus is not none else "NA" }} | {{ entry.comp_net.variant.min_variants_in_consensus if entry.comp_net.variant.min_variants_in_consensus is not none else "NA" }}–{{ entry.comp_net.variant.max_variants_in_consensus if entry.comp_net.variant.max_variants_in_consensus is not none else "NA" }} |
| Variants >=75% AF in VCF | {{ entry.comp_net.variant.median_variants_in_consensus_vcf if entry.comp_net.variant.median_variants_in_consensus_vcf is not none else "NA" }} | {{ entry.comp_net.variant.min_variants_in_consensus_vcf if entry.comp_net.variant.min_variants_in_consensus_vcf is not none else "NA" }}–{{ entry.comp_net.variant.max_variants_in_consensus_vcf if entry.comp_net.variant.max_variants_in_consensus_vcf is not none else "NA" }} |
| Discrepancies between metadata and VCF | {{ entry.comp_net.variant.median_discrepancies_in_reported_variants if entry.comp_net.variant.median_discrepancies_in_reported_variants is not none else "NA" }} | {{ entry.comp_net.variant.min_discrepancies_in_reported_variants if entry.comp_net.variant.min_discrepancies_in_reported_variants is not none else "NA" }}–{{ entry.comp_net.variant.max_discrepancies_in_reported_variants if entry.comp_net.variant.max_discrepancies_in_reported_variants is not none else "NA" }} |
| Total variants in VCF (n={{ entry.comp_net.variant_metadata_reporting.total_variants_in_vcf_reported_n_labs }}) | {{ entry.comp_net.variant.median_variants_in_vcf if entry.comp_net.variant.median_variants_in_vcf is not none else "NA" }} | {{ entry.comp_net.variant.min_variants_in_vcf if entry.comp_net.variant.min_variants_in_vcf is not none else "NA" }}–{{ entry.comp_net.variant.max_variants_in_vcf if entry.comp_net.variant.max_variants_in_vcf is not none else "NA" }} |
{% endif %}
{% endfor %}

{% for entry in classification_appendix_entries.value %}
{% if entry.comp_code == appendix_comp_code %}
#### Lineage, Subtype and Clade Assignment Supplementary Material

**Appendix Table {{ entry.table_num }}. Network-level classification outcomes per sample for {{ entry.comp_code }}.**

| Sample ID | Lineage/Subtype matches (%) | Clade matches (%) |
|---|---:|---:|
{% for s in entry.comp_net.typing.samples %}
| {{ s.collecting_lab_sample_id }} | {{ "%.2f"|format(s.lineage_hit_pct) }} | {{ "%.2f"|format(s.clade_hit_pct) }} |
{% endfor %}
{% endif %}
{% endfor %}

{% for entry in qc_appendix_entries.value %}
{% if entry.comp_code == appendix_comp_code %}
#### Sample Quality Control Assessment Supplementary Material

**Appendix Table {{ entry.table_num }}. Sample-level QC concordance for {{ entry.comp_code }} for reported QC classification.**

| Sample ID | Gold standard QC | % Match | # Matches | # Discrepancies | Total evaluations |
|---|---:|---:|---:|---:|---:|
{% for s in entry.comp_net.qc.samples %}
| {{ s.collecting_lab_sample_id }} | {{ s.gold_standard_qc }} | {{ pct(s.reported_match_rate_pct) }} | {{ s.matches }} | {{ s.discrepancies }} | {{ s.total_evaluations }} |
{% endfor %}
{% endif %}
{% endfor %}

{% set benchmark_heading_ns = namespace(shown=false) %}
{% for entry in benchmark_appendix_entries.value %}
{% if entry.comp_code == appendix_comp_code %}
{% if not benchmark_heading_ns.shown %}
{% set benchmark_heading_ns.shown = true %}
<h4 class="appendix-landscape-heading">Pipeline Benchmarking and Comparative Performance Supplementary Material</h4>
{% endif %}

{% if entry.kind == "dehosting" %}
<h5 class="appendix-landscape-heading">De-hosting</h5>
**Appendix Table {{ entry.table_num }}. Performance summary of declared de-hosting software for {{ entry.comp_code }}.**

| De-hosting software | Version | N labs | Median % host reads |
|---|---:|---:|---:|
{% for p in entry.comp_net.benchmarking.dehosting.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.per_reads_host }} |
{% endfor %}
{% elif entry.kind == "preprocessing" %}
<h5 class="appendix-landscape-heading">Pre-processing</h5>
**Appendix Table {{ entry.table_num }}. Performance summary of declared pre-processing software configurations for {{ entry.comp_code }}.**

| Pre-processing software | Version | N labs | Most common configuration | Median number of reads sequenced | Median reads passing filters |
|---|---:|---:|---:|---:|---:|
{% for p in entry.comp_net.benchmarking.preprocessing.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.params|mdcell }} | {{ p.number_of_reads_sequenced }} | {{ p.pass_reads }} |
{% endfor %}
{% elif entry.kind == "mapping" %}
<h5 class="appendix-landscape-heading">Mapping</h5>
**Appendix Table {{ entry.table_num }}. Performance summary of declared mapping software configurations for {{ entry.comp_code }}.**

| Mapping software | Version | N labs | Most common configuration | Most common depth of coverage threshold | Median % reads virus |
|---|---:|---:|---:|---:|---:|
{% for p in entry.comp_net.benchmarking.mapping.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.params|mdcell }} | {{ p.depth_of_coverage_threshold if p.depth_of_coverage_threshold is not none else "N/A" }} | {{ p.per_reads_virus if p.per_reads_virus is not none else "N/A" }} |
{% endfor %}
{% elif entry.kind == "assembly" %}
<h5 class="appendix-landscape-heading">Assembly</h5>
**Appendix Table {{ entry.table_num }}. Performance summary of declared assembly software configurations for {{ entry.comp_code }}.**

| Assembly software | Version | N labs | Most common configuration | Median consensus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
{% for p in entry.comp_net.benchmarking.assembly.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.params|mdcell }} | {{ p.consensus_genome_length }} | {{ p.median_identity_pct }} | {{ p.median_discrepancies }} |
{% endfor %}
{% elif entry.kind == "consensus_software" %}
<h5 class="appendix-landscape-heading">Consensus software</h5>
**Appendix Table {{ entry.table_num }}. Performance summary of declared consensus software configurations for {{ entry.comp_code }}.**

| Consensus software | Version | N labs | Most common configuration | Median consensus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
{% for p in entry.comp_net.benchmarking.consensus_software.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.params|mdcell }} | {{ p.consensus_genome_length }} | {{ p.median_identity_pct }} | {{ p.median_discrepancies }} |
{% endfor %}
{% elif entry.kind == "variant_calling" %}
<h5 class="appendix-landscape-heading">Variant calling</h5>
**Appendix Table {{ entry.table_num }}. Performance summary of declared variant calling software configurations for {{ entry.comp_code }}.**

{% if entry.comp_code[:3] == "FLU" %}
| Variant calling software | Version | N labs | Most common configuration | Median high and low frequency (%) | Median high frequency only (%) | Median low frequency only (%) | Median variants (AF >=75%) | Median variants in VCF (AF >=75%) | Median variants with effect | Median metadata-VCF discrepancies | Median total variants in VCF |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{% for p in entry.comp_net.benchmarking.variant_calling.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.params|mdcell }} | {{ p.high_and_low_freq_pct }} | {{ p.high_freq_only_pct }} | {{ p.low_freq_only_pct }} | {{ p.number_of_variants_in_consensus }} | {{ p.number_of_variants_in_consensus_vcf }} | {{ p.number_of_variants_with_effect }} | {{ p.discrepancies_in_reported_variants }} | {{ p.number_of_variants_in_vcf }} |
{% endfor %}
{% else %}
| Variant calling software | Version | N labs | {{ "Model used" if entry.comp_code == "SARS2" else "Most common configuration" }} | Median high and low frequency (%) | Median high frequency only (%) | Median low frequency only (%) | Median variants (AF >=75%) | Median variants in VCF (AF >=75%) | Median variants with effect | Median variants with effect in VCF | Median metadata-VCF discrepancies | Median effect discrepancies | Median successful hits | Median total discrepancies |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{% for p in entry.comp_net.benchmarking.variant_calling.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ (p.model if entry.comp_code == "SARS2" else p.params)|mdcell }} | {{ p.high_and_low_freq_pct }} | {{ p.high_freq_only_pct }} | {{ p.low_freq_only_pct }} | {{ p.number_of_variants_in_consensus }} | {{ p.number_of_variants_in_consensus_vcf }} | {{ p.number_of_variants_with_effect }} | {{ p.number_of_variants_with_effect_vcf }} | {{ p.discrepancies_in_reported_variants }} | {{ p.discrepancies_in_reported_variants_effect }} | {{ p.successful_hits }} | {{ p.total_discrepancies }} |
{% endfor %}
{% endif %}
{% elif entry.kind == "clade_assignment" %}
<h5 class="appendix-landscape-heading">Clade assignment</h5>
**Appendix Table {{ entry.table_num }}. Performance summary of declared clade assignment software configurations for {{ entry.comp_code }}.**

| Clade assignment software | Version | N labs | Database version | % of clade match |
|---|---:|---:|---:|---:|
{% for p in entry.comp_net.benchmarking.clade_assignment.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.database_version|mdcell }} | {{ p.clade_hit_pct }} |
{% endfor %}
{% elif entry.kind == "lineage_assignment" %}
<h5 class="appendix-landscape-heading">Lineage assignment</h5>
**Appendix Table {{ entry.table_num }}. Performance summary of declared lineage assignment software configurations for {{ entry.comp_code }}.**

| Lineage Assignment software | Version | N labs | Database version | % of lineage match |
|---|---:|---:|---:|---:|
{% for p in entry.comp_net.benchmarking.lineage_assignment.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.database_version|mdcell }} | {{ p.lineage_hit_pct }} |
{% endfor %}
{% elif entry.kind == "type_assignment" %}
<h5 class="appendix-landscape-heading">Type assignment</h5>
**Appendix Table {{ entry.table_num }}. Performance summary of declared type assignment software configurations for {{ entry.comp_code }}.**

| Type Assignment software | Version | N labs | Database version | % of type match |
|---|---:|---:|---:|---:|
{% for p in entry.comp_net.benchmarking.type_assignment.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.database_version|mdcell }} | {{ p.type_hit_pct }} |
{% endfor %}
{% elif entry.kind == "subtype_assignment" %}
<h5 class="appendix-landscape-heading">Subtype assignment</h5>
**Appendix Table {{ entry.table_num }}. Performance summary of declared subtype assignment software configurations for {{ entry.comp_code }}.**

| Subtype Assignment software | Version | N labs | Database version | % of subtype match |
|---|---:|---:|---:|---:|
{% for p in entry.comp_net.benchmarking.subtype_assignment.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.database_version|mdcell }} | {{ p.subtype_hit_pct }} |
{% endfor %}
{% endif %}
{% endif %}
{% endfor %}

{% endif %}
{% endfor %}

{% if labdata %}
### Individual Laboratory Supplementary Material

{% for appendix_comp_code, appendix_comp_name in [
  ("SARS1", "SARS-CoV-2, Illumina"),
  ("SARS2", "SARS-CoV-2, Oxford Nanopore Technologies"),
  ("FLU1", "Influenza virus, Illumina"),
  ("FLU2", "Influenza virus, Oxford Nanopore Technologies")
] %}
{% set comp_consensus_entries = lab_consensus_appendix_entries.value | selectattr("comp_code", "equalto", appendix_comp_code) | list %}
{% set comp_variant_entries = lab_variant_figure_appendix_entries.value | selectattr("comp_code", "equalto", appendix_comp_code) | list %}
{% set comp_classification_entries = lab_classification_figure_appendix_entries.value | selectattr("comp_code", "equalto", appendix_comp_code) | list %}
{% set comp_qc_entries = lab_qc_figure_appendix_entries.value | selectattr("comp_code", "equalto", appendix_comp_code) | list %}
{% set comp_metadata_entries = lab_metadata_metrics_appendix_entries.value | selectattr("comp_code", "equalto", appendix_comp_code) | list %}
{% set comp_workflow_entries = lab_workflow_figure_appendix_entries.value | selectattr("comp_code", "equalto", appendix_comp_code) | list %}
{% if comp_consensus_entries or comp_variant_entries or comp_classification_entries or comp_qc_entries or comp_metadata_entries or comp_workflow_entries %}
#### {{ appendix_comp_code }} ({{ appendix_comp_name }})

{% if comp_consensus_entries %}
##### Consensus Genome Reconstruction Performance Supplementary Material

{% for entry in comp_consensus_entries %}
**Appendix Table {{ entry.metrics_table_num }}. Per-sample consensus reconstruction metrics for {{ labdata.lab.lab_cod }} ({{ entry.comp_code }}).**

| Sample ID | {{ labdata.lab.lab_cod }} Genome identity (%) | Network Genome Identity Median | {{ labdata.lab.lab_cod }} Total discrepancies | Network total discrepancies median |
|---|---:|---:|---:|---:|
{% for collecting_lab_sample_id, s in entry.comp.samples.items() -%}
| {{ collecting_lab_sample_id }} | {{ pct(s.consensus.genome_identity_pct, 4) }} | {{ general.components[entry.comp_code].consensus.samples[collecting_lab_sample_id].median_identity_pct }} | {{ s.consensus.total_discrepancies }} | {{ general.components[entry.comp_code].consensus.samples[collecting_lab_sample_id].median_discrepancies }} |
{% endfor %}

**Appendix Table {{ entry.breakdown_table_num }}. Discrepancy type breakdown per sample for {{ labdata.lab.lab_cod }} ({{ entry.comp_code }}).**

| Sample ID | Total wrong nucleotides | Total ambiguity instead of nucleotide | Total nucleotide instead of ambiguity | Total stretch of Ns instead of nucleotide stretch | Total nucleotide stretch instead of stretch of Ns | Total insertion relative to gold standard | Total deletion relative to gold standard |
|---|---:|---:|---:|---:|---:|---:|---:|
{% for collecting_lab_sample_id, s in entry.comp.samples.items() -%}
| {{ collecting_lab_sample_id }} | {{ s.consensus.discrepancy_breakdown.wrong_nt }} | {{ s.consensus.discrepancy_breakdown.ambiguity2nt }} | {{ s.consensus.discrepancy_breakdown.nt2ambiguity }} | {{ s.consensus.discrepancy_breakdown.ns2nt }} | {{ s.consensus.discrepancy_breakdown.nt2ns }} | {{ s.consensus.discrepancy_breakdown.insertions }} | {{ s.consensus.discrepancy_breakdown.deletions }} |
{% endfor %}
{% endfor %}
{% endif %}

{% if comp_variant_entries %}
##### Variant Detection Supplementary Material

{% for entry in comp_variant_entries %}
{% set figure_cfg.style = "max-width: 98%;" %}
{{ render_figure(
  entry.variant_metrics_path,
  entry.comp_code ~ ": distribution of variant reporting metrics across the network; black diamond indicates " ~ labdata.lab.lab_cod ~ ".",
  has_panels=True
) }}

{% if entry.comp_code in ["SARS1", "SARS2"] %}
**Appendix Figure {{ entry.figure_num }}. Metadata-reported and VCF-derived variant metrics across participating laboratories ({{ entry.comp_code }}).** Panel A shows reported variants with AF >=75%, Panel B reported variants with effect, Panel C variants in VCF with AF >=75%, Panel D variants with effect in VCF, Panel E metadata-VCF discrepancies for AF >=75% variants, and Panel F metadata-VCF discrepancies for variants with effect across the RELECOV network. Only panels with evaluable data for **{{ labdata.lab.lab_cod }}** are shown. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. The black diamond corresponds to the results obtained by **{{ labdata.lab.lab_cod }}**.
{% else %}
**Appendix Figure {{ entry.figure_num }}. Influenza-specific variant reporting metrics across participating laboratories ({{ entry.comp_code }}).** Panel A shows reported variants with AF >=75%, Panel B VCF-derived variants with AF >=75%, Panel C reported variants with effect, Panel D metadata-VCF discrepancies, and Panel E total variants present in the submitted VCF files across the RELECOV network. Only panels with evaluable data for **{{ labdata.lab.lab_cod }}** are shown. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. The black diamond corresponds to the results obtained by **{{ labdata.lab.lab_cod }}**.
{% endif %}
{% endfor %}
{% endif %}

{% if comp_classification_entries %}
##### Lineage, Subtype and Clade Assignment Supplementary Material

{% for entry in comp_classification_entries %}
{% set figure_cfg.style = "max-width: 98%;" %}
{{ render_figure(
  entry.classification_concordance_path,
  entry.comp_code ~ ": lineage/type and clade classification outcomes across the network; black diamond indicates " ~ labdata.lab.lab_cod ~ ".",
  has_panels=True
) }}

**Appendix Figure {{ entry.figure_num }}. Lineage/type and clade classification outcomes across participating laboratories ({{ entry.comp_code }}).** Panel A shows the proportion of Match, Discrepancy, and Not provided outcomes for lineage/type assignments across participating laboratories for each sample. Panel B shows the corresponding proportions for clade assignments. Stacked bars represent the percentage of laboratories with correct classifications, incorrect classifications, or missing classifications relative to the curated gold standard. The black diamond marks the result reported by **{{ labdata.lab.lab_cod }}**, positioned within the Match, Discrepancy, or Not provided segment for each sample.
{% endfor %}
{% endif %}

{% if comp_qc_entries %}
##### Sample Quality Control Assessment Supplementary Material

{% for entry in comp_qc_entries %}
{% set figure_cfg.style = "max-width: 80%;" %}
{{ render_figure(
  entry.qc_match_rate_path,
  entry.comp_code ~ ": sample-level QC concordance across the network, with " ~ labdata.lab.lab_cod ~ " highlighted."
) }}

**Appendix Figure {{ entry.figure_num }}. Sample-level QC concordance across the network for {{ entry.comp_code }}, with {{ labdata.lab.lab_cod }} highlighted.** Stacked bars represent the network-wide proportions of Match, Discrepancy, and Not provided outcomes relative to the gold standard for each sample. The black diamond indicates whether **{{ labdata.lab.lab_cod }}** reported a Match, a Discrepancy, or did not provide a QC assessment for the corresponding sample. Not provided values are shown separately and are not counted as discrepancies.
{% endfor %}
{% endif %}

{% if comp_metadata_entries %}
##### Metadata-Derived Analytical Metrics Supplementary Material

{% for entry in comp_metadata_entries %}
##### {{ entry.collecting_lab_sample_id }}

**Appendix Table {{ entry.table_num }}. Metadata-derived analytical metrics for {{ labdata.lab.lab_cod }} (component {{ entry.comp_code }}, sample {{ entry.collecting_lab_sample_id }}).**

| Metric | {{ labdata.lab.lab_cod }} | Network median | Network min - max |
|---|---:|---:|---:|
{% for metric_key, metric_label in metadata_metric_labels.items() -%}
| {{ metric_label }} | {{ entry.sample.metadata_metrics[metric_key] if entry.sample.metadata_metrics.get(metric_key) is not none else "NA" }} | {{ entry.network_sample_metrics[metric_key].median if entry.network_sample_metrics and entry.network_sample_metrics.get(metric_key) else "NA" }} | {{ entry.network_sample_metrics[metric_key].min if entry.network_sample_metrics and entry.network_sample_metrics.get(metric_key) else "NA" }} - {{ entry.network_sample_metrics[metric_key].max if entry.network_sample_metrics and entry.network_sample_metrics.get(metric_key) else "NA" }} |
{% endfor %}
{% endfor %}
{% endif %}

{% if comp_workflow_entries %}
##### Workflow Benchmarking Supplementary Material

{% for entry in comp_workflow_entries %}
{% set figure_cfg.style = "max-width: 98%;" %}
{{ render_figure(
  entry.workflow_positioning_path,
  entry.comp_code ~ ": workflow positioning relative to network-wide distributions, with " ~ labdata.lab.lab_cod ~ " highlighted by a black diamond.",
  has_panels=True
) }}

**Appendix Figure {{ entry.figure_num }}. Workflow positioning within the RELECOV network for {{ entry.comp_code }}.** Multi-panel boxplots summarise the laboratory-level distribution across the network for Panel A total consensus discrepancies, Panel B median genome identity, Panel C total classification matches, and Panel D metadata completeness. Only panels with evaluable data for **{{ labdata.lab.lab_cod }}** are shown. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. In Panel B, the y-axis is truncated to highlight differences among high-identity values. The black diamond corresponds to the results obtained by **{{ labdata.lab.lab_cod }}**.
{% endfor %}
{% endif %}

{% endif %}
{% endfor %}
{% endif %}
