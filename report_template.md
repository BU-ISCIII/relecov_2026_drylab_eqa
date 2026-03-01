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
{% macro render_figure(path, caption=None) -%}
  {% if path %}
    <figure>
      <img src="{{ path }}" alt="{{ caption|default('Figure') }}" style="max-width: 100%;"/>
      {% if caption %}<figcaption>{{ caption }}</figcaption>{% endif %}
    </figure>
  {% endif %}
{%- endmacro %}

{% macro network_iqr(iqr_list, decimals=2) -%}
{% set fig_counter = namespace(value=0) %}
{% set table_counter = namespace(value=0) %}

# RELECOV 2.0 - Consolidation of WGS and RT-PCR activities for SARS-CoV-2 in Spain towards sustainable use and integration of enhanced infrastructure and capacities in the RELECOV network</h1>
Authors: Sarai Varona Fernánez
Version: v0.0.1 (25/02/2026)

## Table of Contents

- [Executive Summary](#executive-summary)
- [1. Introduction](#1-introduction)
- [2. Scope of the EQA](#2-scope-of-the-eqa)
- [3. Dataset Design and Sample Selection Criteria](#3-dataset-design-and-sample-selection-criteria)
  - [3.1 Rationale for Dataset Selection](#31-rationale-for-dataset-selection)
  - [3.2 SARS-CoV-2 Dataset Selection](#32-sars-cov-2-dataset-selection)
  - [3.3 Influenza Dataset Selection](#33-influenza-dataset-selection)
    - [In-Silico Influenza Dataset Construction](#in-silico-influenza-dataset-construction)
- [4. Methodology of Evaluation](#4-methodology-of-evaluation)
  - [4.1 Evaluation of Consensus Genome Sequences](#41-evaluation-of-consensus-genome-sequences)
  - [4.2 Evaluation of Variant Detection](#42-evaluation-of-variant-detection)
  - [4.3 Evaluation of Lineage/Type and Clade Assignment](#43-evaluation-of-taxonomic-and-phylogenetic-classification)
  - [4.4 Evaluation of Metadata Completeness and Compliance](#44-evaluation-of-metadata-completeness-and-compliance)
  - [4.5 Pipeline Benchmarking](#45-pipeline-benchmarking)
- [5. General Results](#5-general-results)
  - [5.1 Submission Completeness](#51-submission-completeness)
  - [5.2 Consensus Genome Reconstruction Performance](#52-consensus-genome-reconstruction-performance)
  - [5.3 Variant Detection Accuracy](#53-variant-detection-accuracy)
  - [5.4 Lineage, Type and Clade Assignment](#54-lineage-type-and-clade-assignment)
  - [5.5 Metadata Quality and Interoperability](#55-metadata-quality-and-interoperability)
  - [5.6 Pipeline Benchmarking and Comparative Performance](#56-pipeline-benchmarking-and-comparative-performance)
- [6. Component-specific Results](#6-component-specific-results)
  - [6.1 SARS1 (SARS-CoV-2, Illumina)](#61-sars1-sars-cov-2-illumina)
  - [6.2 SARS2 (SARS-CoV-2, Nanopore)](#62-sars2-sars-cov-2-nanopore)
  - [6.3 FLU1 (Influenza, Illumina)](#63-flu1-influenza-illumina)
  - [6.4 FLU2 (Influenza, Nanopore)](#64-flu2-influenza-nanopore)
- [7. Discussion](#7-discussion)
- [8. Conclusions](#8-conclusions)
{% if labdata %}
- [9. Individual Laboratory Technical Report](#9-individual-laboratory-technical-report)
  - [9.1 Participation Overview](#91-participation-overview)
  - [9.2 SARS1 (SARS-CoV-2, Illumina)](#92-sars1-sars-cov-2-illumina)
  - [9.3 SARS2 (SARS-CoV-2, Nanopore)](#93-sars2-sars-cov-2-nanopore)
  - [9.4 FLU1 (Influenza, Illumina)](#94-flu1-influenza-illumina)
  - [9.5 FLU2 (Influenza, Nanopore)](#95-flu2-influenza-nanopore)
{ endif }

## Executive Summary</h2>
To be completed after final results are consolidated.

This EQA provides the first fully drylab benchmarking of bioinformatic workflows across the RELECOV network, integrating both internationally validated datasets and purpose-designed in-silico scenarios.

## 1. Introduction
The RELECOV Network aims to strengthen genomic surveillance of respiratory viruses by developing and harmonising analytical capacities across the participating laboratories. In this context, it is essential to **assess the consistency, reproducibility and maturity of the bioinformatics workflows implemented within the network**.

To this end, an **external quality assessment (EQA) exercise in dry lab format** will be conducted, inspired by the ECDC’s 2024 dry-lab EQA. The exercise will focus on the bioinformatic characterisation of respiratory viruses, testing key analytical tasks including viral genome reconstruction, variant identification, and lineage and clade assignment.

A central component of this initiative is to evaluate the range of analytical pipelines currently used in Relecov Network, identify their relative performance, and determine which pipeline is best suited for establishing a genomic surveillance standard within the network. This evaluation directly contributes to **Objective 3** of RELECOV 2.0, which focuses on *generating a comprehensive understanding of the analytical and operational workflows currently implemented across Spanish laboratories*. Furthermore, the EQA provides the practical evidence base required for Task T6.1, which aims to identify, compare and prioritise the analytical pipelines available nationally in order to define the workflow that should be integrated into the RELECOV analytical platform.

The exercise is also aligned with **Milestone M3.2**, which pertains to the *establishment of harmonised analytical procedures for respiratory viruses*. It also contributes to **Deliverable D4.1**, focused on the definition of minimum metadata requirements and harmonised reporting formats, by *testing the ability of laboratories to produce interoperable outputs suitable for integration into national and international surveillance systems*. In addition, the exercise informs the development of the RELECOV platform by providing operational insights critical to **Task T5.2**, which addresses *data ingestion, workflow automation and technical specifications for integrating pipeline outputs within the platform*.

Beyond workflow harmonisation, the EQA contributes to capacity-building and performance assessment activities central to the project. By benchmarking analytical performance across laboratories, the exercise provides essential input for **Milestone M6.3**, related to *determining laboratory readiness and identifying areas requiring technical reinforcement or training*. These insights also support **Deliverable D6.2**, which includes *recommendations for strengthening network-wide analytical capacity and ensuring long-term sustainability of genomic surveillance operations*.

The overall objective of the exercise is to **assess the bioinformatic performance of the participating laboratories, identify areas for improvement, and promote the adoption of consistent and comparable analytical practices across the network**. The outcomes will enhance RELECOV’s preparedness and response capacity in routine surveillance and public health emergencies, ensuring the quality and robustness of genomic analyses performed throughout the network and contributing directly to the fulfilment of key project objectives, milestones and deliverables.

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
- A completed harmonised metadata template, documenting analytical tools, software versions, reference genomes used, parameter settings, coverage thresholds, lineage, type or clade assignment tools, and file paths to submitted outputs.

The primary objective of the exercise was to assess the consistency, reproducibility, and comparability of bioinformatic workflows currently implemented across the network. The evaluation focused on core analytical tasks that are essential for routine genomic surveillance and public health response, including:

- **Viral genome reconstruction**: Generation of high-quality consensus genome sequences from raw sequencing reads produced using Illumina and Oxford Nanopore Technologies platforms.
- **Variant identification and reporting**: Detection and annotation of nucleotide variants relative to a chosen reference genome, including evaluation of filtering criteria, allele frequency thresholds, and variant file standardisation.
- **Lineage, type and clade assignment**: Accurate classification of reconstructed genomes using established nomenclature systems and version-controlled databases.
- **Metadata reporting and interoperability**: Completion of a harmonised metadata template capturing software versions, analytical parameters, reference genome selection, and file traceability, ensuring compatibility with automated validation and integration into the RELECOV analytical platform.

## 3. Dataset Design and Sample Selection Criteria
### 3.1 Rationale for Dataset Selection

The 2026 RELECOV Dry-Lab EQA was specifically designed for laboratories operating in a clinical and hospital-based diagnostic context, where routine genomic surveillance primarily involves human respiratory samples.

Sample selection followed three guiding principles:

- Representation of realistic genomic surveillance scenarios.
- Inclusion of defined analytical challenges.
- Ensuring methodological benchmarking robustness.

Datasets were derived from two sources:

- Reused datasets from the 2024 ECDC ESIB Dry-Lab EQA.
- Newly generated in-silico datasets constructed to simulate seasonal human influenza circulation.

The integration of both sources allowed alignment with internationally validated materials while tailoring the exercise to the operational reality of RELECOV clinical laboratories.

### 3.2 SARS-CoV-2 Dataset Selection

SARS-CoV-2 datasets were selected from the 2024 ECDC ESIB EQA to ensure comparability with internationally benchmarked material. Both Illumina and Nanopore panels included samples representing:

- High-quality baseline genomes.
- Low read-depth scenarios.
- Samples with numerous mixed sites.
- Contamination with non-target viral reads.
- Lineages of epidemiological relevance (e.g., recombinant or XBB-related lineages).

Only samples generated using the same ARTIC primer scheme (v4.1) were selected to avoid introducing variability associated with enrichment panel differences. This ensured that observed performance differences reflect analytical workflow characteristics rather than primer design heterogeneity.

{% set table_counter.value = table_counter.value + 1 %}

Table {{ table_counter.value }} summarises the correspondence between RELECOV EQA samples and their original source datasets, including ECDC ESIB references.

_**Table {{ table_counter.value }}**. Overview of SARS-CoV-2 datasets used in the RELECOV 2026 Dry-Lab EQA.
The table details sample origin, sequencing technology (Illumina paired-end or Oxford Nanopore Technologies), amplicon primer scheme version, and specific analytical characteristics intentionally selected to assess workflow robustness under challenging conditions._

| Sample | Source             | Platform | Amplicon primers version | Ref sample | Key Feature                                       | FASTQ files | Read layout |
|--------|--------------------|----------|--------------------------|------------|---------------------------------------------------|-------------|------------|
| SARS1  | ECDC-ESIB EQA 2024 | Illumina | ARTIC v4.1               | SARS2.04   | Influenza virus sample with some SARS-CoV-2 reads | 2           | Paired-end |
| SARS2  | ECDC-ESIB EQA 2024 | Illumina | ARTIC v4.1               | SARS2.01   | High-quality baseline sample                      | 2           | Paired-end |
| SARS3  | ECDC-ESIB EQA 2024 | Illumina | ARTIC v4.1               | SARS2.16   | XBB sample / insertion challenge                  | 2           | Paired-end |
| SARS4  | ECDC-ESIB EQA 2024 | Illumina | ARTIC v4.1               | SARS2.20   | Very low read depth                               | 2           | Paired-end |
| SARS5  | ECDC-ESIB EQA 2024 | Illumina | ARTIC v4.1               | SARS2.13   | >10 mixed sites                                   | 2           | Paired-end |
| SARS6  | ECDC-ESIB EQA 2024 | Nanopore | ARTIC v4.1               | SARS1.01   | High-quality baseline sample                      | 1           | Single-end |
| SARS7  | ECDC-ESIB EQA 2024 | Nanopore | ARTIC v4.1               | SARS1.09   | XBB sample / ambiguity next to a deletion         | 1           | Single-end |
| SARS8  | ECDC-ESIB EQA 2024 | Nanopore | ARTIC v4.1               | SARS1.15   | >10 mixed sites                                   | 1           | Single-end |
| SARS9  | ECDC-ESIB EQA 2024 | Nanopore | ARTIC v4.1               | SARS1.12   | Influenza virus sample with some SARS-CoV-2 reads | 1           | Single-end |
| SARS10 | ECDC-ESIB EQA 2024 | Nanopore | ARTIC v4.1               | SARS1.05   | Very low read depth                               | 1           | Single-end |

### 3.3 Influenza Dataset Selection

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
- Platform-specific read simulation using ART (Illumina) and Badread (Nanopore).

This approach allowed precise control over:

- Variant frequency structure
- Segment coverage distribution
- Contamination levels
- Platform-dependent error profiles

{% set table_counter.value = table_counter.value + 1 %}

Table {{ table_counter.value }} describes the design characteristics of in-silico samples, including virus composition and intended benchmarking challenges.

_**Table {{ table_counter.value }}**. Viral, host and contaminant composition design of in-silico influenza datasets used for benchmarking._

| Sample | Influenza reads | Host reads | Additional Viral reads | Total reads | Analytical Challenge                |
|--------|-----------------|------------------|------------------|-------------|-------------------------------------|
| FLU2   | 1378764         | 462520 | 0                          | 1841284     | Baseline performance assessment     |
| FLU4   | 181626          | 300000 | 200000 SARS-CoV-2 reads    | 681626      | False positive control              |
| FLU5   | 1088000         | 100000 | 0                          | 1188000     | NA segment dropout                  |
| FLU7   | 5677            | 100    | 255 Rhinovirus reads       | 6032        | Cross-virus contamination challenge |
| FLU8   | 5380            | 300    | 0                          | 5680        | Baseline performance assessment     |
| FLU9   | 19989           | 500    | 0                          | 20489       | HA segment dropout                  |

{% set table_counter.value = table_counter.value + 1 %}

Table {{ table_counter.value }} summarises the influenza datasets included in the EQA, detailing enrichment strategy, primer scheme, sequencing technology, and key analytical challenges.

_**Table {{ table_counter.value }}**. Influenza virus samples used in the RELECOV 2026 Dry-Lab EQA, including sequencing platform, enrichment strategy, primer scheme, and key analytical features._

| Sample | Source    | Platform | Enrichment Strategy | Primer Scheme                                   | Read Layout | Ref_sample       | Type   | Clade     | Key Feature                             |
|--------|-----------|----------|---------------------|-------------------------------------------------|-------------|------------------|--------|-----------|-----------------------------------------|
| FLU1   | ESIB 2024 | Illumina | Amplicon            | CommonUni12/13 (Van den Hoecke 2015)            | Paired-end | INFL2.07          | A/H5N1 | 2.3.4.4.b | High-quality baseline sample (zoonotic) |
| FLU2   | In-silico | Illumina | Amplicon            | Zhou 2009 single-reaction genomic amplification | Paired-end | In-silico Sample1 | A/H1N1 | D.3.1.1   | High-quality baseline sample (human)    |
| FLU3   | In-silico | Illumina | No enrichment       | —                                               | Paired-end | INFL2.04          | —      | —         | No influenza (Rhinovirus only)          |
| FLU4   | In-silico | Illumina | Amplicon            | Zhou 2009 single-reaction genomic amplification | Paired-end | In-silico Sample3 | A/H3N2 | K         | Contamination with SARS-CoV-2           |
| FLU5   | In-silico | Illumina | Amplicon            | Zhou 2009 single-reaction genomic amplification | Paired-end | In-silico Sample4 | A/H3N2 | J.2.2     | NA segment dropout                      |
| FLU6   | ESIB 2024 | Nanopore | No enrichment       | —                                               | Single-end | INFL1.02          | A/H5N6 | 2.3.4.4h  | High-quality baseline sample (zoonotic) |
| FLU7   | In-silico | Nanopore | Amplicon            | Zhou 2009 single-reaction genomic amplification | Single-end | In-silico Sample2 | A/H1N1 | C.1.9.3   | Contamination with Rhinovirus           |
| FLU8   | In-silico | Nanopore | Amplicon            | Zhou 2009 single-reaction genomic amplification | Single-end | In-silico Sample3 | A/H3N2 | K         | High-quality baseline sample (human)    |
| FLU9   | In-silico | Nanopore | Amplicon            | Zhou 2009 single-reaction genomic amplification | Single-end | In-silico Sample1 | A/H1N1 | D.3.1.1   | HA segment dropout                      |
| FLU10  | ESIB 2024 | Nanopore | Amplicon            | CommonUni12/13 (Van den Hoecke 2015)            | Single-end | INFL1.08          | A/H5N1 | 2.3.4.4b  | High-quality baseline sample (zoonotic) |

## 4. Methodology of Evaluation

The evaluation framework was designed to ensure objective, reproducible, and comparable assessment of analytical performance across participating laboratories. Submitted outputs were benchmarked against curated gold standard datasets from ECDC ESIB or generated in-silico.

The evaluation was structured into four independent analytical domains:

- Consensus genome reconstruction
- Variant detection
- Lineage/Type and Clade Assignment
- Metadata completeness and compliance

Each domain was assessed using predefined quantitative metrics to allow cross-laboratory comparison and pipeline benchmarking.


### 4.1 Evaluation of Consensus Genome Sequences</h3>

For each sample, a curated gold standard consensus genome was provided by the ECDC or generated internally. These reference sequences served as the ground truth for comparative analysis.

All submitted consensus sequences (.fasta) were:

- Aligned against the corresponding gold standard sequence using Mafft (TODO add version).
- Compared position-by-position relative to the declared reference genome coordinate system.

Differences between submitted sequences and gold standard sequences were categorised into the following classes: (TODO verificar que es verdad)

- Incorrect nucleotide substitution: A nucleotide different from the allowed reference or ambiguity code.
- Excess ambiguity: Ambiguity codes introduced where a defined nucleotide was expected.
- Missing ambiguity: Defined nucleotide provided where an ambiguity code was expected.
- Excess N stretch: Continuous region of Ns where defined bases were expected.
- Missing N stretch: Defined bases provided where Ns were expected.
- Insertion relative to gold standard
- Deletion relative to gold standard
- Each insertion, deletion, or contiguous stretch of Ns was counted as a single event.

For each laboratory and sample, the following metrics were calculated: (TODO verificar que es verdad)

- Total number of nucleotide discrepancies
- Percentage genome identity
- Number of indel events
- Proportion of ambiguous bases (Ns and IUPAC codes)
- Genome completeness (non-N positions / total genome length)

### 4.2 Evaluation of Variant Detection

A curated reference variant set was generated for each sample. Variant positions were standardized relative to a defined coordinate system referred to the references used by Nextclade.

Submitted .vcf files were: (TODO verificar si es verdad)

- Converted to a standardised long table format
- Compared position-by-position with the reference variant set

Each reported variant was categorised as: (TODO verificar si es verdad)

- True Positive (TP): Variant correctly identified.
- False Positive (FP): Variant reported but absent in gold standard.
- False Negative (FN): Variant present in gold standard but not reported.
- Total variant discrepancies were calculated as FP + FN.

For each laboratory and sample we calculated the following performance metrics: (TODO verificar si es verdad)
- Sensitivity = TP / (TP + FN)
- Precision = TP / (TP + FP)

Comparative analyses were performed to assess the influence of: (TODO verificar si es verdad)
- Allele frequency thresholds
- Minimum coverage thresholds
- Variant filtering criteria
- Reference genome selection


### 4.3 Evaluation of Lineage/Type and Clade Assignment

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

For SARS-CoV-2 and Influenza viruses, concordance was assessed as: (TODO verificar si es verdad)

- Exact match, when both lineage/subtype and clade were correct.
- Half match, when one of the classifications is correct but not the other.
- Discordant classification, when both classifications were incorrect.

Discrepancies were investigated to determine whether they resulted from: (TODO verificar si es verdad)

- Use of outdated lineage databases
- Differences in software versioning
- Incomplete consensus sequences
- Excess ambiguous positions

Failure to identify virus presence in positive samples, or misclassification of negative samples, was recorded separately.


### 4.4 Evaluation of Metadata Completeness and Compliance

Metadata assessment focused on analytical transparency and interoperability rather than biological correctness.

For each laboratory, metadata completeness was calculated as: (TODO verificar si es verdad)

$$
\text{Metadata completeness} =
\frac{\text{Number of correctly populated fields}}
{\text{Total number of applicable fields}}
$$

Fields were evaluated for:

 -Completion
- Compliance with controlled vocabularies
- Logical consistency (e.g., software declared vs tool-specific version fields completed)
- Valid file path reporting

Both mandatory and optional analytical fields were included in the completeness assessment, while fields not applicable to a laboratory’s selected components were excluded from scoring.
Metadata entries were considered non-compliant when:

- Controlled vocabulary options were bypassed
- Free-text substitutions replaced defined values
- Inconsistent analytical parameter reporting was observed

This evaluation allowed quantification of metadata standardisation and reproducibility readiness across the network.

### 4.5 Pipeline Benchmarking

The pipeline benchmarking analysis was designed to evaluate analytical performance at the pipeline and software level, rather than solely at the individual laboratory level. The objective was to identify which analytical workflows most consistently generate results that closely match the curated gold standard datasets.

For each declared pipeline or analytical workflow (including software combinations and parameter configurations), performance was aggregated across all laboratories using that approach.

The primary benchmarking criterion was the degree of similarity between submitted consensus genome sequences (.fasta files) and the corresponding gold standard reference sequences, which were derived from:

- The official ECDC reference datasets (for reused samples)
- The known reference genomes used to generate the in-silico simulated datasets.

Consensus-level identity was considered the central indicator of overall analytical robustness, as accurate genome reconstruction underpins downstream taxonomic classification.

For each pipeline/software group, the following metrics were calculated:

- Median number of nucleotide discrepancies relative to the gold standard
- Median percentage genome identity
- Genome completeness

Pipelines were ranked according to their ability to generate consensus sequences with:

- Minimal nucleotide divergence
- Maximum genome identity
- High genome completeness

Although consensus similarity was the primary benchmarking criterion, additional performance dimensions were integrated to contextualise overall workflow robustness:

- Lineage assignment concordance (SARS-CoV-2)
- Type/subtype and clade assignment concordance (influenza)
- Metadata completeness and compliance

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

The median number of components analysed per laboratory was {{ general.median_components_analysed_per_lab }}.

### 5.1 Submission Completeness

Across all components:

- {{ pct(network.submission_rates_pct.fasta) }} of laboratories submitted consensus genome files (.fasta), where applicable.
- {{ pct(network.submission_rates_pct.vcf) }} submitted variant call files (.vcf), where applicable.

The median number of components analysed per participating laboratory was {{ general.median_components_analysed_per_lab }}.

Submission rates were consistent across components, with minor variability reflecting differences in analytical scope and local implementation strategies.

### 5.2 Consensus Genome Reconstruction Performance

Overall performance was high for Illumina-based components (TODO verificar que es verdad), with a median genome identity of {{ pct(general.general_results.consensus.median_identity_illumina_pct, 2) }} for both Illumina components (TODO verificar que es veredad) and median genome identity of {{ pct(general.general_results.consensus.median_identity_nanopore_pct, 2) }} for Nanopore components. For Nanopore-based datasets, greater inter-laboratory variability was observed (TODO verificar si es verdad).

The main sources of variation included: (TODO verificar si es verdad)
- Differences in minimum coverage thresholds
- Handling of homopolymeric regions
- Indel filtering strategies
- Ambiguity and N masking policies

Across all laboratories and components, {{ pct(general.general_results.consensus.pct_genomes_below_discrepancy_threshold) }}
of submitted genomes showed fewer than {{ general.general_results.consensus.discrepancy_threshold }}
nucleotide discrepancies relative to the gold standard.

The most common discrepancies were: (TODO verificar si es verdad)
- Excess Ns in low-coverage regions
- Unfiltered indels in homopolymer stretches.

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} summarises consensus genome reconstruction performance across all components, stratified by sequencing platform.

{{ render_figure(general.figures.consensus_summary, "Network-level consensus reconstruction performance summary.") }}

**_Figure {{ fig_counter.value }}_. Distribution of consensus genome discrepancies relative to the gold standard across components**. Boxplots represent the number of nucleotide discrepancies per genome across participating laboratories. The central line indicates the median, boxes represent the interquartile range, and whiskers denote the full observed range.

### 5.3 Variant Detection Accuracy

Variant detection accuracy was assessed against curated reference variant sets.
Overall, {{ pct(general.general_results.variants.pct_vcfs_below_discrepancy_threshold) }} of submitted VCFs showed fewer than {{ general.general_results.variants.discrepancy_threshold }} discrepancies relative to the reference variant set.

{% set fig_counter.value = fig_counter.value + 1 %}

Illumina-based analyses generally demonstrated higher concordance and lower false-positive rates compared to Nanopore-based analyses (TODO verificar si es verdad). The distribution of variant detection performance across components is presented in Figure {{ fig_counter.value }}. Observed variability in variant detection performance was associated with:

- Allele frequency thresholds used for consensus incorporation
- Filtering of low-frequency variants
- Reference genome selection
- Variant normalization practices

{{ render_figure(general.figures.variant_summary, "Network-level variant detection performance summary.") }}

**_Figure {{ fig_counter.value }}_. Network-level variant detection performance across components**.
Boxplots display the distribution of nucleotide discrepancies between submitted VCF files and the curated reference variant set. The central line represents the median, boxes indicate the interquartile range, and whiskers denote the full observed range across participating laboratories.

### 5.4 Lineage, Type and Clade Assignment

Lineage, type and clade assignments were evaluated for concordance with gold standard classifications. Overall concordance rates were:

- SARS-CoV-2: {{ pct(general.general_results.classification.sars_cov_2_concordance_pct) }}
- Influenza type identification: {{ pct(general.general_results.classification.influenza_type_concordance_pct) }}
- Clade assignment (all viruses): {{ pct(general.general_results.classification.clade_concordance_pct_all_pct) }}

{% set fig_counter.value = fig_counter.value + 1 %}

As shown in Figure {{ fig_counter.value }}, classification concordance was high across components, with limited inter-laboratory variability. Most classification discrepancies were associated with:

- Use of outdated lineage database versions
- Differences in handling ambiguous consensus positions

Across components, median classification concordance exceeded {{ pct(general.general_results.classification.median_assignment_concordance) }}, and classification accuracy remained above {{ pct(general.general_results.classification.median_assignment_accuracy) }}.

{{ render_figure(general.figures.classification_summary, "Network-level classification performance summary.") }}

**_Figure {{ fig_counter.value }}_. Distribution of classification outcomes across participating laboratories**.
Stacked bars represent the proportion of exact matches, minor discrepancies, and incorrect assignments relative to curated gold standard classifications for each component.

### 5.5 Metadata Quality and Interoperability

The evaluation of metadata focused on analytical transparency, reproducibility, and interoperability within the RELECOV network. Completeness and compliance were assessed according to the criteria defined in Section 4.4, including controlled vocabulary adherence, logical consistency, and reporting of analytical parameters.

#### Overall Completeness

{% set fig_counter.value = fig_counter.value + 1 %}

Across all participating laboratories, the metadata template was completed at a mean completeness rate of {{ pct(general.metadata_completeness.mean_pct) }}, with values ranging from {{ pct(general.metadata_completeness.min_pct) }} to {{ pct(general.metadata_completeness.max_pct) }}. As illustrated in Figure {{ fig_counter.value }}, metadata completeness varied across participating laboratories, with a heterogeneous distribution across predefined completeness ranges.

Optional analytical fields contributed disproportionately to incompleteness (TODO comprobar si es verdad), particularly those related to parameter specification and software versioning.

{{ render_figure(general.figures.metadata_completeness_distribution,
  "Distribution of metadata completeness across participating laboratories.") }}

**_Figure {{ fig_counter.value }}_. Distribution of metadata completeness across participating laboratories**.
Boxplots represent the median and interquartile range of metadata completeness percentages. Whiskers denote the full observed range. The distribution reflects variability in reporting of analytical parameters, software versions, and controlled vocabulary adherence.

Additionally, component-level reporting rates for key analytical parameters (e.g., software name, software version, coverage thresholds, allele frequency thresholds, reference genome declaration) were quantified as the proportion of laboratories providing valid entries.

#### Reporting of Analytical Parameters

Although core pipeline tools were generally reported, variability was observed in the level of parameter detail provided.

- {{ pct(general.metadata_completeness.software_names_pct) }} of laboratories reported exact software names.
- {{ pct(general.metadata_completeness.software_version_pct) }} reported exact software versions.
- {{ pct(general.metadata_completeness.coverage_threshold_pct) }} specified minimum coverage thresholds.
- {{ pct(general.metadata_completeness.frequency_threshold_pct) }} declared allele frequency thresholds used for consensus incorporation.
- {{ pct(general.metadata_completeness.reference_genome_pct) }} reported the reference genome accession or identifier.

Incomplete parameter reporting limited the ability to fully reconstruct or reproduce analytical workflows in {{ pct(general.metadata_completeness.incomplete_parameters_pct) }} of submissions.

#### Controlled Vocabulary Compliance

Compliance with predefined controlled vocabularies was evaluated to assess standardisation readiness:

- {{ pct(general.metadata_completeness.fully_compliant_pct) }} of submissions were fully compliant with controlled vocabulary requirements.
- {{ pct(general.metadata_completeness.free_text_predefine_pct) }} contained at least one free-text substitution where a predefined option was required.
- {{ pct(general.metadata_completeness.inconsistent_tool_version) }} presented inconsistencies between declared tools and version fields.

The most common compliance issues included (TODO revisar si es cierto):

- Free-text entry of software names instead of predefined values.
- Inconsistent declaration of lineage assignment tools.
- Ambiguous reporting of reference genome identifiers.

#### Logical Consistency and File Traceability

Logical consistency checks were performed to ensure coherence between declared metadata and submitted outputs.

Instances of inconsistency included:

- Declared tools not matching file output format.
- Missing file paths for submitted consensus or variant files.
- Discrepancies between declared reference genome and variant coordinate system.

Metadata was complete in {{ general.metadata_completeness.fully_compliant_pct }} of submissions, while {{ general.metadata_completeness.clarification_pct }} required clarification or correction during validation.

#### Diversity of Analytical Workflows

The metadata submissions allowed characterisation of the analytical landscape currently implemented across the RELECOV network.

A total of {{ general.metadata_completeness.total_workflows }} distinct analytical workflows were identified across participating laboratories, defined as unique combinations of software tools, versions, and parameter configurations declared in the metadata template.

Substantial diversity was observed in the selection of core analytical tools:

- Consensus reconstruction software ( {{ general.metadata_completeness.total_consensus_softwares }} distinct tool configurations )
- Variant calling tools ( {{ general.metadata_completeness.total_variant_softwares }} distinct tools or configurations )
- Lineage/type and clade assignment software ( {{ general.metadata_completeness.total_lineage_softwares }} tools or database versions )

Variation was not limited to tool choice, but extended to: (TODO revisar si es verdad)

- Reference genome selection
- Coverage thresholds
- Allele frequency cut-offs
- Indel filtering strategies

This diversity provides a descriptive overview of the analytical landscape within the network and forms the basis for the component-specific benchmarking analyses presented in Section 5.

### 5.6 Pipeline Benchmarking and Comparative Performance

The benchmarking framework was designed to assess whether differences in analytical software and parameterisation were associated with measurable variability in performance across participating laboratories.

Substantial heterogeneity was observed in:

- Choice of consensus reconstruction software
- Variant calling strategies
- Lineage and clade assignment tools
- Reference genome selection
- Coverage and allele frequency thresholds


This diversity reflects the decentralised analytical capacity of the RELECOV network (TODO revisar si es verdad).

Comparative performance analyses stratified by component are presented in Section 6, where software-level differences are evaluated within homogeneous analytical contexts (SARS-CoV-2 Illumina, SARS-CoV-2 Nanopore, Influenza Illumina, Influenza Nanopore).

At a network level, no single analytical software solution was universally optimal across all components. (TODO revisar si es verdad) Instead, performance was influenced by the interaction between:

- Software selection
- Parameter configuration
- Sequencing platform characteristics
- Sample complexity

These findings highlight the importance of harmonising minimum analytical criteria while preserving methodological flexibility within the network.

## 6. Component-specific Results

This section presents the analytical results disaggregated by component, allowing a detailed assessment of performance within each dataset and sequencing technology. For each component, results are structured according to participation and submission metrics, consensus genome reconstruction performance, variant detection accuracy, and lineage, type or clade assignment concordance, as applicable.

Component-level analyses enable identification of platform-specific patterns, dataset-dependent challenges, and variability associated with particular sample characteristics. This approach facilitates a more granular interpretation of performance differences observed at the network level and supports targeted harmonisation recommendations.

{% for comp_code in global.components %}
  {% set comp_net = global.components_detail.get(comp_code) %}
  {% if comp_net %}
### 6.{{ loop.index }} {{ comp_code }} ({{ comp_net.name }})

#### 6.{{ loop.index }}.1 Participation and Submissions

A total of {{ comp_net.total_labs }} laboratories submitted results for the {{ comp_code }} component.

- {{ comp_net.total_fasta }} submitted consensus genome sequences (.fasta), where applicable.
- {{ comp_net.total_vcf }} submitted variant call files (.vcf), where applicable.
- The metadata template completeness for {{ comp_code }} submissions had a mean of {{ pct(comp_net.metadata_completeness_mean) }}.

#### 6.{{ loop.index }}.2 Consensus Genome Reconstruction

Consensus sequences were evaluated against the corresponding curated gold standard references for each sample in the {{ comp_code }} component.

Overall, {{ comp_code }} showed a median genome identity of {{ pct(comp_net.consensus.median_identity, 2) }}, with a median of {{ comp_net.consensus.median_discrepancies }} nucleotide discrepancies per sample (range: {{ comp_net.consensus.min_discrepancies }}–{{ comp_net.consensus.max_discrepancies }}).

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level consensus reconstruction metrics per sample for {{ comp_code }}.**

| Sample ID | Median genome identity (%) | Median discrepancies | IQR discrepancies |
|---|---:|---:|---:|
{% for s in comp_net.consensus.samples %}
| {{ s.sample_id }} | {{ "%.2f"|format(s.median_identity_pct) }} | {{ s.median_discrepancies }} | {{ s.iqr_discrepancies[0] }}–{{ s.iqr_discrepancies[1] }} |
{% endfor %}

Figure {{ fig_counter.value + 1 }} presents the distribution of nucleotide discrepancies per sample across participating laboratories for {{ comp_code }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.consensus.fig_discrepancies_boxplot_by_sample,
  "Consensus discrepancies per sample for " ~ comp_code ~ " relative to the curated gold standard."
) }}

**Figure {{ fig_counter.value }}. Distribution of consensus discrepancies per sample for {{ comp_code }}.** Boxplots represent the number of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample. The central line indicates the median, boxes denote the interquartile range, and whiskers represent the full observed range.

Discrepancy type composition (aggregated across all submitted consensus sequences for {{ comp_code }}):

- Incorrect substitutions: {{ pct(comp_net.consensus.incorrect_nt_pct) }} of discrepancies  
- Excess Ns / ambiguous bases: {{ pct(comp_net.consensus.excess_ambiguous_pct) }}  
- Indels: {{ pct(comp_net.consensus.indels_pct) }}

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level discrepancy composition by type for {{ comp_code }}.**

| Discrepancy type | Network median per sample | Network IQR per sample |
|---|---:|---:|
| Incorrect substitutions | {{ comp_net.consensus.by_type.substitutions.median }} | {{ comp_net.consensus.by_type.substitutions.iqr[0] }}–{{ comp_net.consensus.by_type.substitutions.iqr[1] }} |
| Excess Ns / ambiguous bases | {{ comp_net.consensus.by_type.excess_Ns.median }} | {{ comp_net.consensus.by_type.excess_Ns.iqr[0] }}–{{ comp_net.consensus.by_type.excess_Ns.iqr[1] }} |
| Missing Ns (expected Ns not reported) | {{ comp_net.consensus.by_type.missing_Ns.median }} | {{ comp_net.consensus.by_type.missing_Ns.iqr[0] }}–{{ comp_net.consensus.by_type.missing_Ns.iqr[1] }} |
| Insertions | {{ comp_net.consensus.by_type.insertions.median }} | {{ comp_net.consensus.by_type.insertions.iqr[0] }}–{{ comp_net.consensus.by_type.insertions.iqr[1] }} |
| Deletions | {{ comp_net.consensus.by_type.deletions.median }} | {{ comp_net.consensus.by_type.deletions.iqr[0] }}–{{ comp_net.consensus.by_type.deletions.iqr[1] }} |

The most frequent discrepancy pattern observed in {{ comp_code }} was {{ comp_net.consensus.most_frequent_discrepancy_pattern }}.

Figure {{ fig_counter.value + 1 }} summarises the contribution of each discrepancy category observed in {{ comp_code }} relative to the curated gold standard.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.consensus.fig_discrepancy_type_barplot,
  "Composition of consensus discrepancy types for " ~ comp_code ~ " relative to the curated gold standard."
) }}

**Figure {{ fig_counter.value }}. Composition of consensus discrepancy types relative to the curated gold standard for {{ comp_code }}.** Bars represent aggregated discrepancies across all submitted consensus sequences, stratified by discrepancy category (incorrect substitutions, excess ambiguous bases, and indels).

#### 6.{{ loop.index }}.3 Variant Detection Performance

Variant call files (.vcf) submitted for the {{ comp_code }} component were compared against the curated reference variant set corresponding to each sample.

Overall, {{ comp_code }} showed a median sensitivity of {{ pct(comp_net.variant.median_sensitivity, 2) }} and a median precision of {{ pct(comp_net.variant.median_precision, 2) }} across participating laboratories.

The median number of false positives per sample was {{ comp_net.variant.median_false_positives }} (range: {{ comp_net.variant.min_false_positives }}–{{ comp_net.variant.max_false_positives }}), while the median number of missed expected variants (false negatives) was {{ comp_net.variant.median_false_negatives }}.

- {{ pct(comp_net.variant.no_false_positives_pct) }} of submissions reported no false positive variants.
- {{ pct(comp_net.variant.at_least_one_fn_pct) }} showed at least one missed expected variant.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level variant detection performance per sample for {{ comp_code }}.**

| Sample ID | Median sensitivity (%) | Median precision (%) | Median TP | Median FP | Median FN | Median total differences | IQR total differences |
|---|---:|---:|---:|---:|---:|---:|---:|
{% for s in comp_net.variant.samples %}
| {{ s.sample_id }} | {{ "%.2f"|format(s.median_sensitivity_pct) }} | {{ "%.2f"|format(s.median_precision_pct) }} | {{ s.median_tp }} | {{ s.median_fp }} | {{ s.median_fn }} | {{ s.median_total_differences }} | {{ s.iqr_total_differences[0] }}–{{ s.iqr_total_differences[1] }} |
{% endfor %}

Figure {{ fig_counter.value + 1 }} presents the distribution of variant detection discrepancies per sample across participating laboratories for {{ comp_code }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.variant.fig_total_differences_boxplot_by_sample,
  "Variant detection discrepancies per sample for " ~ comp_code ~ " relative to the curated reference variant set."
) }}

**Figure {{ fig_counter.value }}. Distribution of variant detection discrepancies per sample for {{ comp_code }}.**  
Boxplots represent the distribution of total differences (FP + FN) relative to the curated reference variant set across participating laboratories for each sample. The central line indicates the median, boxes denote the interquartile range, and whiskers represent the full observed range.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level summary of variant detection performance metrics for {{ comp_code }}.**

| Metric | Network median | Network IQR |
|---|---:|---:|
| Sensitivity (%) | {{ "%.2f"|format(comp_net.variant.metrics.sensitivity.median_pct) }} | {{ comp_net.variant.metrics.sensitivity.iqr_pct[0] }}–{{ comp_net.variant.metrics.sensitivity.iqr_pct[1] }} |
| Precision (%) | {{ "%.2f"|format(comp_net.variant.metrics.precision.median_pct) }} | {{ comp_net.variant.metrics.precision.iqr_pct[0] }}–{{ comp_net.variant.metrics.precision.iqr_pct[1] }} |
| True positives (TP) | {{ comp_net.variant.metrics.tp.median }} | {{ comp_net.variant.metrics.tp.iqr[0] }}–{{ comp_net.variant.metrics.tp.iqr[1] }} |
| False positives (FP) | {{ comp_net.variant.metrics.fp.median }} | {{ comp_net.variant.metrics.fp.iqr[0] }}–{{ comp_net.variant.metrics.fp.iqr[1] }} |
| False negatives (FN) | {{ comp_net.variant.metrics.fn.median }} | {{ comp_net.variant.metrics.fn.iqr[0] }}–{{ comp_net.variant.metrics.fn.iqr[1] }} |
| Total differences (FP + FN) | {{ comp_net.variant.metrics.total_differences.median }} | {{ comp_net.variant.metrics.total_differences.iqr[0] }}–{{ comp_net.variant.metrics.total_differences.iqr[1] }} |

Figure {{ fig_counter.value + 1 }} summarises the distribution of key variant detection metrics across participating laboratories for {{ comp_code }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.variant.fig_metric_boxplots_by_pipeline_or_overall,
  "Distribution of key variant detection metrics across laboratories for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of key variant detection metrics across laboratories for {{ comp_code }}.**  
The figure summarises laboratory-level distributions of sensitivity, precision, and discrepancy counts relative to the curated reference variant set.

#### 6.{{ loop.index }}.4 Lineage/Type and Clade Assignment

Lineage/Type and clade assignments submitted for the {{ comp_code }} component were evaluated for concordance with the curated gold standard classifications.

Across all participating laboratories:

- Exact concordance (lineage/type and clade correct): {{ pct(comp_net.typing.exact_concordance_pct) }}
- Partial concordance (only lineage/type or clade correct): {{ pct(comp_net.typing.partial_concordance_pct) }}
- Discordant classification: {{ pct(comp_net.typing.discordance_pct) }}

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level classification outcomes per sample for {{ comp_code }}.**

| Sample ID | Exact concordance (%) | Partial concordance (%) | Discordant (%) |
|---|---:|---:|---:|
{% for s in comp_net.typing.samples %}
| {{ s.sample_id }} | {{ "%.2f"|format(s.exact_pct) }} | {{ "%.2f"|format(s.partial_pct) }} | {{ "%.2f"|format(s.discordant_pct) }} |
{% endfor %}

Figure {{ fig_counter.value + 1 }} presents the distribution of classification outcomes per sample across participating laboratories.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.typing.fig_stacked_bar_by_sample,
  "Classification outcome distribution per sample for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Classification outcome distribution per sample for {{ comp_code }}.**  
Stacked bars represent the proportion of Exact, Partial, and Discordant lineage/type and clade assignments across participating laboratories for each sample relative to the curated gold standard classification.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level classification error counts for {{ comp_code }}.**

| Error type | Network median per sample | Network IQR per sample |
|---|---:|---:|
| Number of lineage/type errors | {{ comp_net.typing.errors.lineage.median }} | {{ comp_net.typing.errors.lineage.iqr[0] }}–{{ comp_net.typing.errors.lineage.iqr[1] }} |
| Number of clade errors | {{ comp_net.typing.errors.clade.median }} | {{ comp_net.typing.errors.clade.iqr[0] }}–{{ comp_net.typing.errors.clade.iqr[1] }} |

Figure {{ fig_counter.value + 1 }} compares the distribution of lineage/type and clade assignment errors across participating laboratories.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.typing.fig_error_distribution,
  "Distribution of lineage/type and clade assignment errors for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of lineage/type and clade assignment errors for {{ comp_code }}.**  
Boxplots represent the number of lineage/type and clade errors per sample across participating laboratories. The central line indicates the median, boxes denote the interquartile range, and whiskers represent the full observed range.

#### 6.{{ loop.index }}.5 Declared Software and Pipeline Summary

Based on metadata submissions, {{ comp_net.pipeline.total_number }} distinct pipeline/software configurations were reported for the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Performance summary of declared software/pipeline configurations for {{ comp_code }}.**

| Pipeline / Software | N labs | Median genome identity (%) | Median discrepancies | Median metadata completeness (%) | Exact classification concordance (%) |
|---|---:|---:|---:|---:|---:|
{% for p in comp_net.pipeline.configurations %}
| {{ p.name }} | {{ p.n_labs }} | {{ "%.2f"|format(p.median_identity_pct) }} | {{ p.median_discrepancies }} | {{ "%.1f"|format(p.median_metadata_completeness_pct) }} | {{ "%.1f"|format(p.exact_classification_pct) }} |
{% endfor %}

The comparative positioning of declared workflows within the {{ comp_code }} component is shown in Figure {{ fig_counter.value + 1 }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.pipeline.fig_scatter_positioning,
  "Comparative positioning of declared workflows for " ~ comp_code ~ " based on consensus accuracy and classification concordance."
) }}

**Figure {{ fig_counter.value }}. Positioning of declared workflows for {{ comp_code }}.**  
Each point represents a distinct pipeline/software configuration. The x-axis indicates median genome identity relative to the curated gold standard, and the y-axis indicates exact lineage/type and clade concordance. Point size reflects the number of laboratories reporting each configuration.

Figure {{ fig_counter.value + 1 }} summarises the distribution of key performance metrics stratified by declared pipeline configuration.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.pipeline.fig_metric_boxplots_by_pipeline,
  "Distribution of performance metrics by pipeline configuration for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared pipeline configuration for {{ comp_code }}.**  
Multi-panel boxplots summarise laboratory-level performance stratified by pipeline/software configuration. Panels display genome identity (%), discrepancy counts, metadata completeness (%), and exact classification concordance (%). The central line indicates the median, boxes represent the interquartile range, and whiskers denote the full observed range across laboratories using each configuration.

  {% endif %}
{% endfor %}

## 7. Discussion

(TODO verificar si tiene sentido)

The 2026 RELECOV Dry-Lab EQA represents the first drylab benchmarking exercise focused only on bioinformatic analytical performance across the RELECOV network. By integrating internationally validated ECDC datasets with purpose-designed in-silico scenarios reflecting routine clinical surveillance, this exercise provides a comprehensive overview of the current analytical landscape within RELECOV.

### 7.1. Overall Analytical Robustness

Across components, consensus genome reconstruction showed high overall concordance with curated gold standards, particularly for Illumina-based datasets. This indicates that most laboratories possess mature workflows for routine SARS-CoV-2 genome reconstruction.
However, increased variability observed in Nanopore-based datasets highlights persistent challenges related to:

- Homopolymer-associated indels
- Coverage threshold policies
- Ambiguity handling
- Segment-specific dropout in influenza

These findings are consistent with the known technical characteristics of long-read sequencing technologies and underscore the importance of platform-specific optimisation rather than universal parameter application.

### 7.2. Variant Detection Variability

Variant detection demonstrated generally high sensitivity and Precision in high-quality samples. Nonetheless, variability increased in analytically challenging scenarios, including low read depth and mixed-site samples.

The observed differences were primarily associated with:

- Divergent allele frequency thresholds
- Inconsistent filtering criteria
- Reference genome selection

Although variant detection was not incorporated into cross-pipeline ranking, the variability observed indicates that harmonised minimal variant reporting criteria would enhance inter-laboratory comparability.

### 7.3. Classification Accuracy and Database Versioning

Lineage, type, and clade assignment performance was high overall. Most discrepancies were attributable to outdated lineage databases or incomplete consensus reconstruction rather than fundamental classification errors.

This finding suggests that classification performance is less dependent on core algorithmic differences and more sensitive to:

- Version control practices
- Regular database updates
- Consensus sequence completeness

Thus, governance and update policies may have a greater impact on classification accuracy than software selection alone.

### 7.4. Workflow Diversity and Standardisation Balance
A high diversity of declared analytical workflows was observed across the network. This heterogeneity reflects distributed expertise and technical autonomy within participating laboratories.

However, diversity in:

- Reference genome selection
- Coverage thresholds
- Allele frequency cut-offs
- Indel handling strategies

introduces measurable inter-laboratory variability.

The results indicate that harmonisation efforts should focus on defining minimum performance and parameter standards rather than enforcing a single analytical pipeline.

### 7.5 Implications for RELECOV 2.0

The benchmarking exercise directly informs the development of the RELECOV analytical platform. The findings support:

- Definition of minimum consensus accuracy thresholds
- Establishment of recommended allele frequency cut-offs
- Mandatory version control reporting
- Standardised metadata schema enforcement

Importantly, no single pipeline demonstrated universal superiority across all components. Instead, analytical robustness emerged from the interaction between software choice, parameter configuration, and sequencing platform characteristics.

Therefore, harmonisation within RELECOV 2.0 should prioritise:

- Performance-based criteria
- Metadata interoperability
- Transparent version control
- Platform-aware optimisation

## 8. Conclusions

(TODO verificar si tiene sentido)

The 2026 RELECOV Dry-Lab EQA demonstrates that the network possesses strong bioinformatic capacity for respiratory virus genomic surveillance, with high overall concordance in consensus genome reconstruction and classification tasks.

Illumina-based workflows showed highly consistent performance across laboratories. Nanopore-based analyses exhibited greater variability, particularly in challenging genomic regions, indicating the need for platform-specific harmonisation guidance.

Variant detection performance was generally robust but sensitive to threshold and filtering heterogeneity, supporting the definition of minimal reporting standards.

Metadata evaluation revealed that while core analytical information is routinely reported, variability in parameter documentation and controlled vocabulary compliance remains a limiting factor for full interoperability.

No single analytical workflow was universally optimal. Performance was influenced by the interaction between software selection, parameter configuration, and sequencing technology.

Collectively, these findings provide a technical foundation for:

- Defining minimum analytical performance criteria
- Establishing harmonised metadata standards
- Guiding workflow standardisation within RELECOV 2.0
- Supporting long-term sustainability of national genomic surveillance infrastructure

The EQA therefore provides a robust technical basis for harmonised, performance-driven genomic surveillance within RELECOV 2.0.

{% if labdata %}
# 9. Individual Laboratory Technical Report
## Laboratory: {{ labdata.lab.laboratory_name }} ({{ labdata.lab.lab_cod }})

This section provides a detailed technical assessment of the analytical results submitted by **{{ labdata.lab.lab_cod }}** within the 2026 RELECOV Dry-Lab EQA. Performance metrics are benchmarked against curated gold standards and contextualised relative to aggregated network-wide performance distributions. Network medians and interquartile ranges are provided for comparative interpretation, without disclosure of other laboratories’ identities.

The purpose of this section is to support technical optimisation, parameter harmonisation, and alignment with the analytical standards defined within RELECOV 2.0.

## 9.1 Participation Overview

The laboratory analysed **{{ labdata.lab_overview.components_analysed_count }}** out of 4 components. Network median components analysed per laboratory: **{{ general.median_components_analysed_per_lab }}**.

Analysed components:

{% for c in general.components %}
- {{ c }}: {{ "✔" if c in labdata.lab_overview.components_analysed else "✖" }}
{% endfor %}

Submitted outputs (analysed components):

- `.fasta`: **{{ labdata.lab_overview.submitted_files.fasta_submitted }} / {{ labdata.lab_overview.submitted_files.fasta_expected }}**
- `.vcf`: **{{ labdata.lab_overview.submitted_files.vcf_submitted }} / {{ labdata.lab_overview.submitted_files.vcf_expected }}**

Regarding metadata completeness:

- Metadata completeness for **{{ labdata.lab.lab_cod }}**: **{{ pct(labdata.lab_overview.metadata.completeness_pct) }}**
- Network median metadata completeness: **{{ pct(network.metadata_completeness.median_pct) }}**  
- Network range: **{{ pct(network.metadata_completeness.min_pct) }}–{{ pct(network.metadata_completeness.max_pct) }}**

{% if labdata.lab_overview.metadata.primary_incompleteness_drivers %}
Primary contributors to incompleteness:
{% for d in labdata.lab_overview.metadata.primary_incompleteness_drivers %}
- {{ d }}
{% endfor %}
{% endif %}

{% if labdata.lab_overview.metadata.missing_sample_ids and labdata.lab_overview.metadata.missing_sample_ids|length > 0 %}
Missing submissions were observed in: {{ labdata.lab_overview.metadata.missing_sample_ids|join(", ") }}.
{% else %}
No missing sample-level submissions were observed for the analysed components.
{% endif %}

{% for comp_code, comp in labdata.components.items() %}

# 9.{{ loop.index + 1 }} {{ comp_code }} ({{ comp.display_name }})

The laboratory submitted results for the **{{ comp_code }}** component{% if comp.source %} ({{ comp.source }}){% endif %}.

## 9.{{ loop.index + 1 }}.1 Consensus Genome Reconstruction

Consensus genome sequences (`.fasta`) submitted by **{{ labdata.lab.lab_cod }}** were compared against the curated gold standard reference for each sample included in the {{ comp_code }} component.

### Per-sample summary metrics

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Per-sample consensus reconstruction metrics for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

| Sample ID | Genome identity (%) | Total discrepancies | Indel events | Ambiguous bases (%) | Genome completeness (%) |
|---|---:|---:|---:|---:|---:|
{% for sample_id, s in comp.samples.items() -%}
| {{ sample_id }} | {{ num(s.consensus.genome_identity_pct, 4) }} | {{ s.consensus.total_discrepancies }} | {{ s.consensus.indel_events }} | {{ num(s.consensus.ambiguous_bases_pct, 2) }} | {{ num(s.consensus.genome_completeness_pct, 2) }} |
{% endfor %}

The metrics presented in Table {{ table_counter.value }} summarise overall sequence similarity, discrepancy burden, and completeness relative to the curated gold standard reference.

### Discrepancy type breakdown per sample

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Discrepancy type breakdown per sample for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

| Sample ID | Substitutions | Excess Ns | Missing Ns | Insertions | Deletions |
|---|---:|---:|---:|---:|---:|
{% for sample_id, s in comp.samples.items() -%}
| {{ sample_id }} | {{ s.consensus.discrepancy_breakdown.substitutions }} | {{ s.consensus.discrepancy_breakdown.excess_Ns }} | {{ s.consensus.discrepancy_breakdown.missing_Ns }} | {{ s.consensus.discrepancy_breakdown.insertions }} | {{ s.consensus.discrepancy_breakdown.deletions }} |
{% endfor %}

Table {{ table_counter.value }} provides a detailed characterisation of discrepancy categories contributing to the total differences observed for each sample.

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} illustrates the distribution of mean nucleotide discrepancies per sample across participating laboratories in the network, contextualising the results obtained by **{{ labdata.lab.lab_cod }}**.

{{ render_fig(
  comp.component_figures.consensus_discrepancy_distribution,
  comp_code ~ ": distribution of mean discrepancies per sample across the network; red marker indicates " ~ labdata.lab.lab_cod ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of mean consensus discrepancies per sample across participating laboratories ({{ comp_code }}).**  
Boxplots represent the distribution of nucleotide discrepancies relative to the curated gold standard across the RELECOV network. The central line indicates the median, boxes denote the interquartile range, and whiskers represent the full observed range. The red marker corresponds to the results obtained by **{{ labdata.lab.lab_cod }}**.

## 9.{{ loop.index + 1 }}.2 Variant Detection Performance

Variant call files (`.vcf`) submitted by **{{ labdata.lab.lab_cod }}** were compared against the curated reference variant set for each sample included in the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Per-sample variant detection performance metrics for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

| Sample ID | Sensitivity (%) | Precision (%) | True Positives (TP) | False Positives (FP) | False Negatives (FN) |
|---|---:|---:|---:|---:|---:|
{% for sample_id, s in comp.samples.items() -%}
| {{ sample_id }} | {{ num(s.variants.sensitivity_pct, 2) }} | {{ num(s.variants.precision_pct, 2) }} | {{ s.variants.tp }} | {{ s.variants.fp }} | {{ s.variants.fn }} |
{% endfor %}

The metrics presented in Table {{ table_counter.value }} summarise per-sample variant detection accuracy relative to the curated reference variant set, including sensitivity, precision, and discrepancy counts.

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} illustrates the distribution of variant detection performance metrics across participating laboratories in the network, contextualising the results obtained by **{{ labdata.lab.lab_cod }}**.

{{ render_fig(
  comp.component_figures.variant_metrics_distribution,
  comp_code ~ ": distribution of variant detection metrics across the network; red marker indicates " ~ labdata.lab.lab_cod ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of variant detection performance across participating laboratories ({{ comp_code }}).**  
Boxplots represent the distribution of laboratory-level variant detection metrics relative to the curated reference variant set. The central line indicates the median, boxes denote the interquartile range, and whiskers represent the full observed range across the network. The red marker corresponds to the results obtained by **{{ labdata.lab.lab_cod }}**.

## 9.{{ loop.index + 1 }}.3 Lineage/Type and Clade Assignment

Lineage/type and clade assignments submitted by **{{ labdata.lab.lab_cod }}** were compared against the curated gold standard classifications for each sample included in the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Per-sample lineage/type and clade assignment results for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

| Sample ID | Expected lineage/type | Reported lineage/type | Expected clade | Reported clade | Concordance |
|---|---|---|---|---|---|
{% for sample_id, s in comp.samples.items() -%}
| {{ sample_id }} | {{ s.classification.expected_lineage }} | {{ s.classification.reported_lineage }} | {{ s.classification.expected_clade }} | {{ s.classification.reported_clade }} | {{ s.classification.concordance }} |
{% endfor %}

Table {{ table_counter.value }} summarises the concordance between expected and reported lineage/type and clade classifications for each sample.

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} illustrates the distribution of classification outcomes per sample across participating laboratories, contextualising the results obtained by **{{ labdata.lab.lab_cod }}**.

{{ render_fig(
  comp.component_figures.classification_stackedbars,
  comp_code ~ ": distribution of classification outcomes per sample across the network; red marker indicates " ~ labdata.lab.lab_cod ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of lineage/type and clade classification outcomes across participating laboratories ({{ comp_code }}).**  
Stacked bars represent the proportion of Exact, Partial, and Discordant classifications within the network for each sample relative to the curated gold standard. The red marker indicates the classification outcome reported by **{{ labdata.lab.lab_cod }}**.

{% if comp.get("classification_error_counts") %}

### Classification error counts

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Classification error counts for {{ labdata.lab.lab_cod }} compared to network medians ({{ comp_code }}).**

| Discrepancy type | {{ labdata.lab.lab_cod }} | Network median |
|---|---:|---:|
| Number of lineage/type errors | {{ comp.classification_error_counts.lab.lineage_errors }} | {{ comp.classification_error_counts.network_median.lineage_errors }} |
| Number of clade errors | {{ comp.classification_error_counts.lab.clade_errors }} | {{ comp.classification_error_counts.network_median.clade_errors }} |

Table {{ table_counter.value }} compares the number of lineage/type and clade assignment errors observed for **{{ labdata.lab.lab_cod }}** against the median values calculated across participating laboratories.

{% endif %}

## 9.{{ loop.index + 1 }}.4 Pipeline Performance Positioning within the Network

The analytical workflow declared by **{{ labdata.lab.lab_cod }}** was benchmarked against other workflows implemented across the RELECOV network for the {{ comp_code }} component.

Positioning was evaluated based on two primary performance indicators:

1. Median consensus genome identity relative to the curated gold standard.
2. Exact lineage/type and clade classification concordance.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Workflow performance positioning for {{ labdata.lab.lab_cod }} within the network ({{ comp_code }}).**

| Metric | {{ labdata.lab.lab_cod }} workflow | Network median | Network IQR |
|---|---:|---:|---:|
| Median genome identity (%) | {{ num(comp.workflow_positioning.median_genome_identity_pct, 4) }} | {{ num(comp.workflow_positioning.network_median_identity_pct, 4) }} | {{ iqr_range(comp.workflow_positioning.network_iqr_identity_pct, 4) }} |
| Exact classification concordance (%) | {{ num(comp.workflow_positioning.exact_classification_concordance_pct, 2) }} | {{ num(comp.workflow_positioning.network_median_class_pct, 2) }} | {{ iqr_range(comp.workflow_positioning.network_iqr_class_pct, 2) }} |

Table {{ table_counter.value }} contextualises the performance of the declared workflow relative to aggregated network-level metrics. Network medians and interquartile ranges (IQR) provide a reference distribution against which the positioning of the declared workflow can be interpreted.

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} illustrates the comparative positioning of declared workflows within the network for the {{ comp_code }} component.

{{ render_fig(
  comp.workflow_positioning.scatter_fig,
  comp_code ~ ": workflow positioning across the network (grey), with " ~ labdata.lab.lab_cod ~ " highlighted (red)."
) }}

**Figure {{ fig_counter.value }}. Workflow positioning within the RELECOV network for {{ comp_code }}.**  
Each point represents a declared analytical workflow reported by participating laboratories. The x-axis indicates median consensus genome identity relative to the curated gold standard, and the y-axis indicates exact lineage/type and clade classification concordance. Grey markers represent network workflows, while the red marker corresponds to the workflow declared by **{{ labdata.lab.lab_cod }}**.

## 9.{{ loop.index + 1 }}.5 Metadata-Derived Analytical Metrics (per sample)

This section summarises selected quantitative analytical metrics declared in the metadata submission of **{{ labdata.lab.lab_cod }}**, disaggregated by sample within the {{ comp_code }} component.

Only metrics explicitly provided by the laboratory are included in the comparative assessment. Network-level medians and interquartile ranges (IQR) are shown for contextual interpretation.

{% for sample_id, s in comp.samples.items() %}
{% set m = s.metadata_metrics %}

{% if m %}
{% set table_counter.value = table_counter.value + 1 %}

### {{ sample_id }}

**Table {{ table_counter.value }}. Metadata-derived analytical metrics for {{ labdata.lab.lab_cod }} ({{ comp_code }}, {{ sample_id }}).**

| Metric | {{ labdata.lab.lab_cod }} | Network median | Network IQR |
|---|---:|---:|---:|
{% for metric_key, metric_label in {
  "pct_genome_covered_ge_10x": "% Genome Covered ≥10X",
  "mean_depth": "Mean Depth of Coverage",
  "pct_genome_masked_Ns": "% Genome Masked (Ns)",
  "pct_virus_reads": "% Virus Reads",
  "pct_host_reads": "% Host Reads",
  "total_variants_reported": "Total Variants Reported",
  "variants_with_predicted_effect": "Variants with Predicted Effect"
}.items() %}
{% if m.get(metric_key) is not none %}
| {{ metric_label }} | {{ m[metric_key] }} | {{ s.get("metadata_metrics_network", {}).get(metric_key, {}).get("median", "NA") }} | {{ iqr_range(s.get("metadata_metrics_network", {}).get(metric_key, {}).get("iqr"), 2) if s.get("metadata_metrics_network", {}).get(metric_key) else "NA" }} |
{% endif %}
{% endfor %}

Table {{ table_counter.value }} contextualises laboratory-reported analytical parameters relative to the distribution of values observed across the RELECOV network for the same sample.

{% endif %}
{% endfor %}

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} illustrates the distribution of metadata-derived analytical metrics across participating laboratories for the {{ comp_code }} component.

{{ render_fig(
  comp.component_figures.metadata_metrics_panel,
  comp_code ~ ": distribution of metadata-derived analytical metrics across the network per sample; red marker indicates " ~ labdata.lab.lab_cod ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of metadata-derived analytical metrics across participating laboratories ({{ comp_code }}).**  
Panels summarise the distribution of selected quantitative analytical parameters declared in the metadata template (e.g., genome coverage, depth of coverage, read composition, and variant counts). Grey distributions represent network-level variability, while the red marker corresponds to the values reported by **{{ labdata.lab.lab_cod }}**.

{% endfor %}

# Acknowledgement

We sincerely thank **{{ labdata.lab.lab_cod }}** for its participation in the 2026 RELECOV Dry-Lab EQA. The contribution of each laboratory is fundamental to maintaining analytical comparability, reproducibility, and interoperability across the network.

For any questions, technical clarifications, or follow-up discussions regarding this report, please contact the RELECOV WP.6 coordination team at bioinformatica@isciii.es.
{ endif }
