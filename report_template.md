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
  - [4.3 Evaluation of Taxonomic and Phylogenetic Classification](#43-evaluation-of-taxonomic-and-phylogenetic-classification)
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
- [9. Individual Laboratory Technical Report – [Laboratory Name]](#9-individual-laboratory-technical-report--laboratory-name)
  - [9.1 Participation Overview](#91-participation-overview)
  - [9.2 SARS1 (SARS-CoV-2, Illumina)](#92-sars1-sars-cov-2-illumina)
  - [9.3 SARS2 (SARS-CoV-2, Nanopore)](#93-sars2-sars-cov-2-nanopore)
  - [9.4 FLU1 (Influenza, Illumina)](#94-flu1-influenza-illumina)
  - [9.5 FLU2 (Influenza, Nanopore)](#95-flu2-influenza-nanopore)

## Executive Summary</h2>
To be completed after final results are consolidated.

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

Table 1 summarises the correspondence between RELECOV EQA samples and their original source datasets, including ECDC ESIB references.

_**Table 1**. Overview of SARS-CoV-2 datasets used in the RELECOV 2026 Dry-Lab EQA.
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

Table 2 describes the design characteristics of in-silico samples, including virus composition and intended benchmarking challenges.

_**Table 2**. Viral, host and contaminant composition design of in-silico influenza datasets used for benchmarking._

| Sample | Influenza reads | Host reads | Additional Viral reads | Total reads | Analytical Challenge                |
|--------|-----------------|------------------|------------------|-------------|-------------------------------------|
| FLU2   | 1378764         | 462520 | 0                          | 1841284     | Baseline performance assessment     |
| FLU4   | 181626          | 300000 | 200000 SARS-CoV-2 reads    | 681626      | False positive control              |
| FLU5   | 1088000         | 100000 | 0                          | 1188000     | NA segment dropout                  |
| FLU7   | 5677            | 100    | 255 Rhinovirus reads       | 6032        | Cross-virus contamination challenge |
| FLU8   | 5380            | 300    | 0                          | 5680        | Baseline performance assessment     |
| FLU9   | 19989           | 500    | 0                          | 20489       | HA segment dropout                  |

Table 3 summarises the influenza datasets included in the EQA, detailing enrichment strategy, primer scheme, sequencing technology, and key analytical challenges.

_**Table 3**. Influenza virus samples used in the RELECOV 2026 Dry-Lab EQA, including sequencing platform, enrichment strategy, primer scheme, and key analytical features._

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
- Taxonomic and phylogenetic classification
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

For each laboratory and sample we calculated the following performance metrics: (TODO verificar si es verdad)
- Sensitivity = TP / (TP + FN)
- Precission = TP / (TP + FP)

Comparative analyses were performed to assess the influence of: (TODO verificar si es verdad)
- Allele frequency thresholds
- Minimum coverage thresholds
- Variant filtering criteria
- Reference genome selection


### 4.3 Evaluation of Taxonomic and Phylogenetic Classification

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

- Mean number of nucleotide discrepancies relative to the gold standard
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


