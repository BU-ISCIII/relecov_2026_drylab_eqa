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

## 5. General Results

A total of 52 laboratories comprising all the network were invited to participate. Of these, {{ general.total_participants }} laboratories {{ pct(general.total_participants_pct) }} submitted results for one or more components with the following distribution:

- SARS1 (SARS-CoV-2, Illumina): {{ general.participation_per_component.SARS1 }} laboratories.
- SARS2 (SARS-CoV-2, Oxford Nanopore Technologies): {{ general.participation_per_component.SARS2 }} laboratories.
- FLU1 (Influenza virus, Illumina): {{ general.participation_per_component.FLU1 }} laboratories.
- FLU2 (Influenza virus, Oxford Nanopore Technologies): {{ general.participation_per_component.FLU2 }} laboratories.

The median number of components analysed per laboratory was {{ general.median_components_analysed_per_lab }}.

This EQA provides the first fully drylab benchmarking of bioinformatic workflows across the RELECOV network, integrating both internationally validated datasets and purpose-designed in-silico scenarios.

### 5.1 Submission Completeness

Across all components:

- {{ pct(general.submission_rates_pct.fasta) }} of laboratories submitted consensus genome files (.fasta), where applicable
- {{ pct(general.submission_rates_pct.vcf) }} submitted variant call files (.vcf), where applicable
- The metadata template was completed at a mean completeness rate of {{ pct(general.metadata_completeness.mean_pct) }} per laboratory, with values ranging from {{ pct(general.metadata_completeness.min_pct) }} to {{ pct(general.metadata_completeness.max_pct) }}.

{% set fig_counter.value = fig_counter.value + 1 %}

As illustrated in Figure {{ fig_counter.value }}, metadata completeness varied across laboratories, with a mean completeness of {{ pct(general.metadata_completeness.mean_pct) }}.

{{ render_figure(network.figures.metadata_completeness_distribution,
  "Distribution of metadata completeness across participating laboratories.") }}

**_Figure {{ fig_counter.value }}_. Distribution of metadata completeness across participating laboratories**.
Boxplots represent the median and interquartile range of metadata completeness percentages. Whiskers denote the full observed range. The distribution reflects variability in reporting of analytical parameters, software versions, and controlled vocabulary adherence.

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

{{ render_figure(network.figures.consensus_summary, "Network-level consensus reconstruction performance summary.") }}

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

{{ render_figure(network.figures.variant_summary, "Network-level variant detection performance summary.") }}

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

Across components, median consensus concordance exceeded {{ pct(general.general_results.classification.median_assignment_concordance) }}, and classification accuracy remained above {{ pct(general.general_results.classification.median_assignment_accuracy) }}.

{{ render_figure(network.figures.classification_summary, "Network-level classification performance summary.") }}

**_Figure {{ fig_counter.value }}_. Distribution of classification outcomes across participating laboratories**.
Stacked bars represent the proportion of exact matches, minor discrepancies, and incorrect assignments relative to curated gold standard classifications for each component.

### 5.5 Metadata Quality and Interoperability

The evaluation of metadata focused on analytical transparency, reproducibility, and interoperability within the RELECOV network. Completeness and compliance were assessed according to the criteria defined in Section 6.4, including controlled vocabulary adherence, logical consistency, and reporting of analytical parameters.

#### Overall Completeness

{% set fig_counter.value = fig_counter.value + 1 %}

Across all participating laboratories, the metadata template was completed at a mean completeness rate of {{ pct(general.metadata_completeness.mean_pct) }}, with values ranging from {{ pct(general.metadata_completeness.min_pct) }} to {{ pct(general.metadata_completeness.max_pct) }}. As illustrated in Figure {{ fig_counter.value }}, metadata completeness varied across participating laboratories, with a heterogeneous distribution across predefined completeness ranges.

Optional analytical fields contributed disproportionately to incompleteness (TODO comprobar si es verdad), particularly those related to parameter specification and software versioning.

[Insert figure: Distribution of metadata completeness per laboratory]

**_Figure {{ fig_counter.value }}_. Distribution of metadata completeness across participating laboratories.** Bars represent the proportion of laboratories within different completeness ranges. These ranges reflect overall reporting quality and readiness for automated interoperability.

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
- {{ pct(general.metadata_completeness.free_text_predefine-pct) }} contained at least one free-text substitution where a predefined option was required.
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

Metadata was complete in {{ general.metadata_completeness.fully_compliant_pct }} of submissions, while {{ general.metadata_completeness.carification }} required clarification or correction during validation.

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

Variant detection demonstrated generally high sensitivity and Precission in high-quality samples. Nonetheless, variability increased in analytically challenging scenarios, including low read depth and mixed-site samples.

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

