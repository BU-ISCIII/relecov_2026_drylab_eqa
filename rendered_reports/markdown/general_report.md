# Interlaboratory Comparison Exercise RELECOV 2.0 - Consolidation of WGS and RT-PCR activities for SARS-CoV-2 in Spain towards sustainable use and integration of enhanced infrastructure and capacities in the RELECOV network

##### Sarai Varona, Enrique Sapena, Pablo Mata, Alejandro Bernabéu, Pau Pascual, Magdalena Matito, Juan Ledesma, Emilia Arjona, Victor Lopez, Olga Dolgova, Sara Monzón, Isabel Cuesta

## Table of Contents

- [Executive Summary](#executive-summary)
- [1. Introduction](#1-introduction)
- [2. Scope of the Interlaboratory Comparison Exercise](#2-scope-of-the-interlaboratory-comparison-exercise)
- [3. Dataset Design and Selection Strategy](#3-dataset-design-and-selection-strategy)
    - [3.1. Rationale for Dataset Selection](#31-rationale-for-dataset-selection)
    - [3.2. SARS-CoV-2 Dataset](#32-sars-cov-2-dataset)
    - [3.3. Influenza Dataset](#33-influenza-dataset)
- [4. Methodology of Evaluation](#4-methodology-of-evaluation)
    - [4.1. Submission Completeness](#41-submission-completeness)
    - [4.2. Evaluation of Consensus Genome Reconstruction Performance](#42-evaluation-of-consensus-genome-reconstruction-performance)
    - [4.3. Evaluation of Variant Detection Accuracy](#43-evaluation-of-variant-detection-accuracy)
    - [4.4. Evaluation of Lineage, Subtype and Clade Assignment](#44-evaluation-of-lineage-subtype-and-clade-assignment)
    - [4.5. Evaluation of Metadata Completeness and Compliance](#45-evaluation-of-metadata-completeness-and-compliance)
- [5. General Results](#5-general-results)
    - [5.1. Submission Completeness](#51-submission-completeness)
    - [5.2. Consensus Genome Reconstruction Performance](#52-consensus-genome-reconstruction-performance)
    - [5.3. Variant Detection Accuracy](#53-variant-detection-accuracy)
    - [5.4. Lineage, Subtype and Clade Assignment](#54-lineage-subtype-and-clade-assignment)
    - [5.5. Metadata Completeness and Compliance](#55-metadata-completeness-and-compliance)
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
- [Appendix](#appendix)

## Executive Summary

The 2026 RELECOV Dry-Lab Interlaboratory Comparison Exercise provides the first network-wide dry-lab assessment focused specifically on bioinformatic analytical performance across respiratory virus genomic surveillance workflows. The exercise evaluated 20 distributed datasets grouped into four analytical components, comprising Illumina and Nanopore data for SARS-CoV-2 and influenza virus. Participating laboratories were assessed on consensus genome reconstruction, variant reporting, lineage/type and clade assignment, metadata completeness, and the reproducibility of their declared analytical workflows relative to curated gold standards and network-wide distributions.

Nineteen RELECOV member laboratories participated, corresponding to 36.54% of invited laboratories, with high submission rates for expected analytical outputs: 97.52% for consensus genome files and 89.11% for VCF files.

Across the network, consensus genome reconstruction performed slightly better in the Illumina-based components, with a combined median genome identity of 97.72%, compared with 96.17% in the Nanopore-based components. However, broad identity ranges in SARS2 and FLU2 indicate that outlier submissions remained present, particularly in contexts where masking, coverage thresholds, and consensus-generation choices differed across laboratories.

Variant reporting showed clear methodological heterogeneity. For SARS-CoV-2, the median number of discrepancies relative to the curated reference variant sets was 3 in the Illumina component and 3 in the Nanopore component. Influenza submissions were more heterogeneous structurally: the median number of high-frequency variants reported in metadata by laboratories was 477, whereas the corresponding median obtained from submitted VCF files after filtering to Allele Frequency (AF) >=75% was 188. The median discrepancy between metadata-reported values and VCF-derived values was 377, indicating that metadata declarations and filtered VCF content were not always directly concordant.

Classification performance was consistently higher for lineage/type assignment than for clade assignment. SARS-CoV-2 lineage concordance reached 77.7%, compared with 75.5% for clade assignment, while influenza type/subtype concordance reached 87.8%, compared with 67.3% for clade assignment. Review of submitted files further indicated that part of the excess discordance in clade assignment reflected field completion and nomenclature issues, including missing clade entries and lineage/type-like values entered in the clade field.

Metadata completeness and reporting remain major priorities for harmonisation. The median metadata completeness rate across participating laboratories was 59.1%, with values ranging from 11.3% to 92.0%. Although software names were reported for 71.8% of expected fields, only 60.1% of software-version fields, 45.6% of coverage thresholds, 38.9% of variant-calling parameter fields, and 57.7% of reference genome identifiers were completed. A total of 9 distinct workflows were identified, together with diversity in consensus, variant calling, and classification software.

Overall, the EQA demonstrates substantial analytical capability across the RELECOV network under the evaluated conditions, while also showing that interlaboratory comparability remains limited by heterogeneous thresholds, parameter reporting, reference selection, and uneven completion of metadata and QC fields. These findings support RELECOV 2.0 priorities of establishing minimum performance standards, strengthening metadata requirements, clarifying reporting rules for consensus and variants, and guiding component-aware recommendations for future harmonisation. A separate benchmarking deliverable document complements this report with a workflow-level comparison of analytical pipelines and software combinations.

## 1. Introduction

The RELECOV Network aims to strengthen genomic surveillance of respiratory viruses by developing and harmonising analytical capacities across the participating laboratories. In this context, it was essential to **assess the consistency, reproducibility and maturity of the bioinformatic workflows implemented across the network**.

To this end, an **Interlaboratory Comparison Exercise exercise in dry lab format** was conducted, based on the European Centre for Disease Prevention and Control (ECDC) 2024 dry-lab EQA. The exercise focused on the bioinformatic characterisation of respiratory viruses, covering key analytical tasks including viral genome reconstruction, variant identification, and lineage and clade assignment.

Beyond its role as an external quality assessment of laboratory performance, the exercise was also designed to support the methodological harmonisation objectives of RELECOV 2.0. A central component of this initiative was to characterise the diversity of analytical pipelines implemented across the RELECOV Network, evaluate their performance under the conditions of this exercise, and generate evidence to support future harmonisation activities within the network. This evaluation contributes directly to **Objective 2.1** of RELECOV 2.0, which focuses on _improving deep knowledge of the capacities and methodologies of the laboratories belonging to the network, as well as identifying a common methodology adapted to them and to the needs of the platform_. Furthermore, the exercise provides the practical evidence base required for **Task T6.1**, which aims to  _identify the most suitable bioinformatic analysis method for each sequencing platform, through an intercomparison exercise with simulated data for bioinformaticians_, in order to define the workflow that should be integrated into the RELECOV analytical platform.

The exercise was also aligned with **Milestone M6.3**, which pertains to _define sequencing and analysis protocols for each of the sequencing platforms_. In addition, the exercise provided operational insights relevant to **Task T6.5**, which addresses _the adaptation and improvement of the analysis pipeline for the different sequencing platforms used by the laboratories of the network_. It also contributed to **Task T6.4**, related to _sequence metadata annotation with ontologies, schema generation, parsing and validation_, by highlighting practical issues affecting metadata completeness, controlled-vocabulary use, and the consistency of reported analytical parameters.

The present report focuses on the interlaboratory assessment of analytical performance, reporting quality, and harmonisation needs. A separate benchmarking deliverable document complements this report with a workflow-level comparison of pipelines, software combinations, and parameter configurations.

The overall objective of the exercise was to **assess the bioinformatic performance of the participating laboratories, identify areas for improvement, and promote the adoption of consistent and comparable analytical practices across the network**. The outcomes presented in this report are expected to strengthen RELECOV’s preparedness and response capacity for routine surveillance and public health emergencies, while supporting the harmonisation objectives defined within RELECOV 2.0.

## 2. Scope of the Interlaboratory Comparison Exercise

The 2026 RELECOV Dry-Lab Interlaboratory Comparison Exercise was designed to evaluate the bioinformatic performance of laboratories participating in the RELECOV Network in the context of respiratory virus genomic surveillance.

Participating laboratories were provided with raw sequencing datasets corresponding to four independent analytical components:

- **SARS1**: Five SARS-CoV-2 samples sequenced using paired-end Illumina technology from the 2024 ECDC ESIB EQA
- **SARS2**: Five SARS-CoV-2 samples sequenced using Oxford Nanopore Technologies from the 2024 ECDC ESIB EQA
- **FLU1**: Five influenza virus samples sequenced using paired-end Illumina technology, 3 generated in-silico and 2 from the 2024 ECDC ESIB EQA.
- **FLU2**: Five influenza virus samples sequenced using Oxford Nanopore Technologies, 3 generated in-silico and 2 from the 2024 ECDC ESIB EQA.

Datasets were distributed as raw sequencing reads (.fastq files), and each component could be analysed independently, allowing laboratories to participate according to their technical capacity and routine workflow.

Laboratories were requested to submit the following deliverables:

- For each analysed sample:
    - One consensus genome sequence in `.fasta` format, containing exclusively the target viral genome reconstructed from the provided reads.
    - One or more variant call files in `.vcf` format, listing detected nucleotide variants relative to the reference genome selected by the laboratory.
- A completed harmonised metadata template, documenting analytical tools, software versions, reference genomes used, parameter settings, coverage thresholds, Lineage, Subtype or clade assignment tools, file names of submitted outputs, and the analytical decisions required to interpret and evaluate consensus reconstruction, variant reporting, lineage/type assignment, clade assignment, and quality control results. The values declared in this template were used throughout the evaluation to contextualise laboratory performance and to compare metadata-reported outputs against the submitted files. The template used in this exercise is available here: [Relecov_metadata_template_EQA2026.xlsx](https://github.com/BU-ISCIII/relecov_2026_drylab_eqa/blob/main/Relecov_metadata_template_EQA2026.xlsx).

The evaluation focused on core analytical tasks that are essential for routine genomic surveillance and public health response, including:

- **Viral genome reconstruction**: Generation of high-quality consensus genome sequences from raw sequencing reads produced using Illumina and Oxford Nanopore Technologies platforms.
- **Variant identification and reporting**: Detection and annotation of nucleotide variants relative to a chosen reference genome, including evaluation of filtering criteria, allele frequency thresholds, and variant file standardisation.
- **Lineage, Subtype and clade assignment**: Accurate classification of reconstructed genomes using established nomenclature systems and version-controlled databases.
- **Metadata reporting and interoperability**: Completion of a harmonised metadata template capturing software versions, analytical parameters, reference genome selection, and file traceability, ensuring compatibility with automated validation and integration into the RELECOV analytical platform.
- **Quality control assesment**: Evaluation of laboratory quality-control practices, including the interpretation of sequencing quality, the identification of analytical limitations, the application of quality thresholds, and the completeness and consistency of quality-control information reported throughout the submitted results.

## 3. Dataset Design and Selection Strategy

### 3.1. Rationale for Dataset Selection

The 2026 RELECOV Dry-Lab Interlaboratory Comparison Exercise was specifically designed for microbiology laboratories participating in the RELECOV Network, which perform genomic surveillance of respiratory viruses in Spain.

Sample selection followed three guiding principles:

- Representation of realistic genomic surveillance scenarios.
- Inclusion of predefined analytical challenges.
- Ensuring methodological robustness and representativeness.

Datasets were derived from two sources:

- Reused datasets from the 2024 ECDC ESIB Dry-Lab EQA.
- Newly generated in-silico datasets constructed to simulate seasonal human influenza circulation.

The integration of both sources allowed alignment with internationally validated materials while tailoring the exercise to the operational reality of RELECOV clinical laboratories.

### 3.2. SARS-CoV-2 Dataset

SARS-CoV-2 datasets were selected from the 2024 ECDC ESIB EQA to ensure comparability with internationally benchmarked material. Both Illumina and Nanopore panels included samples representing:

- High-quality baseline genomes.
- Low read-depth scenarios.
- Samples with numerous mixed sites.
- Contamination with non-target viral reads.
- Lineages of epidemiological relevance (e.g., recombinant or XBB-related lineages).

Only samples generated using the same ARTIC primer scheme (v4.1) were selected to avoid introducing variability associated with enrichment panel differences. This ensured that observed performance differences reflect analytical workflow characteristics rather than enrichment strategies heterogeneity.

_**Table 1**. Overview of SARS-CoV-2 datasets used in the RELECOV 2026 Dry-Lab Interlaboratory Comparison Exercise .
The table details sample origin, sequencing technology (Illumina paired-end or Oxford Nanopore Technologies), amplicon primer scheme version, and specific analytical characteristics intentionally selected to assess workflow robustness under challenging conditions._

| Sample | Component | Source             | Platform | Ref sample | Key Feature                                       | FASTQ files | Read layout | Clade Assignment | Lineage Assignment | Quality check |
|--------|-----------|--------------------|----------|------------|---------------------------------------------------|-------------|-------------|------------------|--------------------|---------------|
| SARS1  | SARS1     | ECDC-ESIB EQA 2024 | Illumina | SARS2.04   | Influenza virus sample with some SARS-CoV-2 reads | 2           | Paired-end  | -                | -                  | Bad           |
| SARS2  | SARS1     | ECDC-ESIB EQA 2024 | Illumina | SARS2.01   | High-quality baseline sample                      | 2           | Paired-end  | 21K              | BA.1.13            | Ok            |
| SARS3  | SARS1     | ECDC-ESIB EQA 2024 | Illumina | SARS2.16   | XBB sample / insertion challenge                  | 2           | Paired-end  | 23A              | XBB.1.5            | Ok            |
| SARS4  | SARS1     | ECDC-ESIB EQA 2024 | Illumina | SARS2.20   | Very low read depth                               | 2           | Paired-end  | 22E              | BQ.1.1             | Bad           |
| SARS5  | SARS1     | ECDC-ESIB EQA 2024 | Illumina | SARS2.13   | >10 mixed sites                                   | 2           | Paired-end  | 21K              | BA.1.1             | Bad           |
| SARS6  | SARS2     | ECDC-ESIB EQA 2024 | Nanopore | SARS1.01   | High-quality baseline sample                      | 1           | Single-end  | recombinant      | XCH.1              | Ok            |
| SARS7  | SARS2     | ECDC-ESIB EQA 2024 | Nanopore | SARS1.09   | XBB sample / ambiguity next to a deletion         | 1           | Single-end  | 23A              | XBB.1.5.24         | Ok            |
| SARS8  | SARS2     | ECDC-ESIB EQA 2024 | Nanopore | SARS1.15   | >10 mixed sites                                   | 1           | Single-end  | 23D              | XBB.1.9.1          | Bad           |
| SARS9  | SARS2     | ECDC-ESIB EQA 2024 | Nanopore | SARS1.12   | Influenza virus sample with some SARS-CoV-2 reads | 1           | Single-end  | -                | -                  | Bad           |
| SARS10 | SARS2     | ECDC-ESIB EQA 2024 | Nanopore | SARS1.05   | Very low read depth                               | 1           | Single-end  | -                | -                  | Bad           |

### 3.3. Influenza Dataset

The influenza datasets provided in the 2024 ECDC ESIB EQA predominantly correspond to zoonotic influenza strains of animal origin, including H5N1, H5N6, and reassortant genomes.

While these datasets are valuable for specialised surveillance contexts, they do not represent the routine analytical scenario encountered by most RELECOV laboratories, which primarily process seasonal human Influenza A/H1N1 and A/H3N2.

Given that the objective of this Interlaboratory Comparison Exercise is to benchmark bioinformatic workflows in a clinical hospital environment, it was considered methodologically necessary to include representative seasonal human influenza strains.

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

These design features ensured that individual analytical challenges could be evaluated independently while maintaining biologically plausible sequencing datasets.

_**Table 2**. Viral, host and contaminant composition design of in-silico influenza datasets used for benchmarking._

| Sample | Influenza reads | Host reads | Additional Viral reads  | Total reads | Analytical Challenge                |
|--------|-----------------|------------|-------------------------|-------------|-------------------------------------|
| FLU2   | 1378764         | 462520     | 0                       | 1841284     | Baseline performance assessment     |
| FLU4   | 181626          | 300000     | 200000 SARS-CoV-2 reads | 681626      | Contamination with SARS-CoV-2       |
| FLU5   | 1088000         | 100000     | 0                       | 1188000     | NA segment dropout                  |
| FLU7   | 5677            | 100        | 255 Rhinovirus reads    | 6032        | Cross-virus contamination challenge |
| FLU8   | 5380            | 300        | 0                       | 5680        | Baseline performance assessment     |
| FLU9   | 19989           | 500        | 0                       | 20489       | HA segment dropout                  |

_**Table 3**. Influenza virus samples used in the RELECOV 2026 Dry-Lab Interlaboratory Comparison Exercise, including sequencing platform, enrichment strategy, primer scheme, and key analytical features._

| Sample | Component | Source    | Platform | Enrichment Strategy | Primer Scheme                                   | Read Layout | Ref_sample        | Type   | Clade HA  | Legacy Clade       | Key Feature                             | Quality check |
|--------|-----------|-----------|----------|---------------------|-------------------------------------------------|-------------|-------------------|--------|-----------| ------------------ | ----------------------------------------|---------------|
| FLU1   | FLU1      | ESIB 2024 | Illumina | Amplicon            | CommonUni12/13 (Van den Hoecke 2015)            | Paired-end  | INFL2.07          | A/H5N1 | 2.3.4.4b  | -                  | High-quality baseline sample (zoonotic) | Ok            |
| FLU2   | FLU1      | In-silico | Illumina | Amplicon            | Zhou 2009 single-reaction genomic amplification | Paired-end  | In-silico Sample1 | A/H1N1 | D.3.1.1   | 6B.1A.5a.2a.1      | High-quality baseline sample (human)    | Ok            |
| FLU3   | FLU1      | ESIB 2024 | Illumina | No enrichment       | —                                               | Paired-end  | INFL2.04          | —      | —         | -                  | No influenza (Rhinovirus only)          | Bad           |
| FLU4   | FLU1      | In-silico | Illumina | Amplicon            | Zhou 2009 single-reaction genomic amplification | Paired-end  | In-silico Sample3 | A/H3N2 | K         | 3C.2a1b.2a.2a.3a.1 | Contamination with SARS-CoV-2           | Ok            |
| FLU5   | FLU1      | In-silico | Illumina | Amplicon            | Zhou 2009 single-reaction genomic amplification | Paired-end  | In-silico Sample4 | A/H3N2 (A/H3 or A/H3Nx) | J.2.2     | 3C.2a1b.2a.2a.3a.1 | NA segment dropout                      | Bad           |
| FLU6   | FLU2      | ESIB 2024 | Nanopore | No enrichment       | —                                               | Single-end  | INFL1.02          | A/H5N6 | 2.3.4.4h  | -                  | High-quality baseline sample (zoonotic) | Ok            |
| FLU7   | FLU2      | In-silico | Nanopore | Amplicon            | Zhou 2009 single-reaction genomic amplification | Single-end  | In-silico Sample2 | A/H1N1 | C.1.9.3   | 6B.1A.5a.2a        | Contamination with Rhinovirus           | Ok            |
| FLU8   | FLU2      | In-silico | Nanopore | Amplicon            | Zhou 2009 single-reaction genomic amplification | Single-end  | In-silico Sample3 | A/H3N2 | K         | 3C.2a1b.2a.2a.3a.1 | High-quality baseline sample (human)    | Ok            |
| FLU9   | FLU2      | In-silico | Nanopore | Amplicon            | Zhou 2009 single-reaction genomic amplification | Single-end  | In-silico Sample1 | A/H1N1 (A/N1 or A/HxN1) | D.3.1.1   | 6B.1A.5a.2a.1      | HA segment dropout                      | Bad           |
| FLU10  | FLU2      | ESIB 2024 | Nanopore | Amplicon            | CommonUni12/13 (Van den Hoecke 2015)            | Single-end  | INFL1.08          | A/H5N1 | 2.3.4.4b  | -                  | High-quality baseline sample (zoonotic) | Ok            |

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

<div class="equation"><span class="equation-lhs">File submission rate</span><span class="equation-equals">=</span><span class="equation-fraction"><span class="equation-numerator">Number of submitted files</span><span class="equation-denominator">Total number of expected files</span></span></div>

Missing files were recorded but not penalised beyond descriptive reporting, as laboratories were allowed to participate selectively according to local analytical capacity.

### 4.2. Evaluation of Consensus Genome Reconstruction Performance

For each sample, a curated consensus sequence was provided by the ECDC or generated in silico as previously explained. Hereafter, this consensus sequence will be referred to as the **gold standard**. For the in silico samples, the gold standard corresponds to the original FASTA genome used to simulate the sequencing reads. Influenza gold standard consensus sequences generated for this exercise are available in the project GitHub repository. Gold standard consensus sequences corresponding to ECDC-provided samples are not distributed with the repository, as we are not authorised to provide those datasets.

All submitted consensus sequences (.fasta) were:

- Aligned against the corresponding gold standard genome sequence using [Mafft v7.475](https://mafft.cbrc.jp/alignment/software/).
- Compared position-by-position relative to the declared gold standard coordinate system.

For SARS-CoV-2 samples, evaluated positions were reported relative to the Wuhan reference genome coordinate system. For influenza virus samples, evaluated positions were reported relative to the curated gold standard sequence for each segment. In the influenza gold standards, primer-binding regions were masked to allow evaluation of laboratories that retained primer-derived regions in the submitted consensus sequence. Influenza positions showing two alternative alleles at approximately balanced frequencies (40-60% allele frequency for each allele) were represented using IUPAC ambiguity codes in the gold standard. At these positions, either the corresponding ambiguity code or either of the two represented nucleotides was accepted as concordant with the gold standard.

Differences between submitted sequences and gold standard sequences were categorised into the following classes:

- **Wrong nucleotide**: A nucleotide different from the allowed reference or ambiguity code.
- **Ambiguity instead of nucleotide**: Ambiguity codes introduced where a defined nucleotide was expected.
- **Nucleotide instead of ambiguity**: Defined nucleotide provided where an ambiguity code was expected.
- **Nucleotide stretch instead of stretch of Ns**: Defined bases provided where Ns were expected.
- **Stretch of Ns instead of nucleotide stretch**: Continuous region of Ns where defined bases were expected.
- **Insertion relative to gold standard**: One or more nucleotides present in the submitted sequence at a position where no corresponding bases are present in the gold standard.
- **Deletion relative to gold standard**: One or more nucleotides absent from the submitted sequence at a position where corresponding bases are present in the gold standard.

Each insertion, deletion, or contiguous stretch of Ns was counted as a single event.

For each laboratory and sample, the following summary metrics were compiled:

- Total number of nucleotide discrepancies
- Percentage genome identity relative to the curated gold standard

The proportional contribution of each discrepancy category was calculated relative to the total number of discrepancies observed per component. In the results, "dominant" refers to the discrepancy category with the highest median burden across evaluable observations, rather than the category with the single highest maximum value or the largest cumulative total. This choice reduces sensitivity to a small number of extreme outliers and is intended to reflect the most typical discrepancy pattern observed across laboratories. These patterns should be interpreted in the context of threshold choices, reference-genome selection, and declared pipeline configurations.

### 4.3. Evaluation of Variant Detection Accuracy

#### 4.3.1. SARS-CoV-2

A curated reference variant set was generated for each SARS-CoV-2 sample. Variant positions were standardized relative to a defined coordinate system referred to the references used by Nextclade.

Submitted .vcf files were:

- Converted to a standardised long table format for coordinate comparison
- Compared position-by-position with the reference variant set

Differences between submitted variants and reference variant set were categorised into the following classes:

- **Wrong nucleotide**: A nucleotide different from the allowed reference or ambiguity code.
- **Insertion relative to gold standard**: A variant call indicating inserted nucleotides in the submitted variant set at a position where no insertion is present in the gold standard variant set.
- **Deletion relative to gold standard**: A variant call indicating deleted nucleotides in the submitted variant set at a position where no deletion is present in the gold standard variant set.
- **Missing variant**: Variants present in the reference but missing in the sample.
- **De novo**: Variants present in the sample but missing in the reference set.

Each insertion, deletion, or contiguous stretch of Ns was counted as a single event.

For each laboratory and sample, the total number of nucleotide discrepancies was calculated. In the results, "dominant" refers to the discrepancy category with the highest median burden across evaluable observations, rather than the category with the single highest maximum value or the largest cumulative total. This choice reduces sensitivity to a small number of extreme outliers and is intended to reflect the most typical discrepancy pattern observed across laboratories. These patterns should be interpreted in the context of threshold choices, reference-genome selection, and declared pipeline configurations, all of which can shift the balance between successful hits, missing expected variants, and de novo calls even when laboratories analyse the same raw data.

The proportional contribution of each discrepancy category was calculated relative to the total number of discrepancies observed per component.

Metadata describing the following analytical settings were collected to support result interpretation:

- Allele frequency thresholds
- Minimum coverage thresholds
- Reference genome selection

#### 4.3.2 Influenza

For influenza virus datasets, direct position-by-position comparison of reported variants against the curated reference variant set was not feasible under the same framework applied to SARS-CoV-2.

Unlike SARS-CoV-2, where laboratories predominantly use a globally standardised reference genomes (either MN908947.3 or NC_045512.2), influenza virus analyses exhibited heterogeneity in reference genome selection. As a result:

- Variant coordinates were reported relative to different reference accessions.
- Segment boundaries and numbering schemes varied.
- Insertions and deletions were represented inconsistently across reference backbones.

This heterogeneity prevented robust coordinate harmonisation across submissions without introducing alignment-dependent artifacts and interpretation bias.

#### 4.3.3. Descriptive and Structural Variant Reporting Metrics

In addition to nucleotide-level discrepancy analysis for SARS-CoV-2, both SARS-CoV-2 and influenza submissions were evaluated using descriptive and structural reporting metrics to characterise reporting behaviour and methodological heterogeneity across laboratories.

For both viruses, the following reporting practice metrics were collected:

- Number of laboratories reporting high-frequency variants only.
- Number of laboratories reporting both high- and low-frequency variants.
- Number of laboratories reporting exclusively low-frequency variants.
- Total number of distinct reference genomes employed for variant calling or mapping.

For influenza virus, additional structural summary metrics were calculated because direct coordinate-harmonised comparison of all submitted variants was not methodologically robust across segment-specific references:

- Number of variants with an allele frequency higher than 75% reported in metadata.
- Number of variants with an allele frequency higher than 75% derived from submitted VCF files.
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

### 4.5. Evaluation of Metadata Completeness and Compliance

Metadata assessment focused on analytical transparency and interoperability rather than biological correctness. Before the start of the exercise, a metadata template with controlled-vocabulary dropdowns was distributed among the laboratories to review the available options and incorporation of missing software tools into the schema.

For each submitted sample, metadata completeness was calculated as:

<div class="equation"><span class="equation-lhs">Sample metadata completeness</span><span class="equation-equals">=</span><span class="equation-fraction"><span class="equation-numerator">Number of correctly populated fields</span><span class="equation-denominator">Total number of applicable fields</span></span></div>

Laboratory-level and component-level completeness summaries were then derived from these sample-level values. Fields were evaluated for:

- Completion: Each sample has a list of minimum **recommended** fields, based on the sample characteristics. For each component/sample/lab the total number of completed minimum **recommended** fields was evaluated. Both mandatory and optional analytical fields were included in the completeness assessment, while fields not applicable to a laboratory’s selected components were excluded from scoring.
- Compliance with controlled vocabularies. Metadata entries were considered non-compliant when:
    - Controlled vocabulary options were bypassed
    - Free-text substitutions replaced defined values
- Valid file name reporting

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

<div class="equation"><span class="equation-lhs">QC concordance rate</span><span class="equation-equals">=</span><span class="equation-fraction"><span class="equation-numerator">Number of Matches</span><span class="equation-denominator">Total QC evaluations</span></span></div>

QC evaluations were calculated only for samples analysed by the laboratory.

The QC assessment evaluation was limited to concordance analysis. The exercise did not attempt to infer the internal QC criteria applied by laboratories, but rather assessed agreement with the predefined gold standard QC status to evaluate interpretative consistency across the network.

## 5. General Results

A total of 52 laboratories within the RELECOV network were invited to participate. Of these, 19 laboratories (36.5%) submitted results for one or more components with the following distribution:

- SARS1 (SARS-CoV-2, Illumina): 16 laboratories.
- SARS2 (SARS-CoV-2, Oxford Nanopore Technologies): 10 laboratories.
- FLU1 (Influenza virus, Illumina): 12 laboratories.
- FLU2 (Influenza virus, Oxford Nanopore Technologies): 10 laboratories.

The median number of components analysed per participating laboratory was 2.

The results presented below are interpreted according to the evaluation framework described in [Section 4](#4-methodology-of-evaluation).

### 5.1. Submission Completeness

Across all components:

- 97.5% of laboratories submitted consensus genome files (.fasta), where applicable.
- 89.1% submitted variant call files (.vcf), where applicable.

Component-level submission totals are presented in Section 6 and reflect both the number of participating laboratories and the expected output files for each dataset.

### 5.2. Consensus Genome Reconstruction Performance

Across the two Illumina-based components, the combined median genome identity was 97.72%, compared with 96.17% across the two Nanopore-based components. When grouped by virus, the combined median genome identity was 99.59% across the SARS-CoV-2 components and 95.84% across the influenza components. Notably, the influenza median was also slightly lower than the combined Nanopore-based median, indicating that virus-specific analytical complexity likely contributed in addition to platform-related differences. Nanopore-based datasets also showed broader overall identity ranges, where low-identity outliers were present.

Dominant discrepancy patterns differed by component:

- In SARS1, the most frequent discrepancy category was defined nucleotides in the submitted consensus where stretches of Ns were present in the gold standard (as shown in Appendix Table 3 and Appendix Figure 1).
- In SARS2, the most frequent discrepancy category was stretches of Ns in the submitted consensus where defined nucleotides were present in the gold standard (as shown in Appendix Table 10 and Appendix Figure 3).
- FLU1 and FLU2 were both dominated by deletions relative to the gold standard (as shown in Appendix Table 17 and Appendix Table 23 and Appendix Figure 5 and Appendix Figure 6).

Across components, many discrepancy categories had medians of zero, indicating that errors tended to be concentrated in a smaller number of laboratories or samples rather than being uniformly distributed across the network.

<figure>
<img src="figures/network/consensus_summary.png" alt="Network-level consensus reconstruction performance summary." style="width: 98%; max-width: 98%;"/>
</figure>

**_Figure 1_. Consensus genome reconstruction performance across components**. Panel **A** shows the distribution of nucleotide discrepancies relative to the gold standard across components, and panel **B** shows the corresponding distribution of genome identity values. In both panels, the central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

### 5.3. Variant Detection Accuracy

#### 5.3.1. SARS-CoV-2

For SARS-CoV-2 components (SARS1 and SARS2), variant detection accuracy was assessed against curated reference variant sets. Overall, submitted VCFs showed a median number of 3 discrepancies relative to the reference variant set for both Illumina and Nanopore components.

Variant detection performance differed across components (Figure 2). Contextual factors documented in the metadata that may contribute to these differences included:

- Allele frequency thresholds used for incorporation into vcf files
- Variant normalization practices (variant caller software and params)

<figure>
<img src="figures/network/variant_summary.png" alt="Network-level variant detection performance summary." style="width: 70%; max-width: 70%;"/>
</figure>
**_Figure 2_. SARS-CoV-2 network-level variant detection performance summary**. Boxplots represent the number of variant discrepancies per SARS-CoV-2 component across participating laboratories. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

Variant evaluation included structural reporting characteristics and methodological heterogeneity. As shown in Figure 3:

- 66.7% of laboratories reported both high- and low-frequency variants.
- 33.3% reported only high-frequency variants.
- No laboratories reported exclusively low-frequency variants.

Additionally, a total of 3 distinct reference genomes were employed for variant calling across SARS-CoV-2 components (MN908947.3, NC_045512.2, XBB REFERENCE GENOME).

<figure>
<img src="figures/network/sars_variant_reporting_summary.png" alt="SARS-CoV-2 variant reporting practices across the network." style="width: 70%; max-width: 70%;"/>
</figure>

**_Figure 3_. SARS-CoV-2 variant reporting characteristics across the network**. Summarise the proportion of laboratories reporting high- and/or low-frequency variants.

#### 5.3.2. Influenza virus

For influenza virus components (FLU1 and FLU2), variant evaluation focused on structural reporting characteristics and methodological heterogeneity.

As shown in Figure 4:

- 50.7% of laboratories reported both high- and low-frequency variants.
- 12.3% reported only high-frequency variants.
- 37.0% reported exclusively low-frequency variants.

Additionally, an estimated total of 8 distinct reference genomes were employed for variant calling or mapping across influenza components. This value was rounded to the nearest whole genome by dividing the total number of distinct fragment references (63) by 8 influenza genome segments.

Structural summary metrics derived from submitted influenza VCF files are presented in Table 4. These metrics capture the overall magnitude of reported variants in the metadata file and the discrepancy between reported variants with an allele frequency >= 75% in the metadata file and the VCF file, rather than direct nucleotide-level accuracy against a unified reference coordinate system.

**Table 4. Network-level structural summary of influenza variant reporting.**

| Metric | Network median | Min-max |
|---|---:|---:|
| Variants with AF>=75% | 477 | 377–1797 |
| Variants with AF>=75% in VCF | 188 | 0–1373 |
| Discrepancies in reported variants | 377 | 0–1797 |
| Total variants in VCF | 526 | 0–7903 |

<figure>
<img src="figures/network/influenza_variant_reporting_summary.png" alt="Influenza variant reporting practices across the network." style="width: 96%; max-width: 96%;"/>
</figure>

**_Figure 4_. Influenza variant reporting characteristics across the network**. Summarise the proportion of laboratories reporting high- and/or low-frequency variants.

Together, these results show heterogeneity in influenza variant reporting within the network.

### 5.4. Lineage, Subtype and Clade Assignment

Overall concordance rates were:

- SARS-CoV-2 lineage assignment: **77.7%** concordance.
- Influenza type/subtype identification: **87.8%** concordance.
- SARS-CoV-2 clade assignment: **75.5%** concordance.
- Influenza clade assignment: **67.3%** concordance.

Across components, lineage/type concordance was consistently higher than clade concordance. SARS-CoV-2 lineage assignment reached 77.7%, compared with 75.5% for SARS-CoV-2 clade assignment, while influenza type/subtype identification reached 87.8% compared with 67.3% for influenza clade assignment.

<figure>
<img src="figures/network/classification_summary.png" alt="Distribution of classification outcomes across participating laboratories." style="width: 98%; max-width: 98%;"/>
</figure>

**_Figure 5_. Distribution of classification outcomes across participating laboratories.** Panel **A** shows **lineage/type assignments**, and panel **B** shows **clade assignments**. Stacked bars represent the percentage of all possible sample-level classifications across participating laboratories for each component. Bars are partitioned into **Match** (correct assignments relative to the curated gold standard), **Discrepancy** (incorrect assignments), and **Not provided** (classification not reported).

### 5.5. Metadata completeness and compliance

The evaluation of metadata focused on analytical transparency, reproducibility, and interoperability within the RELECOV network, including controlled vocabulary adherence, logical consistency, and reporting of analytical parameters.

#### Overall Completeness

Across all participating laboratories, the metadata template was completed at a median completeness rate of 59.1%, with values ranging from 11.3% to 92.0%. Component-level median completeness values were similar overall, but the observed ranges remained broad in all components (Figure 6). The leading incompleteness drivers were variant calling, pre-processing, and mapping fields, followed by QC metrics, de-hosting, and consensus analysis fields.

Most frequent incompleteness drivers across the network:
<ul class="compact-list">

<li>Variant calling fields (missing in 18 laboratories)</li>

<li>Pre-processing fields (missing in 18 laboratories)</li>

<li>Mapping fields (missing in 17 laboratories)</li>

<li>QC metrics fields (missing in 17 laboratories)</li>

<li>De-hosting fields (missing in 16 laboratories)</li>

</ul>

<figure>
<img src="figures/network/metadata_completeness_distribution.png" alt="Distribution of metadata completeness across participating laboratories." style="width: 80%; max-width: 80%;"/>
</figure>

**_Figure 6_. Distribution of metadata completeness across participating laboratories**. Boxplots represent the distribution of sample-level metadata completeness percentages across the different components. Completeness was calculated for each submitted sample as the proportion of filled metadata fields relative to the total number of maximum expected metadata fields. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

#### Reporting of Analytical Parameters

Although core pipeline tools were generally reported, variability was observed in the level of parameter detail provided.

- 71.8% of the maximum software-name fields were completed across submitted samples.
- 60.1% of the maximum software-version fields were completed across submitted samples.
- 45.6% specified minimum coverage thresholds.
- 38.9% reported variant calling parameters, containing the potential allele frequency thresholds.
- 57.7% reported the reference genome accession or identifier.

The incomplete reporting of parameters limited the ability to fully reconstruct or reproduce analytical workflows in 60.7% of submissions.

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

Overall, the network achieved 71.1% QC concordance, corresponding to 64 Matches and 26 Discrepancies across 90 evaluated sample-level QC decisions.

QC concordance differed across components, ranging from 62.5% in SARS1 to 100.0% in FLU2, based on reported QC information (Figure 7).

<figure>
<img src="figures/network/qc_match_rate_by_component.png" alt="QC concordance by component (Match, Discrepancy, and Not provided relative to the gold standard)." style="width: 80%; max-width: 80%;"/>
</figure>

**_Figure 7_. QC concordance by component relative to the gold standard.** Stacked bars represent the proportion of sample-level QC outcomes classified as Match, Discrepancy, or Not provided for each component across participating laboratories. Not provided values correspond to missing QC assessments and are shown separately from true discrepancies.

## 6. Component-specific Results

This section presents the analytical results disaggregated by component, allowing a detailed assessment of performance within each dataset and sequencing technology. For each component, results are structured according to participation and submission metrics, consensus genome reconstruction performance, variant detection accuracy, and Lineage, Subtype or clade assignment concordance, as applicable.

Component-level analyses enable identification of platform-specific patterns, dataset-dependent challenges, and variability associated with particular sample characteristics. This approach facilitates a more granular interpretation of performance differences observed at the network level and supports targeted harmonisation recommendations.

All component-level results below are reported using the same evaluation framework described in [Section 4](#4-methodology-of-evaluation).

### 6.1. SARS1 (SARS-CoV-2, Illumina)

#### 6.1.1. Participation and Submissions

A total of 16 laboratories submitted results for the SARS1 component:

- A total of 64 consensus genome sequences (.fasta) were submitted.
- A total of 64 variant call files (.vcf) were submitted.
- The metadata template completeness for SARS1 submissions had a median of 55.8%.

Most frequent incompleteness drivers in SARS1:
<ul class="compact-list">

<li>Lineage assignment fields (missing in 14 laboratories)</li>

<li>SARS-CoV-2 QC metrics (missing in 14 laboratories)</li>

<li>Mapping fields (missing in 12 laboratories)</li>

<li>QC metrics fields (missing in 12 laboratories)</li>

<li>Variant calling fields (missing in 12 laboratories)</li>

</ul>

#### 6.1.2. Consensus Genome Reconstruction Performance

Consensus sequences were evaluated against the corresponding curated gold standard for each sample in the SARS1 component.

Overall, SARS1 showed a median genome identity of 99.59%, with a median of 3 nucleotide discrepancies per sample (range: 1–125) (Figure 8).

<figure>
<img src="figures/SARS1/consensus_discrepancies_boxplot_by_sample.png" alt="Consensus discrepancies per sample for SARS1 relative to the curated gold standard." style="width: 90%; max-width: 90%;"/>
</figure>

**Figure 8. Consensus reconstruction performance by sample for SARS1.** Panel A shows the distribution of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample, and Panel B shows the corresponding distribution of genome identity values. In both panels, the central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. In Panel B, the y-axis is truncated to highlight differences among high-identity values.

Figure 9 presents the distribution of nucleotide discrepancy types per sample across participating laboratories for SARS1.

<figure>
<img src="figures/SARS1/consensus_discrepancies_stacked_by_sample.png" alt="Consensus discrepancy types per sample for SARS1 relative to the curated gold standard." style="width: 80%; max-width: 80%;"/>
</figure>

**Figure 9. Consensus discrepancy type composition per sample for SARS1.** Stacked bars represent the number and type of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample.

The dominant discrepancy pattern observed in SARS1 was Nucleotide stretch instead of stretch of Ns (Figure 9). Sample-level consensus reconstruction summary metrics are provided in Appendix Table 1. A full sample-level breakdown of discrepancy categories is provided in Appendix Table 2, while the aggregated discrepancy composition by type and the corresponding category-wise boxplot can be found in Appendix Table 3 and Appendix Figure 1, respectively.

#### 6.1.3. Variant Detection Accuracy

Variant call files (.vcf) submitted for the SARS1 component were compared against the curated reference variant set corresponding to each sample in the SARS1 component.

Overall, SARS1 showed a median of 3 variant discrepancies per sample (range: 0–518), together with a median of 65 successful hits per sample (Table 5, Figure 10).

<figure>
<img src="figures/SARS1/variant_discrepancies_stacked_by_sample.png" alt="Variant discrepancies per sample for SARS1 relative to the curated gold standard." style="width: 80%; max-width: 80%;"/>
</figure>

**Figure 10. Distribution of variant discrepancies per sample for SARS1.** Stacked bars represent the number of nucleotide discrepancies and discrepancy types relative to the curated gold standard across participating laboratories for each sample.

**Table 5. Network-level SARS-CoV-2 variant reporting metrics per sample for SARS1.**

| Sample ID | Expected hits | Median successful hits (n=15) | Median variants >=75% AF in metadata (n=5) | Median variants >=75% AF in VCF (n=16) | Median variants with effect in metadata (n=5) | Median variants with effect in VCF (n=14) | Median discrepancies metadata vs VCF | Median effect discrepancies metadata vs VCF |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SARS1 | 36 | 32 | 30 | 30 | 22 | 22 | 0 | 2 |
| SARS2 | 65 | 64 | 61 | 61 | 45 | 45 | 0 | 1.5 |
| SARS3 | 95 | 92 | 92 | 92 | 66 | 64 | 0 | 2.5 |
| SARS4 | 75 | 75 | 70 | 73 | 52 | 50.5 | 1 | 1 |
| SARS5 | 65 | 64 | 56 | 56 | 46 | 42 | 0 | 3.5 |

At component level, 5 laboratories reported the number of variants in the metadata, whereas 11 did not report this field for any sample in SARS1. The median number of variants with an allele frequency (AF) >=75% was 61 in the metadata and 64 in the submitted VCF files (Table 5, Figure 11).

Figure 11 summarises the distribution of declared variant reporting modes across submitted sample outputs in SARS1.

<figure>
<img src="figures/SARS1/variant_reporting_practice_by_component.png" alt="Variant reporting practices for SARS1." style="width: 70%; max-width: 70%;"/>
</figure>

**Figure 11. Variant reporting practices for SARS1.** Bars represent the proportion of submitted sample outputs classified as high and low frequency reporting, high frequency only, or low frequency only, according to the metadata declarations associated with the variant outputs for this component.

The dominant discrepancy pattern observed in SARS1 was Missing variant.  The full sample-level variant calling profile is provided in Appendix Table 4, while the aggregated discrepancy composition by type and the corresponding category-wise boxplot can be found in Appendix Table 5 and Appendix Figure 2, respectively.

#### 6.1.4. Lineage, Subtype and Clade Assignment

Lineage, subtype and clade assignments submitted for the SARS1 component were evaluated for concordance with the curated gold standard classifications.

Across all participating laboratories, lineage/subtype concordance reached 85.9%, whereas clade concordance reached 68.8%. The sample-level outcome distribution also shows that part of the observed discordance was associated with missing classifications or inconsistent completion of classification fields rather than with uniform analytical failure across all submissions.

<figure>
<img src="figures/SARS1/typing_outcome_stackedbar_by_sample.png" alt="Classification outcome distribution per sample for SARS1." style="width: 98%; max-width: 98%;"/>
</figure>

**Figure 12. Classification outcome distribution per sample for SARS1.** Panel A shows the proportion of lineage/subtype assignment Match, Discrepancy, and Not provided outcomes across participating laboratories for each sample. Panel B shows the corresponding proportions for clade assignments. Percentages are calculated over all participating laboratories in the component, so the Not provided segment captures samples for which lineage/subtype or clade information was not reported. Detailed sample-level concordance percentages are provided in Appendix Table 6.

#### 6.1.5. Sample Quality Control Assessment

Sample-level QC in SARS1 was evaluated as concordance between the laboratory-reported Pass/Fail classification and the predefined gold standard status. QC concordance was heterogeneous across samples, and some laboratories did not report a formal QC assessment. Network-wide concordance for reported QC decisions was 62.5%.

<figure>
<img src="figures/SARS1/qc_match_by_sample.png" alt="Sample-level QC concordance for SARS1 (Match, Discrepancy, and Not provided relative to the gold standard)." style="width: 90%; max-width: 90%;"/>
</figure>

**_Figure 13_. Sample-level QC concordance for SARS1 relative to the gold standard.** Bars represent the proportion of Match, Discrepancy, and Not provided outcomes per sample across participating laboratories. Higher discrepancy rates indicate samples for which laboratories more frequently diverged from the predefined QC status, whereas the Not provided segment captures missing QC assessments and is not interpreted as analytical disagreement. Detailed sample-level percentages and counts are provided in Appendix Table 7.

### 6.2. SARS2 (SARS-CoV-2, Oxford Nanopore Technologies)

#### 6.2.1. Participation and Submissions

A total of 10 laboratories submitted results for the SARS2 component:

- A total of 37 consensus genome sequences (.fasta) were submitted.
- A total of 35 variant call files (.vcf) were submitted.
- The metadata template completeness for SARS2 submissions had a median of 56.5%.

Most frequent incompleteness drivers in SARS2:
<ul class="compact-list">

<li>Pre-processing fields (missing in 10 laboratories)</li>

<li>Lineage assignment fields (missing in 8 laboratories)</li>

<li>SARS-CoV-2 QC metrics (missing in 8 laboratories)</li>

<li>De-hosting fields (missing in 7 laboratories)</li>

<li>Consensus analysis fields (missing in 7 laboratories)</li>

</ul>

#### 6.2.2. Consensus Genome Reconstruction Performance

Consensus sequences were evaluated against the corresponding curated gold standard for each sample in the SARS2 component.

Overall, SARS2 showed a median genome identity of 99.79%, with a median of 11 nucleotide discrepancies per sample (range: 0–52) (Figure 14).

<figure>
<img src="figures/SARS2/consensus_discrepancies_boxplot_by_sample.png" alt="Consensus discrepancies per sample for SARS2 relative to the curated gold standard." style="width: 90%; max-width: 90%;"/>
</figure>

**Figure 14. Consensus reconstruction performance by sample for SARS2.** Panel A shows the distribution of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample, and Panel B shows the corresponding distribution of genome identity values. In both panels, the central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. In Panel B, the y-axis is truncated to highlight differences among high-identity values.

Figure 15 presents the distribution of nucleotide discrepancy types per sample across participating laboratories for SARS2.

<figure>
<img src="figures/SARS2/consensus_discrepancies_stacked_by_sample.png" alt="Consensus discrepancy types per sample for SARS2 relative to the curated gold standard." style="width: 80%; max-width: 80%;"/>
</figure>

**Figure 15. Consensus discrepancy type composition per sample for SARS2.** Stacked bars represent the number and type of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample.

The dominant discrepancy pattern observed in SARS2 was Stretch of Ns instead of nucleotide stretch (Figure 15). Sample-level consensus reconstruction summary metrics are provided in Appendix Table 8. A full sample-level breakdown of discrepancy categories is provided in Appendix Table 9, while the aggregated discrepancy composition by type and the corresponding category-wise boxplot can be found in Appendix Table 10 and Appendix Figure 3, respectively.

#### 6.2.3. Variant Detection Accuracy

Variant call files (.vcf) submitted for the SARS2 component were compared against the curated reference variant set corresponding to each sample in the SARS2 component.

Overall, SARS2 showed a median of 3 variant discrepancies per sample (range: 0–184), together with a median of 92 successful hits per sample (Table 6, Figure 16).

<figure>
<img src="figures/SARS2/variant_discrepancies_stacked_by_sample.png" alt="Variant discrepancies per sample for SARS2 relative to the curated gold standard." style="width: 80%; max-width: 80%;"/>
</figure>

**Figure 16. Distribution of variant discrepancies per sample for SARS2.** Stacked bars represent the number of nucleotide discrepancies and discrepancy types relative to the curated gold standard across participating laboratories for each sample.

**Table 6. Network-level SARS-CoV-2 variant reporting metrics per sample for SARS2.**

| Sample ID | Expected hits | Median successful hits (n=9) | Median variants >=75% AF in metadata (n=3) | Median variants >=75% AF in VCF (n=8) | Median variants with effect in metadata (n=3) | Median variants with effect in VCF (n=8) | Median discrepancies metadata vs VCF | Median effect discrepancies metadata vs VCF |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SARS6 | 118 | 117 | 117 | 107.5 | 78 | 70 | 0 | 1 |
| SARS7 | 97 | 95 | 94 | 86 | 67 | 59.5 | 0 | 1 |
| SARS8 | 94 | 93 | 93 | 86.5 | 66 | 60.5 | 0 | 1 |
| SARS9 | 16 | 15.5 | 16 | 16 | 8 | 10 | 0 | 1 |
| SARS10 | 11 | 11 | 11 | 10.5 | 5 | 5 | 0 | 0 |

At component level, 3 laboratories reported the number of variants in the metadata, whereas 7 did not report this field for any sample in SARS2. The median number of variants with an allele frequency (AF) >=75% was 73 in the metadata and 73 in the submitted VCF files (Table 6, Figure 17).

Figure 17 summarises the distribution of declared variant reporting modes across submitted sample outputs in SARS2.

<figure>
<img src="figures/SARS2/variant_reporting_practice_by_component.png" alt="Variant reporting practices for SARS2." style="width: 70%; max-width: 70%;"/>
</figure>

**Figure 17. Variant reporting practices for SARS2.** Bars represent the proportion of submitted sample outputs classified as high and low frequency reporting, high frequency only, or low frequency only, according to the metadata declarations associated with the variant outputs for this component.

The dominant discrepancy pattern observed in SARS2 was Missing variant.  The full sample-level variant calling profile is provided in Appendix Table 11, while the aggregated discrepancy composition by type and the corresponding category-wise boxplot can be found in Appendix Table 12 and Appendix Figure 4, respectively.

#### 6.2.4. Lineage, Subtype and Clade Assignment

Lineage, subtype and clade assignments submitted for the SARS2 component were evaluated for concordance with the curated gold standard classifications.

Across all participating laboratories, lineage/subtype concordance reached 60.0%, whereas clade concordance reached 90.0%. The sample-level outcome distribution also shows that part of the observed discordance was associated with missing classifications or inconsistent completion of classification fields rather than with uniform analytical failure across all submissions.

<figure>
<img src="figures/SARS2/typing_outcome_stackedbar_by_sample.png" alt="Classification outcome distribution per sample for SARS2." style="width: 98%; max-width: 98%;"/>
</figure>

**Figure 18. Classification outcome distribution per sample for SARS2.** Panel A shows the proportion of lineage/subtype assignment Match, Discrepancy, and Not provided outcomes across participating laboratories for each sample. Panel B shows the corresponding proportions for clade assignments. Percentages are calculated over all participating laboratories in the component, so the Not provided segment captures samples for which lineage/subtype or clade information was not reported. Detailed sample-level concordance percentages are provided in Appendix Table 13.

#### 6.2.5. Sample Quality Control Assessment

Sample-level QC in SARS2 was evaluated as concordance between the laboratory-reported Pass/Fail classification and the predefined gold standard status. QC concordance was heterogeneous across samples, and some laboratories did not report a formal QC assessment. Network-wide concordance for reported QC decisions was 68.0%.

<figure>
<img src="figures/SARS2/qc_match_by_sample.png" alt="Sample-level QC concordance for SARS2 (Match, Discrepancy, and Not provided relative to the gold standard)." style="width: 90%; max-width: 90%;"/>
</figure>

**_Figure 19_. Sample-level QC concordance for SARS2 relative to the gold standard.** Bars represent the proportion of Match, Discrepancy, and Not provided outcomes per sample across participating laboratories. Higher discrepancy rates indicate samples for which laboratories more frequently diverged from the predefined QC status, whereas the Not provided segment captures missing QC assessments and is not interpreted as analytical disagreement. Detailed sample-level percentages and counts are provided in Appendix Table 14.

### 6.3. FLU1 (Influenza virus, Illumina)

#### 6.3.1. Participation and Submissions

A total of 12 laboratories submitted results for the FLU1 component:

- A total of 47 consensus genome sequences (.fasta) were submitted.
- A total of 42 variant call files (.vcf) were submitted.
- The metadata template completeness for FLU1 submissions had a median of 53.7%.

Most frequent incompleteness drivers in FLU1:
<ul class="compact-list">

<li>Mapping fields (missing in 12 laboratories)</li>

<li>QC metrics fields (missing in 11 laboratories)</li>

<li>Variant calling fields (missing in 11 laboratories)</li>

<li>Pre-processing fields (missing in 10 laboratories)</li>

<li>Consensus analysis fields (missing in 9 laboratories)</li>

</ul>

#### 6.3.2. Consensus Genome Reconstruction Performance

Consensus sequences were evaluated against the corresponding curated gold standard for each sample in the FLU1 component.

Overall, FLU1 showed a median genome identity of 95.96%, with a median of 27 nucleotide discrepancies per sample (range: 8–205) (Figure 20).

<figure>
<img src="figures/FLU1/consensus_discrepancies_boxplot_by_sample.png" alt="Consensus discrepancies per sample for FLU1 relative to the curated gold standard." style="width: 90%; max-width: 90%;"/>
</figure>

**Figure 20. Consensus reconstruction performance by sample for FLU1.** Panel A shows the distribution of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample, and Panel B shows the corresponding distribution of genome identity values. In both panels, the central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. In Panel B, the y-axis is truncated to highlight differences among high-identity values.

Figure 21 presents the distribution of nucleotide discrepancy types per sample across participating laboratories for FLU1.

<figure>
<img src="figures/FLU1/consensus_discrepancies_stacked_by_sample.png" alt="Consensus discrepancy types per sample for FLU1 relative to the curated gold standard." style="width: 80%; max-width: 80%;"/>
</figure>

**Figure 21. Consensus discrepancy type composition per sample for FLU1.** Stacked bars represent the number and type of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample.

The dominant discrepancy pattern observed in FLU1 was Deletion relative to gold standard (Figure 21). Sample-level consensus reconstruction summary metrics are provided in Appendix Table 15. A full sample-level breakdown of discrepancy categories is provided in Appendix Table 16, while the aggregated discrepancy composition by type and the corresponding category-wise boxplot can be found in Appendix Table 17 and Appendix Figure 5, respectively.

#### 6.3.3. Variant Detection Accuracy

For the FLU1 component, variant evaluation focused on the agreement between variants with allele frequency above 75% reported in the metadata template and those represented in the submitted VCF files, together with the overall number of variants present in the VCF output. At component level, 2 laboratories reported the number of variants in the metadata, whereas 10 did not report this field for any sample in FLU1.

Figure 22 summarises the distribution of declared variant reporting modes across submitted sample outputs in FLU1.

<figure>
<img src="figures/FLU1/variant_reporting_practice_by_component.png" alt="Variant reporting practices for FLU1." style="width: 80%; max-width: 80%;"/>
</figure>

**Figure 22. Variant reporting practices for FLU1.** Bars represent the proportion of submitted sample outputs classified as high and low frequency reporting, high frequency only, or low frequency only, according to the metadata declarations associated with the variant outputs for this component.

Overall, FLU1 showed a median of 485.5 variants with allele frequency above 75% reported in the metadata template, compared with 174 corresponding variants represented in the consensus-derived VCF. The median number of discrepancies between both representations was 239.5, while the median total number of variants present in the submitted VCF files was 306 (Table 7, Figure 23).

**Table 7. Network-level influenza variant reporting metrics per sample for FLU1.**

| Sample ID | Median variants >=75% AF in metadata (n=2) | Median variants >=75% AF in VCF (n=11) | Median discrepancies between metadata and VCF | Median total variants in VCF (n=11) |
|---|---:|---:|---:|---:|
| FLU1 | 1587 | 383 | 900.5 | 383 |
| FLU2 | 972 | 168 | 749 | 522 |
| FLU3 | N/A | 203 | N/A | 207.5 |
| FLU4 | 480 | 138.5 | 239.5 | 306.5 |
| FLU5 | 464.5 | 174 | 239.5 | 174 |

These patterns indicate that influenza discrepancies reflect not only analytical differences in variant detection, but also differences in reporting conventions, allele-frequency thresholds, and reference selection. The aggregated structural summary for FLU1 is provided in Appendix Table 18.

<figure>
<img src="figures/FLU1/influenza_variant_reporting_summary_by_sample.png" alt="Influenza variant reporting summary by sample for FLU1." style="width: 90%; max-width: 90%;"/>
</figure>

**Figure 23. Influenza variant reporting summary by sample for FLU1.** Panel A shows, for each sample, the distribution across participating laboratories of the number of variants with allele frequency above 75% reported in the metadata template, the corresponding number represented in the consensus-derived VCF, and the discrepancies between both representations. Panel B shows the distribution across participating laboratories of the total number of variants present in the submitted VCF files for each sample. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range within the plotted scale, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

#### 6.3.4. Lineage, Subtype and Clade Assignment

Lineage, subtype and clade assignments submitted for the FLU1 component were evaluated for concordance with the curated gold standard classifications.

Across all participating laboratories, lineage/subtype concordance reached 95.8%, whereas clade concordance reached 68.8%. The sample-level outcome distribution also shows that part of the observed discordance was associated with missing classifications or inconsistent completion of classification fields rather than with uniform analytical failure across all submissions.

<figure>
<img src="figures/FLU1/typing_outcome_stackedbar_by_sample.png" alt="Classification outcome distribution per sample for FLU1." style="width: 98%; max-width: 98%;"/>
</figure>

**Figure 24. Classification outcome distribution per sample for FLU1.** Panel A shows the proportion of lineage/subtype assignment Match, Discrepancy, and Not provided outcomes across participating laboratories for each sample. Panel B shows the corresponding proportions for clade assignments. Percentages are calculated over all participating laboratories in the component, so the Not provided segment captures samples for which lineage/subtype or clade information was not reported. Detailed sample-level concordance percentages are provided in Appendix Table 19.

#### 6.3.5. Sample Quality Control Assessment

Sample-level QC in FLU1 was evaluated as concordance between the laboratory-reported Pass/Fail classification and the predefined gold standard status. QC concordance was heterogeneous across samples, and some laboratories did not report a formal QC assessment. Network-wide concordance for reported QC decisions was 80.0%.

<figure>
<img src="figures/FLU1/qc_match_by_sample.png" alt="Sample-level QC concordance for FLU1 (Match, Discrepancy, and Not provided relative to the gold standard)." style="width: 90%; max-width: 90%;"/>
</figure>

**_Figure 25_. Sample-level QC concordance for FLU1 relative to the gold standard.** Bars represent the proportion of Match, Discrepancy, and Not provided outcomes per sample across participating laboratories. Higher discrepancy rates indicate samples for which laboratories more frequently diverged from the predefined QC status, whereas the Not provided segment captures missing QC assessments and is not interpreted as analytical disagreement. Detailed sample-level percentages and counts are provided in Appendix Table 20.

### 6.4. FLU2 (Influenza virus, Oxford Nanopore Technologies)

#### 6.4.1. Participation and Submissions

A total of 10 laboratories submitted results for the FLU2 component:

- A total of 49 consensus genome sequences (.fasta) were submitted.
- A total of 39 variant call files (.vcf) were submitted.
- The metadata template completeness for FLU2 submissions had a median of 55.5%.

Most frequent incompleteness drivers in FLU2:
<ul class="compact-list">

<li>Mapping fields (missing in 10 laboratories)</li>

<li>QC metrics fields (missing in 9 laboratories)</li>

<li>Variant calling fields (missing in 9 laboratories)</li>

<li>Pre-processing fields (missing in 9 laboratories)</li>

<li>Clade assignment fields (missing in 8 laboratories)</li>

</ul>

#### 6.4.2. Consensus Genome Reconstruction Performance

Consensus sequences were evaluated against the corresponding curated gold standard for each sample in the FLU2 component.

Overall, FLU2 showed a median genome identity of 95.65%, with a median of 33 nucleotide discrepancies per sample (range: 11–2727) (Figure 26).

<figure>
<img src="figures/FLU2/consensus_discrepancies_boxplot_by_sample.png" alt="Consensus discrepancies per sample for FLU2 relative to the curated gold standard." style="width: 90%; max-width: 90%;"/>
</figure>

**Figure 26. Consensus reconstruction performance by sample for FLU2.** Panel A shows the distribution of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample, and Panel B shows the corresponding distribution of genome identity values. In both panels, the central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. In Panel B, the y-axis is truncated to highlight differences among high-identity values.

Figure 27 presents the distribution of nucleotide discrepancy types per sample across participating laboratories for FLU2.

<figure>
<img src="figures/FLU2/consensus_discrepancies_stacked_by_sample.png" alt="Consensus discrepancy types per sample for FLU2 relative to the curated gold standard." style="width: 80%; max-width: 80%;"/>
</figure>

**Figure 27. Consensus discrepancy type composition per sample for FLU2.** Stacked bars represent the number and type of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample.

The dominant discrepancy pattern observed in FLU2 was Deletion relative to gold standard (Figure 27). Sample-level consensus reconstruction summary metrics are provided in Appendix Table 21. A full sample-level breakdown of discrepancy categories is provided in Appendix Table 22, while the aggregated discrepancy composition by type and the corresponding category-wise boxplot can be found in Appendix Table 23 and Appendix Figure 6, respectively.

#### 6.4.3. Variant Detection Accuracy

For the FLU2 component, variant evaluation focused on the agreement between variants with allele frequency above 75% reported in the metadata template and those represented in the submitted VCF files, together with the overall number of variants present in the VCF output. At component level, 3 laboratories reported the number of variants in the metadata, whereas 7 did not report this field for any sample in FLU2.

Figure 28 summarises the distribution of declared variant reporting modes across submitted sample outputs in FLU2.

<figure>
<img src="figures/FLU2/variant_reporting_practice_by_component.png" alt="Variant reporting practices for FLU2." style="width: 80%; max-width: 80%;"/>
</figure>

**Figure 28. Variant reporting practices for FLU2.** Bars represent the proportion of submitted sample outputs classified as high and low frequency reporting, high frequency only, or low frequency only, according to the metadata declarations associated with the variant outputs for this component.

Overall, FLU2 showed a median of 468.5 variants with allele frequency above 75% reported in the metadata template, compared with 271.5 corresponding variants represented in the consensus-derived VCF. The median number of discrepancies between both representations was 377, while the median total number of variants present in the submitted VCF files was 1085.5 (Table 8, Figure 29).

**Table 8. Network-level influenza variant reporting metrics per sample for FLU2.**

| Sample ID | Median variants >=75% AF in metadata (n=3) | Median variants >=75% AF in VCF (n=6) | Median discrepancies between metadata and VCF | Median total variants in VCF (n=8) |
|---|---:|---:|---:|---:|
| FLU6 | 1542.5 | 408.5 | 1651 | 2833 |
| FLU7 | 430 | 303.5 | 217.5 | 916.5 |
| FLU8 | 475 | 355 | 241.5 | 873.5 |
| FLU9 | 398 | 92.5 | 198.5 | 751.5 |
| FLU10 | 1605 | 505 | 901 | 1321 |

These patterns indicate that influenza discrepancies reflect not only analytical differences in variant detection, but also differences in reporting conventions, allele-frequency thresholds, and reference selection. The aggregated structural summary for FLU2 is provided in Appendix Table 24.

<figure>
<img src="figures/FLU2/influenza_variant_reporting_summary_by_sample.png" alt="Influenza variant reporting summary by sample for FLU2." style="width: 90%; max-width: 90%;"/>
</figure>

**Figure 29. Influenza variant reporting summary by sample for FLU2.** Panel A shows, for each sample, the distribution across participating laboratories of the number of variants with allele frequency above 75% reported in the metadata template, the corresponding number represented in the consensus-derived VCF, and the discrepancies between both representations. Panel B shows the distribution across participating laboratories of the total number of variants present in the submitted VCF files for each sample. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range within the plotted scale, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

#### 6.4.4. Lineage, Subtype and Clade Assignment

Lineage, subtype and clade assignments submitted for the FLU2 component were evaluated for concordance with the curated gold standard classifications.

Across all participating laboratories, lineage/subtype concordance reached 80.0%, whereas clade concordance reached 66.0%. The sample-level outcome distribution also shows that part of the observed discordance was associated with missing classifications or inconsistent completion of classification fields rather than with uniform analytical failure across all submissions.

<figure>
<img src="figures/FLU2/typing_outcome_stackedbar_by_sample.png" alt="Classification outcome distribution per sample for FLU2." style="width: 98%; max-width: 98%;"/>
</figure>

**Figure 30. Classification outcome distribution per sample for FLU2.** Panel A shows the proportion of lineage/subtype assignment Match, Discrepancy, and Not provided outcomes across participating laboratories for each sample. Panel B shows the corresponding proportions for clade assignments. Percentages are calculated over all participating laboratories in the component, so the Not provided segment captures samples for which lineage/subtype or clade information was not reported. Detailed sample-level concordance percentages are provided in Appendix Table 25.

#### 6.4.5. Sample Quality Control Assessment

Sample-level QC in FLU2 was evaluated as concordance between the laboratory-reported Pass/Fail classification and the predefined gold standard status. QC concordance was heterogeneous across samples, and some laboratories did not report a formal QC assessment. Network-wide concordance for reported QC decisions was 100.0%.

<figure>
<img src="figures/FLU2/qc_match_by_sample.png" alt="Sample-level QC concordance for FLU2 (Match, Discrepancy, and Not provided relative to the gold standard)." style="width: 90%; max-width: 90%;"/>
</figure>

**_Figure 31_. Sample-level QC concordance for FLU2 relative to the gold standard.** Bars represent the proportion of Match, Discrepancy, and Not provided outcomes per sample across participating laboratories. Higher discrepancy rates indicate samples for which laboratories more frequently diverged from the predefined QC status, whereas the Not provided segment captures missing QC assessments and is not interpreted as analytical disagreement. Detailed sample-level percentages and counts are provided in Appendix Table 26.

## 7. Discussion

The 2026 RELECOV Dry-Lab Interlaboratory Comparison Exercise provides the first network-wide dry-lab assessment focused specifically on bioinformatic performance across consensus reconstruction, variant reporting, classification, metadata reporting, and QC interpretation. By combining ECDC datasets with in-silico influenza material, the exercise captures both routine-use analytical behaviour and performance under heterogeneous reference and reporting conditions. A separate benchmarking deliverable complements this report with a workflow-level comparison of pipelines and software configurations.

The interpretation of network-wide patterns is also limited by participation. Nineteen of the 52 invited laboratories submitted results, and participation varied between components. The benchmarking therefore describes the analytical practices and performance patterns observed among participating laboratories rather than necessarily representing the full RELECOV network.

### 7.1. Consensus Genome Reconstruction

Consensus reconstruction results were strongest in the Illumina-based components overall, with a combined median genome identity of 97.72% compared with 96.17% across the Nanopore-based components. When grouped by virus, the combined median genome identity was also higher for SARS-CoV-2 components (99.59%) than for influenza components (95.84%). The fact that the influenza median remained slightly below the combined Nanopore-based median suggests that lower identity in influenza cannot be attributed to sequencing platform alone and was also influenced by virus-specific analytical complexity. At component level, SARS1 and SARS2 showed median identities of 99.59% and 99.79%, whereas FLU1 and FLU2 showed lower medians of 95.96% and 95.65%. These results suggest that influenza consensus reconstruction remains less standardised across the network than SARS-CoV-2 reconstruction and would benefit from clearer best-practice guidance on coverage thresholds, masking behaviour, indel handling, and segment-level quality criteria.

At the same time, the ranges observed across laboratories show that high medians did not eliminate outlier behaviour. In particular, the minimum identity values in SARS2 and FLU2 dropped to 5.72% and 42.75%, indicating that a subset of submissions diverged from the curated gold standard.

The dominant discrepancy categories also differed by component. SARS1 component was dominated by defined nucleotides in submitted consensuses where stretches of Ns were present in the gold standard, consistent with more permissive reconstruction in regions that had been masked in the curated reference set. In SARS2 component, by contrast, the most frequent discrepancy category was the presence of a nucleotide when gold standard has an ambiguity, followed by the presence of stretch of Ns when the gold standard has a nucleotide, indicating that disagreement was driven less by a single masking pattern and more by a combination of ambiguity handling and local masking decisions. This was especially clear in samples such as SARS4 in Illumina component and SARS8 in Nanopore component, where low depth or mixed-site complexity increased sensitivity to local reconstruction rules and likely exposed differences between ECDC reconstruction criteria and the thresholds applied within the network.

In influenza, the corrected discrepancy profiles were not uniformly dominated by deletions. FLU1 component showed a mixed pattern with substantial contributions of the presence of ambiguities when the gold standard has a nucleotide and deletions, whereas FLU2 component was dominated primarily by wrong nucleotides followed by the presence of ambiguities when the gold standard has a nucleotide. Although primer-masked regions in the influenza gold standards were not counted as consensus errors, some consensus-generation software still appeared not to reconstruct terminal non-coding regions outside the CDS, which continued to contribute to deletion calls at segment ends. Taken together, these results suggest that influenza inter-laboratory divergence was associated with a combination of ambiguity handling, terminal trimming behaviour, and segment-specific reconstruction rules rather than by end deletions alone.

### 7.2. Variant Detection and Reporting

Variant detection and reporting showed that harmonisation has improved, but it is still incomplete. In SARS-CoV-2, performance was broadly comparable across the two components, suggesting that residual differences were driven less by sequencing platform alone than by sample composition, filtering rules, reporting choices, and software configuration. The remaining discrepancies point above all to the need for stronger standardisation of variant-reporting practice across the network, particularly with respect to allele-frequency thresholds, representation of low-frequency calls, and the use of a common reference framework for VCF-based comparisons.

The SARS-CoV-2 results also show that even when laboratories analyse the same raw data, downstream comparability can still be reduced by heterogeneous output formats and reference choices. In practical terms, comparison becomes more difficult when submissions are generated against different reference backbones, or when some laboratories provide non-VCF tabular outputs that cannot be processed with the same downstream annotation and effect-classification workflow. This suggests that further harmonisation in RELECOV should not focus only on variant calling itself, but also on the expected structure and reference basis of submitted variant files.

The influenza components remain more challenging. Unlike SARS-CoV-2, influenza variant reporting is not yet anchored to a single widely shared comparison framework across the network, and this was reflected in broader heterogeneity in metadata-VCF concordance, total VCF content, and reporting modes. Part of that variability likely reflects analytical differences, but part also reflects differences in how laboratories interpret software outputs and translate them into reportable variant summaries. This is therefore not only a technical issue of sensitivity or specificity, but also a network-level issue of definition: which influenza variants should be reported, relative to which reference framework, in which file format, and under which allele-frequency thresholds.

Taken together, these findings suggest that the next step for RELECOV is not simply to compare laboratory performance against current expectations, but to agree more explicitly on the reporting model itself. For SARS-CoV-2, this means pushing further towards common reference usage and more standardised VCF-like outputs. For influenza, it means opening a broader discussion within the network on how best to define a shared reference and reporting framework that allows variant comparisons to become meaningfully interoperable.

### 7.3. Classification and QC Interpretation

Classification performance was acceptable overall but stronger for lineage/type assignment than for clade assignment. SARS-CoV-2 lineage concordance reached 77.7%, compared with 75.5% for SARS-CoV-2 clade assignment. Influenza type/subtype concordance reached 87.8%, compared with 67.3% for influenza clade assignment.

These results suggest that classification deserves to be considered separately from other downstream outputs. In SARS-CoV-2, part of the excess clade discordance was metadata-driven: among the clade assignments reviewed in the submitted metadata files, some were left empty and others contained values that matched the lineage assignment or had lineage-like syntax rather than a clade designation. At sample level, discordance was concentrated in analytically difficult materials rather than being evenly distributed. In SARS1, the samples showing discordance against the ECDC gold standard were those flagged as low quality, such as SARS4 and SARS5, where low coverage or excess mixed sites likely shifted downstream classification. In SARS2, SARS8 emerged as the most problematic lineage-assignment sample, with several laboratories reporting the parent lineage rather than the exact expected designation, again suggesting that small differences in masking and low-coverage treatment can propagate into classification discordance.

Classification results also need to be interpreted in the context of the databases used by the assigned tools. Concordance is influenced not only by software choice, but also by whether the underlying lineage, clade, type, or subtype database is current and whether the database version can be clearly identified from the submitted metadata. This is particularly relevant for tools such as Nextclade and Pangolin, where software version and database version are distinct and both can affect the final label. Within the network, incomplete or ambiguous reporting of classification database versions makes it harder to determine whether discordance reflects the underlying sequence, the interpretation rules applied by the software, or the use of different database snapshots.

The influenza components showed a related but distinct pattern. A relatively high proportion of laboratories did not report clade values at all, although those that did generally reported concordant assignments; the residual heterogeneity mainly reflected differences between legacy clade labels and more specific subclade labels rather than outright analytical failure. Subtype discordance was more informative, especially in FLU1, where one A/H5N1 sample was reported as A/H5N6 by a laboratory using a BLAST-based assignment strategy, illustrating that subtype performance can also depend on software choice, database interpretation, and the degree to which reference databases are kept up to date.

QC interpretation showed additional between-component differences. Only 9 of the 19 participating laboratories reported at least one sample-level QC assessment in their submitted metadata. Among evaluable QC decisions, network-wide concordance was 71.1%, and component-level concordance ranged from 62.5% in SARS1 to 100.0% in FLU2. QC reporting and concordance varied across components, although interpretation of these results is limited by the small number of laboratories providing sample-level QC assessments.

Because QC assignment depends on how each laboratory interprets coverage, ambiguity, contamination, and other internal acceptance criteria, these results are better understood as reflecting heterogeneity in QC interpretation and reporting rather than a simple measure of analytical correctness alone. In influenza, QC reporting was especially sparse, particularly in FLU2, but among the small number of reported QC decisions the agreement with the expected status was generally high. This suggests that the main challenge in QC was not only interpretation itself, but also the incomplete adoption of formal QC reporting across the network.

### 7.4. Workflow Diversity and Reporting Constraints

The metadata confirms that RELECOV laboratories currently use a diverse analytical landscape. A total of 9 distinct workflows were identified across participating laboratories, together with distinct tools or tool/version combinations for consensus genome generation (16), variant calling (20), SARS-CoV-2 lineage assignment (4), influeza type/subtype assignment (8) and clade assignment (10).

The main incompleteness drivers were variant calling, pre-processing, and mapping fields, followed by QC metrics, de-hosting, consensus analysis, and classification-related metadata. This pattern suggests that laboratories were more consistent in declaring core tool identities than in documenting the exact thresholds and parameter sets that determine analytical behaviour.

### 7.5. Metadata Reporting and Schema Compliance

Metadata quality affected not only interpretability but also interoperability. The exercise showed that laboratories were generally able to declare the core structure of their workflows, yet the level of detail required for reproducibility was often incomplete. This was especially evident for software-version fields, minimum coverage thresholds, variant calling parameters, and reference genome identifiers, all of which are necessary to reconstruct how a consensus or variant set was generated and to determine whether observed differences reflect analytical choice or true performance variation.

The validation process also showed that schema compliance remains an operational issue in its own right. As reported in the metadata results section, only 26.32% of submissions (5 laboratories) were fully compliant with controlled-vocabulary requirements, whereas 73.68% (14 laboratories) required at least one manual correction because non-standard values had been used in controlled fields. The most frequent problems involved free-text software names in dropdown-based fields, inconsistent completion of lineage, clade, type, or subtype assignments, and missing mandatory entries. These issues do not necessarily imply poor analytical practice, but they do reduce the value of the metadata for automated validation, cross-laboratory comparison, and downstream integration into the RELECOV platform.

A recurrent source of non-compliance was the use of free-text entries in fields for which predefined dropdown options were available, particularly for software names. This occurred even though the metadata template, including its controlled-vocabulary dropdowns, had been distributed two weeks before the start of the exercise to give laboratories time to review the available options and identify any missing software tools for possible inclusion in the schema, together with guidance on how mandatory fields should be completed when data were not available. In practice, this means that part of the harmonisation problem lies not only in analytical diversity itself, but in the difficulty of consistently mapping that diversity into a controlled metadata structure.

Another recurrent issue concerned the semantics of version reporting fields. In some submissions, the value entered in the database-version field appears to correspond to the software version rather than the actual database version, particularly for tools such as Nextclade and Pangolin where both identifiers are distinct. This suggests that these metadata fields were not always interpreted consistently and that future versions of the template should include stricter validation rules and clearer examples distinguishing software version from database release. Some laboratories also reported relatively old software versions, including older releases of iVar. The use of different software versions adds another source of workflow heterogeneity and should be considered when interpreting differences between submissions.

### 7.6. Implications for RELECOV 2.0

Taken together, the results support a harmonisation strategy centred on minimum performance and reporting standards rather than on enforcement of a single analytical pipeline. The data do not support a universal workflow ranking that would apply equally across all viruses, platforms, and tasks. Instead, the observed performance patterns were associated with differences in dataset characteristics, reporting conventions, software choice, and parameterisation.

## 8. Conclusions

The 2026 RELECOV Dry-Lab Interlaboratory Comparison Exercise shows that participating laboratories have established bioinformatic capacity for respiratory virus genomic surveillance, while also highlighting differences in analytical performance and reporting practices across the network.

Consensus genome reconstruction showed higher median identity in the Illumina-based components than in the Nanopore-based components. However, because sequencing platform and virus are not independently represented across the four components, these results cannot be used to attribute the observed differences to sequencing technology alone. Within the individual components, differences in masking, coverage thresholds, reference selection, and consensus-generation settings were associated with variation in reconstruction performance.

Variant analysis showed that direct comparison against curated reference sets is feasible for SARS-CoV-2, whereas influenza variant reporting was considerably more heterogeneous. Differences in reference backbones, allele-frequency thresholds, and reporting conventions limited direct comparison between laboratories. A more clearly defined framework for influenza variant reporting would therefore improve the comparability of future exercises.

Classification results were generally more consistent than some of the sequence-level metrics. However, the interpretation of lineage, type/subtype, and clade concordance was affected in some cases by incomplete metadata and inconsistent reporting, particularly for clade assignment. QC reporting also remained incomplete, with only a subset of laboratories providing explicit sample-level assessments.

The metadata analysis identified substantial variation in the way workflows were documented. Differences in software and database versions, reference information, coverage thresholds, and analytical parameters make it difficult to attribute observed performance differences to individual tools or workflow components. More consistent reporting of these elements would improve both the interpretation and reproducibility of future benchmarking exercises.

Overall, the results support RELECOV 2.0 priorities centred on:

- establishing minimum performance criteria for consensus reconstruction and variant reporting;
- defining clearer rules for masking, coverage thresholds, and allele-frequency reporting;
- improving the consistency of software, database, parameter, and reference genome reporting;
- improving the consistency of classification and QC field completion;

The exercise also supports a harmonisation strategy based on common analytical and reporting requirements rather than the adoption of a single pipeline. This approach would allow laboratories to retain flexibility in their choice of tools while providing a more consistent basis for comparing performance across the network.

Taken together, the findings provide a practical basis for improving the comparability and reproducibility of bioinformatic workflows within RELECOV 2.0 and for guiding the development of the RELECOV analytical platform.

## Appendix

This appendix is reserved for supplementary material that may support interpretation of the report but is not essential to the main narrative. Additional figures, extended tables, sensitivity analyses, or other secondary outputs can be included here when relevant.

### SARS1 (SARS-CoV-2, Illumina)

#### Consensus Genome Reconstruction Supplementary Material

**Appendix Table 1. Network-level consensus reconstruction metrics per sample for SARS1.**

| Sample ID | Median genome identity (%) | Median discrepancies | Discrepancies min-max |
|---|---:|---:|---:|
| SARS1 | N/A | N/A | N/A – N/A |
| SARS2 | 99.82 | 2.5 | 1 – 3 |
| SARS3 | 99.59 | 2.5 | 2 – 10 |
| SARS4 | 95.58 | 75.5 | 61 – 125 |
| SARS5 | 99.80 | 5 | 3 – 6 |

**Appendix Table 2. Network-level consensus discrepancy types per sample for SARS1.**

| Sample ID | Median of Wrong nucleotide | Median Nucleotide instead of ambiguity | Median Ambiguity instead of nucleotide | Median Stretch of Ns instead of nucleotide stretch | Median Nucleotide stretch instead of stretch of Ns | Median Insertion relative to gold standard | Median Deletion relative to gold standard |
|---|---:|---:|---:|---:|---:|---:|---:|
| SARS1 | N/A | N/A | N/A | N/A | N/A | N/A | N/A |
| SARS2 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| SARS3 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| SARS4 | 0 | 0 | 0 | 38 | 29 | 0 | 0 |
| SARS5 | 0 | 2 | 0 | 1 | 0 | 0 | 0 |

**Appendix Table 3. Network-level consensus discrepancy composition by type for SARS1.**

| Discrepancy type | Network median per sample | Min-max occurrencies |
|---|---:|---:|
| Incorrect nucleotide | 0 | 0–8 |
| Nucleotide instead of ambiguity | 0 | 0–2 |
| Ambiguity instead of nucleotide | 0 | 0–1 |
| Nucleotide stretch instead of stretch of Ns | 2 | 1–100 |
| Stretch of Ns instead of nucleotide stretch | 0 | 0–99 |
| Insertion relative to gold standard | 0 | 0–3 |
| Deletion relative to gold standard | 0 | 0–1 |

Figure 1 in the appendix summarises the contribution of each discrepancy category observed in SARS1 relative to the curated gold standard.

<figure>
<img src="figures/SARS1/consensus_discrepancy_type_boxplot.png" alt="Composition of consensus discrepancy types for SARS1 relative to the curated gold standard." style="width: 96%; max-width: 96%;"/>
</figure>

**Appendix Figure 1. Composition of consensus discrepancy types relative to the curated gold standard for SARS1.** Boxplots represent aggregated discrepancies across all submitted consensus sequences, stratified by discrepancy category. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

<h4 class="appendix-landscape-heading">Variant Detection Accuracy Supplementary Material</h4>

**Appendix Table 4. Network-level SARS-CoV-2 variant calling profile per sample for SARS1.** The discrepancy-type columns correspond to the median count per sample across participating laboratories.

| Sample ID | Median successful hits | Median discrepancies | Discrepancies min-max | Median wrong nucleotide | Median insertions | Median deletions | Median missing | Median _de novo_ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SARS1 | 32 | 6 | 0 – 53 | 0 | 0 | 0 | 4 | 0 |
| SARS2 | 64 | 3 | 0 – 23 | 0 | 0 | 0 | 1 | 0 |
| SARS3 | 92 | 4 | 0 – 56 | 0 | 0 | 0 | 3 | 0 |
| SARS4 | 75 | 4 | 0 – 518 | 0 | 0 | 1 | 0 | 0 |
| SARS5 | 64 | 2 | 0 – 67 | 0 | 0 | 0 | 1 | 0 |

**Appendix Table 5. Network-level variant discrepancy composition by type for SARS1.** The discrepancy-type columns correspond to the median count per sample across participating laboratories.

| Discrepancy type | Network median per sample | Network min-max per sample |
|---|---:|---:|
| Incorrect nucleotide | 0 | 0–1 |
| Insertion relative to gold standard | 0 | 0–4 |
| Deletions relative to gold standard | 0 | 0–40 |
| Missing expected variants | 1 | 0–65 |
| De novo variants | 0 | 0–474 |

Figure 2 in the appendix summarises the contribution of each discrepancy category observed in SARS1 relative to the curated gold standard.

<figure>
<img src="figures/SARS1/variant_discrepancy_type_boxplot.png" alt="Composition of variant discrepancy types for SARS1 relative to the curated gold standard." style="width: 90%; max-width: 90%;"/>
</figure>

**Appendix Figure 2. Composition of variant discrepancy types relative to the curated gold standard for SARS1.** Boxplots represent aggregated discrepancies across all submitted variant calls, stratified by discrepancy category (incorrect nucleotide, excess ambiguous bases, and indels). Where required, a broken y-axis is used to preserve visual resolution in the lower discrepancy range while still displaying higher values above an empty interval. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

#### Lineage, Subtype and Clade Assignment Supplementary Material

**Appendix Table 6. Network-level classification outcomes per sample for SARS1.**

| Sample ID | Lineage/Subtype matches (%) | Clade matches (%) |
|---|---:|---:|
| SARS1 | N/A | N/A |
| SARS2 | 87.50 | 68.75 |
| SARS3 | 81.25 | 68.75 |
| SARS4 | 87.50 | 68.75 |
| SARS5 | 87.50 | 68.75 |

#### Sample Quality Control Assessment Supplementary Material

**Appendix Table 7. Sample-level QC concordance for SARS1 for reported QC classification.**

| Sample ID | Gold standard QC | % Match | # Matches | # Discrepancies | Total evaluations |
|---|---:|---:|---:|---:|---:|
| SARS1 | Fail | 100.0% | 8 | 0 | 8 |
| SARS2 | Pass | 100.0% | 8 | 0 | 8 |
| SARS3 | Pass | 100.0% | 8 | 0 | 8 |
| SARS4 | Fail | 12.5% | 1 | 7 | 8 |
| SARS5 | Fail | 0.0% | 0 | 8 | 8 |

### SARS2 (SARS-CoV-2, Oxford Nanopore Technologies)

#### Consensus Genome Reconstruction Supplementary Material

**Appendix Table 8. Network-level consensus reconstruction metrics per sample for SARS2.**

| Sample ID | Median genome identity (%) | Median discrepancies | Discrepancies min-max |
|---|---:|---:|---:|
| SARS6 | 99.91 | 2.5 | 0 – 23 |
| SARS7 | 99.88 | 8 | 1 – 23 |
| SARS8 | 98.76 | 32.5 | 29 – 52 |
| SARS9 | N/A | N/A | N/A – N/A |
| SARS10 | 99.04 | 2 | 1 – 42 |

**Appendix Table 9. Network-level consensus discrepancy types per sample for SARS2.**

| Sample ID | Median of Wrong nucleotide | Median Nucleotide instead of ambiguity | Median Ambiguity instead of nucleotide | Median Stretch of Ns instead of nucleotide stretch | Median Nucleotide stretch instead of stretch of Ns | Median Insertion relative to gold standard | Median Deletion relative to gold standard |
|---|---:|---:|---:|---:|---:|---:|---:|
| SARS6 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| SARS7 | 1 | 1 | 0 | 0 | 1 | 1 | 0 |
| SARS8 | 0.5 | 24 | 0 | 0 | 7 | 0 | 0 |
| SARS9 | N/A | N/A | N/A | N/A | N/A | N/A | N/A |
| SARS10 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |

**Appendix Table 10. Network-level consensus discrepancy composition by type for SARS2.**

| Discrepancy type | Network median per sample | Min-max occurrencies |
|---|---:|---:|
| Incorrect nucleotide | 0 | 0–23 |
| Nucleotide instead of ambiguity | 0 | 0–28 |
| Ambiguity instead of nucleotide | 0 | 0–0 |
| Nucleotide stretch instead of stretch of Ns | 0 | 0–42 |
| Stretch of Ns instead of nucleotide stretch | 1 | 0–30 |
| Insertion relative to gold standard | 0 | 0–4 |
| Deletion relative to gold standard | 0 | 0–2 |

Figure 3 in the appendix summarises the contribution of each discrepancy category observed in SARS2 relative to the curated gold standard.

<figure>
<img src="figures/SARS2/consensus_discrepancy_type_boxplot.png" alt="Composition of consensus discrepancy types for SARS2 relative to the curated gold standard." style="width: 96%; max-width: 96%;"/>
</figure>

**Appendix Figure 3. Composition of consensus discrepancy types relative to the curated gold standard for SARS2.** Boxplots represent aggregated discrepancies across all submitted consensus sequences, stratified by discrepancy category. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

<h4 class="appendix-landscape-heading">Variant Detection Accuracy Supplementary Material</h4>

**Appendix Table 11. Network-level SARS-CoV-2 variant calling profile per sample for SARS2.** The discrepancy-type columns correspond to the median count per sample across participating laboratories.

| Sample ID | Median successful hits | Median discrepancies | Discrepancies min-max | Median wrong nucleotide | Median insertions | Median deletions | Median missing | Median _de novo_ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SARS6 | 117 | 1 | 0 – 23 | 0 | 0 | 0 | 1 | 0 |
| SARS7 | 95 | 3 | 0 – 23 | 0 | 0 | 1 | 2 | 0 |
| SARS8 | 93 | 4 | 0 – 29 | 0 | 0 | 0 | 1 | 0 |
| SARS9 | 15.5 | 6 | 0 – 157 | 0 | 0 | 0.5 | 0.5 | 4 |
| SARS10 | 11 | 0 | 0 – 184 | 0 | 0 | 0 | 0 | 0 |

**Appendix Table 12. Network-level variant discrepancy composition by type for SARS2.** The discrepancy-type columns correspond to the median count per sample across participating laboratories.

| Discrepancy type | Network median per sample | Network min-max per sample |
|---|---:|---:|
| Incorrect nucleotide | 0 | 0–0 |
| Insertion relative to gold standard | 0 | 0–27 |
| Deletions relative to gold standard | 0 | 0–25 |
| Missing expected variants | 1 | 0–23 |
| De novo variants | 0 | 0–150 |

Figure 4 in the appendix summarises the contribution of each discrepancy category observed in SARS2 relative to the curated gold standard.

<figure>
<img src="figures/SARS2/variant_discrepancy_type_boxplot.png" alt="Composition of variant discrepancy types for SARS2 relative to the curated gold standard." style="width: 90%; max-width: 90%;"/>
</figure>

**Appendix Figure 4. Composition of variant discrepancy types relative to the curated gold standard for SARS2.** Boxplots represent aggregated discrepancies across all submitted variant calls, stratified by discrepancy category (incorrect nucleotide, excess ambiguous bases, and indels). Where required, a broken y-axis is used to preserve visual resolution in the lower discrepancy range while still displaying higher values above an empty interval. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

#### Lineage, Subtype and Clade Assignment Supplementary Material

**Appendix Table 13. Network-level classification outcomes per sample for SARS2.**

| Sample ID | Lineage/Subtype matches (%) | Clade matches (%) |
|---|---:|---:|
| SARS6 | 90.00 | 90.00 |
| SARS7 | 80.00 | 90.00 |
| SARS8 | 10.00 | 90.00 |
| SARS9 | N/A | N/A |
| SARS10 | N/A | N/A |

#### Sample Quality Control Assessment Supplementary Material

**Appendix Table 14. Sample-level QC concordance for SARS2 for reported QC classification.**

| Sample ID | Gold standard QC | % Match | # Matches | # Discrepancies | Total evaluations |
|---|---:|---:|---:|---:|---:|
| SARS6 | Pass | 80.0% | 4 | 1 | 5 |
| SARS7 | Pass | 80.0% | 4 | 1 | 5 |
| SARS8 | Fail | 20.0% | 1 | 4 | 5 |
| SARS9 | Fail | 80.0% | 4 | 1 | 5 |
| SARS10 | Fail | 80.0% | 4 | 1 | 5 |

### FLU1 (Influenza virus, Illumina)

#### Consensus Genome Reconstruction Supplementary Material

**Appendix Table 15. Network-level consensus reconstruction metrics per sample for FLU1.**

| Sample ID | Median genome identity (%) | Median discrepancies | Discrepancies min-max |
|---|---:|---:|---:|
| FLU1 | 96.07 | 27 | 10 – 178 |
| FLU2 | 95.90 | 17.5 | 8 – 205 |
| FLU3 | N/A | N/A | N/A – N/A |
| FLU4 | 96.06 | 25 | 16 – 187 |
| FLU5 | 95.78 | 25 | 19 – 159 |

**Appendix Table 16. Network-level consensus discrepancy types per sample for FLU1.**

| Sample ID | Median of Wrong nucleotide | Median Nucleotide instead of ambiguity | Median Ambiguity instead of nucleotide | Median Stretch of Ns instead of nucleotide stretch | Median Nucleotide stretch instead of stretch of Ns | Median Insertion relative to gold standard | Median Deletion relative to gold standard |
|---|---:|---:|---:|---:|---:|---:|---:|
| FLU1 | 11 | 0 | 1 | 0 | 0 | 0 | 16 |
| FLU2 | 0 | 0 | 0 | 0 | 0 | 0 | 16 |
| FLU3 | N/A | N/A | N/A | N/A | N/A | N/A | N/A |
| FLU4 | 0 | 0 | 0 | 0 | 0 | 0 | 16 |
| FLU5 | 0 | 0 | 0 | 1 | 6.5 | 0 | 14 |

**Appendix Table 17. Network-level consensus discrepancy composition by type for FLU1.**

| Discrepancy type | Network median per sample | Min-max occurrencies |
|---|---:|---:|
| Incorrect nucleotide | 0 | 0–162 |
| Nucleotide instead of ambiguity | 0 | 0–0 |
| Ambiguity instead of nucleotide | 0 | 0–171 |
| Nucleotide stretch instead of stretch of Ns | 0 | 0–11 |
| Stretch of Ns instead of nucleotide stretch | 0 | 0–120 |
| Insertion relative to gold standard | 0 | 0–1 |
| Deletion relative to gold standard | 16 | 0–18 |

Figure 5 in the appendix summarises the contribution of each discrepancy category observed in FLU1 relative to the curated gold standard.

<figure>
<img src="figures/FLU1/consensus_discrepancy_type_boxplot.png" alt="Composition of consensus discrepancy types for FLU1 relative to the curated gold standard." style="width: 96%; max-width: 96%;"/>
</figure>

**Appendix Figure 5. Composition of consensus discrepancy types relative to the curated gold standard for FLU1.** Boxplots represent aggregated discrepancies across all submitted consensus sequences, stratified by discrepancy category. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

<h4 class="appendix-landscape-heading">Variant Detection Accuracy Supplementary Material</h4>

**Appendix Table 18. Aggregated influenza variant reporting metrics for FLU1.**

| Metric | Network median | Network min-max |
|---|---:|---:|
| Variants >=75% AF in metadata | 485.5 | 449–1797 |
| Variants >=75% AF in VCF | 174 | 0–1373 |
| Discrepancies between metadata and VCF | 239.5 | 0–1797 |
| Total variants in VCF (n=11) | 306 | 0–1382 |

#### Lineage, Subtype and Clade Assignment Supplementary Material

**Appendix Table 19. Network-level classification outcomes per sample for FLU1.**

| Sample ID | Lineage/Subtype matches (%) | Clade matches (%) |
|---|---:|---:|
| FLU1 | 91.67 | 58.33 |
| FLU2 | 100.00 | 75.00 |
| FLU3 | N/A | N/A |
| FLU4 | 91.67 | 66.67 |
| FLU5 | 100.00 | 75.00 |

#### Sample Quality Control Assessment Supplementary Material

**Appendix Table 20. Sample-level QC concordance for FLU1 for reported QC classification.**

| Sample ID | Gold standard QC | % Match | # Matches | # Discrepancies | Total evaluations |
|---|---:|---:|---:|---:|---:|
| FLU1 | Pass | 66.7% | 2 | 1 | 3 |
| FLU2 | Pass | 100.0% | 3 | 0 | 3 |
| FLU3 | Fail | 100.0% | 3 | 0 | 3 |
| FLU4 | Pass | 66.7% | 2 | 1 | 3 |
| FLU5 | Fail | 66.7% | 2 | 1 | 3 |

### FLU2 (Influenza virus, Oxford Nanopore Technologies)

#### Consensus Genome Reconstruction Supplementary Material

**Appendix Table 21. Network-level consensus reconstruction metrics per sample for FLU2.**

| Sample ID | Median genome identity (%) | Median discrepancies | Discrepancies min-max |
|---|---:|---:|---:|
| FLU6 | 95.69 | 30 | 18 – 524 |
| FLU7 | 95.52 | 49.5 | 16 – 344 |
| FLU8 | 95.55 | 40 | 18 – 299 |
| FLU9 | 95.47 | 20.5 | 11 – 2727 |
| FLU10 | 96.05 | 29 | 24 – 1828 |

**Appendix Table 22. Network-level consensus discrepancy types per sample for FLU2.**

| Sample ID | Median of Wrong nucleotide | Median Nucleotide instead of ambiguity | Median Ambiguity instead of nucleotide | Median Stretch of Ns instead of nucleotide stretch | Median Nucleotide stretch instead of stretch of Ns | Median Insertion relative to gold standard | Median Deletion relative to gold standard |
|---|---:|---:|---:|---:|---:|---:|---:|
| FLU6 | 2 | 0 | 0 | 0 | 3.5 | 2 | 16 |
| FLU7 | 0 | 0 | 0 | 7 | 0 | 0 | 16 |
| FLU8 | 0 | 0 | 0 | 0 | 1.5 | 0 | 18 |
| FLU9 | 0 | 0 | 0 | 1 | 0 | 0 | 15 |
| FLU10 | 8 | 0 | 0 | 0 | 0 | 2 | 16 |

**Appendix Table 23. Network-level consensus discrepancy composition by type for FLU2.**

| Discrepancy type | Network median per sample | Min-max occurrencies |
|---|---:|---:|
| Incorrect nucleotide | 0 | 0–2623 |
| Nucleotide instead of ambiguity | 0 | 0–1 |
| Ambiguity instead of nucleotide | 0 | 0–424 |
| Nucleotide stretch instead of stretch of Ns | 0 | 0–16 |
| Stretch of Ns instead of nucleotide stretch | 0 | 0–41 |
| Insertion relative to gold standard | 1 | 0–53 |
| Deletion relative to gold standard | 16 | 1–71 |

Figure 6 in the appendix summarises the contribution of each discrepancy category observed in FLU2 relative to the curated gold standard.

<figure>
<img src="figures/FLU2/consensus_discrepancy_type_boxplot.png" alt="Composition of consensus discrepancy types for FLU2 relative to the curated gold standard." style="width: 96%; max-width: 96%;"/>
</figure>

**Appendix Figure 6. Composition of consensus discrepancy types relative to the curated gold standard for FLU2.** Boxplots represent aggregated discrepancies across all submitted consensus sequences, stratified by discrepancy category. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

<h4 class="appendix-landscape-heading">Variant Detection Accuracy Supplementary Material</h4>

**Appendix Table 24. Aggregated influenza variant reporting metrics for FLU2.**

| Metric | Network median | Network min-max |
|---|---:|---:|
| Variants >=75% AF in metadata | 468.5 | 377–1794 |
| Variants >=75% AF in VCF | 271.5 | 0–1363 |
| Discrepancies between metadata and VCF | 377 | 8–1794 |
| Total variants in VCF (n=8) | 1085.5 | 0–7903 |

#### Lineage, Subtype and Clade Assignment Supplementary Material

**Appendix Table 25. Network-level classification outcomes per sample for FLU2.**

| Sample ID | Lineage/Subtype matches (%) | Clade matches (%) |
|---|---:|---:|
| FLU6 | 70.00 | 80.00 |
| FLU7 | 90.00 | 80.00 |
| FLU8 | 90.00 | 80.00 |
| FLU9 | 80.00 | 20.00 |
| FLU10 | 70.00 | 70.00 |

#### Sample Quality Control Assessment Supplementary Material

**Appendix Table 26. Sample-level QC concordance for FLU2 for reported QC classification.**

| Sample ID | Gold standard QC | % Match | # Matches | # Discrepancies | Total evaluations |
|---|---:|---:|---:|---:|---:|
| FLU6 | Pass | 100.0% | 2 | 0 | 2 |
| FLU7 | Pass | 100.0% | 2 | 0 | 2 |
| FLU8 | Pass | 100.0% | 2 | 0 | 2 |
| FLU9 | Fail | 100.0% | 2 | 0 | 2 |
| FLU10 | Pass | 100.0% | 2 | 0 | 2 |

