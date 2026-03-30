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
- [7. Discussion](#7-discussion)
  - [7.1. Overall Analytical Robustness](#71-overall-analytical-robustness)
  - [7.2. Variant Detection Accuracy](#72-variant-detection-accuracy)
  - [7.3. Classification Accuracy and Database Versioning](#73-classification-accuracy-and-database-versioning)
  - [7.4. Workflow Diversity and Standardisation Balance](#74-workflow-diversity-and-standardisation-balance)
  - [7.5. Implications for RELECOV 2.0](#75-implications-for-relecov-20)
- [8. Conclusions](#8-conclusions)


## Executive Summary

The 2026 RELECOV Dry-Lab EQA provides the first network-wide dry-lab assessment focused specifically on bioinformatic analytical performance across respiratory virus surveillance workflows. Nineteen laboratories participated, corresponding to 36.54% of invited laboratories, with high submission rates for expected analytical outputs: 97.52% for consensus genome files and 87.13% for VCF files.

Across the network, consensus genome reconstruction performed best in the Illumina-based components, with a combined median genome identity of 97.75%, compared with 96.10% in the Nanopore-based components. However, broad identity ranges in SARS2 and FLU2 indicate that outlier submissions remained present, particularly in contexts where masking, coverage thresholds, and consensus-generation choices differed across laboratories.

Variant reporting showed clear methodological heterogeneity. For SARS-CoV-2, the median number of discrepancies relative to the curated reference variant sets was 9.0 in the Illumina component and 5.5 in the Nanopore component. Influenza submissions were more heterogeneous structurally, with a median of 477.0 high-frequency variants reported in metadata, compared with 188.0 derived from submitted VCF files, and a median discrepancy of 377.0 between both representations.

Classification performance was consistently higher for lineage/type assignment than for clade assignment. SARS-CoV-2 lineage concordance reached 77.7%, compared with 73.4% for clade assignment, while influenza type/subtype concordance reached 81.6%, compared with 67.3% for clade assignment. Review of submitted files further indicated that part of the excess discordance in clade assignment reflected field completion and nomenclature issues, including missing clade entries and lineage/type-like values entered in the clade field.

Metadata completeness and reporting remain major priorities for harmonisation. The median metadata completeness rate across participating laboratories was 59.1%, with values ranging from 11.3% to 92.0%. Although software names were reported for 71.8% of expected fields, only 60.1% of software-version fields, 45.6% of coverage thresholds, 38.9% of variant-calling parameter fields, and 57.7% of reference genome identifiers were completed. A total of 10 distinct workflows were identified, together with substantial diversity in consensus, variant calling, and classification software.

Overall, the EQA shows that RELECOV laboratories already have substantial bioinformatic capacity, but also that inter-laboratory comparability is limited by heterogeneous thresholds, parameter reporting, reference selection, and uneven completion of metadata and QC fields. These findings support RELECOV 2.0 priorities centred on minimum performance standards, stronger metadata requirements, clearer reporting rules for consensus and variants, and component-aware benchmarking rather than a single cross-context pipeline ranking.

## 1. Introduction

The RELECOV Network aims to strengthen genomic surveillance of respiratory viruses by developing and harmonising analytical capacities across the participating laboratories. In this context, it was essential to **assess the consistency, reproducibility and maturity of the bioinformatics workflows implemented within the network**.

To this end, an **external quality assessment (EQA) exercise in dry lab format** was conducted, inspired by the ECDC’s 2024 dry-lab EQA. The exercise focused on the bioinformatic characterisation of respiratory viruses, covering key analytical tasks including viral genome reconstruction, variant identification, and lineage and clade assignment.

A central component of this initiative was to evaluate the range of analytical pipelines used across the RELECOV Network, identify their relative performance, and determine which approaches were best suited to support a genomic surveillance standard within the network. This evaluation contributes directly to **Objective 3** of RELECOV 2.0, which focuses on *generating a comprehensive understanding of the analytical and operational workflows currently implemented across Spanish laboratories*. Furthermore, the EQA provides the practical evidence base required for Task T6.1, which aims to identify, compare and prioritise the analytical pipelines available nationally in order to define the workflow that should be integrated into the RELECOV analytical platform.

The exercise was also aligned with **Milestone M3.2**, which pertains to the *establishment of harmonised analytical procedures for respiratory viruses*. It also contributes to **Deliverable D4.1**, focused on the definition of minimum metadata requirements and harmonised reporting formats, by *testing the ability of laboratories to produce interoperable outputs suitable for integration into national and international surveillance systems*. In addition, the exercise informs the development of the RELECOV platform by providing operational insights critical to **Task T5.2**, which addresses *data ingestion, workflow automation and technical specifications for integrating pipeline outputs within the platform*.

Beyond workflow harmonisation, the EQA has contributed to capacity-building and performance assessment activities central to the project. By benchmarking analytical performance across laboratories, the exercise provides essential input for **Milestone M6.3**, related to *determining laboratory readiness and identifying areas requiring technical reinforcement or training*. These insights also support **Deliverable D6.2**, which includes *recommendations for strengthening network-wide analytical capacity and ensuring long-term sustainability of genomic surveillance operations*.

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



Table 1 summarises the correspondence between RELECOV EQA samples and their original source datasets, including ECDC ESIB references.

_**Table 1**. Overview of SARS-CoV-2 datasets used in the RELECOV 2026 Dry-Lab EQA.
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



Table 2 describes the design characteristics of in-silico samples, including virus composition and intended benchmarking challenges.

_**Table 2**. Viral, host and contaminant composition design of in-silico influenza datasets used for benchmarking._

| Sample | Influenza reads | Host reads | Additional Viral reads  | Total reads | Analytical Challenge                |
|--------|-----------------|------------|-------------------------|-------------|-------------------------------------|
| FLU2   | 1378764         | 462520     | 0                       | 1841284     | Baseline performance assessment     |
| FLU4   | 181626          | 300000     | 200000 SARS-CoV-2 reads | 681626      | Contamination with SARS-CoV-2       |
| FLU5   | 1088000         | 100000     | 0                       | 1188000     | NA segment dropout                  |
| FLU7   | 5677            | 100        | 255 Rhinovirus reads    | 6032        | Cross-virus contamination challenge |
| FLU8   | 5380            | 300        | 0                       | 5680        | Baseline performance assessment     |
| FLU9   | 19989           | 500        | 0                       | 20489       | HA segment dropout                  |



Table 3 summarises the influenza datasets included in the EQA, detailing enrichment strategy, primer scheme, sequencing technology, and key analytical challenges.

_**Table 3**. Influenza virus samples used in the RELECOV 2026 Dry-Lab EQA, including sequencing platform, enrichment strategy, primer scheme, and key analytical features._

| Sample | Source    | Platform | Enrichment Strategy | Primer Scheme                                   | Read Layout | Ref_sample        | Type   | Clade HA  | Legacy Clade       | Key Feature                             | Quality check |
|--------|-----------|----------|---------------------|-------------------------------------------------|-------------|-------------------|--------|-----------| ------------------ | ----------------------------------------|---------------|
| FLU1   | ESIB 2024 | Illumina | Amplicon            | CommonUni12/13 (Van den Hoecke 2015)            | Paired-end  | INFL2.07          | A/H5N1 | 2.3.4.4b  | -                  | High-quality baseline sample (zoonotic) | Ok            |
| FLU2   | In-silico | Illumina | Amplicon            | Zhou 2009 single-reaction genomic amplification | Paired-end  | In-silico Sample1 | A/H1N1 | D.3.1.1   | 6B.1A.5a.2a.1      | High-quality baseline sample (human)    | Ok            |
| FLU3   | ESIB 2024 | Illumina | No enrichment       | —                                               | Paired-end  | INFL2.04          | —      | —         | -                  | No influenza (Rhinovirus only)          | Bad           |
| FLU4   | In-silico | Illumina | Amplicon            | Zhou 2009 single-reaction genomic amplification | Paired-end  | In-silico Sample3 | A/H3N2 | K         | 3C.2a1b.2a.2a.3a.1 | Contamination with SARS-CoV-2           | Ok            |
| FLU5   | In-silico | Illumina | Amplicon            | Zhou 2009 single-reaction genomic amplification | Paired-end  | In-silico Sample4 | A/H3N2 | J.2.2     | 3C.2a1b.2a.2a.3a.1 | NA segment dropout                      | Bad           |
| FLU6   | ESIB 2024 | Nanopore | No enrichment       | —                                               | Single-end  | INFL1.02          | A/H5N6 | 2.3.4.4h  | -                  | High-quality baseline sample (zoonotic) | Ok            |
| FLU7   | In-silico | Nanopore | Amplicon            | Zhou 2009 single-reaction genomic amplification | Single-end  | In-silico Sample2 | A/H1N1 | C.1.9.3   | 6B.1A.5a.2a        | Contamination with Rhinovirus           | Ok            |
| FLU8   | In-silico | Nanopore | Amplicon            | Zhou 2009 single-reaction genomic amplification | Single-end  | In-silico Sample3 | A/H3N2 | K         | 3C.2a1b.2a.2a.3a.1 | High-quality baseline sample (human)    | Ok            |
| FLU9   | In-silico | Nanopore | Amplicon            | Zhou 2009 single-reaction genomic amplification | Single-end  | In-silico Sample1 | A/H1N1 | D.3.1.1   | 6B.1A.5a.2a.1      | HA segment dropout                      | Bad           |
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



<div class="equation"><span class="equation-lhs">File submission rate</span><span class="equation-equals">=</span><span class="equation-fraction"><span class="equation-numerator">Number of submitted files</span><span class="equation-denominator">Total number of expected files</span></span></div>



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

Unlike SARS-CoV-2, where laboratories predominantly use a shared and globally standardised reference genomes (either MN908947.3 or NC_045512.2), influenza virus analyses exhibited substantial heterogeneity in reference genome selection. Participating laboratories employed distinct segment-specific reference sequences. As a result:

- Variant coordinates were reported relative to different reference accessions.
- Segment boundaries and numbering schemes varied.
- Insertions and deletions were represented inconsistently across reference backbones.

This heterogeneity prevented robust coordinate harmonisation across submissions without introducing alignment-dependent artefacts and interpretation bias.

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

For influenza virus, additional structural summary metrics were calculated because direct coordinate-harmonised comparison of all submitted variants was not methodologically robust across segment-specific reference backbones:

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

Potential contributors discussed during result interpretation included:

- Use of outdated lineage databases
- Differences in software versioning
- Incomplete consensus sequences
- Excess ambiguous positions

Failure to identify virus presence in positive samples, or misclassification of negative samples, was recorded separately.

### 4.5. Evaluation of Metadata Completeness and Compliance

Metadata assessment focused on analytical transparency and interoperability rather than biological correctness.

For each submitted sample, metadata completeness was calculated as:



<div class="equation"><span class="equation-lhs">Sample metadata completeness</span><span class="equation-equals">=</span><span class="equation-fraction"><span class="equation-numerator">Number of correctly populated fields</span><span class="equation-denominator">Total number of applicable fields</span></span></div>



Laboratory-level and component-level completeness summaries were then derived from these sample-level values. Fields were evaluated for:

- Completion: Each sample has a list of minimum **recomended** fields, based on the sample characteristics. For each component/sample/lab the total number of completed minimum **recomended** fields was evaluated. Both mandatory and optional analytical fields were included in the completeness assessment, while fields not applicable to a laboratory’s selected components were excluded from scoring.
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
- Number of Discrepancies
- QC concordance rate, where:



<div class="equation"><span class="equation-lhs">QC concordance rate</span><span class="equation-equals">=</span><span class="equation-fraction"><span class="equation-numerator">Number of Matches</span><span class="equation-denominator">Total QC evaluations</span></span></div>



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

A total of 52 laboratories within the RELECOV network were invited to participate. Of these, 19 laboratories 36.5% submitted results for one or more components with the following distribution:

- SARS1 (SARS-CoV-2, Illumina): 16 laboratories.
- SARS2 (SARS-CoV-2, Oxford Nanopore Technologies): 10 laboratories.
- FLU1 (Influenza virus, Illumina): 12 laboratories.
- FLU2 (Influenza virus, Oxford Nanopore Technologies): 10 laboratories.

The median number of components analysed per participating laboratory was 2.0.

### 5.1. Submission Completeness

Assessment of submission completeness was conducted in accordance with the criteria outlined in [Section 4.1](#41-submission-completeness). Across all components:

- 97.5% of laboratories submitted consensus genome files (.fasta), where applicable.
- 87.1% submitted variant call files (.vcf), where applicable.

Component-level submission totals are presented in Section 6 and reflect both the number of participating laboratories and the expected output files for each dataset.

### 5.2. Consensus Genome Reconstruction Performance

Consensus genome reconstruction performance was measured using the evaluation criteria detailed in [Section 4.2](#42-evaluation-of-consensus-genome-reconstruction-performance). Across the two Illumina-based components, the combined median genome identity was 97.75%, compared with 96.10% across the two Nanopore-based components. Nanopore-based datasets also showed broader overall identity ranges, where low-identity outliers were present.

Dominant discrepancy patterns differed by component:

- In SARS1, the most frequent discrepancy category was stretches of Ns in the submitted consensus where defined nucleotides were present in the gold standard.
- In SARS2, the most frequent discrepancy category was defined nucleotides in the submitted consensus where stretches of Ns were present in the gold standard.
- FLU1 and FLU2 were both dominated by deletions relative to the gold standard.

These SARS-CoV-2 discrepancy patterns are consistent with differences in masking behaviour and/or minimum coverage thresholds relative to the gold standard reconstruction criteria.

Across components, many discrepancy categories had medians of zero, indicating that errors tended to be concentrated in a smaller number of laboratories or samples rather than being uniformly distributed across the network.



Figure 1 summarises consensus genome reconstruction performance across all components.


<figure>
<img src="figures/network/consensus_summary.png" alt="Network-level consensus reconstruction performance summary." style="max-width: 100%;"/>
<figcaption>Network-level consensus reconstruction performance summary.</figcaption>
</figure>

**_Figure 1_. Distribution of consensus genome discrepancies relative to the gold standard across components**. Boxplots represent the number of nucleotide discrepancies per component across participating laboratories. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

### 5.3. Variant Detection Accuracy

Variant detection accuracy was evaluated following the methodological framework described in [Section 4.3](#43-evaluation-of-variant-detectio-accuracy).

#### 5.3.1. SARS-CoV-2

For SARS-CoV-2 compoents (SARS1 and SARS2), variant detection accuracy was assessed against curated reference variant sets. Overall, submitted VCFs showed a median number of discrepancies of 9.0 for Illumina component and a median number of 5.5 for Nanopore component, discrepancies relative to the reference variant set.



For SARS-CoV-2, the median number of variant discrepancies was 9.0 for the Illumina component and 5.5 for the Nanopore component. The distribution of variant detection performance across components is presented in Figure 2. Contextual factors documented in the metadata that may contribute to these differences included:

- Allele frequency thresholds used for incorporation into vcf files
- Reference genome selection
- Variant normalization practices (variant caller software and params)

<figure>
<img src="figures/network/variant_summary.png" alt="Network-level variant detection performance summary." style="max-width: 100%;"/>
<figcaption>Network-level variant detection performance summary.</figcaption>
</figure>


**_Figure 2_. SARS-CoV-2 network-level variant detection performance summary**. Boxplots represent the number of variant discrepancies per SARS-CoV-2 component across participating laboratories. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

Variant evaluation included structural reporting characteristics and methodological heterogeneity. At network level:

- 66.67 of laboratories reported both high- and low-frequency variants.
- 0.0 reported exclusively low-frequency variants.
- 33.33 reported only high-frequency variants.

Additionally, a total of 3 distinct reference genomes were employed for variant calling across SARS-CoV-2 components.



Figure 3 summarizes the distribution of variant reporting practices across participating laboratories for SARS-CoV-2 components.

<figure>
<img src="figures/network/sars_variant_reporting_summary.png" alt="SARS-CoV-2 variant reporting practices across the network." style="max-width: 100%;"/>
<figcaption>SARS-CoV-2 variant reporting practices across the network.</figcaption>
</figure>


**_Figure 3_. SARS-CoV-2 variant reporting characteristics across the network**. Summarise the proportion of laboratories reporting high- and/or low-frequency variants.

#### 5.3.2. Influenza virus

For influenza virus components (FLU1 and FLU2), variant evaluation focused on structural reporting characteristics and methodological heterogeneity.

At network level:

- 50.68 of laboratories reported both high- and low-frequency variants.
- 36.99 reported exclusively low-frequency variants.
- 12.33 reported only high-frequency variants.

Additionally, a total of 7.875 distinct reference genomes were employed for variant calling or mapping (from a total of 63 distinct fragment references), across influenza components.

Structural summary metrics derived from submitted influenza consensus sequences and VCF files are presented in Table 4. These metrics capture the overall magnitude of reported variats in the metadata file and the discrepancy between reported variants with an allele frequency >= 75% in the metadata file and the VCF file, rather than direct nucleotide-level accuracy against a unified reference coordinate system.


**Table 4. Network-level structural summary of influenza variant reporting.**

| Metric | Network median | Min-max |
|---|---:|---:|
| Variants with AF>=75% | 477.0 | 377.0–1797.0 |
| Variants with AF>=75% in VCF | 188.0 | 0.0–1363.0 |
| Discrepancies in reported variants | 377.0 | 8.0–1797.0 |
| Total variants in VCF | 541.0 | 0.0–39071.0 |



Figure 4 summarizes the distribution of variant reporting practices across participating laboratories for influenza components.

<figure>
<img src="figures/network/influenza_variant_reporting_summary.png" alt="Influenza variant reporting practices across the network." style="max-width: 100%;"/>
<figcaption>Influenza variant reporting practices across the network.</figcaption>
</figure>


**_Figure 4_. Influenza variant reporting characteristics across the network**. Summarise the proportion of laboratories reporting high- and/or low-frequency variants.

Together, Figure 4 and Table 4 show marked heterogeneity in influenza variant reporting within the network. This is reflected by mixed reporting modes across laboratories, the use of multiple reference genomes, and wide ranges in structural summary metrics, including 0.0 to 39071.0 total variants in submitted VCF files and 8.0 to 1797.0 discrepancies between metadata-reported and VCF-derived high-frequency variants.

### 5.4. Lineage, Subtype and Clade Assignment

Lineage, Subtype and clade assignments were evaluated for concordance with gold standard classifications according to [Section 4.4](#44-evaluation-of-lineage-type-and-clade-assignment). Overall concordance rates were:

- SARS-CoV-2 lineage assignment: **77.7%** concordance.
- Influenza type/subtype identification: **81.6%** concordance.
- SARS-CoV-2 clade assignment: **73.4%** concordance.
- Influenza clade assignment: **67.3%** concordance.

Across components, lineage/type concordance was consistently higher than clade concordance. SARS-CoV-2 lineage assignment reached 77.7%, compared with 73.4% for SARS-CoV-2 clade assignment, while influenza type/subtype identification reached 81.6% compared with 67.3% for influenza clade assignment.



<figure>
<img src="figures/network/classification_summary.png" alt="Distribution of classification outcomes across participating laboratories." style="max-width: 100%;"/>
<figcaption>Distribution of classification outcomes across participating laboratories.</figcaption>
</figure>


**_Figure 5_. Distribution of classification outcomes across participating laboratories.** Panel **A** shows **lineage/type assignments**, and panel **B** shows **clade assignments**. Stacked bars represent the total number of classification outcomes across all samples and laboratories for each component. Bars are partitioned into **Hits** (correct assignments relative to the curated gold standard) and **Discrepancies** (incorrect assignments).

### 5.5. Metadata completeness and compliance

The evaluation of metadata focused on analytical transparency, reproducibility, and interoperability within the RELECOV network. Completeness and compliance were assessed according to the criteria defined in [Section 4.5](#45-evaluation-of-metadata-completeness-and-compliance), including controlled vocabulary adherence, logical consistency, and reporting of analytical parameters.

#### Overall Completeness



Across all participating laboratories, the metadata template was completed at a median completeness rate of 59.1%, with values ranging from 11.3% to 92.0%. Component-level median completeness values were similar overall, but the observed ranges remained broad in all components, as shown in Figure 6.

The leading incompleteness drivers were variant calling, pre-processing, and mapping fields, followed by QC metrics, de-hosting, and consensus analysis fields.

<figure>
<img src="figures/network/metadata_completeness_distribution.png" alt="Distribution of metadata completeness across participating laboratories." style="max-width: 100%;"/>
<figcaption>Distribution of metadata completeness across participating laboratories.</figcaption>
</figure>


**_Figure 6_. Distribution of metadata completeness across participating laboratories**. Boxplots represent the distribution of sample-level metadata completeness percentages across the different components. Completeness was calculated for each submitted sample as the proportion of filled metadata fields relative to the total number of maximum expected metadata fields. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

All (100%) laboratories required either clarification through e-mail contact or metadata correction during validation, as reflected by the high proportion of submissions with incomplete parameters or controlled-vocabulary corrections.

#### Reporting of Analytical Parameters

Although core pipeline tools were generally reported, variability was observed in the level of parameter detail provided.

- 71.8% of the maximum software-name fields were completed across submitted samples.
- 60.1% of the maximum software-version fields were completed across submitted samples.
- 45.6% specified minimum coverage thresholds.
- 38.9% reported `variant_calling_params`, containing the potential allele frequency thresholds.
- 57.7% reported the reference genome accession or identifier.

Incomplete parameter reporting limited the ability to fully reconstruct or reproduce analytical workflows in 60.7% of submissions.

#### Controlled Vocabulary Compliance

Compliance with controlled vocabulary requirements was assessed to determine the degree of metadata standardisation achieved across participating laboratories. Only fields expected to contain predefined categorical values were considered in this analysis; free-text fields such as software versions, file names, and parameter descriptions were excluded.

- 26.32% of submissions (5 laboratories) were fully compliant with controlled vocabulary requirements.
- 73.68% (14 laboratories) required at least one manual correction due to the use of a non-standard value in a controlled field.

The most common compliance issues included:

- Use of free-text entries instead of predefined software names in dropdown-based metadata fields, which prevented direct validation against the harmonised template.
- Incorrect or inconsistent completion of lineage, clade, influenza type, or subtype fields.
- Missing mandatory fields requiring subsequent normalisation.

A recurrent source of non-compliance was the use of free-text entries in fields for which predefined dropdown options were available, particularly for software names. This occurred despite the fact that the metadata template, including its controlled-vocabulary dropdowns, had been distributed two weeks before the start of the exercise to allow laboratories to review the available options and request the incorporation of missing software tools into the schema, and despite explicit instructions on how mandatory fields without available data should be completed.

#### Sample Quality Control Assessment

Sample quality control (QC) classifications reported by laboratories (Pass/Fail) were compared against the predefined gold standard QC status for each sample (ECDC or in-silico). QC agreement was evaluated as a binary outcome:

- Match: laboratory QC classification equals the gold standard QC status
- Discrepancy: laboratory QC classification differs from the gold standard QC status

Overall, the network achieved 71.1% QC concordance, corresponding to 64 Matches and 26 Discrepancies across 90 evaluated sample-level QC decisions.



As shown in Figure 7, QC concordance differed across components, ranging from 62.5% in SARS1 to 100.0% in FLU2.

<figure>
<img src="figures/network/qc_match_rate_by_component.png" alt="QC concordance by component (Match vs Discrepancy relative to the gold standard)." style="max-width: 100%;"/>
<figcaption>QC concordance by component (Match vs Discrepancy relative to the gold standard).</figcaption>
</figure>


**_Figure 7_. QC concordance by component relative to the gold standard.** Stacked bars represent the proportion of QC evaluations classified as Match or Discrepancy for each component across participating laboratories.

### 5.6. Pipeline Benchmarking and Comparative Performance

The benchmarking framework as defined in [Section 4.6](#46-pipeline-benchmarking-and-comparative-performance) was designed to assess whether differences in analytical software and parameterisation were associated with measurable variability in performance across participating laboratories.

The submitted metadata documented heterogeneity in:

- Choice of consensus reconstruction software
- Variant calling strategies
- Lineage and clade assignment tools
- Reference genome selection
- Coverage and allele frequency thresholds

#### Diversity of Analytical Workflows

The metadata submissions allowed characterisation of the analytical landscape currently implemented across the RELECOV network.

A total of 10 distinct analytical workflows were identified across participating laboratories, defined as unique combinations of software tools and versions declared in the metadata template.

Substantial diversity was observed in the selection of core analytical tools:

- Consensus reconstruction software ( 16 distinct tools or versions )
- Variant calling tools ( 22 distinct tools or versions )
- SARS-CoV-2 lineage assignment software ( 4 distinct tools or versions )
- Clade assignment software ( 10 distinct tools or versions )
- Influenza type assignment software ( 4 distinct tools or versions )
- Influenza subtype assignment software ( 8 distinct tools or versions )

Comparative performance analyses stratified by component are presented in Section 6, where software-level differences are evaluated within homogeneous analytical contexts (SARS-CoV-2 Illumina, SARS-CoV-2 Nanopore, Influenza Illumina, Influenza Nanopore).

Because performance differed by component and by metric, software-level comparisons are presented in Section 6 within component-specific contexts rather than collapsed into a single cross-component ranking.

This diversity shows that multiple analytical configurations are currently in use across the RELECOV network. These findings highlight the importance of harmonising minimum analytical criteria while preserving methodological flexibility within the network.

## 6. Component-specific Results

This section presents the analytical results disaggregated by component, allowing a detailed assessment of performance within each dataset and sequencing technology. For each component, results are structured according to participation and submission metrics, consensus genome reconstruction performance, variant detection accuracy, and Lineage, Subtype or clade assignment concordance, as applicable.

Component-level analyses enable identification of platform-specific patterns, dataset-dependent challenges, and variability associated with particular sample characteristics. This approach facilitates a more granular interpretation of performance differences observed at the network level and supports targeted harmonisation recommendations.



### 6.1. SARS1 (SARS-CoV-2, Illumina)

#### 6.1.1. Participation and Submissions

A total of 16 laboratories submitted results for the SARS1 component.

- 64 submitted consensus genome sequences (.fasta), where applicable.
- 64 submitted variant call files (.vcf), where applicable.
- The metadata template completeness for SARS1 submissions had a median of 55.8%.

#### 6.1.2. Consensus Genome Reconstruction Performance

Consensus sequences were evaluated against the corresponding curated gold standard references for each sample in the SARS1 component.

Overall, SARS1 showed a median genome identity of 99.59%, with a median of 4.5 nucleotide discrepancies per sample (range: 1.0–125.0).


**Table 5. Network-level consensus reconstruction metrics per sample for SARS1.**

| Sample ID | Median genome identity (%) | Median discrepancies | Discrepancies min-max |
|---|---:|---:|---:|
| SARS1 | N/A | N/A | N/A – N/A |
| SARS2 | 99.82 | 3.0 | 1.0 – 4.0 |
| SARS3 | 99.59 | 3.0 | 2.0 – 11.0 |
| SARS4 | 95.58 | 75.5 | 62.0 – 125.0 |
| SARS5 | 99.80 | 5.5 | 3.0 – 7.0 |


Figure 8 presents the distribution of nucleotide discrepancies per sample across participating laboratories for SARS1.


<figure>
<img src="figures/SARS1/consensus_discrepancies_boxplot_by_sample.png" alt="Consensus discrepancies per sample for SARS1 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Consensus discrepancies per sample for SARS1 relative to the curated gold standard.</figcaption>
</figure>


**Figure 8. Distribution of consensus discrepancies per sample for SARS1.** Boxplots represent the number of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

Considering discrepancy type composition aggregated by sample for SARS1:


**Table 6. Network-level consensus discrepancy types per sample for SARS1.**

| Sample ID | Median of Wrong nucleotide | Median Ambiguity instead of nucleotide | Median Nucleotide instead of ambiguity | Median Stretch of Ns instead of nucleotide stretch | Median Nucleotide stretch instead of stretch of Ns | Median Insertion relative to gold standard | Median Deletion relative to gold standard |
|---|---:|---:|---:|---:|---:|---:|---:|
| SARS1 | N/A | N/A | N/A | N/A | N/A | N/A | N/A
| SARS2 | 0.0 | 0.0 | 0.0 | 1.0 | 0.0 | 0.0 | 0.0
| SARS3 | 0.0 | 0.0 | 0.0 | 2.0 | 0.0 | 0.0 | 0.0
| SARS4 | 0.0 | 0.0 | 0.0 | 38.0 | 29.0 | 0.0 | 0.0
| SARS5 | 2.0 | 0.0 | 0.0 | 1.0 | 0.0 | 0.0 | 1.0


Figure 9 presents the distribution of nucleotide discrepancy types per sample across participating laboratories for SARS1.


<figure>
<img src="figures/SARS1/consensus_discrepancies_stacked_by_sample.png" alt="Consensus discrepancy types per sample for SARS1 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Consensus discrepancy types per sample for SARS1 relative to the curated gold standard.</figcaption>
</figure>


**Figure 9. Distribution of consensus discrepancies per sample for SARS1.** Stacked bars represent the number and type of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample.

Discrepancy type composition aggregated across all submitted consensus sequences for SARS1:


**Table 7. Network-level discrepancy composition by type for SARS1.**

| Discrepancy type | Network median per sample | Min-max occurencies |
|---|---:|---:|
| Incorrect nucleotide | 0.0 | 0.0–8.0 |
| Ambiguity instead of nucleotide | 0.0 | 0.0–0.0 |
| Nucleotide instead of ambiguity | 0.0 | 0.0–1.0 |
| Stretch of Ns instead of nucleotide | 2.0 | 1.0–100.0 |
| Nucleotide stretch instead of stretch of Ns| 0.0 | 0.0–99.0 |
| Insertion relative to gold standard | 0.0 | 0.0–3.0 |
| Deletion relative to gold standard | 0.0 | 0.0–3.0 |

The dominant discrepancy pattern observed in SARS1 was ns2nt.

Figure 10 summarises the contribution of each discrepancy category observed in SARS1 relative to the curated gold standard.


<figure>
<img src="figures/SARS1/consensus_discrepancy_type_boxplot.png" alt="Composition of consensus discrepancy types for SARS1 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Composition of consensus discrepancy types for SARS1 relative to the curated gold standard.</figcaption>
</figure>


**Figure 10. Composition of consensus discrepancy types relative to the curated gold standard for SARS1.** Boxplots represent aggregated discrepancies across all submitted consensus sequences, stratified by discrepancy category. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

#### 6.1.3. Variant Detection Accuracy



Variant call files (.vcf) submitted for the SARS1 component were compared against the curated reference variant set corresponding to each sample in the SARS1 component.

Overall, SARS1 showed a median of 9.0 variant discrepancies per sample (range: 0.0–520.0). The component also showed a median of 65.0 successful hits per sample, with median number of variants with an allele frequency higher than 75% of 61.0 in the metadata and 64.0 in the submitted VCF files. Tables 8 and 9 summarise the descriptive reporting metrics and the qualitative discrepancy profile observed across samples in SARS1. Table 8 summarises the descriptive reporting metrics for SARS1, including successful hits, the number of high-frequency variants reported in the metadata and VCF files, and the concordance between both representations for all variants and effect-annotated variants. Table 9 summarises, for each sample in SARS1, the number of successful reference-variant hits together with the qualitative discrepancy profile, including wrong nucleotide calls, insertions, deletions, missing expected variants, and de novo variants.


**Table 8. Network-level SARS-CoV-2 variant reporting metrics per sample for SARS1.**

| Sample ID | Successful hits | Variants >75% AF in metadata | Variants >75% AF in VCF | Variants with effect in metadata | Variants with effect in VCF | Discrepancies metadata vs VCF | Effect discrepancies metadata vs VCF |
|---|---:|---:|---:|---:|---:|---:|---:|
| SARS1 | 32.0 | 30.0 | 30.0 | 22.0 | 22.0 | 0.0 | 2.0 |
| SARS2 | 62.0 | 61.0 | 61.0 | 45.0 | 45.0 | 0.0 | 1.5 |
| SARS3 | 92.0 | 92.0 | 92.0 | 66.0 | 64.0 | 0.0 | 2.5 |
| SARS4 | 71.0 | 70.0 | 73.0 | 52.0 | 50.5 | 1.0 | 3.0 |
| SARS5 | 63.0 | 56.0 | 56.5 | 46.0 | 42.0 | 0.0 | 7.0 |



**Table 9. Network-level SARS-CoV-2 variant calling profile per sample for SARS1.**

| Sample ID | Successful hits | Median discrepancies | Discrepancies min-max | Wrong nucleotide | Insertions | Deletions | Missing | De novo |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SARS1 | 32.0 | 6.0 | 2.0 – 48.0 | 0.0 | 0.0 | 0.0 | 4.0 | 0.0 |
| SARS2 | 62.0 | 8.0 | 0.0 – 30.0 | 0.0 | 0.0 | 0.0 | 3.0 | 4.0 |
| SARS3 | 92.0 | 20.0 | 0.0 – 70.0 | 0.0 | 0.0 | 0.0 | 3.0 | 7.0 |
| SARS4 | 71.0 | 9.0 | 0.0 – 520.0 | 0.0 | 0.0 | 1.0 | 4.0 | 5.0 |
| SARS5 | 63.0 | 9.0 | 2.0 – 85.0 | 1.0 | 0.0 | 1.0 | 1.0 | 3.0 |


Figure 11 presents the distribution of nucleotide discrepancies per sample across participating laboratories for SARS1.


<figure>
<img src="figures/SARS1/variant_discrepancies_stacked_by_sample.png" alt="Variant discrepancies per sample for SARS1 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Variant discrepancies per sample for SARS1 relative to the curated gold standard.</figcaption>
</figure>

**Figure 11. Distribution of variant discrepancies per sample for SARS1.** Stacked bars represent the number of nucleotide discrepancies and discrepancy types relative to the curated gold standard across participating laboratories for each sample.

Discrepancy type composition (aggregated across all submitted variant calls for SARS1):


**Table 10. Network-level discrepancy composition by type for SARS1.** The discrepancy-type columns correspond to the median count per sample across participating laboratories.

| Discrepancy type | Network median per sample | Network min-max per sample |
|---|---:|---:|
| Incorrect nucleotide | 0.0 | 0.0–3.0 |
| Insertion relative to gold standard | 0.0 | 0.0–6.0 |
| Deletions relative to gold standard | 0.0 | 0.0–46.0 |
| Missing expected variants | 3.0 | 0.0–68.0 |
| De novo variants | 4.0 | 0.0–469.0 |

The dominant discrepancy pattern observed in SARS1 was denovo.

Figure 12 summarises the contribution of each discrepancy category observed in SARS1 relative to the curated gold standard.


<figure>
<img src="figures/SARS1/variant_discrepancy_type_boxplot.png" alt="Composition of variant discrepancy types for SARS1 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Composition of variant discrepancy types for SARS1 relative to the curated gold standard.</figcaption>
</figure>


**Figure 12. Composition of variant discrepancy types relative to the curated gold standard for SARS1.** Boxplots represent aggregated discrepancies across all submitted variant calls, stratified by discrepancy category (incorrect nucleotide, excess ambiguous bases, and indels). The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.



#### 6.1.4. Lineage, Subtype and Clade Assignment

Lineage, subtype and clade assignments submitted for the SARS1 component were evaluated for concordance with the curated gold standard classifications.

Across all participating laboratories:

- Lineage/Subtype matches (lineage/type was correct): 85.9%.
- Clade matches (clade was correct): 68.8%


**Table 11. Network-level classification outcomes per sample for SARS1.**

| Sample ID | Lineage/Subtype matches (%) | Clade matches (%)|
|---|---:|---:|
| SARS1 | N/A | N/A |
| SARS2 | 87.50 | 68.75 |
| SARS3 | 81.25 | 68.75 |
| SARS4 | 87.50 | 68.75 |
| SARS5 | 87.50 | 68.75 |


Table 11 summarises the sample-level lineage/subtype and clade concordance rates for SARS1. Figure 13 presents the distribution of classification outcomes per sample across participating laboratories.


<figure>
<img src="figures/SARS1/typing_outcome_stackedbar_by_sample.png" alt="Classification outcome distribution per sample for SARS1." style="max-width: 100%;"/>
<figcaption>Classification outcome distribution per sample for SARS1.</figcaption>
</figure>


**Figure 13. Classification outcome distribution per sample for SARS1.** Panel A shows the proportion of lineage/subtype assignment matches and discrepancies across participating laboratories for each sample. Panel B shows the corresponding proportion for clade assignments. Stacked bars represent Match and Discrepancy outcomes relative to the curated gold standard classification.

#### 6.1.5. Sample Quality Control Assessment

Laboratory-reported sample QC evaluations (Pass/Fail) for the SARS1 component were compared against the predefined gold standard QC status for each sample. Concordance was assessed as a binary outcome:

- Match: reported QC status equals the gold standard
- Discrepancy: reported QC status differs from the gold standard

Overall, QC concordance for SARS1 was 62.5%, corresponding to 25 Matches and 15 Discrepancies across 40 evaluated QC decisions.


**_Table 12_. Sample-level QC concordance for SARS1.**

| Sample ID | Gold standard QC | % Match | # Matches | # Discrepancies | Total evaluations |
|---|---:|---:|---:|---:|---:|
| SARS1 | Fail | 100.0% | 8 | 0 | 8 |
| SARS2 | Pass | 100.0% | 8 | 0 | 8 |
| SARS3 | Pass | 100.0% | 8 | 0 | 8 |
| SARS4 | Fail | 12.5% | 1 | 7 | 8 |
| SARS5 | Fail | 0.0% | 0 | 8 | 8 |



Table 12 summarises the proportion of laboratories correctly classifying QC status for each sample relative to the gold standard definition, and Figure 14 presents the corresponding sample-level distribution of Match and Discrepancy outcomes within SARS1.

<figure>
<img src="figures/SARS1/qc_match_by_sample.png" alt="Sample-level QC concordance for SARS1 (Match vs Discrepancy relative to the gold standard)." style="max-width: 100%;"/>
<figcaption>Sample-level QC concordance for SARS1 (Match vs Discrepancy relative to the gold standard).</figcaption>
</figure>


**_Figure 14_. Sample-level QC concordance for SARS1 relative to the gold standard.** Bars represent the proportion of Match vs Discrepancy outcomes per sample across participating laboratories. Higher discrepancy rates indicate samples for which laboratories more frequently diverged from the predefined QC status.

#### 6.1.6. Pipeline Benchmarking and Comparative Performance


##### Bioinformatics protocol

Based on metadata submissions, 5 distinct bioinformatics protocols were reported for the SARS1 component.


**Table 13. Performance summary of declared bioinformatics protocols for SARS1.**

| Bioinformatics protocol | Version | N labs | Median genome identity (%) | Median discrepancies | Median metadata completeness (%) | Clade concordance (%) | Lineage/type concordance (%) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Custom pipeline/workflow |  | 5 | 99.63 | 5.0 | 62.7 | 20.0 | 60.0 |
| DRAGEN Targeted Microbial | 1.1.0 | 2 | 99.56 | 4.0 | 44.4 | 100.0 | 87.5 |
| INSaFLU | 2.2.2 | 1 | 99.73 | 2.0 | 44.1 | 100.0 | 100.0 |
| INSaFLU | Web version | 1 | 99.70 | 2.0 | 33.9 | 100.0 | 100.0 |
| nf-core/viralrecon | 3.0.0 | 6 | 99.69 | 6.0 | 91.4 | 83.3 | 100.0 |




Figure 16 summarises the distribution of key performance metrics stratified by declared pipeline configuration.


<figure>
<img src="figures/SARS1/bioinformatics_protocol_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pipeline configuration for SARS1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by pipeline configuration for SARS1.</figcaption>
</figure>


**Figure 16. Distribution of performance metrics by declared pipeline configuration for SARS1.** Multi-panel boxplots summarise sample-level performance stratified by bioinformatics protocols. Panel A displays genome identity (%), Panel B discrepancy counts, Panel C metadata completeness (%), and Panel D exact classification concordance (%). Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### De-hosting software

Based on metadata submissions, 4 distinct de-hosting softwares were reported for the SARS1 component.


**Table 14. Performance summary of declared de-hosting software for SARS1.**

| De-hosting software | Version | N labs | % Host reads |
|---|---:|---:|---:|
| BWA-Mem | V2.2.1 | 1 | N/A |
| Kraken2 | 2.1.6 | 5 | 0.0 |
| Kraken2 | DRAGEN Microbial Amplicon ad hoc version | 1 | N/A |
| bowtie2 | 2.5.1 | 1 | 0.12 |




Figure 18 summarises the percentage of host reads metric stratified by declared dehosting sfotaware version.


<figure>
<img src="figures/SARS1/dehosting_metric_boxplots_by_pipeline.png" alt="Distribution of percentage of host reads metrics by dehosting software version for SARS1." style="max-width: 100%;"/>
<figcaption>Distribution of percentage of host reads metrics by dehosting software version for SARS1.</figcaption>
</figure>

**Figure 18. Distribution of percentage of host reads by declared dehosting software version for SARS1.** Boxplots summarise sample-level percentage of host reads stratified by dehosting software version. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each version, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Preprocessing software

Based on metadata submissions, 6 distinct pre-processing software configurations were reported for the SARS1 component.


**Table 15. Performance summary of declared pre-processing software configurations for SARS1.** The configuration column represents the most frequently reported parameter string among laboratories declaring that software and version.

| Pre-processing software | Version | N labs | Most common configuration | Number of reads sequenced | Reads passing filters |
|---|---:|---:|---:|---:|---:|
| Fastp | 0.20.1 | 1 | --detect_adapter_for_pe --cut_tail --cut_window_size 10 --cut_mean_quality 20 --length_required 35 --json json_file --html html_file --thread 34 | 3024406.0 | 2423968.0 |
| Fastp | 0.24.0 | 1 | N/A | 113415910.0 | N/A |
| Fastp | 1.0.1 | 6 | --detect_adapter_for_pe --cut_front --cut_tail --trim_poly_x --cut_mean_quality 30 --qualified_quality_phred 30 --unqualified_percent_limit 10 --length_required 50 | 3024406.0 | 2115268.0 |
| Trimmomatic | 0.39 | 2 | LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30 | 3024406.0 | 2950068.0 |
| Trimmomatic | 0.39-2 | 1 | N/A | N/A | N/A |
| Trimmomatic | 0.40 | 1 | N/A | 3024406.0 | 2631478.0 |




Figure 20 summarises the distribution of key performance metrics stratified by declared pre-processing software configuration.


<figure>
<img src="figures/SARS1/preprocessing_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pre-processing software configuration for SARS1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by pre-processing software configuration for SARS1.</figcaption>
</figure>

**Figure 20. Distribution of performance metrics by declared pre-processing software configuration for SARS1.** Multi-panel boxplots summarise sample-level performance stratified by pre-processing software. Panel A displays Number of reads sequenced and Panel B Reads passing filters. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Mapping software

Based on metadata submissions, 6 distinct mapping software configurations were reported for the SARS1 component.


**Table 16. Performance summary of declared mapping software configurations for SARS1.**

| Mapping software | Version | N labs | Most common configuration | Depth of coverage threshold | % Reads virus |
|---|---:|---:|---:|---:|---:|
| bowtie2 | 2.5.4 | 6 | --local --very-sensitive-local --seed 1 | 10x | 99.955 |
| bwa mem | 0.7.17-r1188 | 1 | -Y -M -t 34 | 30x | 99.64 |
| bwa mem | 0.7.19 | 1 | N/A | N/A | N/A |
| bwa mem | 0.7.19-r1273 | 1 | -p -Y -v 2 | 10x | N/A |
| bwa mem | V2.2.1 | 1 | -t 4 covid-ba2 | 10X | N/A |
| minimap2 | DRAGEN Microbial Amplicon ad-hoc version | 1 | N/A | N/A | N/A |




Figure 22 summarises the distribution of key performance metrics stratified by declared mapping software configuration.


<figure>
<img src="figures/SARS1/mapping_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by mapping software configuration for SARS1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by mapping software configuration for SARS1.</figcaption>
</figure>

**Figure 22. Distribution of performance metrics by declared mapping software configuration for SARS1.** Boxplots summarise sample-level performance stratified by mapping software. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Assembly software

Based on metadata submissions, 3 distinct assembly software configurations were reported for the SARS1 component.


**Table 17. Performance summary of declared assembly software configurations for SARS1.**

| Assembly software | Version | N labs | Most common configuration | Consnsus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
| Megahit | DRAGEN Microbial Amplicon ad hoc version | 1 | N/A | N/A | 99.5697 |  3.0 |
| SPAdes | 4.0.0 | 1 | --rnaviral --threads 24 | 29808.0 | 99.6049 |  7.0 |
| SPAdes | 4.1.0 | 1 | N/A | 29801.0 | 99.6905 |  6.0 |




Figure 24 summarises the distribution of key performance metrics stratified by declared assembly software configuration.


<figure>
<img src="figures/SARS1/assembly_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by assembly software configuration for SARS1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by assembly software configuration for SARS1.</figcaption>
</figure>

**Figure 24. Distribution of performance metrics by declared assembly software configuration for SARS1.** Multi-panel boxplots summarise sample-level performance stratified by assembly software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Consensus software

Based on metadata submissions, 8 distinct consensus software configurations were reported for the SARS1 component.


**Table 18. Performance summary of declared consensus software configurations for SARS1.**

| Consensus software | Version | N labs | Most common configuration | Consnsus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
| DRAGEN consensus | DRAGEN Microbial Amplicon ad hoc version | 2 | N/A | N/A | 99.5697 |  3.0 |
| Ivar consensus |  | 1 | N/A | N/A | 99.7287 |  2.0 |
| Ivar consensus | 1.4.2 | 1 | -t 0.01 | 29903.0 | 99.6584 |  5.0 |
| Ivar consensus | 1.4.3 | 1 | -q 20 -t 0.8 -m 20 -n N | 29869.0 | 99.7388 |  5.0 |
| Ivar consensus | 1.4.4 | 2 | N/A | N/A | 99.6385 |  6.0 |
| bcftools consensus | 1.2.2 | 1 | -m | 29812.0 | 99.6938 |  6.0 |
| bcftools consensus | 1.21 | 1 | mosdepth -x -Q 1 | awk '$4<=10' | 29.868 | 99.6367 |  5.0 |
| bcftools consensus | 1.22 | 4 | AF > 0.75 | 29885.0 | 99.6938 |  6.0 |




Figure 26 summarises the distribution of key performance metrics stratified by declared consensus software configuration.


<figure>
<img src="figures/SARS1/consensus_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by consensus software configuration for SARS1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by consensus software configuration for SARS1.</figcaption>
</figure>

**Figure 26. Distribution of performance metrics by declared consensus software configuration for SARS1.** Multi-panel boxplots summarise sample-level performance stratified by consensus software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Variant calling software

Based on metadata submissions, 8 distinct variant calling software configurations were reported for the SARS1 component.


**Table 19. Performance summary of declared variant calling software configurations for SARS1.**


| Variant calling software | Version | N labs | Most common configuration | Median high + low freq (%) | Median high freq only (%) | Median low freq only (%) | Median variants (AF >=75%) | Median variants in VCF (AF >=75%) | Median variants with effect | Median variants with effect in VCF | Median metadata-VCF discrepancies | Median effect discrepancies | Median successful hits | Median total discrepancies |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| DRAGEN VariantCaller | DRAGEN Microbial Amplicon ad hoc version | 2 | N/A | 0.0 | 100.0 | 0.0 | N/A | 75.0 | N/A | 52.0 | N/A | N/A | 68.5 | 9.0 |
| Ivar | 1.22 | 1 | N/A | 100.0 | 0.0 | 0.0 | N/A | 61.0 | N/A | 45.0 | N/A | N/A | 65.0 | 2.0 |
| Ivar | 1.4.2 | 1 | -t 0.03 -m 10 | 100.0 | 0.0 | 0.0 | 65.0 | 51.0 | N/A | 37.0 | 3.0 | N/A | 55.0 | 22.0 |
| Ivar | 1.4.3 | 1 | -q 15 -t 0.01 -m 1 -r NC_045512.2.fasta -g NC_045512.2.gff3 | 100.0 | 0.0 | 0.0 | 61.0 | 76.0 | 47.0 | N/A | 4.0 | N/A | 65.0 | 63.0 |
| Ivar | 1.4.4 | 6 | -t 0.25 -q 20 -m 10 | 100.0 | 0.0 | 0.0 | 61.0 | 61.0 | 46.5 | 45.0 | 0.0 | 2.5 | 65.0 | 2.0 |
| Octopus | 0.7.4 | 1 | -P 1 | 0.0 | 100.0 | 0.0 | N/A | 65.0 | N/A | 46.0 | N/A | N/A | 62.0 | 14.0 |
| freebayes | 1.1.0 | 1 | N/A | 100.0 | 0.0 | 0.0 | N/A | 57.0 | N/A | 40.0 | N/A | N/A | 49.0 | 69.0 |
| lofreq | 2.1.5 | 1 | indelqual --dindel // call-parallel --pp-threads 64 --call-indels | 100.0 | 0.0 | 0.0 | N/A | 51.0 | N/A | N/A | N/A | N/A | N/A | N/A |





Figure 28 summarises the distribution of key performance metrics stratified by declared variant calling software configuration.


<figure>
<img src="figures/SARS1/variant_calling_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by variant calling software configuration for SARS1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by variant calling software configuration for SARS1.</figcaption>
</figure>


**Figure 28. Distribution of performance metrics by declared variant calling software configuration for SARS1.** Panel A is a stacked bar chart showing the number of evaluable samples assigned to each allele frequency reporting pattern for each software configuration. Boxplot Panel B displays discrepancies in reported variants with AF >=75% in the submitted VCF, Panel C discrepancies in reported variants with effect, Panel D successful hits, and Panel E total discrepancies. Only panels with evaluable data are shown. In the boxplots, the central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.





##### Clade Assignment Software

Based on metadata submissions, 7 distinct clade assignment software configurations were reported for the SARS1 component.


**Table 20. Performance summary of declared clade assignment software configurations for SARS1.**

| Clade assignment software | Version | N labs | Database version | % of clade match | % of clade discrepancy |
|---|---:|---:|---:|---:|---:|
| Nextclade | 2.12.0 | 1 | 2024-10-17--16-48-48Z | 100.0 | 0.0 |
| Nextclade | 3.11.0 | 1 | N/A | 0.0 | 100.0 |
| Nextclade | 3.13.1 | 1 | 2026-01-06--14-59-32Z | 0.0 | 100.0 |
| Nextclade | 3.18.0 | 1 | 2026-01-14–19-24-43Z | 100.0 | 0.0 |
| Nextclade | 3.18.1 | 7 | 2026-01-06--14-59-32Z | 88.89 | 11.11 |
| Nextclade | 3.8.1 | 1 | 2026-01-06T14:59:32Z | 100.0 | 0.0 |
| Nextclade | 3.9.0 | 1 | 2026-01-06T14:59:32Z | 0.0 | 100.0 |




Figure 30 summarises the distribution of key performance metrics stratified by declared clade assignment software configuration.


<figure>
<img src="figures/SARS1/clade_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by clade assignment software configuration for SARS1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by clade assignment software configuration for SARS1.</figcaption>
</figure>

**Figure 30. Distribution of performance metrics by declared clade assignment software configuration for SARS1.** Multi-panel boxplots summarise sample-level performance stratified by clade assignment software. Panel A displays the % of clade matches and Panel B the % of clade discrepancies. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Lineage Assignment Software Name

Based on metadata submissions, 2 distinct lineage assignment software configurations were reported for the SARS1 component.


**Table 21. Performance summary of declared lineage assignment software configurations for SARS1.**

| Lineage Assignment software | Version | N labs | Database version | % of lineage match | % of lineage discrepancy |
|---|---:|---:|---:|---:|---:|
| Pangolin | 4.3.1 | 6 | 1.37 | 100.0 | 0.0 |
| Pangolin | 4.3.4 | 7 | 4.3.4 | 84.38 | 15.62 |




Figure 32 summarises the distribution of key performance metrics stratified by declared lineage assignment software configuration.


<figure>
<img src="figures/SARS1/lineage_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by lineage assignment software configuration for SARS1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by lineage assignment software configuration for SARS1.</figcaption>
</figure>

**Figure 32. Distribution of performance metrics by declared lineage assignment software configuration for SARS1.** Multi-panel boxplots summarise sample-level performance stratified by lineage assignment software. Panel A displays the % of lineage matches and Panel B the % of lineage discrepancies. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.









### 6.2. SARS2 (SARS-CoV-2, Oxford Nanopore Technologies)

#### 6.2.1. Participation and Submissions

A total of 10 laboratories submitted results for the SARS2 component.

- 37 submitted consensus genome sequences (.fasta), where applicable.
- 35 submitted variant call files (.vcf), where applicable.
- The metadata template completeness for SARS2 submissions had a median of 56.5%.

#### 6.2.2. Consensus Genome Reconstruction Performance

Consensus sequences were evaluated against the corresponding curated gold standard references for each sample in the SARS2 component.

Overall, SARS2 showed a median genome identity of 99.79%, with a median of 8.0 nucleotide discrepancies per sample (range: 0.0–45.0).


**Table 22. Network-level consensus reconstruction metrics per sample for SARS2.**

| Sample ID | Median genome identity (%) | Median discrepancies | Discrepancies min-max |
|---|---:|---:|---:|
| SARS6 | 99.91 | 2.5 | 0.0 – 23.0 |
| SARS7 | 99.88 | 8.5 | 0.0 – 23.0 |
| SARS8 | 98.76 | 14.5 | 2.0 – 30.0 |
| SARS9 | N/A | N/A | N/A – N/A |
| SARS10 | 99.04 | 2.0 | 1.0 – 45.0 |


Figure 33 presents the distribution of nucleotide discrepancies per sample across participating laboratories for SARS2.


<figure>
<img src="figures/SARS2/consensus_discrepancies_boxplot_by_sample.png" alt="Consensus discrepancies per sample for SARS2 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Consensus discrepancies per sample for SARS2 relative to the curated gold standard.</figcaption>
</figure>


**Figure 33. Distribution of consensus discrepancies per sample for SARS2.** Boxplots represent the number of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

Considering discrepancy type composition aggregated by sample for SARS2:


**Table 23. Network-level consensus discrepancy types per sample for SARS2.**

| Sample ID | Median of Wrong nucleotide | Median Ambiguity instead of nucleotide | Median Nucleotide instead of ambiguity | Median Stretch of Ns instead of nucleotide stretch | Median Nucleotide stretch instead of stretch of Ns | Median Insertion relative to gold standard | Median Deletion relative to gold standard |
|---|---:|---:|---:|---:|---:|---:|---:|
| SARS6 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0
| SARS7 | 1.0 | 0.0 | 0.0 | 0.0 | 1.0 | 1.0 | 0.0
| SARS8 | 0.5 | 0.0 | 0.0 | 0.0 | 7.0 | 0.0 | 0.0
| SARS9 | N/A | N/A | N/A | N/A | N/A | N/A | N/A
| SARS10 | 0.0 | 0.0 | 0.0 | 0.0 | 1.0 | 0.0 | 0.0


Figure 34 presents the distribution of nucleotide discrepancy types per sample across participating laboratories for SARS2.


<figure>
<img src="figures/SARS2/consensus_discrepancies_stacked_by_sample.png" alt="Consensus discrepancy types per sample for SARS2 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Consensus discrepancy types per sample for SARS2 relative to the curated gold standard.</figcaption>
</figure>


**Figure 34. Distribution of consensus discrepancies per sample for SARS2.** Stacked bars represent the number and type of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample.

Discrepancy type composition aggregated across all submitted consensus sequences for SARS2:


**Table 24. Network-level discrepancy composition by type for SARS2.**

| Discrepancy type | Network median per sample | Min-max occurencies |
|---|---:|---:|
| Incorrect nucleotide | 0.0 | 0.0–23.0 |
| Ambiguity instead of nucleotide | 0.0 | 0.0–0.0 |
| Nucleotide instead of ambiguity | 0.0 | 0.0–0.0 |
| Stretch of Ns instead of nucleotide | 0.0 | 0.0–42.0 |
| Nucleotide stretch instead of stretch of Ns| 1.0 | 0.0–30.0 |
| Insertion relative to gold standard | 0.0 | 0.0–4.0 |
| Deletion relative to gold standard | 0.0 | 0.0–2.0 |

The dominant discrepancy pattern observed in SARS2 was nt2ns.

Figure 35 summarises the contribution of each discrepancy category observed in SARS2 relative to the curated gold standard.


<figure>
<img src="figures/SARS2/consensus_discrepancy_type_boxplot.png" alt="Composition of consensus discrepancy types for SARS2 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Composition of consensus discrepancy types for SARS2 relative to the curated gold standard.</figcaption>
</figure>


**Figure 35. Composition of consensus discrepancy types relative to the curated gold standard for SARS2.** Boxplots represent aggregated discrepancies across all submitted consensus sequences, stratified by discrepancy category. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

#### 6.2.3. Variant Detection Accuracy



Variant call files (.vcf) submitted for the SARS2 component were compared against the curated reference variant set corresponding to each sample in the SARS2 component.

Overall, SARS2 showed a median of 5.5 variant discrepancies per sample (range: 0.0–203.0). The component also showed a median of 81.5 successful hits per sample, with median number of variants with an allele frequency higher than 75% of 73.0 in the metadata and 73.0 in the submitted VCF files. Tables 25 and 26 summarise the descriptive reporting metrics and the qualitative discrepancy profile observed across samples in SARS2. Table 25 summarises the descriptive reporting metrics for SARS2, including successful hits, the number of high-frequency variants reported in the metadata and VCF files, and the concordance between both representations for all variants and effect-annotated variants. Table 26 summarises, for each sample in SARS2, the number of successful reference-variant hits together with the qualitative discrepancy profile, including wrong nucleotide calls, insertions, deletions, missing expected variants, and de novo variants.


**Table 25. Network-level SARS-CoV-2 variant reporting metrics per sample for SARS2.**

| Sample ID | Successful hits | Variants >75% AF in metadata | Variants >75% AF in VCF | Variants with effect in metadata | Variants with effect in VCF | Discrepancies metadata vs VCF | Effect discrepancies metadata vs VCF |
|---|---:|---:|---:|---:|---:|---:|---:|
| SARS6 | 117.0 | 117.0 | 107.5 | 78.0 | 70.0 | 0.0 | 1.0 |
| SARS7 | 93.0 | 94.0 | 86.0 | 67.0 | 59.5 | 0.0 | 1.0 |
| SARS8 | 87.0 | 93.0 | 86.5 | 66.0 | 60.5 | 0.0 | 1.0 |
| SARS9 | 15.5 | 16.0 | 16.0 | 8.0 | 10.0 | 0.0 | 1.0 |
| SARS10 | 11.0 | 11.0 | 10.5 | 5.0 | 5.0 | 0.0 | 0.0 |



**Table 26. Network-level SARS-CoV-2 variant calling profile per sample for SARS2.**

| Sample ID | Successful hits | Median discrepancies | Discrepancies min-max | Wrong nucleotide | Insertions | Deletions | Missing | De novo |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SARS6 | 117.0 | 2.0 | 0.0 – 32.0 | 0.0 | 0.0 | 0.0 | 1.0 | 0.0 |
| SARS7 | 93.0 | 19.0 | 0.0 – 26.0 | 0.0 | 0.0 | 2.0 | 4.0 | 0.0 |
| SARS8 | 87.0 | 21.0 | 0.0 – 46.0 | 0.0 | 0.0 | 0.0 | 7.0 | 0.0 |
| SARS9 | 15.5 | 6.0 | 0.0 – 157.0 | 0.0 | 0.0 | 0.5 | 0.5 | 4.0 |
| SARS10 | 11.0 | 1.0 | 0.0 – 203.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |


Figure 36 presents the distribution of nucleotide discrepancies per sample across participating laboratories for SARS2.


<figure>
<img src="figures/SARS2/variant_discrepancies_stacked_by_sample.png" alt="Variant discrepancies per sample for SARS2 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Variant discrepancies per sample for SARS2 relative to the curated gold standard.</figcaption>
</figure>

**Figure 36. Distribution of variant discrepancies per sample for SARS2.** Stacked bars represent the number of nucleotide discrepancies and discrepancy types relative to the curated gold standard across participating laboratories for each sample.

Discrepancy type composition (aggregated across all submitted variant calls for SARS2):


**Table 27. Network-level discrepancy composition by type for SARS2.** The discrepancy-type columns correspond to the median count per sample across participating laboratories.

| Discrepancy type | Network median per sample | Network min-max per sample |
|---|---:|---:|
| Incorrect nucleotide | 0.0 | 0.0–0.0 |
| Insertion relative to gold standard | 0.0 | 0.0–27.0 |
| Deletions relative to gold standard | 0.0 | 0.0–25.0 |
| Missing expected variants | 1.0 | 0.0–23.0 |
| De novo variants | 1.0 | 0.0–169.0 |

The dominant discrepancy pattern observed in SARS2 was denovo.

Figure 37 summarises the contribution of each discrepancy category observed in SARS2 relative to the curated gold standard.


<figure>
<img src="figures/SARS2/variant_discrepancy_type_boxplot.png" alt="Composition of variant discrepancy types for SARS2 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Composition of variant discrepancy types for SARS2 relative to the curated gold standard.</figcaption>
</figure>


**Figure 37. Composition of variant discrepancy types relative to the curated gold standard for SARS2.** Boxplots represent aggregated discrepancies across all submitted variant calls, stratified by discrepancy category (incorrect nucleotide, excess ambiguous bases, and indels). The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.



#### 6.2.4. Lineage, Subtype and Clade Assignment

Lineage, subtype and clade assignments submitted for the SARS2 component were evaluated for concordance with the curated gold standard classifications.

Across all participating laboratories:

- Lineage/Subtype matches (lineage/type was correct): 60.0%.
- Clade matches (clade was correct): 83.3%


**Table 28. Network-level classification outcomes per sample for SARS2.**

| Sample ID | Lineage/Subtype matches (%) | Clade matches (%)|
|---|---:|---:|
| SARS6 | 90.00 | 70.00 |
| SARS7 | 80.00 | 90.00 |
| SARS8 | 10.00 | 90.00 |
| SARS9 | N/A | N/A |
| SARS10 | N/A | N/A |


Table 28 summarises the sample-level lineage/subtype and clade concordance rates for SARS2. Figure 38 presents the distribution of classification outcomes per sample across participating laboratories.


<figure>
<img src="figures/SARS2/typing_outcome_stackedbar_by_sample.png" alt="Classification outcome distribution per sample for SARS2." style="max-width: 100%;"/>
<figcaption>Classification outcome distribution per sample for SARS2.</figcaption>
</figure>


**Figure 38. Classification outcome distribution per sample for SARS2.** Panel A shows the proportion of lineage/subtype assignment matches and discrepancies across participating laboratories for each sample. Panel B shows the corresponding proportion for clade assignments. Stacked bars represent Match and Discrepancy outcomes relative to the curated gold standard classification.

#### 6.2.5. Sample Quality Control Assessment

Laboratory-reported sample QC evaluations (Pass/Fail) for the SARS2 component were compared against the predefined gold standard QC status for each sample. Concordance was assessed as a binary outcome:

- Match: reported QC status equals the gold standard
- Discrepancy: reported QC status differs from the gold standard

Overall, QC concordance for SARS2 was 68.0%, corresponding to 17 Matches and 8 Discrepancies across 25 evaluated QC decisions.


**_Table 29_. Sample-level QC concordance for SARS2.**

| Sample ID | Gold standard QC | % Match | # Matches | # Discrepancies | Total evaluations |
|---|---:|---:|---:|---:|---:|
| SARS6 | Pass | 80.0% | 4 | 1 | 5 |
| SARS7 | Pass | 80.0% | 4 | 1 | 5 |
| SARS8 | Fail | 20.0% | 1 | 4 | 5 |
| SARS9 | Fail | 80.0% | 4 | 1 | 5 |
| SARS10 | Fail | 80.0% | 4 | 1 | 5 |



Table 29 summarises the proportion of laboratories correctly classifying QC status for each sample relative to the gold standard definition, and Figure 39 presents the corresponding sample-level distribution of Match and Discrepancy outcomes within SARS2.

<figure>
<img src="figures/SARS2/qc_match_by_sample.png" alt="Sample-level QC concordance for SARS2 (Match vs Discrepancy relative to the gold standard)." style="max-width: 100%;"/>
<figcaption>Sample-level QC concordance for SARS2 (Match vs Discrepancy relative to the gold standard).</figcaption>
</figure>


**_Figure 39_. Sample-level QC concordance for SARS2 relative to the gold standard.** Bars represent the proportion of Match vs Discrepancy outcomes per sample across participating laboratories. Higher discrepancy rates indicate samples for which laboratories more frequently diverged from the predefined QC status.

#### 6.2.6. Pipeline Benchmarking and Comparative Performance


##### Bioinformatics protocol

Based on metadata submissions, 5 distinct bioinformatics protocols were reported for the SARS2 component.


**Table 30. Performance summary of declared bioinformatics protocols for SARS2.**

| Bioinformatics protocol | Version | N labs | Median genome identity (%) | Median discrepancies | Median metadata completeness (%) | Clade concordance (%) | Lineage/type concordance (%) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Artic pipeline | 1.8.5 | 1 | 99.48 | 7.0 | 79.7 | 100.0 | 66.7 |
| Custom pipeline/workflow |  | 2 | 99.51 | 15.0 | 56.8 | 100.0 | 33.3 |
| INSaFLU | 2.2.2 | 1 | 99.88 | 4.0 | 37.1 | 100.0 | 66.7 |
| INSaFLU | Web version | 1 | 99.79 | 8.0 | 30.5 | 66.7 | 66.7 |
| nf-core/viralrecon | 3.0.0 | 4 | 99.88 | 2.5 | 91.5 | 75.0 | 58.3 |




Figure 41 summarises the distribution of key performance metrics stratified by declared pipeline configuration.


<figure>
<img src="figures/SARS2/bioinformatics_protocol_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pipeline configuration for SARS2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by pipeline configuration for SARS2.</figcaption>
</figure>


**Figure 41. Distribution of performance metrics by declared pipeline configuration for SARS2.** Multi-panel boxplots summarise sample-level performance stratified by bioinformatics protocols. Panel A displays genome identity (%), Panel B discrepancy counts, Panel C metadata completeness (%), and Panel D exact classification concordance (%). Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### De-hosting software

Based on metadata submissions, 1 distinct de-hosting softwares were reported for the SARS2 component.


**Table 31. Performance summary of declared de-hosting software for SARS2.**

| De-hosting software | Version | N labs | % Host reads |
|---|---:|---:|---:|
| Kraken2 | 2.1.6 | 4 | 0.01 |




Figure 43 summarises the percentage of host reads metric stratified by declared dehosting sfotaware version.


<figure>
<img src="figures/SARS2/dehosting_metric_boxplots_by_pipeline.png" alt="Distribution of percentage of host reads metrics by dehosting software version for SARS2." style="max-width: 100%;"/>
<figcaption>Distribution of percentage of host reads metrics by dehosting software version for SARS2.</figcaption>
</figure>

**Figure 43. Distribution of percentage of host reads by declared dehosting software version for SARS2.** Boxplots summarise sample-level percentage of host reads stratified by dehosting software version. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each version, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Preprocessing software

Based on metadata submissions, 3 distinct pre-processing software configurations were reported for the SARS2 component.


**Table 32. Performance summary of declared pre-processing software configurations for SARS2.** The configuration column represents the most frequently reported parameter string among laboratories declaring that software and version.

| Pre-processing software | Version | N labs | Most common configuration | Number of reads sequenced | Reads passing filters |
|---|---:|---:|---:|---:|---:|
| Fastp | 1.0.1 | 1 | N/A | N/A | N/A |
| NanoFilt | 2.6.0 | 1 | -q 8 -l 50 --headcrop 30 --tailcrop 30 --maxlength 50000 | N/A | N/A |
| Porechop | 0.3.2pre | 1 | --discard_middle | 71109.0 | 56514.0 |




Figure 45 summarises the distribution of key performance metrics stratified by declared pre-processing software configuration.


<figure>
<img src="figures/SARS2/preprocessing_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pre-processing software configuration for SARS2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by pre-processing software configuration for SARS2.</figcaption>
</figure>

**Figure 45. Distribution of performance metrics by declared pre-processing software configuration for SARS2.** Multi-panel boxplots summarise sample-level performance stratified by pre-processing software. Panel A displays Number of reads sequenced and Panel B Reads passing filters. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Mapping software

Based on metadata submissions, 3 distinct mapping software configurations were reported for the SARS2 component.


**Table 33. Performance summary of declared mapping software configurations for SARS2.**

| Mapping software | Version | N labs | Most common configuration | Depth of coverage threshold | % Reads virus |
|---|---:|---:|---:|---:|---:|
| minimap2 | 2.28-r1209 | 3 | -a -x map-ont -t 12 | 20x | 92.16 |
| minimap2 | 2.30 | 1 | N/A | N/A | N/A |
| minimap2 | 2.30-r1287 | 3 | -x map-ont | 10x | 47.6089 |




Figure 47 summarises the distribution of key performance metrics stratified by declared mapping software configuration.


<figure>
<img src="figures/SARS2/mapping_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by mapping software configuration for SARS2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by mapping software configuration for SARS2.</figcaption>
</figure>

**Figure 47. Distribution of performance metrics by declared mapping software configuration for SARS2.** Boxplots summarise sample-level performance stratified by mapping software. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.






##### Consensus software

Based on metadata submissions, 6 distinct consensus software configurations were reported for the SARS2 component.


**Table 34. Performance summary of declared consensus software configurations for SARS2.**

| Consensus software | Version | N labs | Most common configuration | Consnsus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
| Medaka consensus | 2.2.0 | 1 | --model r941_min_sup_g507 --chunk_len 800 --chunk_ovlp 400 | N/A | 99.4591 |  12.0 |
| Medaka consensus | 2.2.1 | 1 | -m r941_min_hac_variant_g507 -r N | N/A | 99.7823 |  17.0 |
| bcftools consensus | 1.17 | 1 | --min-depth 20 | 29856.0 | 99.5049 |  1.0 |
| bcftools consensus | 1.2.2 | 1 | -m | 29260.0 | 99.8896 |  3.0 |
| bcftools consensus | 1.22 | 3 | --min-depth 20 | 29898.5 | 99.4647 |  7.0 |
| bcftools consensus | 2.5 | 1 | N/A | N/A | 99.8794 |  4.0 |




Figure 49 summarises the distribution of key performance metrics stratified by declared consensus software configuration.


<figure>
<img src="figures/SARS2/consensus_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by consensus software configuration for SARS2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by consensus software configuration for SARS2.</figcaption>
</figure>

**Figure 49. Distribution of performance metrics by declared consensus software configuration for SARS2.** Multi-panel boxplots summarise sample-level performance stratified by consensus software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Variant calling software

Based on metadata submissions, 8 distinct variant calling software configurations were reported for the SARS2 component.


**Table 35. Performance summary of declared variant calling software configurations for SARS2.**


| Variant calling software | Version | N labs | Most common configuration | Median high + low freq (%) | Median high freq only (%) | Median low freq only (%) | Median variants (AF >=75%) | Median variants in VCF (AF >=75%) | Median variants with effect | Median variants with effect in VCF | Median metadata-VCF discrepancies | Median effect discrepancies | Median successful hits | Median total discrepancies |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Clair3 | 1.0.10 | 1 | --enable_long_indel --chunk_size=10000 --haploid_precise --no_phasing_for_fa --threads='12' --platform='ont' --model_path='/opt/conda/bin/models/r941_prom_hac_g360+g422' --include_all_ctgs | 0.0 | 100.0 | 0.0 | 94.0 | 94.0 | 66.0 | 66.0 | 0.0 | 0.0 | 94.0 | 0.0 |
| Clair3 | 1.0.11 | 1 | --enable_long_indel --chunk_size=10000 --haploid_precise –no_phasing_for_fa --model_path='/opt/conda/bin/models/r941_prom_hac_g360+g422' --include_all_ctgs | 0.0 | 100.0 | 0.0 | 72.0 | 72.0 | 46.0 | 47.0 | 0.0 | 1.0 | 72.0 | 22.0 |
| Clair3 | 1.2.0 | 1 | run_clair3.sh --enable_long_indel --chunk_size 10000 --haploid_sensitive --no_phasing_for_fa --bam_fn ${bam} --ref_fn ${reference} --platform ont --model_path “r1041_e82_400bps_sup_g615" --include_all_ctgs --enable_variant_calling_at_sequence_head_and_tail | 0.0 | 100.0 | 0.0 | N/A | 156.0 | N/A | 101.0 | N/A | N/A | 87.0 | 1.0 |
| Clair3 | 1.6.2 | 1 | N/A | 0.0 | 100.0 | 0.0 | N/A | 176.0 | N/A | 123.0 | N/A | N/A | 94.0 | 0.0 |
| Clair3 | artic_1.8.5 | 1 | --platform='ont' --model_path='r1041_e82_400bps_hac_v430' --haploid_sensitive --enable_long_indel --no_phasing_for_fa --enable_variant_calling_at_sequence_head_and_tail | 0.0 | 100.0 | 0.0 | 93.0 | 93.0 | 71.0 | 66.0 | 0.0 | 5.0 | 93.0 | 1.0 |
| Medaka | 1.2.1 | 1 | N/A | 100.0 | 0.0 | 0.0 | N/A | 75.0 | N/A | 52.0 | N/A | N/A | 82.0 | 20.5 |
| Medaka | 2.2.0 | 1 | --model r941_min_sup_variant_g507 --chunk_len 800 --chunk_ovlp 400 | 100.0 | 0.0 | 0.0 | N/A | N/A | N/A | N/A | N/A | N/A | 92.0 | 22.0 |
| Medaka | 2.2.1 | 1 | -m r941_min_hac_variant_g507 | N/A | N/A | N/A | N/A | 0.0 | N/A | 0.0 | N/A | N/A | 79.0 | 46.0 |





Figure 51 summarises the distribution of key performance metrics stratified by declared variant calling software configuration.


<figure>
<img src="figures/SARS2/variant_calling_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by variant calling software configuration for SARS2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by variant calling software configuration for SARS2.</figcaption>
</figure>


**Figure 51. Distribution of performance metrics by declared variant calling software configuration for SARS2.** Panel A is a stacked bar chart showing the number of evaluable samples assigned to each allele frequency reporting pattern for each software configuration. Boxplot Panel B displays discrepancies in reported variants with AF >=75% in the submitted VCF, Panel C discrepancies in reported variants with effect, Panel D successful hits, and Panel E total discrepancies. Only panels with evaluable data are shown. In the boxplots, the central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.





##### Clade Assignment Software

Based on metadata submissions, 6 distinct clade assignment software configurations were reported for the SARS2 component.


**Table 36. Performance summary of declared clade assignment software configurations for SARS2.**

| Clade assignment software | Version | N labs | Database version | % of clade match | % of clade discrepancy |
|---|---:|---:|---:|---:|---:|
| Nextclade | 3.11.0 | 1 | N/A | 0.0 | 100.0 |
| Nextclade | 3.18.0 | 1 | 2026-01-14–19-24-43Z | 100.0 | 0.0 |
| Nextclade | 3.18.1 | 6 | 2026-01-06--14-59-32Z | 90.0 | 10.0 |
| Nextclade | 3.18.2 | 1 | 3.0.0 | 100.0 | 0.0 |
| Nextclade | 3.18.3 | 1 | 3.0.0 | N/A | N/A |
| Nextclade | 3.8.1 | 1 | 2026-01-06T14:59:32Z | 100.0 | 0.0 |




Figure 53 summarises the distribution of key performance metrics stratified by declared clade assignment software configuration.


<figure>
<img src="figures/SARS2/clade_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by clade assignment software configuration for SARS2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by clade assignment software configuration for SARS2.</figcaption>
</figure>

**Figure 53. Distribution of performance metrics by declared clade assignment software configuration for SARS2.** Multi-panel boxplots summarise sample-level performance stratified by clade assignment software. Panel A displays the % of clade matches and Panel B the % of clade discrepancies. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Lineage Assignment Software Name

Based on metadata submissions, 4 distinct lineage assignment software configurations were reported for the SARS2 component.


**Table 37. Performance summary of declared lineage assignment software configurations for SARS2.**

| Lineage Assignment software | Version | N labs | Database version | % of lineage match | % of lineage discrepancy |
|---|---:|---:|---:|---:|---:|
| INSaFLU | 2.2.2 | 1 | 4.3.4 | 66.67 | 33.33 |
| Nextclade | 3.18.1 | 2 | 4.3.4 | 83.33 | 16.67 |
| Pangolin | 4.3.1 | 2 | 1.37 | 50.0 | 50.0 |
| Pangolin | 4.3.4 | 4 | v1.37 | 50.0 | 50.0 |




Figure 55 summarises the distribution of key performance metrics stratified by declared lineage assignment software configuration.


<figure>
<img src="figures/SARS2/lineage_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by lineage assignment software configuration for SARS2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by lineage assignment software configuration for SARS2.</figcaption>
</figure>

**Figure 55. Distribution of performance metrics by declared lineage assignment software configuration for SARS2.** Multi-panel boxplots summarise sample-level performance stratified by lineage assignment software. Panel A displays the % of lineage matches and Panel B the % of lineage discrepancies. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.









### 6.3. FLU1 (Influenza virus, Illumina)

#### 6.3.1. Participation and Submissions

A total of 12 laboratories submitted results for the FLU1 component.

- 47 submitted consensus genome sequences (.fasta), where applicable.
- 38 submitted variant call files (.vcf), where applicable.
- The metadata template completeness for FLU1 submissions had a median of 53.7%.

#### 6.3.2. Consensus Genome Reconstruction Performance

Consensus sequences were evaluated against the corresponding curated gold standard references for each sample in the FLU1 component.

Overall, FLU1 showed a median genome identity of 96.07%, with a median of 21.0 nucleotide discrepancies per sample (range: 1.0–246.0).


**Table 38. Network-level consensus reconstruction metrics per sample for FLU1.**

| Sample ID | Median genome identity (%) | Median discrepancies | Discrepancies min-max |
|---|---:|---:|---:|
| FLU1 | 96.07 | 11.0 | 7.0 – 162.0 |
| FLU2 | 96.15 | 6.5 | 1.0 – 246.0 |
| FLU3 | N/A | N/A | N/A – N/A |
| FLU4 | 96.21 | 23.0 | 15.0 – 188.0 |
| FLU5 | 95.95 | 27.5 | 15.0 – 186.0 |


Figure 56 presents the distribution of nucleotide discrepancies per sample across participating laboratories for FLU1.


<figure>
<img src="figures/FLU1/consensus_discrepancies_boxplot_by_sample.png" alt="Consensus discrepancies per sample for FLU1 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Consensus discrepancies per sample for FLU1 relative to the curated gold standard.</figcaption>
</figure>


**Figure 56. Distribution of consensus discrepancies per sample for FLU1.** Boxplots represent the number of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

Considering discrepancy type composition aggregated by sample for FLU1:


**Table 39. Network-level consensus discrepancy types per sample for FLU1.**

| Sample ID | Median of Wrong nucleotide | Median Ambiguity instead of nucleotide | Median Nucleotide instead of ambiguity | Median Stretch of Ns instead of nucleotide stretch | Median Nucleotide stretch instead of stretch of Ns | Median Insertion relative to gold standard | Median Deletion relative to gold standard |
|---|---:|---:|---:|---:|---:|---:|---:|
| FLU1 | 11.0 | 0.0 | 0.5 | 0.0 | 0.0 | 0.0 | 0.0
| FLU2 | 1.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 5.0
| FLU3 | N/A | N/A | N/A | N/A | N/A | N/A | N/A
| FLU4 | 3.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 15.0
| FLU5 | 5.5 | 0.0 | 0.0 | 0.0 | 2.5 | 0.0 | 15.0


Figure 57 presents the distribution of nucleotide discrepancy types per sample across participating laboratories for FLU1.


<figure>
<img src="figures/FLU1/consensus_discrepancies_stacked_by_sample.png" alt="Consensus discrepancy types per sample for FLU1 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Consensus discrepancy types per sample for FLU1 relative to the curated gold standard.</figcaption>
</figure>


**Figure 57. Distribution of consensus discrepancies per sample for FLU1.** Stacked bars represent the number and type of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample.

Discrepancy type composition aggregated across all submitted consensus sequences for FLU1:


**Table 40. Network-level discrepancy composition by type for FLU1.**

| Discrepancy type | Network median per sample | Min-max occurencies |
|---|---:|---:|
| Incorrect nucleotide | 5.0 | 0.0–162.0 |
| Ambiguity instead of nucleotide | 0.0 | 0.0–0.0 |
| Nucleotide instead of ambiguity | 0.0 | 0.0–208.0 |
| Stretch of Ns instead of nucleotide | 0.0 | 0.0–0.0 |
| Nucleotide stretch instead of stretch of Ns| 0.0 | 0.0–121.0 |
| Insertion relative to gold standard | 0.0 | 0.0–1.0 |
| Deletion relative to gold standard | 12.0 | 0.0–18.0 |

The dominant discrepancy pattern observed in FLU1 was deletions.

Figure 58 summarises the contribution of each discrepancy category observed in FLU1 relative to the curated gold standard.


<figure>
<img src="figures/FLU1/consensus_discrepancy_type_boxplot.png" alt="Composition of consensus discrepancy types for FLU1 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Composition of consensus discrepancy types for FLU1 relative to the curated gold standard.</figcaption>
</figure>


**Figure 58. Composition of consensus discrepancy types relative to the curated gold standard for FLU1.** Boxplots represent aggregated discrepancies across all submitted consensus sequences, stratified by discrepancy category. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

#### 6.3.3. Variant Detection Accuracy



For the FLU1 component, variant evaluation focused on the agreement between variants with allele frequency above 75% reported in the metadata template and those represented in the submitted VCF files, together with the overall number of variants present in the VCF output. Overall, FLU1 showed a median of 485.5 variants with allele frequency above 75% reported in the metadata template, compared with 174.0 corresponding variants represented in the consensus-derived VCF. The median number of discrepancies between both representations was 247.5, while the median total number of variants present in the submitted VCF files was 306.0.


**Table 41. Network-level influenza variant reporting metrics per sample for FLU1.**

| Sample ID | Variants >75% AF in metadata | Variants >75% AF in consensus VCF | Discrepancies between metadata and VCF | Total variants in VCF |
|---|---:|---:|---:|---:|
| FLU1 | 1587.0 | 383.0 | 909.0 | 383.0 |
| FLU2 | 972.0 | 168.0 | 751.5 | 522.0 |
| FLU3 | N/A | 203.0 | N/A | 207.5 |
| FLU4 | 480.0 | 138.5 | 247.5 | 306.5 |
| FLU5 | 467.0 | 174.0 | 247.0 | 174.0 |


Table 41 summarises, for each sample in FLU1, the median number of high-frequency variants reported in the metadata template, the corresponding number represented in the submitted VCF, the discrepancies between both representations, and the total number of variants observed in the submitted VCF files.


**Table 42. Aggregated influenza variant reporting metrics for FLU1.**

| Metric | Network median | Network min-max |
|---|---:|---:|
| Variants >75% AF in metadata | 485.5 | 449.0–1797.0 |
| Variants >75% AF in consensus VCF | 174.0 | 0.0–1356.0 |
| Discrepancies between metadata and VCF | 247.5 | 8.0–1797.0 |
| Total variants in VCF | 306.0 | 0.0–39071.0 |

Figure 59 summarises influenza variant reporting patterns across samples in FLU1, showing the distribution of laboratory-level observations for agreement between metadata-reported and VCF-derived high-frequency variants and for the overall number of variants observed in the submitted VCF files.


<figure>
<img src="figures/FLU1/influenza_variant_reporting_summary_by_sample.png" alt="Influenza variant reporting summary by sample for FLU1." style="max-width: 100%;"/>
<figcaption>Influenza variant reporting summary by sample for FLU1.</figcaption>
</figure>


**Figure 59. Influenza variant reporting summary by sample for FLU1.** Panel A shows, for each sample, the distribution across participating laboratories of the number of variants with allele frequency above 75% reported in the metadata template, the corresponding number represented in the consensus-derived VCF, and the discrepancies between both representations. Panel B shows the distribution across participating laboratories of the total number of variants present in the submitted VCF files for each sample. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range within the plotted scale, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.



#### 6.3.4. Lineage, Subtype and Clade Assignment

Lineage, subtype and clade assignments submitted for the FLU1 component were evaluated for concordance with the curated gold standard classifications.

Across all participating laboratories:

- Lineage/Subtype matches (lineage/type was correct): 89.6%.
- Clade matches (clade was correct): 68.8%


**Table 43. Network-level classification outcomes per sample for FLU1.**

| Sample ID | Lineage/Subtype matches (%) | Clade matches (%)|
|---|---:|---:|
| FLU1 | 91.67 | 58.33 |
| FLU2 | 100.00 | 75.00 |
| FLU3 | N/A | N/A |
| FLU4 | 91.67 | 66.67 |
| FLU5 | 75.00 | 75.00 |


Table 43 summarises the sample-level lineage/subtype and clade concordance rates for FLU1. Figure 60 presents the distribution of classification outcomes per sample across participating laboratories.


<figure>
<img src="figures/FLU1/typing_outcome_stackedbar_by_sample.png" alt="Classification outcome distribution per sample for FLU1." style="max-width: 100%;"/>
<figcaption>Classification outcome distribution per sample for FLU1.</figcaption>
</figure>


**Figure 60. Classification outcome distribution per sample for FLU1.** Panel A shows the proportion of lineage/subtype assignment matches and discrepancies across participating laboratories for each sample. Panel B shows the corresponding proportion for clade assignments. Stacked bars represent Match and Discrepancy outcomes relative to the curated gold standard classification.

#### 6.3.5. Sample Quality Control Assessment

Laboratory-reported sample QC evaluations (Pass/Fail) for the FLU1 component were compared against the predefined gold standard QC status for each sample. Concordance was assessed as a binary outcome:

- Match: reported QC status equals the gold standard
- Discrepancy: reported QC status differs from the gold standard

Overall, QC concordance for FLU1 was 80.0%, corresponding to 12 Matches and 3 Discrepancies across 15 evaluated QC decisions.


**_Table 44_. Sample-level QC concordance for FLU1.**

| Sample ID | Gold standard QC | % Match | # Matches | # Discrepancies | Total evaluations |
|---|---:|---:|---:|---:|---:|
| FLU1 | Pass | 66.7% | 2 | 1 | 3 |
| FLU2 | Pass | 100.0% | 3 | 0 | 3 |
| FLU3 | Fail | 100.0% | 3 | 0 | 3 |
| FLU4 | Pass | 66.7% | 2 | 1 | 3 |
| FLU5 | Fail | 66.7% | 2 | 1 | 3 |



Table 44 summarises the proportion of laboratories correctly classifying QC status for each sample relative to the gold standard definition, and Figure 61 presents the corresponding sample-level distribution of Match and Discrepancy outcomes within FLU1.

<figure>
<img src="figures/FLU1/qc_match_by_sample.png" alt="Sample-level QC concordance for FLU1 (Match vs Discrepancy relative to the gold standard)." style="max-width: 100%;"/>
<figcaption>Sample-level QC concordance for FLU1 (Match vs Discrepancy relative to the gold standard).</figcaption>
</figure>


**_Figure 61_. Sample-level QC concordance for FLU1 relative to the gold standard.** Bars represent the proportion of Match vs Discrepancy outcomes per sample across participating laboratories. Higher discrepancy rates indicate samples for which laboratories more frequently diverged from the predefined QC status.

#### 6.3.6. Pipeline Benchmarking and Comparative Performance


##### Bioinformatics protocol

Based on metadata submissions, 6 distinct bioinformatics protocols were reported for the FLU1 component.


**Table 45. Performance summary of declared bioinformatics protocols for FLU1.**

| Bioinformatics protocol | Version | N labs | Median genome identity (%) | Median discrepancies | Median metadata completeness (%) | Clade concordance (%) | Lineage/type concordance (%) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Custom pipeline/workflow |  | 4 | 96.18 | 11.0 | 60.7 | 43.8 | 93.8 |
| DRAGEN Targeted Microbial | 1.1.0 | 1 | 98.21 | 14.0 | 39.3 | 100.0 | 100.0 |
| INSaFLU | 2.2.2 | 1 | 92.91 | 162.0 | 41.1 | 75.0 | 50.0 |
| IRMA |  | 1 | 96.11 | 11.0 | 32.1 | 100.0 | 100.0 |
| IRMA | 1.2.0 | 3 | 95.89 | 18.5 | 39.3 | 58.3 | 91.7 |
| IRMA | 1.3.1 | 2 | 96.01 | 13.5 | 63.0 | 100.0 | 87.5 |




Figure 63 summarises the distribution of key performance metrics stratified by declared pipeline configuration.


<figure>
<img src="figures/FLU1/bioinformatics_protocol_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pipeline configuration for FLU1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by pipeline configuration for FLU1.</figcaption>
</figure>


**Figure 63. Distribution of performance metrics by declared pipeline configuration for FLU1.** Multi-panel boxplots summarise sample-level performance stratified by bioinformatics protocols. Panel A displays genome identity (%), Panel B discrepancy counts, Panel C metadata completeness (%), and Panel D exact classification concordance (%). Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### De-hosting software

Based on metadata submissions, 5 distinct de-hosting softwares were reported for the FLU1 component.


**Table 46. Performance summary of declared de-hosting software for FLU1.**

| De-hosting software | Version | N labs | % Host reads |
|---|---:|---:|---:|
| Kraken2 | 2.1.3 | 1 | 25.04 |
| Kraken2 | 2.1.6 | 1 | N/A |
| Kraken2 | DRAGEN Microbial Amplicon | 1 | N/A |
| NCBI Human Read Scrubber | 2.0.0 | 1 | 23.69 |
| bowtie2 | 2.5.1 | 1 | 10.24 |




Figure 65 summarises the percentage of host reads metric stratified by declared dehosting sfotaware version.


<figure>
<img src="figures/FLU1/dehosting_metric_boxplots_by_pipeline.png" alt="Distribution of percentage of host reads metrics by dehosting software version for FLU1." style="max-width: 100%;"/>
<figcaption>Distribution of percentage of host reads metrics by dehosting software version for FLU1.</figcaption>
</figure>

**Figure 65. Distribution of percentage of host reads by declared dehosting software version for FLU1.** Boxplots summarise sample-level percentage of host reads stratified by dehosting software version. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each version, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Preprocessing software

Based on metadata submissions, 6 distinct pre-processing software configurations were reported for the FLU1 component.


**Table 47. Performance summary of declared pre-processing software configurations for FLU1.** The configuration column represents the most frequently reported parameter string among laboratories declaring that software and version.

| Pre-processing software | Version | N labs | Most common configuration | Number of reads sequenced | Reads passing filters |
|---|---:|---:|---:|---:|---:|
| Custom preprocessing script |  | 1 | FaQCs --ascii 33 --min_L 50 --avg_q 30 | 1841284.0 | 1060538.0 |
| Fastp | 0.20.0 | 1 | --cut_front --cut_tail --cut_mean_quality 15 --qualified_quality_phred 15 --trim_poly_x --length_required 50 --detect_adapter_for_pe | 1841284.0 | 1841260.0 |
| Fastp | 0.24.0 | 3 | --detect_adapter_for_pe --correction --length_required 45 --qualified_quality_phred 25 --average_qual 30 --cut_front --cut_tail --trim_front1 5 --trim_tail1 5 --trim_front2 5 --trim_tail2 5 | 1841359.0 | 1295030.0 |
| Fastp | 1.0.1 | 1 | N/A | N/A | N/A |
| IRMA custom script | 1.2.0 | 1 | N/A | N/A | N/A |
| Trimmomatic | 0.39 | 1 | LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30 | 1841284.0 | 535168.0 |




Figure 67 summarises the distribution of key performance metrics stratified by declared pre-processing software configuration.


<figure>
<img src="figures/FLU1/preprocessing_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pre-processing software configuration for FLU1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by pre-processing software configuration for FLU1.</figcaption>
</figure>

**Figure 67. Distribution of performance metrics by declared pre-processing software configuration for FLU1.** Multi-panel boxplots summarise sample-level performance stratified by pre-processing software. Panel A displays Number of reads sequenced and Panel B Reads passing filters. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Mapping software

Based on metadata submissions, 6 distinct mapping software configurations were reported for the FLU1 component.


**Table 48. Performance summary of declared mapping software configurations for FLU1.**

| Mapping software | Version | N labs | Most common configuration | Depth of coverage threshold | % Reads virus |
|---|---:|---:|---:|---:|---:|
| BLAT | IRMA 1.3.0-rc1 | 1 | Default (IRMA module FLU_AD) | 1x | 25.97 |
| BLAT | v.36 | 1 | N/A | N/A | N/A |
| BLAT | v35 | 2 | N/A | 10x | 63.72 |
| bwa mem | 0.7.19-r1273 | 2 | -p -Y -v 2 | 10x | N/A |
| minimap2 | 2.26-r1175 | 1 | -ax sr | 20x | N/A |
| minimap2 | DRAGEN Microbial Amplicon ad-hoc version | 1 | N/A | N/A | N/A |




Figure 69 summarises the distribution of key performance metrics stratified by declared mapping software configuration.


<figure>
<img src="figures/FLU1/mapping_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by mapping software configuration for FLU1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by mapping software configuration for FLU1.</figcaption>
</figure>

**Figure 69. Distribution of performance metrics by declared mapping software configuration for FLU1.** Boxplots summarise sample-level performance stratified by mapping software. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Assembly software

Based on metadata submissions, 4 distinct assembly software configurations were reported for the FLU1 component.


**Table 49. Performance summary of declared assembly software configurations for FLU1.**

| Assembly software | Version | N labs | Most common configuration | Consnsus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
| LABEL | IRMA 1.3.0-rc1 | 1 | Default (IRMA module FLU_AD) | 13139.0 | 96.11 |  11.0 |
| LABEL | v0.6.5 | 1 | N/A | 13134.5 | 95.7597 |  41.0 |
| Megahit | DRAGEN Microbial Amplicon ad hoc version | 1 | N/A | N/A | 98.2127 |  14.0 |
| SPAdes |  | 1 | N/A | N/A | 96.11 |  11.0 |




Figure 71 summarises the distribution of key performance metrics stratified by declared assembly software configuration.


<figure>
<img src="figures/FLU1/assembly_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by assembly software configuration for FLU1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by assembly software configuration for FLU1.</figcaption>
</figure>

**Figure 71. Distribution of performance metrics by declared assembly software configuration for FLU1.** Multi-panel boxplots summarise sample-level performance stratified by assembly software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Consensus software

Based on metadata submissions, 7 distinct consensus software configurations were reported for the FLU1 component.


**Table 50. Performance summary of declared consensus software configurations for FLU1.**

| Consensus software | Version | N labs | Most common configuration | Consnsus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
| IRMA consensus |  | 1 | N/A | N/A | 95.7184 |  47.0 |
| IRMA consensus | 1.2.0 | 2 | MIN_AMBIG=0.25 MIN_CONS_SUPPORT=10 | 13134.5 | 96.009 |  13.5 |
| IRMA consensus | 1.3.1 | 1 | MIN_AMBIG=0.75; MIN_CONS_SUPPORT=9 | 13139.0 | 96.11 |  13.5 |
| Ivar consensus |  | 1 | N/A | N/A | 92.912 |  162.0 |
| Ivar consensus | 1.4.3 | 1 | N/A | N/A | 97.7533 |  8.0 |
| bcftools consensus | 1.17 | 1 | AF > 0.50 | 13136.0 | 96.535 |  48.0 |
| bcftools consensus | 1.21 | 1 | mosdepth -x -Q 1 | awk '$4<=10' | 13.1835 | 96.6422 |  1.0 |




Figure 73 summarises the distribution of key performance metrics stratified by declared consensus software configuration.


<figure>
<img src="figures/FLU1/consensus_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by consensus software configuration for FLU1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by consensus software configuration for FLU1.</figcaption>
</figure>

**Figure 73. Distribution of performance metrics by declared consensus software configuration for FLU1.** Multi-panel boxplots summarise sample-level performance stratified by consensus software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Variant calling software

Based on metadata submissions, 6 distinct variant calling software configurations were reported for the FLU1 component.


**Table 51. Performance summary of declared variant calling software configurations for FLU1.**


| Variant calling software | Version | N labs | Most common configuration | Median high + low freq (%) | Median high freq only (%) | Median low freq only (%) | Median variants (AF >=75%) | Median variants in VCF (AF >=75%) | Median variants with effect | Median metadata-VCF discrepancies | Median total variants in VCF |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Custom variant calling script |  | 1 | -f 0.25 -d 10 | 100.0 | 0.0 | 0.0 | 473.5 | 453.0 | 104.5 | 18.0 | 36259.5 |
| IRMA custom VariantCaller |  | 1 | N/A | 0.0 | 0.0 | 100.0 | N/A | 0.0 | N/A | N/A | 187.0 |
| IRMA custom VariantCaller | 1.2.0 | 2 | N/A | 0.0 | 0.0 | 100.0 | 987.0 | 0.0 | 165.0 | 987.0 | 166.0 |
| Ivar | 1.4.3 | 1 | N/A | 100.0 | 0.0 | 0.0 | N/A | 591.0 | N/A | N/A | 623.0 |
| Octopus | 0.7.4 | 1 | -P 1 | 0.0 | 100.0 | 0.0 | N/A | 233.5 | N/A | N/A | 233.5 |
| lofreq | 2.1.5 | 1 | indelqual --dindel // call-parallel --pp-threads 64 --call-indels | 100.0 | 0.0 | 0.0 | N/A | 233.0 | N/A | N/A | 528.5 |





Figure 75 summarises the distribution of key performance metrics stratified by declared variant calling software configuration.


<figure>
<img src="figures/FLU1/variant_calling_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by variant calling software configuration for FLU1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by variant calling software configuration for FLU1.</figcaption>
</figure>


**Figure 75. Distribution of performance metrics by declared variant calling software configuration for FLU1.** Panel A is a stacked bar chart showing the number of evaluable samples assigned to each allele frequency reporting pattern for each software configuration. Boxplot Panel B displays the number of reported variants with AF >=75%, Panel C the number of variants with AF >=75% in the submitted VCF, Panel D the number of variants with effect, Panel E metadata-VCF discrepancies, and Panel F the total number of variants present in the submitted VCF files. Only panels with evaluable data are shown. In the boxplots, the central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.





##### Clade Assignment Software

Based on metadata submissions, 3 distinct clade assignment software configurations were reported for the FLU1 component.


**Table 52. Performance summary of declared clade assignment software configurations for FLU1.**

| Clade assignment software | Version | N labs | Database version | % of clade match | % of clade discrepancy |
|---|---:|---:|---:|---:|---:|
| Nextclade | 3.13.1 | 1 | 2026-01-14--08-53-00Z | 75.0 | 25.0 |
| Nextclade | 3.18.1 | 5 | 2026-01-14--19-24-43Z | 96.3 | 3.7 |
| Nextclade | 3.9.1 | 1 | 2026-01-14--19-24-43Z | 100.0 | 0.0 |




Figure 77 summarises the distribution of key performance metrics stratified by declared clade assignment software configuration.


<figure>
<img src="figures/FLU1/clade_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by clade assignment software configuration for FLU1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by clade assignment software configuration for FLU1.</figcaption>
</figure>

**Figure 77. Distribution of performance metrics by declared clade assignment software configuration for FLU1.** Multi-panel boxplots summarise sample-level performance stratified by clade assignment software. Panel A displays the % of clade matches and Panel B the % of clade discrepancies. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.






##### Type Assignment Software Name

Based on metadata submissions, 3 distinct type assignment software configurations were reported for the FLU1 component.


**Table 53. Performance summary of declared type assignment software configurations for FLU1.**

| Type Assignment software | Version | N labs | Database version | % of type match | % of type discrepancy |
|---|---:|---:|---:|---:|---:|
| Custom typing script |  | 4 | v2 | 87.5 | 12.5 |
| INSaFLU |  | 1 | 2.2.2 | 50.0 | 50.0 |
| IRMA typing |  | 5 | 1.2.0 | 95.0 | 5.0 |




Figure 79 summarises the distribution of key performance metrics stratified by declared xxx software configuration.


<figure>
<img src="figures/FLU1/type_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by type assignment software configuration for FLU1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by type assignment software configuration for FLU1.</figcaption>
</figure>

**Figure 79. Distribution of performance metrics by declared type assignment software configuration for FLU1.** Multi-panel boxplots summarise sample-level performance stratified by type assignment software. Panel A displays the % of type matches and Panel B the % of type discrepancies. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.





##### Subtype Assignment Software Name

Based on metadata submissions, 7 distinct subtype assignment software configurations were reported for the FLU1 component.


**Table 54. Performance summary of declared subtype assignment software configurations for FLU1.**

| Subtype Assignment software | Version | N labs | Database version | % of subtype match | % of subtype discrepancy |
|---|---:|---:|---:|---:|---:|
| Custom subtyping script |  | 3 | v2 | 83.33 | 16.67 |
| INSaFLU | 2.2.2 | 1 | 2.2.2 | 50.0 | 50.0 |
| IRMA subtyping |  | 1 | influenza A & B, genome sets version 2 | 75.0 | 25.0 |
| IRMA subtyping | 1.2.0 | 2 | 1 | 100.0 | 0.0 |
| IRMA subtyping | 1.3.0-rc1 | 1 | 1.3.0-rc1 | 100.0 | 0.0 |
| IRMA subtyping | 1.3.1 | 1 | N/A | 100.0 | 0.0 |
| MASH | 2.3 | 1 | 2.3 | 100.0 | 0.0 |




Figure 81 summarises the distribution of key performance metrics stratified by declared subtype assignment software configuration.


<figure>
<img src="figures/FLU1/subtype_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by subtype assignment software configuration for FLU1." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by subtype assignment software configuration for FLU1.</figcaption>
</figure>

**Figure 81. Distribution of performance metrics by declared subtype assignment software configuration for FLU1.** Multi-panel boxplots summarise sample-level performance stratified by subtype assignment software. Panel A displays the % of subtype matches and Panel B the % of subtype discrepancies. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.





### 6.4. FLU2 (Influenza virus, Oxford Nanopore Technologies)

#### 6.4.1. Participation and Submissions

A total of 10 laboratories submitted results for the FLU2 component.

- 49 submitted consensus genome sequences (.fasta), where applicable.
- 39 submitted variant call files (.vcf), where applicable.
- The metadata template completeness for FLU2 submissions had a median of 55.5%.

#### 6.4.2. Consensus Genome Reconstruction Performance

Consensus sequences were evaluated against the corresponding curated gold standard references for each sample in the FLU2 component.

Overall, FLU2 showed a median genome identity of 95.73%, with a median of 20.0 nucleotide discrepancies per sample (range: 0.0–3051.0).


**Table 55. Network-level consensus reconstruction metrics per sample for FLU2.**

| Sample ID | Median genome identity (%) | Median discrepancies | Discrepancies min-max |
|---|---:|---:|---:|
| FLU6 | 95.69 | 0.0 | 0.0 – 378.0 |
| FLU7 | 95.55 | 41.0 | 7.0 – 348.0 |
| FLU8 | 95.80 | 40.0 | 20.0 – 305.0 |
| FLU9 | 95.47 | 68.0 | 4.0 – 3051.0 |
| FLU10 | 96.05 | 0.0 | 0.0 – 403.0 |


Figure 82 presents the distribution of nucleotide discrepancies per sample across participating laboratories for FLU2.


<figure>
<img src="figures/FLU2/consensus_discrepancies_boxplot_by_sample.png" alt="Consensus discrepancies per sample for FLU2 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Consensus discrepancies per sample for FLU2 relative to the curated gold standard.</figcaption>
</figure>


**Figure 82. Distribution of consensus discrepancies per sample for FLU2.** Boxplots represent the number of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

Considering discrepancy type composition aggregated by sample for FLU2:


**Table 56. Network-level consensus discrepancy types per sample for FLU2.**

| Sample ID | Median of Wrong nucleotide | Median Ambiguity instead of nucleotide | Median Nucleotide instead of ambiguity | Median Stretch of Ns instead of nucleotide stretch | Median Nucleotide stretch instead of stretch of Ns | Median Insertion relative to gold standard | Median Deletion relative to gold standard |
|---|---:|---:|---:|---:|---:|---:|---:|
| FLU6 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0
| FLU7 | 5.5 | 0.0 | 0.0 | 6.5 | 0.0 | 0.0 | 5.0
| FLU8 | 9.0 | 0.0 | 0.0 | 0.0 | 3.0 | 0.0 | 17.0
| FLU9 | 19.0 | 0.0 | 0.0 | 0.0 | 1.0 | 0.0 | 5.5
| FLU10 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0


Figure 83 presents the distribution of nucleotide discrepancy types per sample across participating laboratories for FLU2.


<figure>
<img src="figures/FLU2/consensus_discrepancies_stacked_by_sample.png" alt="Consensus discrepancy types per sample for FLU2 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Consensus discrepancy types per sample for FLU2 relative to the curated gold standard.</figcaption>
</figure>


**Figure 83. Distribution of consensus discrepancies per sample for FLU2.** Stacked bars represent the number and type of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample.

Discrepancy type composition aggregated across all submitted consensus sequences for FLU2:


**Table 57. Network-level discrepancy composition by type for FLU2.**

| Discrepancy type | Network median per sample | Min-max occurencies |
|---|---:|---:|
| Incorrect nucleotide | 3.0 | 0.0–2967.0 |
| Ambiguity instead of nucleotide | 0.0 | 0.0–0.0 |
| Nucleotide instead of ambiguity | 0.0 | 0.0–345.0 |
| Stretch of Ns instead of nucleotide | 0.0 | 0.0–16.0 |
| Nucleotide stretch instead of stretch of Ns| 0.0 | 0.0–64.0 |
| Insertion relative to gold standard | 0.0 | 0.0–45.0 |
| Deletion relative to gold standard | 5.0 | 0.0–70.0 |

The dominant discrepancy pattern observed in FLU2 was deletions.

Figure 84 summarises the contribution of each discrepancy category observed in FLU2 relative to the curated gold standard.


<figure>
<img src="figures/FLU2/consensus_discrepancy_type_boxplot.png" alt="Composition of consensus discrepancy types for FLU2 relative to the curated gold standard." style="max-width: 100%;"/>
<figcaption>Composition of consensus discrepancy types for FLU2 relative to the curated gold standard.</figcaption>
</figure>


**Figure 84. Composition of consensus discrepancy types relative to the curated gold standard for FLU2.** Boxplots represent aggregated discrepancies across all submitted consensus sequences, stratified by discrepancy category. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

#### 6.4.3. Variant Detection Accuracy



For the FLU2 component, variant evaluation focused on the agreement between variants with allele frequency above 75% reported in the metadata template and those represented in the submitted VCF files, together with the overall number of variants present in the VCF output. Overall, FLU2 showed a median of 468.5 variants with allele frequency above 75% reported in the metadata template, compared with 271.5 corresponding variants represented in the consensus-derived VCF. The median number of discrepancies between both representations was 377.0, while the median total number of variants present in the submitted VCF files was 1143.0.


**Table 58. Network-level influenza variant reporting metrics per sample for FLU2.**

| Sample ID | Variants >75% AF in metadata | Variants >75% AF in consensus VCF | Discrepancies between metadata and VCF | Total variants in VCF |
|---|---:|---:|---:|---:|
| FLU6 | 1542.5 | 408.5 | 1651.0 | 2833.0 |
| FLU7 | 430.0 | 303.5 | 217.5 | 970.5 |
| FLU8 | 475.0 | 355.0 | 241.5 | 1000.0 |
| FLU9 | 399.0 | 92.5 | 198.5 | 751.5 |
| FLU10 | 1605.0 | 505.0 | 901.0 | 4541.0 |


Table 58 summarises, for each sample in FLU2, the median number of high-frequency variants reported in the metadata template, the corresponding number represented in the submitted VCF, the discrepancies between both representations, and the total number of variants observed in the submitted VCF files.


**Table 59. Aggregated influenza variant reporting metrics for FLU2.**

| Metric | Network median | Network min-max |
|---|---:|---:|
| Variants >75% AF in metadata | 468.5 | 377.0–1794.0 |
| Variants >75% AF in consensus VCF | 271.5 | 0.0–1363.0 |
| Discrepancies between metadata and VCF | 377.0 | 8.0–1794.0 |
| Total variants in VCF | 1143.0 | 0.0–21455.0 |

Figure 85 summarises influenza variant reporting patterns across samples in FLU2, showing the distribution of laboratory-level observations for agreement between metadata-reported and VCF-derived high-frequency variants and for the overall number of variants observed in the submitted VCF files.


<figure>
<img src="figures/FLU2/influenza_variant_reporting_summary_by_sample.png" alt="Influenza variant reporting summary by sample for FLU2." style="max-width: 100%;"/>
<figcaption>Influenza variant reporting summary by sample for FLU2.</figcaption>
</figure>


**Figure 85. Influenza variant reporting summary by sample for FLU2.** Panel A shows, for each sample, the distribution across participating laboratories of the number of variants with allele frequency above 75% reported in the metadata template, the corresponding number represented in the consensus-derived VCF, and the discrepancies between both representations. Panel B shows the distribution across participating laboratories of the total number of variants present in the submitted VCF files for each sample. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range within the plotted scale, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.



#### 6.4.4. Lineage, Subtype and Clade Assignment

Lineage, subtype and clade assignments submitted for the FLU2 component were evaluated for concordance with the curated gold standard classifications.

Across all participating laboratories:

- Lineage/Subtype matches (lineage/type was correct): 74.0%.
- Clade matches (clade was correct): 66.0%


**Table 60. Network-level classification outcomes per sample for FLU2.**

| Sample ID | Lineage/Subtype matches (%) | Clade matches (%)|
|---|---:|---:|
| FLU6 | 70.00 | 80.00 |
| FLU7 | 90.00 | 80.00 |
| FLU8 | 90.00 | 80.00 |
| FLU9 | 50.00 | 20.00 |
| FLU10 | 70.00 | 70.00 |


Table 60 summarises the sample-level lineage/subtype and clade concordance rates for FLU2. Figure 86 presents the distribution of classification outcomes per sample across participating laboratories.


<figure>
<img src="figures/FLU2/typing_outcome_stackedbar_by_sample.png" alt="Classification outcome distribution per sample for FLU2." style="max-width: 100%;"/>
<figcaption>Classification outcome distribution per sample for FLU2.</figcaption>
</figure>


**Figure 86. Classification outcome distribution per sample for FLU2.** Panel A shows the proportion of lineage/subtype assignment matches and discrepancies across participating laboratories for each sample. Panel B shows the corresponding proportion for clade assignments. Stacked bars represent Match and Discrepancy outcomes relative to the curated gold standard classification.

#### 6.4.5. Sample Quality Control Assessment

Laboratory-reported sample QC evaluations (Pass/Fail) for the FLU2 component were compared against the predefined gold standard QC status for each sample. Concordance was assessed as a binary outcome:

- Match: reported QC status equals the gold standard
- Discrepancy: reported QC status differs from the gold standard

Overall, QC concordance for FLU2 was 100.0%, corresponding to 10 Matches and 0 Discrepancies across 10 evaluated QC decisions.


**_Table 61_. Sample-level QC concordance for FLU2.**

| Sample ID | Gold standard QC | % Match | # Matches | # Discrepancies | Total evaluations |
|---|---:|---:|---:|---:|---:|
| FLU6 | Pass | 100.0% | 2 | 0 | 2 |
| FLU7 | Pass | 100.0% | 2 | 0 | 2 |
| FLU8 | Pass | 100.0% | 2 | 0 | 2 |
| FLU9 | Fail | 100.0% | 2 | 0 | 2 |
| FLU10 | Pass | 100.0% | 2 | 0 | 2 |



Table 61 summarises the proportion of laboratories correctly classifying QC status for each sample relative to the gold standard definition, and Figure 87 presents the corresponding sample-level distribution of Match and Discrepancy outcomes within FLU2.

<figure>
<img src="figures/FLU2/qc_match_by_sample.png" alt="Sample-level QC concordance for FLU2 (Match vs Discrepancy relative to the gold standard)." style="max-width: 100%;"/>
<figcaption>Sample-level QC concordance for FLU2 (Match vs Discrepancy relative to the gold standard).</figcaption>
</figure>


**_Figure 87_. Sample-level QC concordance for FLU2 relative to the gold standard.** Bars represent the proportion of Match vs Discrepancy outcomes per sample across participating laboratories. Higher discrepancy rates indicate samples for which laboratories more frequently diverged from the predefined QC status.

#### 6.4.6. Pipeline Benchmarking and Comparative Performance


##### Bioinformatics protocol

Based on metadata submissions, 6 distinct bioinformatics protocols were reported for the FLU2 component.


**Table 62. Performance summary of declared bioinformatics protocols for FLU2.**

| Bioinformatics protocol | Version | N labs | Median genome identity (%) | Median discrepancies | Median metadata completeness (%) | Clade concordance (%) | Lineage/type concordance (%) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Custom pipeline/workflow |  | 3 | 95.99 | 20.5 | 60.7 | 57.1 | 78.6 |
| INSaFLU | 2.2.2 | 1 | 95.70 | 40.0 | 39.3 | 80.0 | 60.0 |
| INSaFLU | Web version | 1 | 96.54 | 4.0 | 32.1 | 60.0 | 60.0 |
| IRMA | 1.2.0 | 1 | 94.89 | 268.0 | 87.5 | 80.0 | 100.0 |
| IRMA | 1.3.1 | 3 | 95.64 | 6.0 | 64.3 | 81.8 | 90.9 |
| IRMA | v1.3.1 | 1 | 95.90 | 7.0 | 67.9 | 100.0 | 100.0 |




Figure 89 summarises the distribution of key performance metrics stratified by declared pipeline configuration.


<figure>
<img src="figures/FLU2/bioinformatics_protocol_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pipeline configuration for FLU2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by pipeline configuration for FLU2.</figcaption>
</figure>


**Figure 89. Distribution of performance metrics by declared pipeline configuration for FLU2.** Multi-panel boxplots summarise sample-level performance stratified by bioinformatics protocols. Panel A displays genome identity (%), Panel B discrepancy counts, Panel C metadata completeness (%), and Panel D exact classification concordance (%). Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### De-hosting software

Based on metadata submissions, 3 distinct de-hosting softwares were reported for the FLU2 component.


**Table 63. Performance summary of declared de-hosting software for FLU2.**

| De-hosting software | Version | N labs | % Host reads |
|---|---:|---:|---:|
| Kraken2 | 2.1.2 | 1 | 5.8 |
| Kraken2 | 2.1.3 | 1 | 2.7289 |
| Kraken2 | 2.1.6 | 1 | N/A |




Figure 91 summarises the percentage of host reads metric stratified by declared dehosting sfotaware version.


<figure>
<img src="figures/FLU2/dehosting_metric_boxplots_by_pipeline.png" alt="Distribution of percentage of host reads metrics by dehosting software version for FLU2." style="max-width: 100%;"/>
<figcaption>Distribution of percentage of host reads metrics by dehosting software version for FLU2.</figcaption>
</figure>

**Figure 91. Distribution of percentage of host reads by declared dehosting software version for FLU2.** Boxplots summarise sample-level percentage of host reads stratified by dehosting software version. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each version, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Preprocessing software

Based on metadata submissions, 5 distinct pre-processing software configurations were reported for the FLU2 component.


**Table 64. Performance summary of declared pre-processing software configurations for FLU2.** The configuration column represents the most frequently reported parameter string among laboratories declaring that software and version.

| Pre-processing software | Version | N labs | Most common configuration | Number of reads sequenced | Reads passing filters |
|---|---:|---:|---:|---:|---:|
| FiltLong | 0.2.1 | 1 | --min_length 100 --keep_percent 98 | 20489.0 | 17217.0 |
| IRMA custom script |  | 1 | Default IRMA | 20489.0 | 17942.0 |
| IRMA custom script | 1.3.1 | 1 | N/A | N/A | N/A |
| NanoFilt | 2.6.0 | 1 | -q 8 -l 50 --headcrop 30 --tailcrop 30 --maxlength 50000 | N/A | N/A |
| Porechop | v0.2.4 | 2 | -q 7 -l 300 --maxlength 2500 | 20489.0 | 14129.5 |




Figure 93 summarises the distribution of key performance metrics stratified by declared pre-processing software configuration.


<figure>
<img src="figures/FLU2/preprocessing_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pre-processing software configuration for FLU2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by pre-processing software configuration for FLU2.</figcaption>
</figure>

**Figure 93. Distribution of performance metrics by declared pre-processing software configuration for FLU2.** Multi-panel boxplots summarise sample-level performance stratified by pre-processing software. Panel A displays Number of reads sequenced and Panel B Reads passing filters. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Mapping software

Based on metadata submissions, 4 distinct mapping software configurations were reported for the FLU2 component.


**Table 65. Performance summary of declared mapping software configurations for FLU2.**

| Mapping software | Version | N labs | Most common configuration | Depth of coverage threshold | % Reads virus |
|---|---:|---:|---:|---:|---:|
| BLAT |  | 1 | Default IRMA | Default IRMA | 77.6 |
| BLAT | IRMA 1.3.0-rc1 | 1 | Default (IRMA module FLU_AD) | 1x | 33.57 |
| BLAT | v35 | 2 | N/A | 50x | 78.9 |
| minimap2 | 2.30-r1287 | 2 | -ax map-ont | N/A | N/A |




Figure 95 summarises the distribution of key performance metrics stratified by declared mapping software configuration.


<figure>
<img src="figures/FLU2/mapping_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by mapping software configuration for FLU2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by mapping software configuration for FLU2.</figcaption>
</figure>

**Figure 95. Distribution of performance metrics by declared mapping software configuration for FLU2.** Boxplots summarise sample-level performance stratified by mapping software. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Assembly software

Based on metadata submissions, 3 distinct assembly software configurations were reported for the FLU2 component.


**Table 66. Performance summary of declared assembly software configurations for FLU2.**

| Assembly software | Version | N labs | Most common configuration | Consnsus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
| Flye | 2.9.6-b1802 | 1 | --nano-raw --iterations 2 | N/A | 96.0597 |  0.0 |
| LABEL | IRMA 1.3.0-rc1 | 1 | Default (IRMA module FLU_AD) | 13132.0 | 95.9184 |  7.0 |
| LABEL | v0.6.5 | 1 | N/A | 13123.0 | 94.889 |  268.0 |




Figure 97 summarises the distribution of key performance metrics stratified by declared assembly software configuration.


<figure>
<img src="figures/FLU2/assembly_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by assembly software configuration for FLU2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by assembly software configuration for FLU2.</figcaption>
</figure>

**Figure 97. Distribution of performance metrics by declared assembly software configuration for FLU2.** Multi-panel boxplots summarise sample-level performance stratified by assembly software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Consensus software

Based on metadata submissions, 5 distinct consensus software configurations were reported for the FLU2 component.


**Table 67. Performance summary of declared consensus software configurations for FLU2.**

| Consensus software | Version | N labs | Most common configuration | Consnsus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
| IRMA consensus |  | 2 | Default IRMA | N/A | 95.605 |  13.5 |
| IRMA consensus | 1.2.0 | 1 | MIN_AMBIG=0.25 MIN_CONS_QUALITY=10 MIN_CONS_SUPPORT=20 QUAL_THRESHOLD=10 | 13123.0 | 94.889 |  268.0 |
| IRMA consensus | 1.3.1 | 2 | MIN_AMBIG=0.75; MIN_CONS_SUPPORT=9 | 13123.0 | 95.9114 |  3.0 |
| Medaka consensus | 2.2.1 | 1 | -m r941_min_hac_variant_g507 -r N | N/A | 97.1551 |  37.5 |
| bcftools consensus | 1.23 | 1 | N/A | N/A | 95.1296 |  248.0 |




Figure 99 summarises the distribution of key performance metrics stratified by declared consensus software configuration.


<figure>
<img src="figures/FLU2/consensus_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by consensus software configuration for FLU2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by consensus software configuration for FLU2.</figcaption>
</figure>

**Figure 99. Distribution of performance metrics by declared consensus software configuration for FLU2.** Multi-panel boxplots summarise sample-level performance stratified by consensus software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.




##### Variant calling software

Based on metadata submissions, 7 distinct variant calling software configurations were reported for the FLU2 component.


**Table 68. Performance summary of declared variant calling software configurations for FLU2.**


| Variant calling software | Version | N labs | Most common configuration | Median high + low freq (%) | Median high freq only (%) | Median low freq only (%) | Median variants (AF >=75%) | Median variants in VCF (AF >=75%) | Median variants with effect | Median metadata-VCF discrepancies | Median total variants in VCF |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Custom variant calling script |  | 1 | -f 0.25 -d 10 | 100.0 | 0.0 | 0.0 | 463.5 | 419.0 | 95.0 | 18.0 | 4531.0 |
| IRMA custom VariantCaller |  | 1 | N/A | 0.0 | 0.0 | 100.0 | N/A | 0.0 | N/A | N/A | 4840.0 |
| IRMA custom VariantCaller | 1.2.0 | 1 | N/A | 0.0 | 0.0 | 100.0 | 462.0 | 0.0 | 107.0 | 462.0 | 4588.0 |
| IRMA custom VariantCaller | 1.3.1 | 1 | N/A | N/A | N/A | N/A | N/A | N/A | N/A | N/A | 4872.0 |
| Medaka |  | 1 | default MEDAKA | N/A | N/A | N/A | 475.0 | N/A | 130.0 | N/A | 1135.0 |
| Medaka | 2.2.1 | 1 | -m r941_min_hac_variant_g507 | N/A | N/A | N/A | N/A | N/A | N/A | N/A | 564.0 |
| bcftools | 1.23 | 1 | N/A | 100.0 | 0.0 | 0.0 | N/A | 481.0 | N/A | N/A | 804.5 |





Figure 101 summarises the distribution of key performance metrics stratified by declared variant calling software configuration.


<figure>
<img src="figures/FLU2/variant_calling_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by variant calling software configuration for FLU2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by variant calling software configuration for FLU2.</figcaption>
</figure>


**Figure 101. Distribution of performance metrics by declared variant calling software configuration for FLU2.** Panel A is a stacked bar chart showing the number of evaluable samples assigned to each allele frequency reporting pattern for each software configuration. Boxplot Panel B displays the number of reported variants with AF >=75%, Panel C the number of variants with AF >=75% in the submitted VCF, Panel D the number of variants with effect, Panel E metadata-VCF discrepancies, and Panel F the total number of variants present in the submitted VCF files. Only panels with evaluable data are shown. In the boxplots, the central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.





##### Clade Assignment Software

Based on metadata submissions, 2 distinct clade assignment software configurations were reported for the FLU2 component.


**Table 69. Performance summary of declared clade assignment software configurations for FLU2.**

| Clade assignment software | Version | N labs | Database version | % of clade match | % of clade discrepancy |
|---|---:|---:|---:|---:|---:|
| Nextclade | 3.18.1 | 6 | 2026-01-14--19-24-43Z | 87.88 | 12.12 |
| Nextclade | 3.9.1 | 1 | 2025-09-09--12-13-13Z | 100.0 | 0.0 |




Figure 103 summarises the distribution of key performance metrics stratified by declared clade assignment software configuration.


<figure>
<img src="figures/FLU2/clade_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by clade assignment software configuration for FLU2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by clade assignment software configuration for FLU2.</figcaption>
</figure>

**Figure 103. Distribution of performance metrics by declared clade assignment software configuration for FLU2.** Multi-panel boxplots summarise sample-level performance stratified by clade assignment software. Panel A displays the % of clade matches and Panel B the % of clade discrepancies. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.






##### Type Assignment Software Name

Based on metadata submissions, 4 distinct type assignment software configurations were reported for the FLU2 component.


**Table 70. Performance summary of declared type assignment software configurations for FLU2.**

| Type Assignment software | Version | N labs | Database version | % of type match | % of type discrepancy |
|---|---:|---:|---:|---:|---:|
| ABRicate |  | 1 | 2.2.2 | 100.0 | 0.0 |
| Custom typing script |  | 1 | N/A | 40.0 | 60.0 |
| INSaFLU |  | 1 | 2.2.2 | 60.0 | 40.0 |
| IRMA typing |  | 4 | influenza A & B, genome sets version 2 | 95.0 | 5.0 |




Figure 105 summarises the distribution of key performance metrics stratified by declared xxx software configuration.


<figure>
<img src="figures/FLU2/type_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by type assignment software configuration for FLU2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by type assignment software configuration for FLU2.</figcaption>
</figure>

**Figure 105. Distribution of performance metrics by declared type assignment software configuration for FLU2.** Multi-panel boxplots summarise sample-level performance stratified by type assignment software. Panel A displays the % of type matches and Panel B the % of type discrepancies. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.





##### Subtype Assignment Software Name

Based on metadata submissions, 7 distinct subtype assignment software configurations were reported for the FLU2 component.


**Table 71. Performance summary of declared subtype assignment software configurations for FLU2.**

| Subtype Assignment software | Version | N labs | Database version | % of subtype match | % of subtype discrepancy |
|---|---:|---:|---:|---:|---:|
| ABRicate | 1.0.1 | 1 | 2.2.2 | 100.0 | 0.0 |
| Custom subtyping script |  | 1 | N/A | 40.0 | 60.0 |
| INSaFLU | 2.2.2 | 1 | 2.2.2 | 60.0 | 40.0 |
| IRMA subtyping |  | 1 | influenza A & B, genome sets version 2 | 80.0 | 20.0 |
| IRMA subtyping | 1.2.0 | 1 | 1 | 100.0 | 0.0 |
| IRMA subtyping | 1.3.0-rc1 | 1 | 1.3.0-rc1 | 100.0 | 0.0 |
| IRMA subtyping | 1.3.1 | 2 | v1.3.1 | 100.0 | 0.0 |




Figure 107 summarises the distribution of key performance metrics stratified by declared subtype assignment software configuration.


<figure>
<img src="figures/FLU2/subtype_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by subtype assignment software configuration for FLU2." style="max-width: 100%;"/>
<figcaption>Distribution of performance metrics by subtype assignment software configuration for FLU2.</figcaption>
</figure>

**Figure 107. Distribution of performance metrics by declared subtype assignment software configuration for FLU2.** Multi-panel boxplots summarise sample-level performance stratified by subtype assignment software. Panel A displays the % of subtype matches and Panel B the % of subtype discrepancies. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.





## 7. Discussion

The 2026 RELECOV Dry-Lab EQA provides the first network-wide dry-lab assessment focused specifically on bioinformatic performance across consensus reconstruction, variant reporting, classification, metadata reporting, and QC interpretation. By combining ECDC datasets with in-silico influenza material, the exercise captures both routine-use analytical behaviour and performance under heterogeneous reference and reporting conditions.

### 7.1. Consensus Genome Reconstruction

Consensus reconstruction results were strongest in the Illumina-based components overall, with a combined median genome identity of 97.75% compared with 96.10% across the Nanopore-based components. At component level, SARS1 and SARS2 showed median identities of 99.59% and 99.79%, whereas FLU1 and FLU2 showed lower medians of 96.07% and 95.73%.

At the same time, the ranges observed across laboratories show that high medians did not eliminate outlier behaviour. In particular, the minimum identity values in SARS2 and FLU2 dropped to 5.72% and 0.06%, indicating that a subset of submissions diverged markedly from the curated gold standard.

The dominant discrepancy categories also differed by component. SARS1 was dominated by stretches of Ns in submitted consensuses where defined nucleotides were present in the gold standard (`ns2nt`), whereas SARS2 was dominated by defined nucleotides where stretches of Ns were present in the gold standard (`nt2ns`). These patterns are consistent with differences in masking behaviour and minimum coverage policies relative to the gold standard reconstruction criteria. In influenza, both FLU1 and FLU2 were dominated by deletions relative to the gold standard, suggesting that consensus generation parameters remain important contributors to inter-laboratory divergence.

### 7.2. Variant Detection and Reporting

For SARS-CoV-2, variant detection performance did not follow a simple platform ranking. The median number of discrepancies relative to the curated variant set was 9.0 in the Illumina component and 5.5 in the Nanopore component. This indicates that platform effects were present, but that they interacted with sample composition, reporting choices, and software configuration rather than determining performance on their own.

The submitted metadata documented substantial diversity in variant reporting behaviour. Across SARS-CoV-2 submissions, 66.67 of laboratories reported both high- and low-frequency variants, whereas 33.33 reported high-frequency variants only. Influenza reporting was more heterogeneous: 50.68 of laboratories reported both high- and low-frequency variants, 36.99 reported low-frequency variants only, and 12.33 reported high-frequency variants only.

Influenza results especially highlight the consequences of heterogeneous structural reporting. The network-level median number of variants with AF >=75% reported in metadata was 477.0, whereas the corresponding median derived from submitted VCF files was 188.0. The median discrepancy between these two representations was 377.0, and the total number of variants present in submitted VCF files ranged from 0.0 to 39071.0. Together, these values indicate that influenza variant outputs were not directly comparable under a single harmonised coordinate framework and that reporting conventions differed markedly across laboratories.

### 7.3. Classification and QC Interpretation

Classification performance was acceptable overall but clearly stronger for lineage/type assignment than for clade assignment. SARS-CoV-2 lineage concordance reached 77.7%, compared with 73.4% for SARS-CoV-2 clade assignment. Influenza type/subtype concordance reached 81.6%, compared with 67.3% for influenza clade assignment. In the SARS-CoV-2 submissions, this difference is consistent with metadata-level reporting problems in the clade field itself: among the clade assignments reviewed in the submitted JSON files, some were left empty and others contained values that matched the lineage assignment or had lineage-like syntax rather than a clade designation. This suggests that part of the excess clade discordance reflects field completion and nomenclature issues in addition to true analytical misclassification.

QC interpretation showed additional between-component differences. Only 9 of the 19 participating laboratories reported at least one sample-level QC assessment in their submitted metadata. Among evaluable QC decisions, network-wide concordance was 71.1%, and component-level concordance ranged from 62.5% in SARS1 to 100.0% in FLU2. This indicates that QC interpretation was not equally stable across all datasets, while also showing that many laboratories either did not apply or did not report a formal sample-level QC decision in the metadata template.

### 7.4. Workflow Diversity and Reporting Constraints

The metadata confirm that RELECOV laboratories currently use a diverse analytical landscape. A total of 10 distinct workflows were identified across participating laboratories, together with 16 distinct consensus tools or tool/version combinations, 22 distinct variant tools, 4 SARS-CoV-2 lineage assignment tools, and 10 clade assignment tools.

This diversity is analytically valuable, but its interpretation is constrained by incomplete metadata reporting. Only 60.1% of software-version fields were completed, 45.6% of submitted samples specified a minimum coverage threshold, 38.9% reported variant calling parameters, and 57.7% reported a reference genome accession or identifier. For that reason, some plausible explanations for performance differences can only be discussed as contributing context rather than demonstrated causal effects.

The main incompleteness drivers were variant calling, pre-processing, and mapping fields, followed by QC metrics, de-hosting, consensus analysis, and classification-related metadata. This pattern suggests that laboratories were more consistent in declaring core tool identities than in documenting the exact thresholds and parameter sets that determine analytical behaviour.

The submitted metadata also support the view that parameter heterogeneity contributed to consensus and variant calling variability. Across all submitted samples, laboratories reported at least 8 different conventions for minimum coverage thresholds, 13 distinct consensus parameter strings, 13 distinct mapping parameter strings, and 14 distinct variant calling parameter strings. These differences do not prove causality for any individual discrepancy, but they do show that laboratories were not applying a uniform set of masking, filtering, or coverage rules.

### 7.5. Implications for RELECOV 2.0

Taken together, the results support a harmonisation strategy centred on minimum performance and reporting standards rather than on enforcement of a single analytical pipeline. The data do not support a universal workflow ranking that would apply equally across all viruses, platforms, and tasks. Instead, they show that performance depends on the interaction between dataset characteristics, reporting conventions, software choice, and parameterisation.

For RELECOV 2.0, the clearest priorities emerging from this exercise are:

- clearer minimum metadata requirements for reference genomes, software versions, coverage thresholds, and variant calling parameters
- more explicit consensus masking and ambiguity-handling criteria
- harmonised rules for variant reporting thresholds and metadata-versus-VCF consistency
- regular maintenance of classification databases and version tracking
- component-aware benchmarking rather than cross-context software ranking

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

