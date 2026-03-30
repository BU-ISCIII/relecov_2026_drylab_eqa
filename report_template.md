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
{% if path -%}
<figure>
<img src="{{ path }}" alt="{{ caption|default('Figure') }}" style="max-width: 100%;"/>
{% if caption %}<figcaption>{{ caption }}</figcaption>{% endif %}
</figure>
{%- endif %}
{%- endmacro %}
{% macro software_label(name, version=None, db_version=None) -%}
  {% if name %}
    {{ name }}{% if version %} ({{ version }}){% endif %}{% if db_version %}; DB {{ db_version }}{% endif %}
  {% else %}
    NA
  {% endif %}
{%- endmacro %}

{% set fig_counter = namespace(value=0) %}
{% set table_counter = namespace(value=0) %}

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
{% if labdata %}
- [9. Individual Laboratory Technical Report](#9-individual-laboratory-technical-report)
{% endif %}

## Executive Summary

The 2026 RELECOV Dry-Lab EQA provides the first network-wide dry-lab assessment focused specifically on bioinformatic analytical performance across respiratory virus surveillance workflows. ADD two or three lines summarizing samples sended and eqa scope.

Nineteen laboratories participated, corresponding to {{ pct(general.total_participants_pct, 2) }} of invited laboratories, with high submission rates for expected analytical outputs: {{ pct(general.submission_rates_pct.fasta, 2) }} for consensus genome files and {{ pct(general.submission_rates_pct.vcf, 2) }} for VCF files.

Across the network, consensus genome reconstruction performed best in the Illumina-based components, with a combined median genome identity of {{ pct(general.general_results.consensus.median_identity_illumina_pct, 2) }}, compared with {{ pct(general.general_results.consensus.median_identity_nanopore_pct, 2) }} in the Nanopore-based components. However, broad identity ranges in SARS2 and FLU2 indicate that outlier submissions remained present, particularly in contexts where masking, coverage thresholds, and consensus-generation choices differed across laboratories.

Variant reporting showed clear methodological heterogeneity. For SARS-CoV-2, the median number of discrepancies relative to the curated reference variant sets was {{ general.general_results.sars_variants.median_discrepancy_illumina }} in the Illumina component and {{ general.general_results.sars_variants.median_discrepancy_nanopore }} in the Nanopore component. Influenza submissions were more heterogeneous structurally, with a median of {{ general.general_results.influenza_variants.median_variants_in_consensus }} high-frequency variants reported in metadata, compared with {{ general.general_results.influenza_variants.median_variants_in_consensus_vcf }} derived from submitted VCF files, and a median discrepancy of {{ general.general_results.influenza_variants.median_discrepancies_in_reported_variants }} between both representations.

Classification performance was consistently higher for lineage/type assignment than for clade assignment. SARS-CoV-2 lineage concordance reached {{ pct(general.general_results.classification.sars_lineage_concordance_pct) }}, compared with {{ pct(general.general_results.classification.sars_clade_concordance_pct) }} for clade assignment, while influenza type/subtype concordance reached {{ pct(general.general_results.classification.influenza_type_concordance_pct) }}, compared with {{ pct(general.general_results.classification.influenza_clade_concordance_pct) }} for clade assignment. Review of submitted files further indicated that part of the excess discordance in clade assignment reflected field completion and nomenclature issues, including missing clade entries and lineage/type-like values entered in the clade field.

Metadata completeness and reporting remain major priorities for harmonisation. The median metadata completeness rate across participating laboratories was {{ pct(general.metadata_completeness.median_pct) }}, with values ranging from {{ pct(general.metadata_completeness.min_pct) }} to {{ pct(general.metadata_completeness.max_pct) }}. Although software names were reported for {{ pct(general.metadata_completeness.software_names_pct) }} of expected fields, only {{ pct(general.metadata_completeness.software_version_pct) }} of software-version fields, {{ pct(general.metadata_completeness.coverage_threshold_pct) }} of coverage thresholds, {{ pct(general.metadata_completeness.variant_calling_params_pct) }} of variant-calling parameter fields, and {{ pct(general.metadata_completeness.reference_genome_pct) }} of reference genome identifiers were completed. A total of {{ general.metadata_completeness.total_workflows }} distinct workflows were identified, together with substantial diversity in consensus, variant calling, and classification software.

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

{% set table_counter.value = table_counter.value + 1 %}

Table {{ table_counter.value }} summarises the correspondence between RELECOV EQA samples and their original source datasets, including ECDC ESIB references.

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

Table {{ table_counter.value }} describes the design characteristics of in-silico samples, including virus composition and intended benchmarking challenges.

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

Table {{ table_counter.value }} summarises the influenza datasets included in the EQA, detailing enrichment strategy, primer scheme, sequencing technology, and key analytical challenges.

_**Table {{ table_counter.value }}**. Influenza virus samples used in the RELECOV 2026 Dry-Lab EQA, including sequencing platform, enrichment strategy, primer scheme, and key analytical features._

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

Potential contributors discussed during result interpretation included:

- Use of outdated lineage databases
- Differences in software versioning
- Incomplete consensus sequences
- Excess ambiguous positions

Failure to identify virus presence in positive samples, or misclassification of negative samples, was recorded separately.

### 4.5. Evaluation of Metadata Completeness and Compliance

Metadata assessment focused on analytical transparency and interoperability rather than biological correctness.

For each submitted sample, metadata completeness was calculated as:

$$
\text{Sample metadata completeness} =
\frac{\text{Number of correctly populated fields}}
{\text{Total number of applicable fields}}
$$

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

$$
\text{QC concordance rate} =
\frac{\text{Number of Matches}}
{\text{Total QC evaluations}}
$$

QC evaluations were calculated only for samples analysed by the laboratory.

The QC assessment evaluation was limited to concordance analysis. The exercise did not attempt to infer the internal QC criteria applied by laboratories, but rather assessed agreement with the predefined gold standard QC status to evaluate interpretative consistency across the network.

INCLUDE WHAT WE EVALUATE AS GOOD QUALITY?

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

The median number of components analysed per participating laboratory was {{ general.median_components_analysed_per_lab }}.

### 5.1. Submission Completeness

Assessment of submission completeness was conducted in accordance with the criteria outlined in [Section 4.1](#41-submission-completeness). Across all components:

- {{ pct(general.submission_rates_pct.fasta) }} of laboratories submitted consensus genome files (.fasta), where applicable.
- {{ pct(general.submission_rates_pct.vcf) }} submitted variant call files (.vcf), where applicable.

Component-level submission totals are presented in Section 6 and reflect both the number of participating laboratories and the expected output files for each dataset.

### 5.2. Consensus Genome Reconstruction Performance

Consensus genome reconstruction performance was measured using the evaluation criteria detailed in [Section 4.2](#42-evaluation-of-consensus-genome-reconstruction-performance). Across the two Illumina-based components, the combined median genome identity was {{ pct(general.general_results.consensus.median_identity_illumina_pct, 2) }}, compared with {{ pct(general.general_results.consensus.median_identity_nanopore_pct, 2) }} across the two Nanopore-based components. Nanopore-based datasets also showed broader overall identity ranges, where low-identity outliers were present.

Dominant discrepancy patterns differed by component:

- In SARS1, the most frequent discrepancy category was stretches of Ns in the submitted consensus where defined nucleotides were present in the gold standard.
- In SARS2, the most frequent discrepancy category was defined nucleotides in the submitted consensus where stretches of Ns were present in the gold standard.
- FLU1 and FLU2 were both dominated by deletions relative to the gold standard.

These SARS-CoV-2 discrepancy patterns are consistent with differences in masking behaviour and/or minimum coverage thresholds relative to the gold standard reconstruction criteria. (IF YOU SEPARATED DISCUSSION THIS WOULD BE DISCUSSION)

Across components, many discrepancy categories had medians of zero, indicating that errors tended to be concentrated in a smaller number of laboratories or samples rather than being uniformly distributed across the network.

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} summarises consensus genome reconstruction performance across all components.


{{ render_figure(general.figures.consensus_summary, "Network-level consensus reconstruction performance summary.") }}

**_Figure {{ fig_counter.value }}_. Distribution of consensus genome discrepancies relative to the gold standard across components**. Boxplots represent the number of nucleotide discrepancies per component across participating laboratories. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

### 5.3. Variant Detection Accuracy

Variant detection accuracy was evaluated following the methodological framework described in [Section 4.3](#43-evaluation-of-variant-detectio-accuracy).

#### 5.3.1. SARS-CoV-2

For SARS-CoV-2 components (SARS1 and SARS2), variant detection accuracy was assessed against curated reference variant sets. Overall, submitted VCFs showed a median number of discrepancies of {{ general.general_results.sars_variants.median_discrepancy_illumina }} for Illumina component and a median number of {{ general.general_results.sars_variants.median_discrepancy_nanopore }} for Nanopore component, discrepancies relative to the reference variant set.

{% set fig_counter.value = fig_counter.value + 1 %}

The distribution of variant detection performance across components is presented in Figure {{ fig_counter.value }}. Contextual factors documented in the metadata that may contribute to these differences included:

- Allele frequency thresholds used for incorporation into vcf files
- Reference genome selection
- Variant normalization practices (variant caller software and params)

{{ render_figure(general.figures.variant_summary, "Network-level variant detection performance summary.") }}


**_Figure {{ fig_counter.value }}_. SARS-CoV-2 network-level variant detection performance summary**. Boxplots represent the number of variant discrepancies per SARS-CoV-2 component across participating laboratories. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

Variant evaluation included structural reporting characteristics and methodological heterogeneity. At network level:

- {{ general.general_results.sars_variants.high_and_low_freq_pct }} of laboratories reported both high- and low-frequency variants.
- {{ general.general_results.sars_variants.low_freq_only_pct }} reported exclusively low-frequency variants.
- {{ general.general_results.sars_variants.high_freq_only_pct }} reported only high-frequency variants.

Additionally, a total of {{ general.general_results.sars_variants.total_distinct_references }} distinct reference genomes were employed for variant calling across SARS-CoV-2 components.

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} summarizes the distribution of variant reporting practices across participating laboratories for SARS-CoV-2 components.

{{ render_figure(
general.figures.sars_variant_reporting_summary,
"SARS-CoV-2 variant reporting practices across the network."
) }}


**_Figure {{ fig_counter.value }}_. SARS-CoV-2 variant reporting characteristics across the network**. Summarise the proportion of laboratories reporting high- and/or low-frequency variants.

#### 5.3.2. Influenza virus

For influenza virus components (FLU1 and FLU2), variant evaluation focused on structural reporting characteristics and methodological heterogeneity.

At network level:

- {{ general.general_results.influenza_variants.high_and_low_freq_pct }} of laboratories reported both high- and low-frequency variants.
- {{ general.general_results.influenza_variants.low_freq_only_pct }} reported exclusively low-frequency variants.
- {{ general.general_results.influenza_variants.high_freq_only_pct }} reported only high-frequency variants.

Additionally, a total of {{ general.general_results.influenza_variants.total_distinct_references }} distinct reference genomes were employed for variant calling or mapping (from a total of {{ general.general_results.influenza_variants.total_distinct_fragments }} distinct fragment references), across influenza components.

Structural summary metrics derived from submitted influenza consensus sequences and VCF files are presented in Table {{ table_counter.value + 1 }}. These metrics capture the overall magnitude of reported variats in the metadata file and the discrepancy between reported variants with an allele frequency >= 75% in the metadata file and the VCF file, rather than direct nucleotide-level accuracy against a unified reference coordinate system.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level structural summary of influenza variant reporting.**

| Metric | Network median | Min-max |
|---|---:|---:|
| Variants with AF>=75% | {{ general.general_results.influenza_variants.median_variants_in_consensus }} | {{ general.general_results.influenza_variants.min_variants_in_consensus }}–{{ general.general_results.influenza_variants.max_variants_in_consensus }} |
| Variants with AF>=75% in VCF | {{ general.general_results.influenza_variants.median_variants_in_consensus_vcf }} | {{ general.general_results.influenza_variants.min_variants_in_consensus_vcf }}–{{ general.general_results.influenza_variants.max_variants_in_consensus_vcf }} |
| Discrepancies in reported variants | {{ general.general_results.influenza_variants.median_discrepancies_in_reported_variants }} | {{ general.general_results.influenza_variants.min_discrepancies_in_reported_variants }}–{{ general.general_results.influenza_variants.max_discrepancies_in_reported_variants }} |
| Total variants in VCF | {{ general.general_results.influenza_variants.median_variants_in_vcf }} | {{ general.general_results.influenza_variants.min_variants_in_vcf }}–{{ general.general_results.influenza_variants.max_variants_in_vcf }} |

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} summarizes the distribution of variant reporting practices across participating laboratories for influenza components.

{{ render_figure(
general.figures.influenza_variant_reporting_summary,
"Influenza variant reporting practices across the network."
) }}


**_Figure {{ fig_counter.value }}_. Influenza variant reporting characteristics across the network**. Summarise the proportion of laboratories reporting high- and/or low-frequency variants.

Together, Figure {{ fig_counter.value }} and Table {{ table_counter.value }} show marked heterogeneity in influenza variant reporting within the network. This is reflected by mixed reporting modes across laboratories, the use of multiple reference genomes, and wide ranges in structural summary metrics, including {{ general.general_results.influenza_variants.min_variants_in_vcf }} to {{ general.general_results.influenza_variants.max_variants_in_vcf }} total variants in submitted VCF files and {{ general.general_results.influenza_variants.min_discrepancies_in_reported_variants }} to {{ general.general_results.influenza_variants.max_discrepancies_in_reported_variants }} discrepancies between metadata-reported and VCF-derived high-frequency variants.

### 5.4. Lineage, Subtype and Clade Assignment

Lineage, Subtype and clade assignments were evaluated for concordance with gold standard classifications according to [Section 4.4](#44-evaluation-of-lineage-type-and-clade-assignment). Overall concordance rates were:

- SARS-CoV-2 lineage assignment: **{{ pct(general.general_results.classification.sars_lineage_concordance_pct) }}** concordance.
- Influenza type/subtype identification: **{{ pct(general.general_results.classification.influenza_type_concordance_pct) }}** concordance.
- SARS-CoV-2 clade assignment: **{{ pct(general.general_results.classification.sars_clade_concordance_pct) }}** concordance.
- Influenza clade assignment: **{{ pct(general.general_results.classification.influenza_clade_concordance_pct) }}** concordance.

Across components, lineage/type concordance was consistently higher than clade concordance. SARS-CoV-2 lineage assignment reached {{ pct(general.general_results.classification.sars_lineage_concordance_pct) }}, compared with {{ pct(general.general_results.classification.sars_clade_concordance_pct) }} for SARS-CoV-2 clade assignment, while influenza type/subtype identification reached {{ pct(general.general_results.classification.influenza_type_concordance_pct) }} compared with {{ pct(general.general_results.classification.influenza_clade_concordance_pct) }} for influenza clade assignment.

{% set fig_counter.value = fig_counter.value + 1 %}

{{ render_figure(general.figures.classification_summary, "Distribution of classification outcomes across participating laboratories.") }}


**_Figure {{ fig_counter.value }}_. Distribution of classification outcomes across participating laboratories.** Panel **A** shows **lineage/type assignments**, and panel **B** shows **clade assignments**. Stacked bars represent the percentage of all possible sample-level classifications across participating laboratories for each component. Bars are partitioned into **Match** (correct assignments relative to the curated gold standard), **Discrepancy** (incorrect assignments), and **Not provided** (classification not reported).

### 5.5. Metadata completeness and compliance

The evaluation of metadata focused on analytical transparency, reproducibility, and interoperability within the RELECOV network. Completeness and compliance were assessed according to the criteria defined in [Section 4.5](#45-evaluation-of-metadata-completeness-and-compliance), including controlled vocabulary adherence, logical consistency, and reporting of analytical parameters.

#### Overall Completeness

{% set fig_counter.value = fig_counter.value + 1 %}

Across all participating laboratories, the metadata template was completed at a median completeness rate of {{ pct(general.metadata_completeness.median_pct) }}, with values ranging from {{ pct(general.metadata_completeness.min_pct) }} to {{ pct(general.metadata_completeness.max_pct) }}. Component-level median completeness values were similar overall, but the observed ranges remained broad in all components, as shown in Figure {{ fig_counter.value }}.

The leading incompleteness drivers were variant calling, pre-processing, and mapping fields, followed by QC metrics, de-hosting, and consensus analysis fields.

{{ render_figure(general.figures.metadata_completeness_distribution,
  "Distribution of metadata completeness across participating laboratories.") }}


**_Figure {{ fig_counter.value }}_. Distribution of metadata completeness across participating laboratories**. Boxplots represent the distribution of sample-level metadata completeness percentages across the different components. Completeness was calculated for each submitted sample as the proportion of filled metadata fields relative to the total number of maximum expected metadata fields. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

All (100%) laboratories required either clarification through e-mail contact or metadata correction during validation, as reflected by the high proportion of submissions with incomplete parameters or controlled-vocabulary corrections.

#### Reporting of Analytical Parameters

Although core pipeline tools were generally reported, variability was observed in the level of parameter detail provided.

- {{ pct(general.metadata_completeness.software_names_pct) }} of the maximum software-name fields were completed across submitted samples.
- {{ pct(general.metadata_completeness.software_version_pct) }} of the maximum software-version fields were completed across submitted samples.
- {{ pct(general.metadata_completeness.coverage_threshold_pct) }} specified minimum coverage thresholds.
- {{ pct(general.metadata_completeness.variant_calling_params_pct) }} reported variant calling parameters, containing the potential allele frequency thresholds.
- {{ pct(general.metadata_completeness.reference_genome_pct) }} reported the reference genome accession or identifier.

Incomplete parameter reporting limited the ability to fully reconstruct or reproduce analytical workflows in {{ pct(general.metadata_completeness.incomplete_parameters_pct) }} of submissions.

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

Overall, the network achieved {{ pct(general.qc.match_rate_pct) }} QC concordance, corresponding to {{ general.qc.matches }} Matches and {{ general.qc.discrepancies }} Discrepancies across {{ general.qc.total_evaluations }} evaluated sample-level QC decisions.

{% set fig_counter.value = fig_counter.value + 1 %}

As shown in Figure {{ fig_counter.value }}, QC concordance differed across components, ranging from {{ pct(general.components.SARS1.qc.match_rate_pct) }} in SARS1 to {{ pct(general.components.FLU2.qc.match_rate_pct) }} in FLU2.

{{ render_figure(
general.figures.qc_match_rate_by_component,
"QC concordance by component (Match vs Discrepancy relative to the gold standard)."
) }}


**_Figure {{ fig_counter.value }}_. QC concordance by component relative to the gold standard.** Stacked bars represent the proportion of QC evaluations classified as Match or Discrepancy for each component across participating laboratories.

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

A total of {{ general.metadata_completeness.total_workflows }} distinct analytical workflows were identified across participating laboratories, defined as unique combinations of software tools and versions declared in the metadata template.

Substantial diversity was observed in the selection of core analytical tools:

- Consensus reconstruction software ( {{ general.metadata_completeness.total_consensus_softwares }} distinct tools or versions )
- Variant calling tools ( {{ general.metadata_completeness.total_variant_softwares }} distinct tools or versions )
- SARS-CoV-2 lineage assignment software ( {{ general.metadata_completeness.total_lineage_assignment_softwares }} distinct tools or versions )
- Clade assignment software ( {{ general.metadata_completeness.total_clade_assignment_softwares }} distinct tools or versions )
- Influenza type assignment software ( {{ general.metadata_completeness.total_type_assignment_softwares }} distinct tools or versions )
- Influenza subtype assignment software ( {{ general.metadata_completeness.total_subtype_assignment_softwares }} distinct tools or versions )

Comparative performance analyses stratified by component are presented in Section 6, where software-level differences are evaluated within homogeneous analytical contexts (SARS-CoV-2 Illumina, SARS-CoV-2 Nanopore, Influenza Illumina, Influenza Nanopore).

Because performance differed by component and by metric, software-level comparisons are presented in Section 6 within component-specific contexts rather than collapsed into a single cross-component ranking.

This diversity shows that multiple analytical configurations are currently in use across the RELECOV network. These findings highlight the importance of harmonising minimum analytical criteria while preserving methodological flexibility within the network.

## 6. Component-specific Results

This section presents the analytical results disaggregated by component, allowing a detailed assessment of performance within each dataset and sequencing technology. For each component, results are structured according to participation and submission metrics, consensus genome reconstruction performance, variant detection accuracy, and Lineage, Subtype or clade assignment concordance, as applicable.

Component-level analyses enable identification of platform-specific patterns, dataset-dependent challenges, and variability associated with particular sample characteristics. This approach facilitates a more granular interpretation of performance differences observed at the network level and supports targeted harmonisation recommendations.

{% for comp_code, comp_net in general.components.items() %}

### 6.{{ loop.index }}. {{ comp_code }} ({{ comp_net.name }})

#### 6.{{ loop.index }}.1. Participation and Submissions

A total of {{ comp_net.total_labs }} laboratories submitted results for the {{ comp_code }} component.

- {{ comp_net.total_fasta }} submitted consensus genome sequences (.fasta), where applicable.
- {{ comp_net.total_vcf }} submitted variant call files (.vcf), where applicable.
- The metadata template completeness for {{ comp_code }} submissions had a median of {{ pct(comp_net.metadata_completeness_median) }}.

#### 6.{{ loop.index }}.2. Consensus Genome Reconstruction Performance

Consensus sequences were evaluated against the corresponding curated gold standard references for each sample in the {{ comp_code }} component.

Overall, {{ comp_code }} showed a median genome identity of {{ pct(comp_net.consensus.median_identity, 2) }}, with a median of {{ comp_net.consensus.median_discrepancies }} nucleotide discrepancies per sample (range: {{ comp_net.consensus.min_discrepancies }}–{{ comp_net.consensus.max_discrepancies }}).

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level consensus reconstruction metrics per sample for {{ comp_code }}.**

| Sample ID | Median genome identity (%) | Median discrepancies | Discrepancies min-max |
|---|---:|---:|---:|
{% for s in comp_net.consensus.samples %}
| {{ s.collecting_lab_sample_id }} | {{ "%.2f"|format(s.median_identity_pct) }} | {{ s.median_discrepancies }} | {{ s.min }} – {{ s.max }} |
{% endfor %}

Figure {{ fig_counter.value + 1 }} presents the distribution of nucleotide discrepancies per sample across participating laboratories for {{ comp_code }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.consensus.fig_discrepancies_boxplot_by_sample,
  "Consensus discrepancies per sample for " ~ comp_code ~ " relative to the curated gold standard."
) }}


**Figure {{ fig_counter.value }}. Distribution of consensus discrepancies per sample for {{ comp_code }}.** Boxplots represent the number of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

Considering discrepancy type composition aggregated by sample for {{ comp_code }}:

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level consensus discrepancy types per sample for {{ comp_code }}.**

| Sample ID | Median of Wrong nucleotide | Median Ambiguity instead of nucleotide | Median Nucleotide instead of ambiguity | Median Stretch of Ns instead of nucleotide stretch | Median Nucleotide stretch instead of stretch of Ns | Median Insertion relative to gold standard | Median Deletion relative to gold standard |
|---|---:|---:|---:|---:|---:|---:|---:|
{% for s in comp_net.consensus.samples %}
| {{ s.collecting_lab_sample_id }} | {{ s.wrong_nt }} | {{ s.ambiguity2nt }} | {{ s.nt2ambiguity }} | {{ s.ns2nt }} | {{ s.nt2ns }} | {{ s.insertions }} | {{ s.deletions }}
{% endfor %}

Figure {{ fig_counter.value + 1 }} presents the distribution of nucleotide discrepancy types per sample across participating laboratories for {{ comp_code }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.consensus.fig_discrepancies_stacked_by_sample,
  "Consensus discrepancy types per sample for " ~ comp_code ~ " relative to the curated gold standard."
) }}


**Figure {{ fig_counter.value }}. Distribution of consensus discrepancies per sample for {{ comp_code }}.** Stacked bars represent the number and type of nucleotide discrepancies relative to the curated gold standard across participating laboratories for each sample.

Discrepancy type composition aggregated across all submitted consensus sequences for {{ comp_code }}:

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level discrepancy composition by type for {{ comp_code }}.**

| Discrepancy type | Network median per sample | Min-max occurencies |
|---|---:|---:|
| Incorrect nucleotide | {{ comp_net.consensus.discrepancy_breakdown.wrong_nt.median }} | {{ comp_net.consensus.discrepancy_breakdown.wrong_nt.min }}–{{ comp_net.consensus.discrepancy_breakdown.wrong_nt.max }} |
| Ambiguity instead of nucleotide | {{ comp_net.consensus.discrepancy_breakdown.ambiguity2nt.median }} | {{ comp_net.consensus.discrepancy_breakdown.ambiguity2nt.min }}–{{ comp_net.consensus.discrepancy_breakdown.ambiguity2nt.max }} |
| Nucleotide instead of ambiguity | {{ comp_net.consensus.discrepancy_breakdown.nt2ambiguity.median }} | {{ comp_net.consensus.discrepancy_breakdown.nt2ambiguity.min }}–{{ comp_net.consensus.discrepancy_breakdown.nt2ambiguity.max }} |
| Stretch of Ns instead of nucleotide | {{ comp_net.consensus.discrepancy_breakdown.ns2nt.median }} | {{ comp_net.consensus.discrepancy_breakdown.ns2nt.min }}–{{ comp_net.consensus.discrepancy_breakdown.ns2nt.max }} |
| Nucleotide stretch instead of stretch of Ns| {{ comp_net.consensus.discrepancy_breakdown.nt2ns.median }} | {{ comp_net.consensus.discrepancy_breakdown.nt2ns.min }}–{{ comp_net.consensus.discrepancy_breakdown.nt2ns.max }} |
| Insertion relative to gold standard | {{ comp_net.consensus.discrepancy_breakdown.insertions.median }} | {{ comp_net.consensus.discrepancy_breakdown.insertions.min }}–{{ comp_net.consensus.discrepancy_breakdown.insertions.max }} |
| Deletion relative to gold standard | {{ comp_net.consensus.discrepancy_breakdown.deletions.median }} | {{ comp_net.consensus.discrepancy_breakdown.deletions.min }}–{{ comp_net.consensus.discrepancy_breakdown.deletions.max }} |

The dominant discrepancy pattern observed in {{ comp_code }} was {{ comp_net.consensus.dominant_discrepancy_pattern }}.

Figure {{ fig_counter.value + 1 }} summarises the contribution of each discrepancy category observed in {{ comp_code }} relative to the curated gold standard.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.consensus.fig_discrepancy_type_boxplot,
  "Composition of consensus discrepancy types for " ~ comp_code ~ " relative to the curated gold standard."
) }}


**Figure {{ fig_counter.value }}. Composition of consensus discrepancy types relative to the curated gold standard for {{ comp_code }}.** Boxplots represent aggregated discrepancies across all submitted consensus sequences, stratified by discrepancy category. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

#### 6.{{ loop.index }}.3. Variant Detection Accuracy

{% if comp_net.variant.median_discrepancies is defined %}

Variant call files (.vcf) submitted for the {{ comp_code }} component were compared against the curated reference variant set corresponding to each sample in the {{ comp_code }} component.

Overall, {{ comp_code }} showed a median of {{ comp_net.variant.median_discrepancies }} variant discrepancies per sample (range: {{ comp_net.variant.min_discrepancies }}–{{ comp_net.variant.max_discrepancies }}). The component also showed a median of {{ comp_net.variant.median_successful_hits if comp_net.variant.median_successful_hits is not none else "NA" }} successful hits per sample, with median number of variants with an allele frequency higher than 75% of {{ comp_net.variant.median_variants_in_consensus if comp_net.variant.median_variants_in_consensus is not none else "NA" }} in the metadata and {{ comp_net.variant.median_variants_in_consensus_vcf if comp_net.variant.median_variants_in_consensus_vcf is not none else "NA" }} in the submitted VCF files. Tables {{ table_counter.value + 1 }} and {{ table_counter.value + 2 }} summarise the descriptive reporting metrics and the qualitative discrepancy profile observed across samples in {{ comp_code }}. Table {{ table_counter.value +1 }} summarises the descriptive reporting metrics for {{ comp_code }}, including successful hits, the number of high-frequency variants reported in the metadata and VCF files, and the concordance between both representations for all variants and effect-annotated variants. Table {{ table_counter.value + 2 }} summarises, for each sample in {{ comp_code }}, the number of successful reference-variant hits together with the qualitative discrepancy profile, including wrong nucleotide calls, insertions, deletions, missing expected variants, and de novo variants.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level SARS-CoV-2 variant reporting metrics per sample for {{ comp_code }}.**

| Sample ID | Successful hits | Variants >75% AF in metadata | Variants >75% AF in VCF | Variants with effect in metadata | Variants with effect in VCF | Discrepancies metadata vs VCF | Effect discrepancies metadata vs VCF |
|---|---:|---:|---:|---:|---:|---:|---:|
{% for s in comp_net.variant.samples %}
| {{ s.collecting_lab_sample_id }} | {{ s.median_successful_hits if s.median_successful_hits is not none else "NA" }} | {{ s.variants_in_consensus.median if s.variants_in_consensus and s.variants_in_consensus.median is not none else "NA" }} | {{ s.variants_in_consensus_vcf.median if s.variants_in_consensus_vcf and s.variants_in_consensus_vcf.median is not none else "NA" }} | {{ s.variants_with_effect.median if s.variants_with_effect and s.variants_with_effect.median is not none else "NA" }} | {{ s.variants_with_effect_vcf.median if s.variants_with_effect_vcf and s.variants_with_effect_vcf.median is not none else "NA" }} | {{ s.discrepancies_in_reported_variants.median if s.discrepancies_in_reported_variants and s.discrepancies_in_reported_variants.median is not none else "NA" }} | {{ s.discrepancies_in_reported_variants_effect.median if s.discrepancies_in_reported_variants_effect and s.discrepancies_in_reported_variants_effect.median is not none else "NA" }} |
{% endfor %}

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level SARS-CoV-2 variant calling profile per sample for {{ comp_code }}.**

| Sample ID | Successful hits | Median discrepancies | Discrepancies min-max | Wrong nucleotide | Insertions | Deletions | Missing | De novo |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
{% for s in comp_net.variant.samples %}
| {{ s.collecting_lab_sample_id }} | {{ s.median_successful_hits if s.median_successful_hits is not none else "NA" }} | {{ s.median_discrepancies if s.median_discrepancies is not none else "NA" }} | {{ s.min if s.min is not none else "NA" }} – {{ s.max if s.max is not none else "NA" }} | {{ s.wrong_nt if s.wrong_nt is not none else "NA" }} | {{ s.insertions if s.insertions is not none else "NA" }} | {{ s.deletions if s.deletions is not none else "NA" }} | {{ s.missing if s.missing is not none else "NA" }} | {{ s.denovo if s.denovo is not none else "NA" }} |
{% endfor %}

Figure {{ fig_counter.value + 1 }} presents the distribution of nucleotide discrepancies per sample across participating laboratories for {{ comp_code }}.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.variant.fig_discrepancies_stacked_by_sample,
  "Variant discrepancies per sample for " ~ comp_code ~ " relative to the curated gold standard."
) }}

**Figure {{ fig_counter.value }}. Distribution of variant discrepancies per sample for {{ comp_code }}.** Stacked bars represent the number of nucleotide discrepancies and discrepancy types relative to the curated gold standard across participating laboratories for each sample.

Discrepancy type composition (aggregated across all submitted variant calls for {{ comp_code }}):

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level discrepancy composition by type for {{ comp_code }}.** The discrepancy-type columns correspond to the median count per sample across participating laboratories.

| Discrepancy type | Network median per sample | Network min-max per sample |
|---|---:|---:|
| Incorrect nucleotide | {{ comp_net.variant.discrepancy_breakdown.wrong_nt.median }} | {{ comp_net.variant.discrepancy_breakdown.wrong_nt.min }}–{{ comp_net.variant.discrepancy_breakdown.wrong_nt.max }} |
| Insertion relative to gold standard | {{ comp_net.variant.discrepancy_breakdown.insertions.median }} | {{ comp_net.variant.discrepancy_breakdown.insertions.min }}–{{ comp_net.variant.discrepancy_breakdown.insertions.max }} |
| Deletions relative to gold standard | {{ comp_net.variant.discrepancy_breakdown.deletions.median }} | {{ comp_net.variant.discrepancy_breakdown.deletions.min }}–{{ comp_net.variant.discrepancy_breakdown.deletions.max }} |
| Missing expected variants | {{ comp_net.variant.discrepancy_breakdown.missing.median }} | {{ comp_net.variant.discrepancy_breakdown.missing.min }}–{{ comp_net.variant.discrepancy_breakdown.missing.max }} |
| De novo variants | {{ comp_net.variant.discrepancy_breakdown.denovo.median }} | {{ comp_net.variant.discrepancy_breakdown.denovo.min }}–{{ comp_net.variant.discrepancy_breakdown.denovo.max }} |

The dominant discrepancy pattern observed in {{ comp_code }} was {{ comp_net.variant.dominant_discrepancy_pattern }}.

Figure {{ fig_counter.value + 1 }} summarises the contribution of each discrepancy category observed in {{ comp_code }} relative to the curated gold standard.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.variant.fig_discrepancy_type_boxplot,
  "Composition of variant discrepancy types for " ~ comp_code ~ " relative to the curated gold standard."
) }}


**Figure {{ fig_counter.value }}. Composition of variant discrepancy types relative to the curated gold standard for {{ comp_code }}.** Boxplots represent aggregated discrepancies across all submitted variant calls, stratified by discrepancy category (incorrect nucleotide, excess ambiguous bases, and indels). The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

{% else %}

For the {{ comp_code }} component, variant evaluation focused on the agreement between variants with allele frequency above 75% reported in the metadata template and those represented in the submitted VCF files, together with the overall number of variants present in the VCF output. Overall, {{ comp_code }} showed a median of {{ comp_net.variant.median_variants_in_consensus if comp_net.variant.median_variants_in_consensus is not none else "NA" }} variants with allele frequency above 75% reported in the metadata template, compared with {{ comp_net.variant.median_variants_in_consensus_vcf if comp_net.variant.median_variants_in_consensus_vcf is not none else "NA" }} corresponding variants represented in the consensus-derived VCF. The median number of discrepancies between both representations was {{ comp_net.variant.median_discrepancies_in_reported_variants if comp_net.variant.median_discrepancies_in_reported_variants is not none else "NA" }}, while the median total number of variants present in the submitted VCF files was {{ comp_net.variant.median_variants_in_vcf if comp_net.variant.median_variants_in_vcf is not none else "NA" }}.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level influenza variant reporting metrics per sample for {{ comp_code }}.**

| Sample ID | Variants >75% AF in metadata | Variants >75% AF in consensus VCF | Discrepancies between metadata and VCF | Total variants in VCF |
|---|---:|---:|---:|---:|
{% for s in comp_net.variant.samples %}
| {{ s.collecting_lab_sample_id }} | {{ s.variants_in_consensus.median if s.variants_in_consensus and s.variants_in_consensus.median is not none else "NA" }} | {{ s.variants_in_consensus_vcf.median if s.variants_in_consensus_vcf and s.variants_in_consensus_vcf.median is not none else "NA" }} | {{ s.discrepancies_in_reported_variants.median if s.discrepancies_in_reported_variants and s.discrepancies_in_reported_variants.median is not none else "NA" }} | {{ s.variants_in_vcf.median if s.variants_in_vcf and s.variants_in_vcf.median is not none else "NA" }} |
{% endfor %}

Table {{ table_counter.value }} summarises, for each sample in {{ comp_code }}, the median number of high-frequency variants reported in the metadata template, the corresponding number represented in the submitted VCF, the discrepancies between both representations, and the total number of variants observed in the submitted VCF files.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Aggregated influenza variant reporting metrics for {{ comp_code }}.**

| Metric | Network median | Network min-max |
|---|---:|---:|
| Variants >75% AF in metadata | {{ comp_net.variant.median_variants_in_consensus if comp_net.variant.median_variants_in_consensus is not none else "NA" }} | {{ comp_net.variant.min_variants_in_consensus if comp_net.variant.min_variants_in_consensus is not none else "NA" }}–{{ comp_net.variant.max_variants_in_consensus if comp_net.variant.max_variants_in_consensus is not none else "NA" }} |
| Variants >75% AF in consensus VCF | {{ comp_net.variant.median_variants_in_consensus_vcf if comp_net.variant.median_variants_in_consensus_vcf is not none else "NA" }} | {{ comp_net.variant.min_variants_in_consensus_vcf if comp_net.variant.min_variants_in_consensus_vcf is not none else "NA" }}–{{ comp_net.variant.max_variants_in_consensus_vcf if comp_net.variant.max_variants_in_consensus_vcf is not none else "NA" }} |
| Discrepancies between metadata and VCF | {{ comp_net.variant.median_discrepancies_in_reported_variants if comp_net.variant.median_discrepancies_in_reported_variants is not none else "NA" }} | {{ comp_net.variant.min_discrepancies_in_reported_variants if comp_net.variant.min_discrepancies_in_reported_variants is not none else "NA" }}–{{ comp_net.variant.max_discrepancies_in_reported_variants if comp_net.variant.max_discrepancies_in_reported_variants is not none else "NA" }} |
| Total variants in VCF | {{ comp_net.variant.median_variants_in_vcf if comp_net.variant.median_variants_in_vcf is not none else "NA" }} | {{ comp_net.variant.min_variants_in_vcf if comp_net.variant.min_variants_in_vcf is not none else "NA" }}–{{ comp_net.variant.max_variants_in_vcf if comp_net.variant.max_variants_in_vcf is not none else "NA" }} |

Figure {{ fig_counter.value + 1 }} summarises influenza variant reporting patterns across samples in {{ comp_code }}, showing the distribution of laboratory-level observations for agreement between metadata-reported and VCF-derived high-frequency variants and for the overall number of variants observed in the submitted VCF files.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.variant.fig_reporting_summary_by_sample,
  "Influenza variant reporting summary by sample for " ~ comp_code ~ "."
) }}


**Figure {{ fig_counter.value }}. Influenza variant reporting summary by sample for {{ comp_code }}.** Panel A shows, for each sample, the distribution across participating laboratories of the number of variants with allele frequency above 75% reported in the metadata template, the corresponding number represented in the consensus-derived VCF, and the discrepancies between both representations. Panel B shows the distribution across participating laboratories of the total number of variants present in the submitted VCF files for each sample. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range within the plotted scale, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

#### 6.{{ loop.index }}.4. Lineage, Subtype and Clade Assignment

Lineage, subtype and clade assignments submitted for the {{ comp_code }} component were evaluated for concordance with the curated gold standard classifications.

Across all participating laboratories:

- Lineage/Subtype matches (lineage/type was correct): {{ pct(comp_net.typing.lineage_hit_pct) }}.
- Clade matches (clade was correct): {{ pct(comp_net.typing.clade_hit_pct) }}

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Network-level classification outcomes per sample for {{ comp_code }}.**

| Sample ID | Lineage/Subtype matches (%) | Clade matches (%)|
|---|---:|---:|
{% for s in comp_net.typing.samples %}
| {{ s.collecting_lab_sample_id }} | {{ "%.2f"|format(s.lineage_hit_pct) }} | {{ "%.2f"|format(s.clade_hit_pct) }} |
{% endfor %}

Table {{ table_counter.value }} summarises the sample-level lineage/subtype and clade concordance rates for {{ comp_code }}. Figure {{ fig_counter.value + 1 }} presents the distribution of classification outcomes per sample across participating laboratories.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.typing.fig_stacked_bar_by_sample,
  "Classification outcome distribution per sample for " ~ comp_code ~ "."
) }}


**Figure {{ fig_counter.value }}. Classification outcome distribution per sample for {{ comp_code }}.** Panel A shows the proportion of lineage/subtype assignment Match, Discrepancy, and Not provided outcomes across participating laboratories for each sample. Panel B shows the corresponding proportions for clade assignments. Percentages are calculated over all participating laboratories in the component, so the Not provided segment captures samples for which lineage/subtype or clade information was not reported.

#### 6.{{ loop.index }}.5. Sample Quality Control Assessment

Laboratory-reported sample QC evaluations (Pass/Fail) for the {{ comp_code }} component were compared against the predefined gold standard QC status for each sample. Concordance was assessed as a binary outcome:

- Match: reported QC status equals the gold standard
- Discrepancy: reported QC status differs from the gold standard

Overall, QC concordance for {{ comp_code }} was {{ pct(comp_net.qc.match_rate_pct) }}, corresponding to {{ comp_net.qc.matches }} Matches and {{ comp_net.qc.discrepancies }} Discrepancies across {{ comp_net.qc.total_evaluations }} evaluated QC decisions.

{% set table_counter.value = table_counter.value + 1 %}
**_Table {{ table_counter.value }}_. Sample-level QC concordance for {{ comp_code }}.**

| Sample ID | Gold standard QC | % Match | # Matches | # Discrepancies | Total evaluations |
|---|---:|---:|---:|---:|---:|
{% for s in comp_net.qc.samples %}
| {{ s.collecting_lab_sample_id }} | {{ s.gold_standard_qc }} | {{ pct(s.match_rate_pct) }} | {{ s.matches }} | {{ s.discrepancies }} | {{ s.total_evaluations }} |
{% endfor %}

{% set fig_counter.value = fig_counter.value + 1 %}
Table {{ table_counter.value }} summarises the proportion of laboratories correctly classifying QC status for each sample relative to the gold standard definition, and Figure {{ fig_counter.value }} presents the corresponding sample-level distribution of Match and Discrepancy outcomes within {{ comp_code }}.

{{ render_figure(
comp_net.qc.fig_qc_match_by_sample,
"Sample-level QC concordance for " ~ comp_code ~ " (Match vs Discrepancy relative to the gold standard)."
) }}


**_Figure {{ fig_counter.value }}_. Sample-level QC concordance for {{ comp_code }} relative to the gold standard.** Bars represent the proportion of Match vs Discrepancy outcomes per sample across participating laboratories. Higher discrepancy rates indicate samples for which laboratories more frequently diverged from the predefined QC status.

#### 6.{{ loop.index }}.6. Pipeline Benchmarking and Comparative Performance

{% if comp_net.benchmarking.bioinformatics_protocol %}
##### Bioinformatics protocol

Based on metadata submissions, {{ comp_net.benchmarking.bioinformatics_protocol.total_number }} distinct bioinformatics protocols were reported for the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Performance summary of declared bioinformatics protocols for {{ comp_code }}.**

| Bioinformatics protocol | Version | N labs | Median genome identity (%) | Median discrepancies | Median metadata completeness (%) | Clade concordance (%) | Lineage/type concordance (%) |
|---|---:|---:|---:|---:|---:|---:|---:|
{% for p in comp_net.benchmarking.bioinformatics_protocol.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ "%.2f"|format(p.median_identity_pct) }} | {{ p.median_discrepancies }} | {{ "%.1f"|format(p.median_metadata_completeness_pct) }} | {{ "%.1f"|format(p.clade_hit_pct) }} | {{ "%.1f"|format(p.lineage_hit_pct) }} |
{% endfor %}

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value + 1 }} summarises the distribution of key performance metrics stratified by declared pipeline configuration.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.benchmarking.bioinformatics_protocol.fig_metric_boxplots,
  "Distribution of performance metrics by pipeline configuration for " ~ comp_code ~ "."
) }}


**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared pipeline configuration for {{ comp_code }}.** Multi-panel boxplots summarise sample-level performance stratified by bioinformatics protocols. Panel A displays genome identity (%), Panel B discrepancy counts, Panel C metadata completeness (%), and Panel D exact classification concordance (%). Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.dehosting %}
##### De-hosting software

Based on metadata submissions, {{ comp_net.benchmarking.dehosting.total_number }} distinct de-hosting softwares were reported for the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Performance summary of declared de-hosting software for {{ comp_code }}.**

| De-hosting software | Version | N labs | % Host reads |
|---|---:|---:|---:|
{% for p in comp_net.benchmarking.dehosting.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.per_reads_host }} |
{% endfor %}

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value + 1 }} summarises the percentage of host reads metric stratified by declared dehosting sfotaware version.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.benchmarking.dehosting.fig_metric_boxplots,
  "Distribution of percentage of host reads metrics by dehosting software version for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of percentage of host reads by declared dehosting software version for {{ comp_code }}.** Boxplots summarise sample-level percentage of host reads stratified by dehosting software version. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each version, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.preprocessing %}
##### Preprocessing software

Based on metadata submissions, {{ comp_net.benchmarking.preprocessing.total_number }} distinct pre-processing software configurations were reported for the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Performance summary of declared pre-processing software configurations for {{ comp_code }}.** The configuration column represents the most frequently reported parameter string among laboratories declaring that software and version.

| Pre-processing software | Version | N labs | Most common configuration | Number of reads sequenced | Reads passing filters |
|---|---:|---:|---:|---:|---:|
{% for p in comp_net.benchmarking.preprocessing.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.params }} | {{ p.number_of_reads_sequenced }} | {{ p.pass_reads }} |
{% endfor %}

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value + 1 }} summarises the distribution of key performance metrics stratified by declared pre-processing software configuration.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.benchmarking.preprocessing.fig_metric_boxplots,
  "Distribution of performance metrics by pre-processing software configuration for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared pre-processing software configuration for {{ comp_code }}.** Multi-panel boxplots summarise sample-level performance stratified by pre-processing software. Panel A displays Number of reads sequenced and Panel B Reads passing filters. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.mapping %}
##### Mapping software

Based on metadata submissions, {{ comp_net.benchmarking.mapping.total_number }} distinct mapping software configurations were reported for the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Performance summary of declared mapping software configurations for {{ comp_code }}.**

| Mapping software | Version | N labs | Most common configuration | Depth of coverage threshold | % Reads virus |
|---|---:|---:|---:|---:|---:|
{% for p in comp_net.benchmarking.mapping.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.params }} | {{ p.depth_of_coverage_threshold if p.depth_of_coverage_threshold is not none else "N/A" }} | {{ p.per_reads_virus if p.per_reads_virus is not none else "N/A" }} |
{% endfor %}

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value + 1 }} summarises the distribution of key performance metrics stratified by declared mapping software configuration.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.benchmarking.mapping.fig_metric_boxplots,
  "Distribution of performance metrics by mapping software configuration for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared mapping software configuration for {{ comp_code }}.** Boxplots summarise sample-level performance stratified by mapping software. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.assembly %}
##### Assembly software

Based on metadata submissions, {{ comp_net.benchmarking.assembly.total_number }} distinct assembly software configurations were reported for the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Performance summary of declared assembly software configurations for {{ comp_code }}.**

| Assembly software | Version | N labs | Most common configuration | Consnsus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
{% for p in comp_net.benchmarking.assembly.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.params }} | {{ p.consensus_genome_length }} | {{ p.median_identity_pct }} |  {{ p.median_discrepancies }} |
{% endfor %}

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value + 1 }} summarises the distribution of key performance metrics stratified by declared assembly software configuration.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.benchmarking.assembly.fig_metric_boxplots,
  "Distribution of performance metrics by assembly software configuration for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared assembly software configuration for {{ comp_code }}.** Multi-panel boxplots summarise sample-level performance stratified by assembly software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.consensus_software %}
##### Consensus software

Based on metadata submissions, {{ comp_net.benchmarking.consensus_software.total_number }} distinct consensus software configurations were reported for the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Performance summary of declared consensus software configurations for {{ comp_code }}.**

| Consensus software | Version | N labs | Most common configuration | Consnsus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
{% for p in comp_net.benchmarking.consensus_software.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.params }} | {{ p.consensus_genome_length }} | {{ p.median_identity_pct }} |  {{ p.median_discrepancies }} |
{% endfor %}

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value + 1 }} summarises the distribution of key performance metrics stratified by declared consensus software configuration.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.benchmarking.consensus_software.fig_metric_boxplots,
  "Distribution of performance metrics by consensus software configuration for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared consensus software configuration for {{ comp_code }}.** Multi-panel boxplots summarise sample-level performance stratified by consensus software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.variant_calling %}
##### Variant calling software

Based on metadata submissions, {{ comp_net.benchmarking.variant_calling.total_number }} distinct variant calling software configurations were reported for the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Performance summary of declared variant calling software configurations for {{ comp_code }}.**

{% if comp_code[:3] == "FLU" %}
| Variant calling software | Version | N labs | Most common configuration | Median high + low freq (%) | Median high freq only (%) | Median low freq only (%) | Median variants (AF >=75%) | Median variants in VCF (AF >=75%) | Median variants with effect | Median metadata-VCF discrepancies | Median total variants in VCF |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{% for p in comp_net.benchmarking.variant_calling.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.params }} | {{ p.high_and_low_freq_pct }} | {{ p.high_freq_only_pct }} | {{ p.low_freq_only_pct }} | {{ p.number_of_variants_in_consensus }} | {{ p.number_of_variants_in_consensus_vcf }} | {{ p.number_of_variants_with_effect }} | {{ p.discrepancies_in_reported_variants }} | {{ p.number_of_variants_in_vcf }} |
{% endfor %}
{% else %}
| Variant calling software | Version | N labs | Most common configuration | Median high + low freq (%) | Median high freq only (%) | Median low freq only (%) | Median variants (AF >=75%) | Median variants in VCF (AF >=75%) | Median variants with effect | Median variants with effect in VCF | Median metadata-VCF discrepancies | Median effect discrepancies | Median successful hits | Median total discrepancies |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{% for p in comp_net.benchmarking.variant_calling.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.params }} | {{ p.high_and_low_freq_pct }} | {{ p.high_freq_only_pct }} | {{ p.low_freq_only_pct }} | {{ p.number_of_variants_in_consensus }} | {{ p.number_of_variants_in_consensus_vcf }} | {{ p.number_of_variants_with_effect }} | {{ p.number_of_variants_with_effect_vcf }} | {{ p.discrepancies_in_reported_variants }} | {{ p.discrepancies_in_reported_variants_effect }} | {{ p.successful_hits }} | {{ p.total_discrepancies }} |
{% endfor %}
{% endif %}

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value + 1 }} summarises the distribution of key performance metrics stratified by declared variant calling software configuration.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.benchmarking.variant_calling.fig_metric_boxplots,
  "Distribution of performance metrics by variant calling software configuration for " ~ comp_code ~ "."
) }}

{% if comp_code[:3] == "FLU" %}
**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared variant calling software configuration for {{ comp_code }}.** Panel A is a stacked bar chart showing the number of evaluable samples assigned to each allele frequency reporting pattern for each software configuration. Boxplot Panel B displays the number of reported variants with AF >=75%, Panel C the number of variants with AF >=75% in the submitted VCF, Panel D the number of variants with effect, Panel E metadata-VCF discrepancies, and Panel F the total number of variants present in the submitted VCF files. Only panels with evaluable data are shown. In the boxplots, the central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.
{% else %}
**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared variant calling software configuration for {{ comp_code }}.** Panel A is a stacked bar chart showing the number of evaluable samples assigned to each allele frequency reporting pattern for each software configuration. Boxplot Panel B displays discrepancies in reported variants with AF >=75% in the submitted VCF, Panel C discrepancies in reported variants with effect, Panel D successful hits, and Panel E total discrepancies. Only panels with evaluable data are shown. In the boxplots, the central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.
{% endif %}

{% endif %}

{% if comp_net.benchmarking.clade_assignment %}
##### Clade Assignment Software

Based on metadata submissions, {{ comp_net.benchmarking.clade_assignment.total_number }} distinct clade assignment software configurations were reported for the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Performance summary of declared clade assignment software configurations for {{ comp_code }}.**

| Clade assignment software | Version | N labs | Database version | % of clade match | % of clade discrepancy |
|---|---:|---:|---:|---:|---:|
{% for p in comp_net.benchmarking.clade_assignment.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.database_version }} | {{ p.clade_hit_pct }} | {{ p.clade_discordance_pct }} |
{% endfor %}

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value + 1 }} summarises the distribution of key performance metrics stratified by declared clade assignment software configuration.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.benchmarking.clade_assignment.fig_metric_boxplots,
  "Distribution of performance metrics by clade assignment software configuration for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared clade assignment software configuration for {{ comp_code }}.** Multi-panel boxplots summarise sample-level performance stratified by clade assignment software. Panel A displays the % of clade matches and Panel B the % of clade discrepancies. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.lineage_assignment %}
##### Lineage Assignment Software Name

Based on metadata submissions, {{ comp_net.benchmarking.lineage_assignment.total_number }} distinct lineage assignment software configurations were reported for the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Performance summary of declared lineage assignment software configurations for {{ comp_code }}.**

| Lineage Assignment software | Version | N labs | Database version | % of lineage match | % of lineage discrepancy |
|---|---:|---:|---:|---:|---:|
{% for p in comp_net.benchmarking.lineage_assignment.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.database_version }} | {{ p.lineage_hit_pct }} | {{ p.lineage_discordance_pct }} |
{% endfor %}

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value + 1 }} summarises the distribution of key performance metrics stratified by declared lineage assignment software configuration.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.benchmarking.lineage_assignment.fig_metric_boxplots,
  "Distribution of performance metrics by lineage assignment software configuration for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared lineage assignment software configuration for {{ comp_code }}.** Multi-panel boxplots summarise sample-level performance stratified by lineage assignment software. Panel A displays the % of lineage matches and Panel B the % of lineage discrepancies. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.type_assignment %}
##### Type Assignment Software Name

Based on metadata submissions, {{ comp_net.benchmarking.type_assignment.total_number }} distinct type assignment software configurations were reported for the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Performance summary of declared type assignment software configurations for {{ comp_code }}.**

| Type Assignment software | Version | N labs | Database version | % of type match | % of type discrepancy |
|---|---:|---:|---:|---:|---:|
{% for p in comp_net.benchmarking.type_assignment.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.database_version }} | {{ p.type_hit_pct }} | {{ p.type_discordance_pct }} |
{% endfor %}

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value + 1 }} summarises the distribution of key performance metrics stratified by declared xxx software configuration.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.benchmarking.type_assignment.fig_metric_boxplots,
  "Distribution of performance metrics by type assignment software configuration for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared type assignment software configuration for {{ comp_code }}.** Multi-panel boxplots summarise sample-level performance stratified by type assignment software. Panel A displays the % of type matches and Panel B the % of type discrepancies. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% if comp_net.benchmarking.subtype_assignment %}

##### Subtype Assignment Software Name

Based on metadata submissions, {{ comp_net.benchmarking.subtype_assignment.total_number }} distinct subtype assignment software configurations were reported for the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Performance summary of declared subtype assignment software configurations for {{ comp_code }}.**

| Subtype Assignment software | Version | N labs | Database version | % of subtype match | % of subtype discrepancy |
|---|---:|---:|---:|---:|---:|
{% for p in comp_net.benchmarking.subtype_assignment.softwares %}
| {{ p.name }} | {{ p.version }} | {{ p.n_labs }} | {{ p.database_version }} | {{ p.subtype_hit_pct }} | {{ p.subtype_discordance_pct }} |
{% endfor %}

{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value + 1 }} summarises the distribution of key performance metrics stratified by declared subtype assignment software configuration.

{% set fig_counter.value = fig_counter.value + 1 %}
{{ render_figure(
  comp_net.benchmarking.subtype_assignment.fig_metric_boxplots,
  "Distribution of performance metrics by subtype assignment software configuration for " ~ comp_code ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of performance metrics by declared subtype assignment software configuration for {{ comp_code }}.** Multi-panel boxplots summarise sample-level performance stratified by subtype assignment software. Panel A displays the % of subtype matches and Panel B the % of subtype discrepancies. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

{% endif %}

{% endfor %}

## 7. Discussion

The 2026 RELECOV Dry-Lab EQA provides the first network-wide dry-lab assessment focused specifically on bioinformatic performance across consensus reconstruction, variant reporting, classification, metadata reporting, and QC interpretation. By combining ECDC datasets with in-silico influenza material, the exercise captures both routine-use analytical behaviour and performance under heterogeneous reference and reporting conditions.

### 7.1. Consensus Genome Reconstruction

Consensus reconstruction results were strongest in the Illumina-based components overall, with a combined median genome identity of {{ pct(general.general_results.consensus.median_identity_illumina_pct, 2) }} compared with {{ pct(general.general_results.consensus.median_identity_nanopore_pct, 2) }} across the Nanopore-based components. At component level, SARS1 and SARS2 showed median identities of {{ pct(general.components.SARS1.consensus.median_identity_pct, 2) }} and {{ pct(general.components.SARS2.consensus.median_identity_pct, 2) }}, whereas FLU1 and FLU2 showed lower medians of {{ pct(general.components.FLU1.consensus.median_identity_pct, 2) }} and {{ pct(general.components.FLU2.consensus.median_identity_pct, 2) }}.

At the same time, the ranges observed across laboratories show that high medians did not eliminate outlier behaviour. In particular, the minimum identity values in SARS2 and FLU2 dropped to {{ pct(general.components.SARS2.consensus.identity_pct_min, 2) }} and {{ pct(general.components.FLU2.consensus.identity_pct_min, 2) }}, indicating that a subset of submissions diverged markedly from the curated gold standard.

The dominant discrepancy categories also differed by component. SARS1 was dominated by stretches of Ns in submitted consensuses where defined nucleotides were present in the gold standard (`ns2nt`), whereas SARS2 was dominated by defined nucleotides where stretches of Ns were present in the gold standard (`nt2ns`). These patterns are consistent with differences in masking behaviour and minimum coverage policies relative to the gold standard reconstruction criteria. In influenza, both FLU1 and FLU2 were dominated by deletions relative to the gold standard, suggesting that consensus generation parameters remain important contributors to inter-laboratory divergence.

### 7.2. Variant Detection and Reporting

For SARS-CoV-2, variant detection performance did not follow a simple platform ranking. The median number of discrepancies relative to the curated variant set was {{ general.general_results.sars_variants.median_discrepancy_illumina }} in the Illumina component and {{ general.general_results.sars_variants.median_discrepancy_nanopore }} in the Nanopore component. This indicates that platform effects were present, but that they interacted with sample composition, reporting choices, and software configuration rather than determining performance on their own.

The submitted metadata documented substantial diversity in variant reporting behaviour. Across SARS-CoV-2 submissions, {{ general.general_results.sars_variants.high_and_low_freq_pct }} of laboratories reported both high- and low-frequency variants, whereas {{ general.general_results.sars_variants.high_freq_only_pct }} reported high-frequency variants only. Influenza reporting was more heterogeneous: {{ general.general_results.influenza_variants.high_and_low_freq_pct }} of laboratories reported both high- and low-frequency variants, {{ general.general_results.influenza_variants.low_freq_only_pct }} reported low-frequency variants only, and {{ general.general_results.influenza_variants.high_freq_only_pct }} reported high-frequency variants only.

Influenza results especially highlight the consequences of heterogeneous structural reporting. The network-level median number of variants with AF >=75% reported in metadata was {{ general.general_results.influenza_variants.median_variants_in_consensus }}, whereas the corresponding median derived from submitted VCF files was {{ general.general_results.influenza_variants.median_variants_in_consensus_vcf }}. The median discrepancy between these two representations was {{ general.general_results.influenza_variants.median_discrepancies_in_reported_variants }}, and the total number of variants present in submitted VCF files ranged from {{ general.general_results.influenza_variants.min_variants_in_vcf }} to {{ general.general_results.influenza_variants.max_variants_in_vcf }}. Together, these values indicate that influenza variant outputs were not directly comparable under a single harmonised coordinate framework and that reporting conventions differed markedly across laboratories.

### 7.3. Classification and QC Interpretation

Classification performance was acceptable overall but clearly stronger for lineage/type assignment than for clade assignment. SARS-CoV-2 lineage concordance reached {{ pct(general.general_results.classification.sars_lineage_concordance_pct) }}, compared with {{ pct(general.general_results.classification.sars_clade_concordance_pct) }} for SARS-CoV-2 clade assignment. Influenza type/subtype concordance reached {{ pct(general.general_results.classification.influenza_type_concordance_pct) }}, compared with {{ pct(general.general_results.classification.influenza_clade_concordance_pct) }} for influenza clade assignment. In the SARS-CoV-2 submissions, this difference is consistent with metadata-level reporting problems in the clade field itself: among the clade assignments reviewed in the submitted JSON files, some were left empty and others contained values that matched the lineage assignment or had lineage-like syntax rather than a clade designation. This suggests that part of the excess clade discordance reflects field completion and nomenclature issues in addition to true analytical misclassification.

QC interpretation showed additional between-component differences. Only 9 of the 19 participating laboratories reported at least one sample-level QC assessment in their submitted metadata. Among evaluable QC decisions, network-wide concordance was {{ pct(general.qc.match_rate_pct) }}, and component-level concordance ranged from {{ pct(general.components.SARS1.qc.match_rate_pct) }} in SARS1 to {{ pct(general.components.FLU2.qc.match_rate_pct) }} in FLU2. This indicates that QC interpretation was not equally stable across all datasets, while also showing that many laboratories either did not apply or did not report a formal sample-level QC decision in the metadata template.

### 7.4. Workflow Diversity and Reporting Constraints

The metadata confirm that RELECOV laboratories currently use a diverse analytical landscape. A total of {{ general.metadata_completeness.total_workflows }} distinct workflows were identified across participating laboratories, together with {{ general.metadata_completeness.total_consensus_softwares }} distinct consensus tools or tool/version combinations, {{ general.metadata_completeness.total_variant_softwares }} distinct variant tools, {{ general.metadata_completeness.total_lineage_assignment_softwares }} SARS-CoV-2 lineage assignment tools, and {{ general.metadata_completeness.total_clade_assignment_softwares }} clade assignment tools.

This diversity is analytically valuable, but its interpretation is constrained by incomplete metadata reporting. Only {{ pct(general.metadata_completeness.software_version_pct) }} of software-version fields were completed, {{ pct(general.metadata_completeness.coverage_threshold_pct) }} of submitted samples specified a minimum coverage threshold, {{ pct(general.metadata_completeness.variant_calling_params_pct) }} reported variant calling parameters, and {{ pct(general.metadata_completeness.reference_genome_pct) }} reported a reference genome accession or identifier. For that reason, some plausible explanations for performance differences can only be discussed as contributing context rather than demonstrated causal effects.

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

{% if labdata %}
{% set lab_code = labdata.lab.lab_cod | default(labdata.lab.submitting_institution_id) %}
# 9. Individual Laboratory Technical Report

## Laboratory: {{ labdata.lab.laboratory_name }} ({{ labdata.lab.lab_cod }})

This section provides a detailed technical assessment of the analytical results submitted by **{{ labdata.lab.lab_cod }}** within the 2026 RELECOV Dry-Lab EQA. Performance metrics are benchmarked against curated gold standards and contextualised relative to aggregated network-wide performance distributions. Network medians and interquartile ranges are provided for comparative interpretation, without disclosure of other laboratories’ identities.

The purpose of this section is to support technical optimisation, parameter harmonisation, and alignment with the analytical standards defined within RELECOV 2.0.

Only files, metadata fields, and derived analytical metrics actually provided by the laboratory are displayed in this individual report. If a file was not submitted, or a metadata field was not provided, the corresponding table entries, panels, or figures are omitted for that laboratory.

## 9.1. Participation Overview

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
{% for d in labdata.metadata.primary_incompleteness_drivers %}
- {{ d }}
{% endfor %}
{% endif %}

{% for comp_code, comp in labdata.components.items() %}

# 9.{{ loop.index + 1 }}. {{ comp_code }} ({{ comp.display_name }})

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
{% for d in comp.metadata.primary_incompleteness_drivers %}
- {{ d }}
{% endfor %}
{% endif %}

## 9.{{ loop.index + 1 }}.1. Consensus Genome Reconstruction Performance

Consensus genome sequences (`.fasta`) submitted by **{{ labdata.lab.lab_cod }}** were compared against the curated gold standard reference for each sample included in the {{ comp_code }} component.

### Per-sample summary metrics

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Per-sample consensus reconstruction metrics for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

| Sample ID | {{ labdata.lab.lab_cod }} Genome identity (%) | Network Genome Identity Median | {{ labdata.lab.lab_cod }} Total discrepancies | Network total discrepancies median |
|---|---:|---:|---:|---:|
{% for collecting_lab_sample_id, s in comp.samples.items() -%}
| {{ collecting_lab_sample_id }} | {{ pct(s.consensus.genome_identity_pct, 4) }} | {{ general.components[comp_code].consensus.samples[collecting_lab_sample_id].median_identity_pct }} | {{ s.consensus.total_discrepancies }} | {{ general.components[comp_code].consensus.samples[collecting_lab_sample_id].median_discrepancies }} |
{% endfor %}

The metrics presented in Table {{ table_counter.value }} summarise overall sequence similarity and discrepancy burden relative to the curated gold standard reference for {{ labdata.lab.lab_cod }} compared to the Network's median.

### Discrepancy type breakdown per sample

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Discrepancy type breakdown per sample for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

| Sample ID | Total wrong nucleotides | Total ambiguity instead of nucleotide | Total nucleotide instead of ambiguity | Total stretch of Ns instead of nucleotide stretch | Total sucleotide stretch instead of stretch of Ns | Total insertion relative to gold standard | Total deletion relative to gold standard |
|---|---:|---:|---:|---:|---:|---:|---:|
{% for collecting_lab_sample_id, s in comp.samples.items() -%}
| {{ collecting_lab_sample_id }} | {{ s.consensus.discrepancy_breakdown.wrong_nt }} | {{ s.consensus.discrepancy_breakdown.ambiguity2nt }} | {{ s.consensus.discrepancy_breakdown.nt2ambiguity }} | {{ s.consensus.discrepancy_breakdown.ns2nt }} | {{ s.consensus.discrepancy_breakdown.nt2ns }} | {{ s.consensus.discrepancy_breakdown.insertions }} | {{ s.consensus.discrepancy_breakdown.deletions }} |
{% endfor %}

Table {{ table_counter.value }} provides a detailed characterisation of discrepancy categories contributing to the total differences observed for each sample.

{% set consensus_distribution_panel_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/consensus_distribution_panel.png" %}
{% if path_exists(consensus_distribution_panel_path) %}
{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} illustrates the distribution of consensus discrepancies and genome identity per sample across participating laboratories in the network, contextualising the results obtained by **{{ labdata.lab.lab_cod }}**.

{{ render_figure(
  consensus_distribution_panel_path,
  comp_code ~ ": distribution of consensus discrepancies and genome identity per sample across the network; black diamond indicates " ~ labdata.lab.lab_cod ~ "."
) }}

**Figure {{ fig_counter.value }}. Consensus reconstruction performance across participating laboratories ({{ comp_code }}).** Panel A shows the distribution of total consensus discrepancies per sample relative to the curated gold standard across the RELECOV network. Panel B shows the corresponding distribution of genome identity values per sample. In both panels, the central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. The black diamond corresponds to the results obtained by **{{ labdata.lab.lab_cod }}**.
{% endif %}

{% set consensus_breakdown_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/consensus_discrepancy_breakdown_by_sample.png" %}
{% if path_exists(consensus_breakdown_path) %}
{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} summarises the discrepancy profile reported by **{{ labdata.lab.lab_cod }}** across samples in the {{ comp_code }} component.

{{ render_figure(
  consensus_breakdown_path,
  comp_code ~ ": discrepancy type breakdown by sample for " ~ labdata.lab.lab_cod ~ "."
) }}

**Figure {{ fig_counter.value }}. Discrepancy type breakdown by sample for {{ labdata.lab.lab_cod }} ({{ comp_code }}).** Stacked bars show the contribution of each discrepancy category to the total consensus differences observed for each sample submitted by **{{ labdata.lab.lab_cod }}**.
{% endif %}

{% if comp.metadata.vcf_submitted >=1 %}

## 9.{{ loop.index + 1 }}.2. Variant Detection Performance

{% if comp_code in ["SARS1", "SARS2"] %}

For SARS-CoV-2, variant call files (`.vcf`) submitted by **{{ labdata.lab.lab_cod }}** were compared against the curated reference variant set for each sample included in the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Per-sample variant detection performance metrics for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

| Sample ID | Reporting mode | {{ labdata.lab.lab_cod }} total discrepancies | Network median total discrepancies | {{ labdata.lab.lab_cod }} successful hits | Network median successful hits | Wrong variants | Insertions | Deletions | Missing expected variants | De novo variants |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{% for collecting_lab_sample_id, s in comp.samples.items() -%}
{% set ns = (general.components[comp_code].variant.samples | selectattr("collecting_lab_sample_id","equalto",collecting_lab_sample_id) | list | first) %}
| {{ collecting_lab_sample_id }} | {{ "High + low freq" if s.variants.high_and_low_freq else ("High freq only" if s.variants.high_freq_only else ("Low freq only" if s.variants.low_freq_only else "NA")) }} | {{ s.variants.total_discrepancies if s.variants.total_discrepancies is not none else "NA" }} | {{ ns.median_discrepancies if ns else "NA" }} | {{ s.variants.successful_hits if s.variants.successful_hits is not none else "NA" }} | {{ ns.median_successful_hits if ns else "NA" }} | {{ s.variants.wrong_nt if s.variants.wrong_nt is not none else "NA" }} | {{ s.variants.insertions if s.variants.insertions is not none else "NA" }} | {{ s.variants.deletions if s.variants.deletions is not none else "NA" }} | {{ s.variants.missing if s.variants.missing is not none else "NA" }} | {{ s.variants.denovo if s.variants.denovo is not none else "NA" }} |
{% endfor %}

The metrics presented in Table {{ table_counter.value }} summarise per-sample variant detection accuracy relative to the curated reference variant set and benchmark the laboratory’s results against the network median for the same sample.

{% set variant_detection_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/variant_metadata_vs_vcf_distribution.png" %}
{% if path_exists(variant_detection_path) %}
{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} illustrates the distribution of variant detection performance metrics across participating laboratories in the network, contextualising the results obtained by **{{ labdata.lab.lab_cod }}**.

{{ render_figure(
  variant_detection_path,
  comp_code ~ ": distribution of variant detection metrics across the network; black diamond indicates " ~ labdata.lab.lab_cod ~ "."
) }}

**Figure {{ fig_counter.value }}. Variant detection performance across participating laboratories ({{ comp_code }}).** Panel A shows the distribution of total variant discrepancies per sample across the RELECOV network. Panel B shows the corresponding distribution of successful hits per sample. In both panels, the central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range across the network, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. The black diamond corresponds to the results obtained by **{{ labdata.lab.lab_cod }}**.
{% endif %}

{% endif %}

{% if comp_code in ["FLU1", "FLU2"] %}
For influenza components, evaluation focused on structural reporting metrics and concordance between metadata-reported and VCF-derived variant counts for each sample.
{% else %}
The laboratory-reported variant counts declared in the metadata were also compared against the values derived directly from the submitted VCF files for each sample.
{% endif %}

{% set table_counter.value = table_counter.value + 1 %}
{% if comp_code in ["SARS1", "SARS2"] %}
**Table {{ table_counter.value }}. Metadata-reported and VCF-derived variant metrics for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

| Sample ID | Reporting mode | Metadata: variants (AF >=75%) | VCF-derived variants (AF >=75%) | Metadata: variants with effect | VCF-derived variants with effect | Metadata-VCF discrepancies | Effect discrepancies |
|---|---|---:|---:|---:|---:|---:|---:|
{% for collecting_lab_sample_id, s in comp.samples.items() -%}
| {{ collecting_lab_sample_id }} | {{ "High + low freq" if s.variants.high_and_low_freq else ("High freq only" if s.variants.high_freq_only else ("Low freq only" if s.variants.low_freq_only else "NA")) }} | {{ s.variants.number_of_variants_in_consensus if s.variants and s.variants.number_of_variants_in_consensus is not none else "NA" }} | {{ s.variants.number_of_variants_in_consensus_vcf if s.variants and s.variants.number_of_variants_in_consensus_vcf is not none else "NA" }} | {{ s.variants.number_of_variants_with_effect if s.variants and s.variants.number_of_variants_with_effect is not none else "NA" }} | {{ s.variants.number_of_variants_with_effect_vcf if s.variants and s.variants.number_of_variants_with_effect_vcf is not none else "NA" }} | {{ s.variants.discrepancies_in_reported_variants if s.variants and s.variants.discrepancies_in_reported_variants is not none else "NA" }} | {{ s.variants.discrepancies_in_reported_variants_effect if s.variants and s.variants.discrepancies_in_reported_variants_effect is not none else "NA" }} |
{% endfor %}
{% else %}
**Table {{ table_counter.value }}. Metadata-reported and VCF-derived variant metrics for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

| Sample ID | Reporting mode | Metadata: variants (AF >=75%) | VCF-derived variants (AF >=75%) | Metadata: variants with effect | Metadata-VCF discrepancies | Total variants in VCF |
|---|---|---:|---:|---:|---:|---:|
{% for collecting_lab_sample_id, s in comp.samples.items() -%}
| {{ collecting_lab_sample_id }} | {{ "High + low freq" if s.variants.high_and_low_freq else ("High freq only" if s.variants.high_freq_only else ("Low freq only" if s.variants.low_freq_only else "NA")) }} | {{ s.variants.number_of_variants_in_consensus if s.variants and s.variants.number_of_variants_in_consensus is not none else "NA" }} | {{ s.variants.number_of_variants_in_consensus_vcf if s.variants and s.variants.number_of_variants_in_consensus_vcf is not none else "NA" }} | {{ s.variants.number_of_variants_with_effect if s.variants and s.variants.number_of_variants_with_effect is not none else "NA" }} | {{ s.variants.discrepancies_in_reported_variants if s.variants and s.variants.discrepancies_in_reported_variants is not none else "NA" }} | {{ s.variants.number_of_variants_in_vcf if s.variants and s.variants.number_of_variants_in_vcf is not none else "NA" }} |
{% endfor %}
{% endif %}
{% set variant_metrics_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/variant_metrics_distribution.png" %}
{% if path_exists(variant_metrics_path) %}
{% set fig_counter.value = fig_counter.value + 1 %}

{% if comp_code in ["SARS1", "SARS2"] %}
Figure {{ fig_counter.value }} illustrates the distribution of metadata-reported and VCF-derived variant metrics across participating laboratories in the network, contextualising the results obtained by **{{ labdata.lab.lab_cod }}**.
{% else %}
Figure {{ fig_counter.value }} illustrates the distribution of influenza-specific reporting metrics across participating laboratories in the network, contextualising the results obtained by **{{ labdata.lab.lab_cod }}**.
{% endif %}

{{ render_figure(
  variant_metrics_path,
  comp_code ~ ": distribution of variant reporting metrics across the network; black diamond indicates " ~ labdata.lab.lab_cod ~ "."
) }}

{% if comp_code in ["SARS1", "SARS2"] %}
**Figure {{ fig_counter.value }}. Metadata-reported and VCF-derived variant metrics across participating laboratories ({{ comp_code }}).** Panel A shows reported variants with AF >=75%, Panel B VCF-derived variants with AF >=75%, Panel C reported variants with effect, Panel D VCF-derived variants with effect, Panel E metadata-VCF discrepancies, and Panel F effect discrepancies across the RELECOV network. Only panels with evaluable data for **{{ labdata.lab.lab_cod }}** are shown. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. The black diamond corresponds to the results obtained by **{{ labdata.lab.lab_cod }}**.
{% else %}
**Figure {{ fig_counter.value }}. Influenza-specific variant reporting metrics across participating laboratories ({{ comp_code }}).** Panel A shows reported variants with AF >=75%, Panel B VCF-derived variants with AF >=75%, Panel C reported variants with effect, Panel D metadata-VCF discrepancies, and Panel E total variants present in the submitted VCF files across the RELECOV network. Only panels with evaluable data for **{{ labdata.lab.lab_cod }}** are shown. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. The black diamond corresponds to the results obtained by **{{ labdata.lab.lab_cod }}**.
{% endif %}
{% endif %}
{% endif %}

## 9.{{ loop.index + 1 }}.3. Lineage, Subtype and Clade Assignment

Lineage/type and clade assignments submitted by **{{ labdata.lab.lab_cod }}** were compared against the curated gold standard classifications for each sample included in the {{ comp_code }} component.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Per-sample lineage/type and clade assignment results for {{ labdata.lab.lab_cod }} ({{ comp_code }}).**

| Sample ID | Expected lineage/type | Reported lineage/type | Expected clade | Reported clade | Number of matches | Number of discrepancies |
|---|---|---|---|---|---|---|
{% for collecting_lab_sample_id, s in comp.samples.items() -%}
| {{ collecting_lab_sample_id }} | {{ s.classification.expected_lineage }} | {{ s.classification.lineage_assignment }} | {{ s.classification.expected_clade }} | {{ s.classification.clade_assignment }} | {{ s.classification.number_matches }} | {{ s.classification.number_discrepancies }} |
{% endfor %}

Table {{ table_counter.value }} summarises the concordance between expected and reported lineage/type and clade classifications for each sample.

{% set classification_concordance_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/classification_dimension_concordance.png" %}
{% if path_exists(classification_concordance_path) %}
{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} presents the distribution of classification outcomes across participating laboratories for each sample included in the {{ comp_code }} component, while highlighting the specific result reported by **{{ labdata.lab.lab_cod }}**.

{{ render_figure(
  classification_concordance_path,
  comp_code ~ ": lineage/type and clade classification outcomes across the network; black diamond indicates " ~ labdata.lab.lab_cod ~ "."
) }}

**_Figure {{ fig_counter.value }}_. Lineage/type and clade classification outcomes across participating laboratories ({{ comp_code }}).** Panel A shows the proportion of Match, Discrepancy, and Not provided outcomes for lineage/type assignments across participating laboratories for each sample. Panel B shows the corresponding proportions for clade assignments. Stacked bars represent the percentage of laboratories with correct classifications, incorrect classifications, or missing classifications relative to the curated gold standard. The black diamond marks the result reported by **{{ labdata.lab.lab_cod }}**, positioned within the Match, Discrepancy, or Not provided segment for each sample.
{% endif %}

## 9.{{ loop.index + 1 }}.4. Pipeline Benchmarking and Comparative Performance

The analytical workflow declared by **{{ labdata.lab.lab_cod }}** was benchmarked against other workflows implemented across the RELECOV network for the {{ comp_code }} component.

Positioning was evaluated based on four primary performance indicators:

1. Total number of discrepancies
2. Median consensus genome identity relative to the curated gold standard.
3. Total number of lineage/type and clade classification matches.
4. Metadata completeness

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

Table {{ table_counter.value }} summarises the software configuration declared by **{{ labdata.lab.lab_cod }}** for each analysed sample in {{ comp_code }}.

{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Workflow performance positioning for {{ labdata.lab.lab_cod }} within the network ({{ comp_code }}).**

| Metric | {{ labdata.lab.lab_cod }} workflow | Network median | Network min - max |
|---|---:|---:|---:|
| Total number of discrepancies | {{ comp.total_number_discrepancies_consensus if comp.total_number_discrepancies_consensus is not none else "NA" }} | {{ general.components[comp_code].consensus.total_median_discrepancies if general.components[comp_code].consensus.total_median_discrepancies is not none else "NA" }} | {{ general.components[comp_code].consensus.min_discrepancies if general.components[comp_code].consensus.min_discrepancies is not none else "NA" }} - {{ general.components[comp_code].consensus.max_discrepancies if general.components[comp_code].consensus.max_discrepancies is not none else "NA" }} |
| Median genome identity (%) | {{ pct(comp.median_genome_identity_pct, 4) if comp.median_genome_identity_pct is not none else "NA" }} | {{ pct(general.components[comp_code].consensus.median_identity_pct, 4) if general.components[comp_code].consensus.median_identity_pct is not none else "NA" }} | {{ general.components[comp_code].consensus.identity_pct_min if general.components[comp_code].consensus.identity_pct_min is not none else "NA" }} - {{ general.components[comp_code].consensus.identity_pct_max if general.components[comp_code].consensus.identity_pct_max is not none else "NA" }} |
| Total classification matches | {{ comp.total_classification_matches if comp.total_classification_matches is not none else "NA" }} | {{ general.components[comp_code].typing.total_classification_matches_median if general.components[comp_code].typing.total_classification_matches_median is not none else "NA" }} | {{ general.components[comp_code].typing.total_classification_matches_min if general.components[comp_code].typing.total_classification_matches_min is not none else "NA" }} - {{ general.components[comp_code].typing.total_classification_matches_max if general.components[comp_code].typing.total_classification_matches_max is not none else "NA" }}|
| Metadata completeness (%) | {{ pct(comp.metadata.completeness_pct, 2) if comp.metadata.completeness_pct is not none else "NA" }} | {{ pct(general.components[comp_code].metadata_completeness_median, 2) if general.components[comp_code].metadata_completeness_median is not none else "NA" }} | {{ general.components[comp_code].metadata_completeness_min_pct if general.components[comp_code].metadata_completeness_min_pct is not none else "NA" }} - {{ general.components[comp_code].metadata_completeness_max_pct if general.components[comp_code].metadata_completeness_max_pct is not none else "NA" }} |

Table {{ table_counter.value }} contextualises the performance of the declared workflow relative to aggregated network-level metrics. Network medians and (min-max) ranges provide a reference distribution against which the positioning of the declared workflow can be interpreted.

{% set workflow_positioning_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/workflow_positioning_boxplots.png" %}
{% if path_exists(workflow_positioning_path) %}
{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} illustrates the position of the workflow declared by **{{ labdata.lab.lab_cod }}** within the network-wide distribution of key performance indicators for the {{ comp_code }} component.

{{ render_figure(
  workflow_positioning_path,
  comp_code ~ ": workflow positioning relative to network-wide distributions, with " ~ labdata.lab.lab_cod ~ " highlighted by a black diamond."
) }}

**Figure {{ fig_counter.value }}. Workflow positioning within the RELECOV network for {{ comp_code }}.** Multi-panel boxplots summarise the laboratory-level distribution across the network for Panel A total consensus discrepancies, Panel B median genome identity, Panel C total classification matches, and Panel D metadata completeness. Only panels with evaluable data for **{{ labdata.lab.lab_cod }}** are shown. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. The black diamond corresponds to the results obtained by **{{ labdata.lab.lab_cod }}**.
{% endif %}

## 9.{{ loop.index + 1 }}.5. Metadata-Derived Analytical Metrics (per sample)

This section summarises selected quantitative analytical metrics declared in the metadata submission of **{{ labdata.lab.lab_cod }}**, disaggregated by sample within the {{ comp_code }} component.

Only metrics explicitly provided by the laboratory are included in the comparative assessment. Because laboratories may not complete all quantitative metadata fields for every sample, tables and panels below include only those metrics that were actually reported by **{{ labdata.lab.lab_cod }}**. Network-level medians and (min-max) ranges are shown for contextual interpretation.

#### Sample Quality Control Assessment

{{ labdata.lab.lab_cod }} QC evaluations (Pass/Fail) were compared against the predefined gold standard QC status for each sample within {{ comp_code }}. Samples without a laboratory-reported QC assessment are shown as `NA` in the table and are omitted from the comparative figure.
{% set table_counter.value = table_counter.value + 1 %}
**Table {{ table_counter.value }}. Sample-level QC assessment for {{ labdata.lab.lab_cod }} ({{ comp_code }}), benchmarked against network-level QC concordance.**

| Sample ID | Reported QC | Gold standard QC | Network % Match |
|---|---|---|---:|
{% for collecting_lab_sample_id, s in comp.samples.items() %}
{% set ns = (general.components[comp_code].qc.samples | selectattr("collecting_lab_sample_id","equalto",collecting_lab_sample_id) | list | first) %}
| {{ collecting_lab_sample_id }} | {{ s.qc_test if s.qc_test is not none else "NA" }} | {{ ns.gold_standard_qc if ns else "NA" }} | {{ pct(ns.match_rate_pct) if ns and ns.match_rate_pct is not none else "NA" }} |
{% endfor %}

Table {{ table_counter.value }} summarises the QC decision reported by **{{ labdata.lab.lab_cod }}** for each sample and benchmarks it against the network-level QC concordance for the same sample.

{% set qc_tests_reported = (comp.samples.values() | selectattr("qc_test", "ne", none) | list | length) > 0 %}
{% set qc_match_rate_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/qc_match_rate.png" %}
{% if qc_tests_reported and path_exists(qc_match_rate_path) %}
{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} contextualises the laboratory’s QC decision performance reported by **{{ labdata.lab.lab_cod }}** relative to network-wide QC concordance for each evaluable sample.

{{ render_figure(
  qc_match_rate_path,
  comp_code ~ ": sample-level QC concordance across the network, with " ~ labdata.lab.lab_cod ~ " highlighted."
) }}

**_Figure {{ fig_counter.value }}_. Sample-level QC concordance across the network for {{ comp_code }}, with {{ labdata.lab.lab_cod }} highlighted.** Stacked bars represent the network-wide proportions of Match and Discrepancy outcomes relative to the gold standard for each sample. The black diamond indicates whether **{{ labdata.lab.lab_cod }}** reported a Match or a Discrepancy for the corresponding evaluable sample.
{% else %}

No comparative QC concordance figure is shown for {{ comp_code }} because **{{ labdata.lab.lab_cod }}** did not report any sample-level QC assessment for this component.
{% endif %}

#### Other metrics

{% set metadata_metric_labels = {
  "per_genome_greater_10x": "% Genome > 10x",
  "depth_of_coverage_value": "Depth of coverage mean value",
  "per_Ns": "% Ns",
  "per_reads_virus": "% Reads virus",
  "per_reads_host": "% Reads host"
} %}
{% set metadata_metrics_reported = namespace(count=0) %}
{% for collecting_lab_sample_id, s in comp.samples.items() %}
{% set m = s.metadata_metrics %}
{% if m %}
{% set ns = (general.components[comp_code].metadata_metrics.samples | selectattr("sample_id","equalto",collecting_lab_sample_id) | list | first) %}
{% set sample_metrics = namespace(count=0) %}
{% for metric_key in metadata_metric_labels.keys() %}
{% if m.get(metric_key) is not none %}
{% set sample_metrics.count = sample_metrics.count + 1 %}
{% endif %}
{% endfor %}
{% if sample_metrics.count > 0 %}
{% set metadata_metrics_reported.count = metadata_metrics_reported.count + sample_metrics.count %}

{% set table_counter.value = table_counter.value + 1 %}

### {{ collecting_lab_sample_id }}

**Table {{ table_counter.value }}. Metadata-derived analytical metrics for {{ labdata.lab.lab_cod }} ({{ comp_code }}, {{ collecting_lab_sample_id }}).**

| Metric | {{ labdata.lab.lab_cod }} | Network median | Network min - max |
|---|---:|---:|---:|
{% for metric_key, metric_label in metadata_metric_labels.items() %}
{% if m.get(metric_key) is not none %}

| {{ metric_label }} | {{ m[metric_key] }} | {{ ns[metric_key].median if ns and ns.get(metric_key) else "NA" }} | {{ ns[metric_key].min if ns and ns.get(metric_key) else "NA" }} - {{ ns[metric_key].max if ns and ns.get(metric_key) else "NA" }} |
{% endif %}
{% endfor %}

Table {{ table_counter.value }} contextualises laboratory-reported analytical parameters relative to the distribution of values observed across the RELECOV network for the same sample.

{% endif %}
{% endif %}
{% endfor %}

{% set metadata_metrics_panel_path = "figures/labs/" ~ lab_code ~ "/" ~ comp_code ~ "/metadata_metrics_panel.png" %}
{% if metadata_metrics_reported.count > 0 and path_exists(metadata_metrics_panel_path) %}
{% set fig_counter.value = fig_counter.value + 1 %}

Figure {{ fig_counter.value }} illustrates the distribution of metadata-derived analytical metrics across participating laboratories for the {{ comp_code }} component.

{{ render_figure(
  metadata_metrics_panel_path,
  comp_code ~ ": distribution of metadata-derived analytical metrics across the network per sample; black diamond indicates " ~ labdata.lab.lab_cod ~ "."
) }}

**Figure {{ fig_counter.value }}. Distribution of metadata-derived analytical metrics across participating laboratories ({{ comp_code }}).**
Panel A shows genome coverage above 10x, Panel B depth of coverage, Panel C proportion of Ns, Panel D viral reads, and Panel E host reads. Only metrics actually reported by **{{ labdata.lab.lab_cod }}** are shown, so only panels with evaluable data are displayed. The central line indicates the median, boxes denote the interquartile range, whiskers represent the full observed range, translucent points correspond to individual laboratory observations, and hollow circles beyond the whiskers indicate outliers. The black diamond corresponds to the values reported by **{{ labdata.lab.lab_cod }}**.
{% else %}

No comparative metadata-derived analytical metrics figure is shown for {{ comp_code }} because **{{ labdata.lab.lab_cod }}** did not report any evaluable quantitative metadata metrics for this component.
{% endif %}

{% endfor %}

## Acknowledgement

We sincerely thank **{{ labdata.lab.lab_cod }}** for its participation in the 2026 RELECOV Dry-Lab EQA. The contribution of each laboratory is fundamental to maintaining analytical comparability, reproducibility, and interoperability across the network.

For any questions, technical clarifications, or follow-up discussions regarding this report, please contact the RELECOV WP.6 coordination team at [bioinformatica@isciii.es](mailto:bioinformatica@isciii.es).
{% endif %}
