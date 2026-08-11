{# =========================
  JINJA2 TEMPLATE: RELECOV Interlaboratory Comparison Exercise REPORT
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
  "ambiguity2nt": "Nucleotide instead of ambiguity",
  "nt2ambiguity": "Ambiguity instead of nucleotide",
  "ns2nt": "Nucleotide stretch instead of stretch of Ns",
  "nt2ns": "Stretch of Ns instead of nucleotide stretch",
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

# Benchmarking of the Interlaboratory Comparison Exercise RELECOV 2.0

##### Sarai Varona, Enrique Sapena, Pablo Mata, Alejandro Bernabéu, Pau Pascual, Magdalena Matito, Juan Ledesma, Emilia Arjona, Victor Lopez, Olga Dolgova, Sara Monzón, Isabel Cuesta

## Table of Contents

- [1. Introduction](#1-introduction)
- [2. Pipeline Benchmarking and Comparative Performance](#2-pipeline-benchmarking-and-comparative-performance)
- [3. Component-Specific Results](#3-component-specific-results)
    - [3.1. SARS1 (SARS-CoV-2, Illumina)](#31-sars1-sars-cov-2-illumina)
    - [3.2. SARS2 (SARS-CoV-2, Oxford Nanopore Technologies)](#32-sars2-sars-cov-2-oxford-nanopore-technologies)
    - [3.3. FLU1 (Influenza virus, Illumina)](#33-flu1-influenza-virus-illumina)
    - [3.4. FLU2 (Influenza virus, Oxford Nanopore Technologies)](#34-flu2-influenza-virus-oxford-nanopore-technologies)
- [4. Discussion](#7-discussion)
- [5. Conclusions](#8-conclusions)
- [Appendix](#appendix)

## 1. Introduction

The RELECOV Network aims to strengthen genomic surveillance of respiratory viruses by developing and harmonising analytical capacities across the participating laboratories. In this context, it was essential to **assess the consistency, reproducibility and maturity of the bioinformatic workflows implemented across the network**.

To this end, an **Interlaboratory Comparison Exercise exercise in dry lab format** was conducted, based on the European Centre for Disease Prevention and Control (ECDC) 2024 dry-lab EQA. The exercise focused on the bioinformatic characterisation of respiratory viruses, covering key analytical tasks including viral genome reconstruction, variant identification, and lineage and clade assignment.

Beyond its role as an external quality assessment of laboratory performance, the exercise was also designed to support the methodological harmonisation objectives of RELECOV 2.0. A central component of this initiative was to characterise the diversity of analytical pipelines implemented across the RELECOV Network, evaluate their performance under the conditions of this exercise, and generate evidence to support future harmonisation activities within the network. This evaluation contributes directly to **Objective 2.1** of RELECOV 2.0, which focuses on _improving deep knowledge of the capacities and methodologies of the laboratories belonging to the network, as well as identifying a common methodology adapted to them and to the needs of the platform_. Furthermore, the exercise provides the practical evidence base required for **Task T6.1**, which aims to  _identify the most suitable bioinformatic analysis method for each sequencing platform, through an intercomparison exercise with simulated data for bioinformaticians_, in order to define the workflow that should be integrated into the RELECOV analytical platform.

The exercise was also aligned with **Milestone M6.3**, which pertains to _define sequencing and analysis protocols for each of the sequencing platforms_. In addition, the exercise provided operational insights relevant to **Task T6.5**, which addresses _the adaptation and improvement of the analysis pipeline for the different sequencing platforms used by the laboratories of the network_. It also contributed to **Task T6.4**, related to _sequence metadata annotation with ontologies, schema generation, parsing and validation_, by highlighting practical issues affecting metadata completeness, controlled-vocabulary use, and the consistency of reported analytical parameters.

The overall objective of the exercise was to **assess the bioinformatic performance of the participating laboratories, identify areas for improvement, and promote the adoption of consistent and comparable analytical practices across the network**. The outcomes presented in this report are expected to strengthen RELECOV’s preparedness and response capacity for routine surveillance and public health emergencies, while supporting the harmonisation objectives defined within RELECOV 2.0.

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

### 2. Pipeline Benchmarking and Comparative Performance

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

## 3. Component-specific Results

This section presents the analytical results disaggregated by component, allowing a detailed assessment of performance within each dataset and sequencing technology. For each component, results are structured according to participation and submission metrics, consensus genome reconstruction performance, variant detection accuracy, and Lineage, Subtype or clade assignment concordance, as applicable.

Component-level analyses enable identification of platform-specific patterns, dataset-dependent challenges, and variability associated with particular sample characteristics. This approach facilitates a more granular interpretation of performance differences observed at the network level and supports targeted harmonisation recommendations.

All component-level results below are reported using the same evaluation framework described in [Section 4](#4-methodology-of-evaluation).

{% for comp_code, comp_net in general.components.items() %}

### 3.{{ loop.index }}. {{ comp_code }} ({{ comp_net.name }})


This section presents an exploratory comparative analysis of declared workflow configurations within {{ comp_code }}. Because laboratories differed in reference selection, software versions, parameterisation, reporting detail, and internal decision criteria, the results below should be interpreted as descriptive comparisons of observed performance patterns rather than as a controlled ranking of pipelines.

{% if comp_net.benchmarking.bioinformatics_protocol %}
#### 3.1. Bioinformatics protocol

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
#### 3.2. De-hosting software

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
#### 3.3. Preprocessing software

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
#### 3.4. Mapping software

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
#### 3.5. Assembly software

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
#### 3.6. Consensus software

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
#### 3.7. Variant calling software

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
#### 3.8. Clade Assignment Software

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
#### 3.9. Lineage Assignment Software Name

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
#### 3.10. Type Assignment Software Name

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

#### 3.11. Subtype Assignment Software Name

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

## 4. Discussion

The 2026 RELECOV Dry-Lab Interlaboratory Comparison Exercise provides the first network-wide dry-lab assessment focused specifically on bioinformatic performance across consensus reconstruction, variant reporting, classification, metadata reporting, and QC interpretation. By combining ECDC datasets with in-silico influenza material, the exercise captures both routine-use analytical behaviour and performance under heterogeneous reference and reporting conditions.

The benchmarking results also suggest that apparent top-performing configurations should be interpreted against the number of laboratories supporting them. In SARS1 and SARS2, several of the most favourable raw performance values were associated with single-laboratory configurations, whereas nf-core/viralrecon was supported by multiple laboratories and combined consistently high genome identity with low discrepancy burdens and strong metadata completeness. That pattern makes it a more informative indicator of reproducible network performance than a nominally better single-observation configuration. At the same time, the SARS results do not show a simple linear relationship between consensus discrepancy counts and classification concordance: most pipelines retained very high lineage or clade performance despite modest sequence-level variation, suggesting that SARS classification is robust to moderate reconstruction differences but still sensitive to reporting quality.

The influenza benchmarking profiles were more heterogeneous and showed stronger configuration-dependent effects. In FLU1, DRAGEN achieved the highest identity but only on a single observation, whereas custom pipelines and different IRMA versions showed broader performance ranges. INSaFLU performed worse in both identity and discrepancy burden in that component, suggesting that workflow choice can have more impact in influenza than in SARS-CoV-2 under these datasets. In FLU2, IRMA v1.3.1 provided the most balanced profile across discrepancy burden and classification performance, whereas other versions, particularly IRMA v1.2.0, showed poorer consensus reconstruction behaviour despite acceptable classification fields. Taken together, these results indicate that influenza benchmarking is more sensitive to software versioning, parameterisation, and component-specific sample properties, and therefore requires more explicit best-practice recommendations rather than simple transfer of SARS-CoV-2 assumptions.

Taken together, the results support a harmonisation strategy centred on minimum performance and reporting standards rather than on enforcement of a single analytical pipeline. The data do not support a universal workflow ranking that would apply equally across all viruses, platforms, and tasks. Instead, they show that performance depends on the interaction between dataset characteristics, reporting conventions, software choice, and parameterisation.

## 5. Conclusions

The 2026 RELECOV Dry-Lab Interlaboratory Comparison Exercise shows that participating laboratories already have bioinformatic capacity for respiratory virus genomic surveillance, but that performance and comparability still depend strongly on the analytical context in which each task is performed.

Consensus genome reconstruction was generally strongest in the Illumina-based components, while broader performance ranges in SARS2 and FLU2 indicate that a subset of submissions remained highly sensitive to masking behaviour, coverage thresholds, and consensus-generation choices. Variant analysis showed that direct SARS-CoV-2 comparison against curated reference sets is feasible, whereas influenza reporting remained much more heterogeneous because of mixed allele frequency reporting strategies, multiple reference backbones, and large discrepancies between metadata-reported and VCF-derived summaries.

Classification and QC interpretation further showed that harmonisation challenges are not limited to core sequence processing. Lineage/type assignment was more concordant than clade assignment, and part of the excess clade discordance in SARS-CoV-2 appears to reflect metadata completion and nomenclature problems in the clade field itself. QC interpretation was also unevenly reported, with only a subset of laboratories providing explicit sample-level QC assessments in the metadata template.

Overall, the results support RELECOV 2.0 priorities centred on:

- minimum performance standards for consensus reconstruction and variant reporting
- clearer rules for masking, coverage thresholds, and allele frequency reporting
- stronger metadata requirements for software versions, parameters, and reference genomes
- improved consistency in classification and QC field completion
- component-aware benchmarking rather than a single cross-context workflow ranking

Taken together, these findings provide a practical basis for harmonising analytical expectations across the network while preserving the methodological flexibility needed for different pathogens, sequencing platforms, and surveillance scenarios.

The Interlaboratory Comparison Exercise therefore provides a technical basis for harmonised, performance-driven genomic surveillance within RELECOV 2.0.

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
