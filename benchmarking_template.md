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

- [1. Introduction and scope](#1-introduction-and-scope)
- [2. Benchmarking approach](#2-benchmarking-approach)
- [3. Analytical workflow diversity across the RELECOV network](#3-analytical-workflow-diversity-across-the-relecov-network)
- [4. Component-Specific Results](#4-component-specific-results)
    - [4.1. SARS1 (SARS-CoV-2, Illumina)](#41-sars1-sars-cov-2-illumina)
    - [4.2. SARS2 (SARS-CoV-2, Oxford Nanopore Technologies)](#42-sars2-sars-cov-2-oxford-nanopore-technologies)
    - [4.3. FLU1 (Influenza virus, Illumina)](#43-flu1-influenza-virus-illumina)
    - [4.4. FLU2 (Influenza virus, Oxford Nanopore Technologies)](#44-flu2-influenza-virus-oxford-nanopore-technologies)
- [5. Discussion](#5-discussion)
- [6. Conclusions](#6-conclusions)
- [Appendix](#appendix)

## 1. Introduction and scope

This document presents the benchmarking component of the 2026 RELECOV Dry-Lab Interlaboratory Comparison Exercise. The exercise contributes to **Objective 2.1** of RELECOV 2.0, which focuses on *improving deep knowledge of the capacities and methodologies of the laboratories belonging to the network, as well as identifying a common methodology adapted to them and to the needs of the platform*. It also supports **Task T6.1**, which aims to *identify the most suitable bioinformatic analysis method for each sequencing platform, through an intercomparison exercise with simulated data for bioinformaticians*. The results will contribute to defining the workflow to be integrated into the RELECOV analytical platform.

The exercise is also related to **Milestone M6.3**, which concerns the *definition of sequencing and analysis protocols for each of the sequencing platforms*, and to **Task T6.5**, focused on *the adaptation and improvement of the analysis pipeline for the different sequencing platforms used by the laboratories of the network*.

The benchmarking exercise compares the analytical workflows reported by participating laboratories using SARS-CoV-2 and influenza datasets. The comparison considers factors that may affect analytical performance, including software selection, software and database versions, parameter settings, reference strategies, and reporting practices. These aspects are considered together with the metadata needed to interpret and reproduce the reported results.

The aim of the benchmarking is to assess analytical performance at the workflow and pipeline level rather than solely at the individual laboratory level. In particular, the analysis examines whether differences in workflow configuration are associated with differences in performance when laboratories analyse the same datasets under comparable conditions.

Further information on the exercise, including its design and the selection of samples, is provided in the main Interlaboratory Comparison Exercise document.

## 2. Benchmarking approach

For each declared pipeline or analytical workflow, including the software combinations and parameter configurations reported by the laboratories, performance was assessed across all laboratories using that approach. This approach allows workflow characteristics to be compared while taking into account the results obtained across the network.

The primary benchmarking criteria were based on the following performance indicators:

- Median consensus genome identity relative to the curated gold standard.
- Median number of discrepancies relative to the curated gold standard.
- Exact lineage/type and clade classification concordance.
- Median metadata completeness.

These metrics were used to assess whether workflows producing high consensus similarity also showed consistent downstream analytical performance. The analysis also considered the consistency of declared workflows across laboratories and the potential effect of software versions, reference genome selection, and parameter settings.

The benchmarking results were used to identify workflow configurations showing more consistent performance, as well as parameter settings or other workflow characteristics associated with systematic discrepancies. They also provide information on the potential effects of software versioning and reference genome selection on the results.

Because laboratories used different reference strategies, parameter settings, software versions, and output formats, the results should be interpreted as a descriptive comparison of the performance patterns observed across the network rather than as a strict ranking of tools or pipelines.

The benchmarking results will support the further development of RELECOV 2.0 by helping to identify workflow configurations associated with more consistent performance, define minimum performance criteria for network harmonisation, clarify the metadata and reporting requirements needed to interpret differences in analytical performance, and inform recommendations for the standardisation and further development of the RELECOV analytical platform.

## 3. Analytical workflow diversity across the RELECOV network

The metadata submissions provide an overview of the analytical workflows currently used across the RELECOV network. A total of {{ general.metadata_completeness.total_workflows }} distinct analytical workflows were identified, based on unique combinations of the software tools and versions reported in the metadata template.

The submitted metadata show considerable diversity in the software used for the main analytical steps. Based on the declared software name and version, where available, the following numbers of distinct software identities were identified:

- Consensus reconstruction software ({{ general.metadata_completeness.total_consensus_softwares }} distinct declared software identities)
- Variant calling tools ({{ general.metadata_completeness.total_variant_softwares }} distinct declared software identities)
- SARS-CoV-2 lineage assignment software ({{ general.metadata_completeness.total_lineage_assignment_softwares }} distinct declared software identities)
- Clade assignment software ({{ general.metadata_completeness.total_clade_assignment_softwares }} distinct declared software identities)
- Influenza type assignment software ({{ general.metadata_completeness.total_type_assignment_softwares }} distinct declared software identities)
- Influenza subtype assignment software ({{ general.metadata_completeness.total_subtype_assignment_softwares }} distinct declared software identities)

For the lineage, clade, type, and subtype benchmarking presented in [Section 4](#4-component-specific-results), these software identities are further separated by database version when this information was available. As a result, the categories used for benchmarking may be more detailed than the overall software diversity counts presented above.

The performance of individual software components is assessed in [Section 4](#4-component-specific-results) within the relevant analytical context: SARS-CoV-2 Illumina, SARS-CoV-2 Nanopore, Influenza Illumina, and Influenza Nanopore. This component-specific approach is important because performance varied depending on both the analytical component and the metric considered. Software comparisons were therefore not combined into a single ranking across the different components.

## 4. Component-specific Results

This section presents the analytical results stratified by component, allowing a detailed assessment of performance for each dataset and sequencing technology. For each component, the results are structured according to participation and submission metrics, consensus genome reconstruction performance, variant detection accuracy, and lineage, type/subtype, or clade assignment concordance, as applicable. Component-level analyses allow the identification of platform-specific patterns, differences associated with specific datasets, and variability between workflow configurations. This approach provides a more granular view of the performance differences observed across the network and can inform targeted harmonisation recommendations.

{% for comp_code, comp_net in general.components.items() %}

### 4.{{ loop.index }}. {{ comp_code }} ({{ comp_net.name }})

This section presents an exploratory comparative analysis of declared workflow configurations within {{ comp_code }}. Because laboratories differed in reference selection, software versions, parameterisation, reporting detail, and internal decision criteria, the results below should be interpreted as descriptive comparisons of observed performance patterns rather than as a controlled ranking of pipelines.

{% if comp_net.benchmarking.bioinformatics_protocol %}
#### 4.1. Bioinformatics protocol

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
#### 4.2. De-hosting software

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
#### 4.3. Preprocessing software

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
#### 4.4. Mapping software

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
#### 4.5. Assembly software

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
#### 4.6. Consensus software

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
#### 4.7. Variant calling software

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
#### 4.8. Clade Assignment Software

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
#### 4.9. Lineage Assignment Software Name

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
#### 4.10. Type Assignment Software Name

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

#### 4.11. Subtype Assignment Software Name

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

## 5. Discussion

The benchmarking results show that observed workflow performance varied across analytical components and was associated with differences in software configuration, reference strategies, parameterisation, and reporting quality. The results do not point to a single universally best-performing pipeline. Instead, they highlight the importance of considering performance together with the consistency of the results across laboratories and the information available to describe the workflow.

The metadata confirm that a diverse range of analytical configurations is currently in use across the RELECOV network. This diversity is analytically valuable, as it allows different approaches to be compared within the same benchmarking exercise. At the same time, interpretation is limited by incomplete reporting. Only {{ pct(general.metadata_completeness.software_version_pct) }} of software-version fields were completed, while minimum coverage thresholds, variant calling parameters, and reference genome identifiers were reported for {{ pct(general.metadata_completeness.coverage_threshold_pct) }}, {{ pct(general.metadata_completeness.variant_calling_params_pct) }}, and {{ pct(general.metadata_completeness.reference_genome_pct) }} of submitted samples, respectively. This means that some plausible explanations for observed performance differences can only be considered as contributing factors rather than demonstrated causal effects.

The metadata also show substantial heterogeneity in the parameters used by laboratories. At least 8 different conventions for minimum coverage thresholds, 13 distinct consensus parameter strings, 13 mapping parameter strings, and 14 variant calling parameter strings were reported. These differences do not demonstrate that a particular parameter caused an observed discrepancy, but they show that laboratories were applying different masking, filtering, and coverage criteria. This variation is therefore an important part of the context in which workflow performance should be interpreted.

The SARS-CoV-2 components showed relatively consistent performance across the main workflow configurations. In particular, nf-core/viralrecon was represented by several laboratories in both SARS1 and SARS2 and showed high median genome identity, relatively low discrepancy counts, and high metadata completeness. This makes it a more informative configuration for assessing reproducibility across the network than configurations represented by a single laboratory. However, these results should not be interpreted as evidence that nf-core/viralrecon is universally superior, since the different workflows were not evaluated under fully standardised laboratory conditions.

The SARS-CoV-2 results also indicate that sequence-level discrepancies and downstream classification concordance do not always vary in parallel. Several configurations retained high lineage or clade concordance despite differences in consensus reconstruction. This suggests that modest sequence-level differences did not necessarily translate into differences in classification for the datasets evaluated here. The interpretation of this pattern should nevertheless remain specific to the datasets and workflows included in the exercise.

The influenza components showed greater heterogeneity between declared workflow configurations. In FLU1, some of the most favourable raw performance values were associated with configurations represented by a single laboratory, whereas custom pipelines and different IRMA showed broader performance ranges. In FLU2, differences between workflow configurations and software versions were also apparent across consensus reconstruction and classification metrics. In FLU2, IRMA v1.3.1 provided the most balanced profile across discrepancy burden and classificationperformance, whereas other versions, particularly IRMA v1.2.0, showed poorer consensus reconstruction behaviour despite acceptable classification fields. These patterns suggest that influenza benchmarking may be more sensitive to workflow configuration, parameterisation, and sample characteristics in the datasets evaluated here. However, the available data do not allow the individual contribution of these factors to be separated.

An important limitation of the benchmarking is therefore the unequal representation of workflow configurations across laboratories. Some configurations were reported by several laboratories, whereas others were represented by only one laboratory. In addition, software version, reference selection, parameterisation, and other workflow characteristics were not always independent. The observed differences should consequently be interpreted as descriptive performance patterns rather than evidence of causal effects attributable to individual software tools or parameters.

Another important finding is that analytical performance and reporting quality should be considered as related but separate dimensions. Good analytical performance does not necessarily imply complete or reproducible workflow documentation, while complete metadata do not by themselves guarantee high analytical performance. More consistent reporting of software and database versions, reference sequences, and key parameters would therefore improve the interpretation of future benchmarking results.

Overall, the results support a harmonisation strategy based on minimum analytical and reporting requirements rather than the adoption of a single pipeline. Establishing common criteria for workflow documentation, reference selection, parameter reporting, and minimum performance would improve comparability across laboratories while retaining flexibility in the choice of analytical tools.

## 6. Conclusions

The benchmarking exercise identified differences in workflow performance across the four analytical components, but these differences were context-dependent and could not be reduced to a single ranking of pipelines.

The main conclusions are:

- Workflow performance should be evaluated within the relevant analytical context, taking into account the pathogen, sequencing platform, and analytical component.
- Workflow configurations represented by multiple laboratories provide a stronger basis for assessing reproducibility than configurations represented by a single laboratory.
- Analytical performance and reporting quality should be assessed as complementary but separate dimensions.
- Software versions, database versions, reference sequences, and key parameters should be reported consistently to improve the interpretation and reproducibility of benchmarking results.
- Harmonisation should focus on minimum analytical and reporting requirements rather than enforcing a single pipeline across the network.
- Future interlaboratory exercises would benefit from more standardised workflow definitions and reporting requirements, allowing the effects of individual workflow components to be assessed more reliably.

Taken together, these findings provide a practical basis for improving the comparability of bioinformatic workflows within RELECOV and for defining recommendations to support the further harmonisation and development of the RELECOV analytical platform.

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
