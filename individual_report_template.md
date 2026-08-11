{# =========================
  JINJA2 TEMPLATE: Individual Laboratory Technical Report
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

# Individual Laboratory Technical Report - Interlaboratory Comparison Exercise RELECOV 2.0

##### Sarai Varona, Enrique Sapena, Pablo Mata, Alejandro Bernabéu, Pau Pascual, Magdalena Matito, Juan Ledesma, Emilia Arjona, Victor Lopez, Olga Dolgova, Sara Monzón, Isabel Cuesta

<!-- TEMPLATE_TOC -->
## Table of Contents

- [1. Participation Overview](#1-participation-overview)
{% for comp_code, comp in labdata.components.items() %}
- [{{ loop.index + 1 }}. {{ comp_code }} ({{ comp.display_name }})](#component-{{ loop.index }}-{{ comp_code|lower|replace(' ', '-') }})
    - [{{ loop.index + 1 }}.1. Consensus Genome Reconstruction Performance](#component-{{ loop.index }}-{{ comp_code|lower|replace(' ', '-') }}-consensus)
    - [{{ loop.index + 1 }}.2. Variant Detection Performance](#component-{{ loop.index }}-{{ comp_code|lower|replace(' ', '-') }}-variant)
    - [{{ loop.index + 1 }}.3. Lineage, Subtype and Clade Assignment](#component-{{ loop.index }}-{{ comp_code|lower|replace(' ', '-') }}-classification)
    - [{{ loop.index + 1 }}.4. Pipeline Benchmarking and Comparative Performance](#component-{{ loop.index }}-{{ comp_code|lower|replace(' ', '-') }}-workflow)
    - [{{ loop.index + 1 }}.5. Metadata-Derived Analytical Metrics (per sample)](#component-{{ loop.index }}-{{ comp_code|lower|replace(' ', '-') }}-metadata)
{% endfor %}
- [Appendix](#appendix)

{% set lab_code = labdata.lab.lab_cod | default(labdata.lab.submitting_institution_id) %}

This section provides a detailed technical assessment of the analytical results submitted by **{{ labdata.lab.lab_cod }}** within the 2026 RELECOV Dry-Lab Interlaboratory Comparison Exercise. Performance metrics are benchmarked against curated gold standards and contextualised relative to aggregated network-wide performance distributions. Network medians and interquartile ranges are provided for comparative interpretation, without disclosure of other laboratories’ identities.

The purpose of this section is to support technical optimisation, parameter harmonisation, and alignment with the analytical standards defined within RELECOV 2.0.

Only files, metadata fields, and derived analytical metrics actually provided by the laboratory are displayed in this individual report. If a file was not submitted, or a metadata field was not provided, the corresponding table entries, panels, or figures are omitted for that laboratory.

<h2 id="1-participation-overview" class="no-page-break">1. Participation Overview</h2>

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

## {{ loop.index + 1 }}. {{ comp_code }} ({{ comp.display_name }})

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

### {{ loop.index + 1 }}.1. Consensus Genome Reconstruction Performance

Consensus genome sequences (`.fasta`) submitted by **{{ labdata.lab.lab_cod }}** were compared against the curated gold standard for each sample included in the {{ comp_code }} component.

#### Per-sample summary metrics

{% set appendix_table_counter.value = appendix_table_counter.value + 1 %}
{% set lab_consensus_metrics_table_num = appendix_table_counter.value %}
The detailed per-sample consensus reconstruction metrics for **{{ labdata.lab.lab_cod }}** are provided in Appendix Table {{ lab_consensus_metrics_table_num }}. The figure below summarises overall sequence similarity and discrepancy burden relative to the curated gold standard for {{ labdata.lab.lab_cod }} compared to the network.


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

#### Discrepancy type breakdown per sample
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

### {{ loop.index + 1 }}.2. Variant Detection Performance

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

### {{ loop.index + 1 }}.3. Lineage, Subtype and Clade Assignment

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

### {{ loop.index + 1 }}.4. Pipeline Benchmarking and Comparative Performance

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
| Total number of discrepancies in consensus | {{ comp.total_number_discrepancies_consensus if comp.total_number_discrepancies_consensus is not none else "NA" }} | {{ general.components[comp_code].workflow_total_discrepancies_median if general.components[comp_code].workflow_total_discrepancies_median is not none else "NA" }} | {{ general.components[comp_code].workflow_total_discrepancies_min if general.components[comp_code].workflow_total_discrepancies_min is not none else "NA" }} - {{ general.components[comp_code].workflow_total_discrepancies_max if general.components[comp_code].workflow_total_discrepancies_max is not none else "NA" }} |
| Median genome identity (%) | {{ pct(comp.median_genome_identity_pct) if comp.median_genome_identity_pct is not none else "NA" }} | {{ pct(general.components[comp_code].workflow_median_identity_pct_median) if general.components[comp_code].workflow_median_identity_pct_median is not none else "NA" }} | {{ pct(general.components[comp_code].workflow_median_identity_pct_min) if general.components[comp_code].workflow_median_identity_pct_min is not none else "NA" }} - {{ pct(general.components[comp_code].workflow_median_identity_pct_max) if general.components[comp_code].workflow_median_identity_pct_max is not none else "NA" }} |
| Total classification matches | {{ comp.total_classification_matches if comp.total_classification_matches is not none else "NA" }} | {{ general.components[comp_code].typing.total_classification_matches_median if general.components[comp_code].typing.total_classification_matches_median is not none else "NA" }} | {{ general.components[comp_code].typing.total_classification_matches_min if general.components[comp_code].typing.total_classification_matches_min is not none else "NA" }} - {{ general.components[comp_code].typing.total_classification_matches_max if general.components[comp_code].typing.total_classification_matches_max is not none else "NA" }}|
| Metadata completeness (%) | {{ pct(comp.metadata.completeness_pct, 2) if comp.metadata.completeness_pct is not none else "NA" }} | {{ pct(general.components[comp_code].metadata_completeness_median, 2) if general.components[comp_code].metadata_completeness_median is not none else "NA" }} | {{ general.components[comp_code].metadata_completeness_min_pct if general.components[comp_code].metadata_completeness_min_pct is not none else "NA" }} - {{ general.components[comp_code].metadata_completeness_max_pct if general.components[comp_code].metadata_completeness_max_pct is not none else "NA" }} |

### {{ loop.index + 1 }}.5. Metadata-Derived Analytical Metrics (per sample)

This section summarises selected quantitative analytical metrics declared in the metadata submission of **{{ labdata.lab.lab_cod }}**, disaggregated by sample within the {{ comp_code }} component.

Only metrics explicitly provided by the laboratory are included in the comparative assessment. Because laboratories may not complete all quantitative metadata fields for every sample, tables and panels below include only those metrics that were actually reported by **{{ labdata.lab.lab_cod }}**. Network-level medians and (min-max) ranges are shown for contextual interpretation.

#### Sample Quality Control Assessment

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

#### Other metrics

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

## Acknowledgement

We sincerely thank **{{ labdata.lab.lab_cod }}** for its participation in the 2026 RELECOV Dry-Lab Interlaboratory Comparison Exercise. The contribution of each laboratory is fundamental to maintaining analytical comparability, reproducibility, and interoperability across the network.

For any questions, technical clarifications, or follow-up discussions regarding this report, please contact the RELECOV WP.6 coordination team at [bioinformatica@isciii.es](mailto:bioinformatica@isciii.es).

## Appendix

This appendix for {{ labdata.lab.lab_cod }} is reserved for supplementary material that may support interpretation of the report but is not essential to the main narrative. Additional figures, extended tables, sensitivity analyses, or other secondary outputs can be included here when relevant.

{# Use `appendix_fig_counter` and `appendix_table_counter` for supplementary material moved here.
   Refer to them from the main text as "Appendix Figure X" and "Appendix Table X". #}

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
