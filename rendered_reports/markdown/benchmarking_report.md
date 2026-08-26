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

The metadata submissions provide an overview of the analytical workflows currently used across the RELECOV network. A total of 9 distinct analytical workflows were identified, based on unique combinations of the software tools and versions reported in the metadata template.

The submitted metadata show considerable diversity in the software used for the main analytical steps. Based on the declared software name and version, where available, the following numbers of distinct software identities were identified:

- Consensus reconstruction software (16 distinct declared software identities)
- Variant calling tools (20 distinct declared software identities)
- SARS-CoV-2 lineage assignment software (4 distinct declared software identities)
- Clade assignment software (10 distinct declared software identities)
- Influenza type assignment software (4 distinct declared software identities)
- Influenza subtype assignment software (8 distinct declared software identities)

For the lineage, clade, type, and subtype benchmarking presented in [Section 4](#4-component-specific-results), these software identities are further separated by database version when this information was available. As a result, the categories used for benchmarking may be more detailed than the overall software diversity counts presented above.

The performance of individual software components is assessed in [Section 4](#4-component-specific-results) within the relevant analytical context: SARS-CoV-2 Illumina, SARS-CoV-2 Nanopore, Influenza Illumina, and Influenza Nanopore. This component-specific approach is important because performance varied depending on both the analytical component and the metric considered. Software comparisons were therefore not combined into a single ranking across the different components.

## 4. Component-specific Results

This section presents the analytical results stratified by component, allowing a detailed assessment of performance for each dataset and sequencing technology. For each component, the results are structured according to participation and submission metrics, consensus genome reconstruction performance, variant detection accuracy, and lineage, type/subtype, or clade assignment concordance, as applicable. Component-level analyses allow the identification of platform-specific patterns, differences associated with specific datasets, and variability between workflow configurations. This approach provides a more granular view of the performance differences observed across the network and can inform targeted harmonisation recommendations.

### 4.1. SARS1 (SARS-CoV-2, Illumina)

This section presents an exploratory comparative analysis of declared workflow configurations within SARS1. Because laboratories differed in reference selection, software versions, parameterisation, reporting detail, and internal decision criteria, the results below should be interpreted as descriptive comparisons of observed performance patterns rather than as a controlled ranking of pipelines.

#### 4.1. Bioinformatics protocol

Based on metadata submissions, 5 distinct bioinformatics protocols were reported for the SARS1 component. These summaries compare declared workflow configurations as they were used in practice across participating laboratories.

<figure>
<img src="figures/SARS1/bioinformatics_protocol_discrepancies_boxplot_by_pipeline.png" alt="Distribution of consensus discrepancies by pipeline configuration for SARS1." style="width: 80%; max-width: 80%;"/>
</figure>

**Figure 1. Distribution of consensus discrepancies by declared pipeline configuration for SARS1.** This boxplot summarises sample-level consensus discrepancies stratified by bioinformatics protocol. The left y-axis shows discrepancy counts, while the right y-axis overlays lineage/type and clade classification accuracy for the same software configuration. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

**Table 1. Performance summary of declared bioinformatics protocols for SARS1.**

| Bioinformatics protocol | Version | N labs | Median genome identity (%) | Median discrepancies | Median metadata completeness (%) | Clade concordance (%) | Lineage/type concordance (%) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Custom pipeline/workflow |  | 5 | 99.63 | 5 | 62.7 | 20.0 | 60.0 |
| DRAGEN Targeted Microbial | 1.1.0 | 2 | 99.56 | 3 | 44.4 | 100.0 | 87.5 |
| INSaFLU | 2.2.2 | 1 | 99.73 | 2.5 | 44.1 | 100.0 | 100.0 |
| INSaFLU | Web version | 1 | 99.70 | 3 | 33.9 | 100.0 | 100.0 |
| nf-core/viralrecon | 3.0.0 | 6 | 99.69 | 6 | 91.4 | 83.3 | 100.0 |

The observed differences across configurations should be read in the context of heterogeneous laboratory practices, including differences in reference choice, parameterisation, and thresholding. The table and figures therefore help identify recurrent performance patterns within SARS1, but they do not support a strong cross-laboratory ranking of bioinformatics protocols.

<figure class="benchmark-figure">
<img src="figures/SARS1/bioinformatics_protocol_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pipeline configuration for SARS1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 2. Distribution of performance metrics by declared pipeline configuration for SARS1.** Multi-panel boxplots summarise sample-level performance stratified by bioinformatics protocols. Panel A displays genome identity (%), Panel B metadata completeness (%), and Panel C exact classification concordance (%). X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Where required, Panel A uses a truncated y-axis to highlight differences among high-identity values. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.2. De-hosting software

4 distinct de-hosting software declarations were reported for the SARS1 component.

The distribution below reflects only configurations with evaluable percentage of host reads values in the reported metadata, so the number of boxplots may be lower than the total number of declared de-hosting configurations. The full list of declared configurations and associated summary values is provided in Appendix Table 1.

<figure class="benchmark-figure">
<img src="figures/SARS1/dehosting_metric_boxplots_by_pipeline.png" alt="Distribution of percentage of host reads metrics by dehosting software version for SARS1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 3. Distribution of percentage of host reads by declared dehosting software version for SARS1.** Boxplots summarise sample-level percentage of host reads stratified by dehosting software version. Only configurations with evaluable percentage of host reads values are displayed, so some declared software categories may be absent from the plot. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each version, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.3. Preprocessing software

6 distinct pre-processing software configurations were reported for the SARS1 component.

Only pre-processing configurations with evaluable observations for the displayed metrics contribute to the figure, so some declared categories may not appear in the plot. The complete list of declared configurations and their summary values is provided in Appendix Table 2.

<figure class="benchmark-figure">
<img src="figures/SARS1/preprocessing_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pre-processing software configuration for SARS1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 4. Distribution of performance metrics by declared pre-processing software configuration for SARS1.** Multi-panel boxplots summarise sample-level performance stratified by pre-processing software. Panel A displays Number of reads sequenced and Panel B Reads passing filters. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.4. Mapping software

6 distinct mapping software configurations were reported for the SARS1 component.

The mapping boxplots include only configurations for which the relevant performance metrics were available, which means that fewer categories may be plotted than were originally declared. Full configuration-level summaries are reported in Appendix Table 3.

<figure class="benchmark-figure">
<img src="figures/SARS1/mapping_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by mapping software configuration for SARS1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 5. Distribution of performance metrics by declared mapping software configuration for SARS1.** Boxplots summarise sample-level performance stratified by mapping software. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.5. Assembly software

3 distinct assembly software configurations were reported for the SARS1 component.

The assembly figures are restricted to configurations with evaluable values for the displayed metrics. As a result, some declared assembly categories may be absent from the plots; the full set of declared configurations and summary values is provided in Appendix Table 4.

<figure class="benchmark-figure">
<img src="figures/SARS1/assembly_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by assembly software configuration for SARS1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 6. Distribution of performance metrics by declared assembly software configuration for SARS1.** Multi-panel boxplots summarise sample-level performance stratified by assembly software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers. Panel B uses a truncated y-axis to highlight differences among high-identity values.

#### 4.6. Consensus software

8 distinct consensus software configurations were reported for the SARS1 component.

Only consensus software configurations with sufficient evaluable data are visualised in the figure below, so the plotted set may be smaller than the total set of declarations. All declared configurations and their associated summary values can be reviewed in Appendix Table 5.

<figure class="benchmark-figure">
<img src="figures/SARS1/consensus_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by consensus software configuration for SARS1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 7. Distribution of performance metrics by declared consensus software configuration for SARS1.** Multi-panel boxplots summarise sample-level performance stratified by consensus software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers. For SARS1 and FLU1, Panel B uses a truncated y-axis to highlight differences among high-identity values.

#### 4.7. Variant calling software

8 distinct variant calling software configurations were reported for the SARS1 component.

The plotted variant calling categories correspond only to configurations with evaluable observations for the displayed metrics. Consequently, the figure may show fewer configurations than were declared overall; the complete summaries are listed in Appendix Table 6.

<figure class="benchmark-figure">
<img src="figures/SARS1/variant_calling_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by variant calling software configuration for SARS1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 8. Distribution of performance metrics by declared variant calling software configuration for SARS1.** Panel A is a stacked bar chart showing the number of evaluable samples assigned to each allele frequency reporting pattern for each software configuration. Boxplot Panel B displays discrepancies in reported variants with AF >=75% in the submitted VCF, Panel C discrepancies in reported variants with effect, Panel D successful hits, and Panel E total discrepancies. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. In the boxplots, the central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.8. Clade Assignment Software

10 distinct clade assignment software configurations were reported for the SARS1 component. For this category, configurations were counted as unique combinations of software name, software version, and clade assignment database version when available.

Because clade concordance could not be evaluated for every declared configuration, the boxplot includes only categories with usable observations. The full list of declared configurations and their summary values is available in Appendix Table 7.

<figure class="benchmark-figure landscape-benchmark-figure">
<img src="figures/SARS1/clade_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by clade assignment software configuration for SARS1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 9. Distribution of clade concordance by declared clade assignment software configuration for SARS1.** This boxplot summarises sample-level clade concordance stratified by clade assignment software configuration, where each configuration corresponds to a unique combination of software name, software version, and clade assignment database version when available. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.9. Lineage Assignment Software Name

6 distinct lineage assignment software configurations were reported for the SARS1 component. For this category, configurations were counted as unique combinations of software name, software version, and lineage assignment database version when available.

Lineage assignment configurations are shown only when concordance values were evaluable for the submitted observations, so the plotted categories may represent only a subset of the declarations. The complete configuration-level summary is provided in Appendix Table 8.

<figure class="benchmark-figure landscape-benchmark-figure">
<img src="figures/SARS1/lineage_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by lineage assignment software configuration for SARS1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 10. Distribution of lineage concordance by declared lineage assignment software configuration for SARS1.** This boxplot summarises sample-level lineage concordance stratified by lineage assignment software configuration, where each configuration corresponds to a unique combination of software name, software version, and lineage assignment database version when available. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

### 4.2. SARS2 (SARS-CoV-2, Oxford Nanopore Technologies)

This section presents an exploratory comparative analysis of declared workflow configurations within SARS2. Because laboratories differed in reference selection, software versions, parameterisation, reporting detail, and internal decision criteria, the results below should be interpreted as descriptive comparisons of observed performance patterns rather than as a controlled ranking of pipelines.

#### 4.1. Bioinformatics protocol

Based on metadata submissions, 5 distinct bioinformatics protocols were reported for the SARS2 component. These summaries compare declared workflow configurations as they were used in practice across participating laboratories.

<figure>
<img src="figures/SARS2/bioinformatics_protocol_discrepancies_boxplot_by_pipeline.png" alt="Distribution of consensus discrepancies by pipeline configuration for SARS2." style="width: 80%; max-width: 80%;"/>
</figure>

**Figure 11. Distribution of consensus discrepancies by declared pipeline configuration for SARS2.** This boxplot summarises sample-level consensus discrepancies stratified by bioinformatics protocol. The left y-axis shows discrepancy counts, while the right y-axis overlays lineage/type and clade classification accuracy for the same software configuration. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

**Table 2. Performance summary of declared bioinformatics protocols for SARS2.**

| Bioinformatics protocol | Version | N labs | Median genome identity (%) | Median discrepancies | Median metadata completeness (%) | Clade concordance (%) | Lineage/type concordance (%) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Artic pipeline | 1.8.5 | 1 | 99.48 | 5 | 79.7 | 100.0 | 66.7 |
| Custom pipeline/workflow |  | 2 | 99.51 | 15 | 56.8 | 100.0 | 33.3 |
| INSaFLU | 2.2.2 | 1 | 99.88 | 11 | 37.1 | 100.0 | 66.7 |
| INSaFLU | Web version | 1 | 99.79 | 11 | 30.5 | 100.0 | 66.7 |
| nf-core/viralrecon | 3.0.0 | 4 | 99.88 | 9.5 | 91.5 | 75.0 | 58.3 |

The observed differences across configurations should be read in the context of heterogeneous laboratory practices, including differences in reference choice, parameterisation, and thresholding. The table and figures therefore help identify recurrent performance patterns within SARS2, but they do not support a strong cross-laboratory ranking of bioinformatics protocols.

<figure class="benchmark-figure">
<img src="figures/SARS2/bioinformatics_protocol_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pipeline configuration for SARS2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 12. Distribution of performance metrics by declared pipeline configuration for SARS2.** Multi-panel boxplots summarise sample-level performance stratified by bioinformatics protocols. Panel A displays genome identity (%), Panel B metadata completeness (%), and Panel C exact classification concordance (%). X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Where required, Panel A uses a truncated y-axis to highlight differences among high-identity values. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.2. De-hosting software

1 distinct de-hosting software declarations were reported for the SARS2 component.

The distribution below reflects only configurations with evaluable percentage of host reads values in the reported metadata, so the number of boxplots may be lower than the total number of declared de-hosting configurations. The full list of declared configurations and associated summary values is provided in Appendix Table 9.

<figure class="benchmark-figure">
<img src="figures/SARS2/dehosting_metric_boxplots_by_pipeline.png" alt="Distribution of percentage of host reads metrics by dehosting software version for SARS2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 13. Distribution of percentage of host reads by declared dehosting software version for SARS2.** Boxplots summarise sample-level percentage of host reads stratified by dehosting software version. Only configurations with evaluable percentage of host reads values are displayed, so some declared software categories may be absent from the plot. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each version, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.3. Preprocessing software

3 distinct pre-processing software configurations were reported for the SARS2 component.

Only pre-processing configurations with evaluable observations for the displayed metrics contribute to the figure, so some declared categories may not appear in the plot. The complete list of declared configurations and their summary values is provided in Appendix Table 10.

<figure class="benchmark-figure">
<img src="figures/SARS2/preprocessing_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pre-processing software configuration for SARS2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 14. Distribution of performance metrics by declared pre-processing software configuration for SARS2.** Multi-panel boxplots summarise sample-level performance stratified by pre-processing software. Panel A displays Number of reads sequenced and Panel B Reads passing filters. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.4. Mapping software

3 distinct mapping software configurations were reported for the SARS2 component.

The mapping boxplots include only configurations for which the relevant performance metrics were available, which means that fewer categories may be plotted than were originally declared. Full configuration-level summaries are reported in Appendix Table 11.

<figure class="benchmark-figure">
<img src="figures/SARS2/mapping_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by mapping software configuration for SARS2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 15. Distribution of performance metrics by declared mapping software configuration for SARS2.** Boxplots summarise sample-level performance stratified by mapping software. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.6. Consensus software

6 distinct consensus software configurations were reported for the SARS2 component.

Only consensus software configurations with sufficient evaluable data are visualised in the figure below, so the plotted set may be smaller than the total set of declarations. All declared configurations and their associated summary values can be reviewed in Appendix Table 12.

<figure class="benchmark-figure">
<img src="figures/SARS2/consensus_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by consensus software configuration for SARS2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 16. Distribution of performance metrics by declared consensus software configuration for SARS2.** Multi-panel boxplots summarise sample-level performance stratified by consensus software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers. For SARS1 and FLU1, Panel B uses a truncated y-axis to highlight differences among high-identity values.

#### 4.7. Variant calling software

8 distinct variant calling software configurations were reported for the SARS2 component.

The plotted variant calling categories correspond only to configurations with evaluable observations for the displayed metrics. Consequently, the figure may show fewer configurations than were declared overall; the complete summaries are listed in Appendix Table 13.

<figure class="benchmark-figure">
<img src="figures/SARS2/variant_calling_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by variant calling software configuration for SARS2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 17. Distribution of performance metrics by declared variant calling software configuration for SARS2.** Panel A is a stacked bar chart showing the number of evaluable samples assigned to each allele frequency reporting pattern for each software configuration. Boxplot Panel B displays discrepancies in reported variants with AF >=75% in the submitted VCF, Panel C discrepancies in reported variants with effect, Panel D successful hits, and Panel E total discrepancies. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. In the boxplots, the central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.8. Clade Assignment Software

8 distinct clade assignment software configurations were reported for the SARS2 component. For this category, configurations were counted as unique combinations of software name, software version, and clade assignment database version when available.

Because clade concordance could not be evaluated for every declared configuration, the boxplot includes only categories with usable observations. The full list of declared configurations and their summary values is available in Appendix Table 14.

<figure class="benchmark-figure landscape-benchmark-figure">
<img src="figures/SARS2/clade_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by clade assignment software configuration for SARS2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 18. Distribution of clade concordance by declared clade assignment software configuration for SARS2.** This boxplot summarises sample-level clade concordance stratified by clade assignment software configuration, where each configuration corresponds to a unique combination of software name, software version, and clade assignment database version when available. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.9. Lineage Assignment Software Name

8 distinct lineage assignment software configurations were reported for the SARS2 component. For this category, configurations were counted as unique combinations of software name, software version, and lineage assignment database version when available.

Lineage assignment configurations are shown only when concordance values were evaluable for the submitted observations, so the plotted categories may represent only a subset of the declarations. The complete configuration-level summary is provided in Appendix Table 15.

<figure class="benchmark-figure landscape-benchmark-figure">
<img src="figures/SARS2/lineage_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by lineage assignment software configuration for SARS2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 19. Distribution of lineage concordance by declared lineage assignment software configuration for SARS2.** This boxplot summarises sample-level lineage concordance stratified by lineage assignment software configuration, where each configuration corresponds to a unique combination of software name, software version, and lineage assignment database version when available. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

### 4.3. FLU1 (Influenza virus, Illumina)

This section presents an exploratory comparative analysis of declared workflow configurations within FLU1. Because laboratories differed in reference selection, software versions, parameterisation, reporting detail, and internal decision criteria, the results below should be interpreted as descriptive comparisons of observed performance patterns rather than as a controlled ranking of pipelines.

#### 4.1. Bioinformatics protocol

Based on metadata submissions, 6 distinct bioinformatics protocols were reported for the FLU1 component. These summaries compare declared workflow configurations as they were used in practice across participating laboratories.

<figure>
<img src="figures/FLU1/bioinformatics_protocol_discrepancies_boxplot_by_pipeline.png" alt="Distribution of consensus discrepancies by pipeline configuration for FLU1." style="width: 80%; max-width: 80%;"/>
</figure>

**Figure 20. Distribution of consensus discrepancies by declared pipeline configuration for FLU1.** This boxplot summarises sample-level consensus discrepancies stratified by bioinformatics protocol. The left y-axis shows discrepancy counts, while the right y-axis overlays lineage/type and clade classification accuracy for the same software configuration. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

**Table 3. Performance summary of declared bioinformatics protocols for FLU1.**

| Bioinformatics protocol | Version | N labs | Median genome identity (%) | Median discrepancies | Median metadata completeness (%) | Clade concordance (%) | Lineage/type concordance (%) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Custom pipeline/workflow |  | 5 | 96.07 | 21 | 60.7 | 43.8 | 100.0 |
| DRAGEN Targeted Microbial | 1.1.0 | 1 | 98.07 | 16.5 | 39.3 | 100.0 | 100.0 |
| INSaFLU | 2.2.2 | 1 | 92.94 | 178 | 41.1 | 75.0 | 75.0 |
| IRMA |  | 1 | 96.01 | 18 | 32.1 | 100.0 | 100.0 |
| IRMA | 1.2.0 | 3 | 95.70 | 33.5 | 39.3 | 58.3 | 91.7 |
| IRMA | 1.3.1 | 2 | 95.93 | 28 | 63.0 | 100.0 | 100.0 |

The observed differences across configurations should be read in the context of heterogeneous laboratory practices, including differences in reference choice, parameterisation, and thresholding. The table and figures therefore help identify recurrent performance patterns within FLU1, but they do not support a strong cross-laboratory ranking of bioinformatics protocols.

<figure class="benchmark-figure">
<img src="figures/FLU1/bioinformatics_protocol_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pipeline configuration for FLU1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 21. Distribution of performance metrics by declared pipeline configuration for FLU1.** Multi-panel boxplots summarise sample-level performance stratified by bioinformatics protocols. Panel A displays genome identity (%), Panel B metadata completeness (%), and Panel C exact classification concordance (%). X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Where required, Panel A uses a truncated y-axis to highlight differences among high-identity values. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.2. De-hosting software

5 distinct de-hosting software declarations were reported for the FLU1 component.

The distribution below reflects only configurations with evaluable percentage of host reads values in the reported metadata, so the number of boxplots may be lower than the total number of declared de-hosting configurations. The full list of declared configurations and associated summary values is provided in Appendix Table 16.

<figure class="benchmark-figure">
<img src="figures/FLU1/dehosting_metric_boxplots_by_pipeline.png" alt="Distribution of percentage of host reads metrics by dehosting software version for FLU1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 22. Distribution of percentage of host reads by declared dehosting software version for FLU1.** Boxplots summarise sample-level percentage of host reads stratified by dehosting software version. Only configurations with evaluable percentage of host reads values are displayed, so some declared software categories may be absent from the plot. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each version, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.3. Preprocessing software

6 distinct pre-processing software configurations were reported for the FLU1 component.

Only pre-processing configurations with evaluable observations for the displayed metrics contribute to the figure, so some declared categories may not appear in the plot. The complete list of declared configurations and their summary values is provided in Appendix Table 17.

<figure class="benchmark-figure">
<img src="figures/FLU1/preprocessing_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pre-processing software configuration for FLU1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 23. Distribution of performance metrics by declared pre-processing software configuration for FLU1.** Multi-panel boxplots summarise sample-level performance stratified by pre-processing software. Panel A displays Number of reads sequenced and Panel B Reads passing filters. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.4. Mapping software

6 distinct mapping software configurations were reported for the FLU1 component.

The mapping boxplots include only configurations for which the relevant performance metrics were available, which means that fewer categories may be plotted than were originally declared. Full configuration-level summaries are reported in Appendix Table 18.

<figure class="benchmark-figure">
<img src="figures/FLU1/mapping_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by mapping software configuration for FLU1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 24. Distribution of performance metrics by declared mapping software configuration for FLU1.** Boxplots summarise sample-level performance stratified by mapping software. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.5. Assembly software

4 distinct assembly software configurations were reported for the FLU1 component.

The assembly figures are restricted to configurations with evaluable values for the displayed metrics. As a result, some declared assembly categories may be absent from the plots; the full set of declared configurations and summary values is provided in Appendix Table 19.

<figure class="benchmark-figure">
<img src="figures/FLU1/assembly_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by assembly software configuration for FLU1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 25. Distribution of performance metrics by declared assembly software configuration for FLU1.** Multi-panel boxplots summarise sample-level performance stratified by assembly software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers. Panel B uses a truncated y-axis to highlight differences among high-identity values.

#### 4.6. Consensus software

7 distinct consensus software configurations were reported for the FLU1 component.

Only consensus software configurations with sufficient evaluable data are visualised in the figure below, so the plotted set may be smaller than the total set of declarations. All declared configurations and their associated summary values can be reviewed in Appendix Table 20.

<figure class="benchmark-figure">
<img src="figures/FLU1/consensus_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by consensus software configuration for FLU1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 26. Distribution of performance metrics by declared consensus software configuration for FLU1.** Multi-panel boxplots summarise sample-level performance stratified by consensus software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers. For SARS1 and FLU1, Panel B uses a truncated y-axis to highlight differences among high-identity values.

#### 4.7. Variant calling software

5 distinct variant calling software configurations were reported for the FLU1 component.

The plotted variant calling categories correspond only to configurations with evaluable observations for the displayed metrics. Consequently, the figure may show fewer configurations than were declared overall; the complete summaries are listed in Appendix Table 21.

<figure class="benchmark-figure">
<img src="figures/FLU1/variant_calling_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by variant calling software configuration for FLU1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 27. Distribution of performance metrics by declared variant calling software configuration for FLU1.** Panel A is a stacked bar chart showing the number of evaluable samples assigned to each allele frequency reporting pattern for each software configuration. Boxplot Panel B displays the number of reported variants with AF >=75%, Panel C the number of variants with AF >=75% in the submitted VCF, Panel D the number of variants with effect, Panel E metadata-VCF discrepancies, and Panel F the total number of variants present in the submitted VCF files. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. In the boxplots, the central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.8. Clade Assignment Software

9 distinct clade assignment software configurations were reported for the FLU1 component. For this category, configurations were counted as unique combinations of software name, software version, and clade assignment database version when available.

Because clade concordance could not be evaluated for every declared configuration, the boxplot includes only categories with usable observations. The full list of declared configurations and their summary values is available in Appendix Table 22.

<figure class="benchmark-figure landscape-benchmark-figure">
<img src="figures/FLU1/clade_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by clade assignment software configuration for FLU1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 28. Distribution of clade concordance by declared clade assignment software configuration for FLU1.** This boxplot summarises sample-level clade concordance stratified by clade assignment software configuration, where each configuration corresponds to a unique combination of software name, software version, and clade assignment database version when available. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.10. Type Assignment Software Name

6 distinct type assignment software configurations were reported for the FLU1 component. For this category, configurations were counted as unique combinations of software name, software version, and type assignment database version when available.

The type assignment plot is limited to configurations with evaluable concordance results, and therefore may contain fewer categories than the total number declared in metadata. The complete list of declared configurations and summary values is provided in Appendix Table 23.

<figure class="benchmark-figure landscape-benchmark-figure">
<img src="figures/FLU1/type_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by type assignment software configuration for FLU1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 29. Distribution of type concordance by declared type assignment software configuration for FLU1.** This boxplot summarises sample-level type concordance stratified by type assignment software configuration, where each configuration corresponds to a unique combination of software name, software version, and type assignment database version when available. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.11. Subtype Assignment Software Name

8 distinct subtype assignment software configurations were reported for the FLU1 component. For this category, configurations were counted as unique combinations of software name, software version, and subtype assignment database version when available.

Subtype assignment configurations are plotted only when evaluable concordance data were available, so some declared categories may not be represented in the figure. Appendix Table 24 contains the full list of declarations and their summary values.

<figure class="benchmark-figure landscape-benchmark-figure">
<img src="figures/FLU1/subtype_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by subtype assignment software configuration for FLU1." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 30. Distribution of subtype concordance by declared subtype assignment software configuration for FLU1.** This boxplot summarises sample-level subtype concordance stratified by subtype assignment software configuration, where each configuration corresponds to a unique combination of software name, software version, and subtype assignment database version when available. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

### 4.4. FLU2 (Influenza virus, Oxford Nanopore Technologies)

This section presents an exploratory comparative analysis of declared workflow configurations within FLU2. Because laboratories differed in reference selection, software versions, parameterisation, reporting detail, and internal decision criteria, the results below should be interpreted as descriptive comparisons of observed performance patterns rather than as a controlled ranking of pipelines.

#### 4.1. Bioinformatics protocol

Based on metadata submissions, 5 distinct bioinformatics protocols were reported for the FLU2 component. These summaries compare declared workflow configurations as they were used in practice across participating laboratories.

<figure>
<img src="figures/FLU2/bioinformatics_protocol_discrepancies_boxplot_by_pipeline.png" alt="Distribution of consensus discrepancies by pipeline configuration for FLU2." style="width: 80%; max-width: 80%;"/>
</figure>

**Figure 31. Distribution of consensus discrepancies by declared pipeline configuration for FLU2.** This boxplot summarises sample-level consensus discrepancies stratified by bioinformatics protocol. The left y-axis shows discrepancy counts, while the right y-axis overlays lineage/type and clade classification accuracy for the same software configuration. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

**Table 4. Performance summary of declared bioinformatics protocols for FLU2.**

| Bioinformatics protocol | Version | N labs | Median genome identity (%) | Median discrepancies | Median metadata completeness (%) | Clade concordance (%) | Lineage/type concordance (%) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Custom pipeline/workflow |  | 3 | 95.95 | 33.5 | 60.7 | 57.1 | 78.6 |
| INSaFLU | 2.2.2 | 1 | 95.45 | 52 | 39.3 | 80.0 | 80.0 |
| INSaFLU | Web version | 1 | 96.48 | 39 | 32.1 | 60.0 | 80.0 |
| IRMA | 1.2.0 | 1 | 95.21 | 262 | 87.5 | 80.0 | 100.0 |
| IRMA | 1.3.1 | 4 | 95.70 | 21 | 67.9 | 87.5 | 100.0 |

The observed differences across configurations should be read in the context of heterogeneous laboratory practices, including differences in reference choice, parameterisation, and thresholding. The table and figures therefore help identify recurrent performance patterns within FLU2, but they do not support a strong cross-laboratory ranking of bioinformatics protocols.

<figure class="benchmark-figure">
<img src="figures/FLU2/bioinformatics_protocol_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pipeline configuration for FLU2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 32. Distribution of performance metrics by declared pipeline configuration for FLU2.** Multi-panel boxplots summarise sample-level performance stratified by bioinformatics protocols. Panel A displays genome identity (%), Panel B metadata completeness (%), and Panel C exact classification concordance (%). X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Where required, Panel A uses a truncated y-axis to highlight differences among high-identity values. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.2. De-hosting software

3 distinct de-hosting software declarations were reported for the FLU2 component.

The distribution below reflects only configurations with evaluable percentage of host reads values in the reported metadata, so the number of boxplots may be lower than the total number of declared de-hosting configurations. The full list of declared configurations and associated summary values is provided in Appendix Table 25.

<figure class="benchmark-figure">
<img src="figures/FLU2/dehosting_metric_boxplots_by_pipeline.png" alt="Distribution of percentage of host reads metrics by dehosting software version for FLU2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 33. Distribution of percentage of host reads by declared dehosting software version for FLU2.** Boxplots summarise sample-level percentage of host reads stratified by dehosting software version. Only configurations with evaluable percentage of host reads values are displayed, so some declared software categories may be absent from the plot. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each version, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.3. Preprocessing software

4 distinct pre-processing software configurations were reported for the FLU2 component.

Only pre-processing configurations with evaluable observations for the displayed metrics contribute to the figure, so some declared categories may not appear in the plot. The complete list of declared configurations and their summary values is provided in Appendix Table 26.

<figure class="benchmark-figure">
<img src="figures/FLU2/preprocessing_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by pre-processing software configuration for FLU2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 34. Distribution of performance metrics by declared pre-processing software configuration for FLU2.** Multi-panel boxplots summarise sample-level performance stratified by pre-processing software. Panel A displays Number of reads sequenced and Panel B Reads passing filters. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.4. Mapping software

4 distinct mapping software configurations were reported for the FLU2 component.

The mapping boxplots include only configurations for which the relevant performance metrics were available, which means that fewer categories may be plotted than were originally declared. Full configuration-level summaries are reported in Appendix Table 27.

<figure class="benchmark-figure">
<img src="figures/FLU2/mapping_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by mapping software configuration for FLU2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 35. Distribution of performance metrics by declared mapping software configuration for FLU2.** Boxplots summarise sample-level performance stratified by mapping software. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.5. Assembly software

3 distinct assembly software configurations were reported for the FLU2 component.

The assembly figures are restricted to configurations with evaluable values for the displayed metrics. As a result, some declared assembly categories may be absent from the plots; the full set of declared configurations and summary values is provided in Appendix Table 28.

<figure class="benchmark-figure">
<img src="figures/FLU2/assembly_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by assembly software configuration for FLU2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 36. Distribution of performance metrics by declared assembly software configuration for FLU2.** Multi-panel boxplots summarise sample-level performance stratified by assembly software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers. Panel B uses a truncated y-axis to highlight differences among high-identity values.

#### 4.6. Consensus software

5 distinct consensus software configurations were reported for the FLU2 component.

Only consensus software configurations with sufficient evaluable data are visualised in the figure below, so the plotted set may be smaller than the total set of declarations. All declared configurations and their associated summary values can be reviewed in Appendix Table 29.

<figure class="benchmark-figure">
<img src="figures/FLU2/consensus_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by consensus software configuration for FLU2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 37. Distribution of performance metrics by declared consensus software configuration for FLU2.** Multi-panel boxplots summarise sample-level performance stratified by consensus software. Panel A displays consensus genome length, Panel B genome identity, and Panel C discrepancy counts. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers. For SARS1 and FLU1, Panel B uses a truncated y-axis to highlight differences among high-identity values.

#### 4.7. Variant calling software

5 distinct variant calling software configurations were reported for the FLU2 component.

The plotted variant calling categories correspond only to configurations with evaluable observations for the displayed metrics. Consequently, the figure may show fewer configurations than were declared overall; the complete summaries are listed in Appendix Table 30.

<figure class="benchmark-figure">
<img src="figures/FLU2/variant_calling_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by variant calling software configuration for FLU2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 38. Distribution of performance metrics by declared variant calling software configuration for FLU2.** Panel A is a stacked bar chart showing the number of evaluable samples assigned to each allele frequency reporting pattern for each software configuration. Boxplot Panel B displays the number of reported variants with AF >=75%, Panel C the number of variants with AF >=75% in the submitted VCF, Panel D the number of variants with effect, Panel E metadata-VCF discrepancies, and Panel F the total number of variants present in the submitted VCF files. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. Only panels with evaluable data are shown. In the boxplots, the central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.8. Clade Assignment Software

6 distinct clade assignment software configurations were reported for the FLU2 component. For this category, configurations were counted as unique combinations of software name, software version, and clade assignment database version when available.

Because clade concordance could not be evaluated for every declared configuration, the boxplot includes only categories with usable observations. The full list of declared configurations and their summary values is available in Appendix Table 31.

<figure class="benchmark-figure landscape-benchmark-figure">
<img src="figures/FLU2/clade_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by clade assignment software configuration for FLU2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 39. Distribution of clade concordance by declared clade assignment software configuration for FLU2.** This boxplot summarises sample-level clade concordance stratified by clade assignment software configuration, where each configuration corresponds to a unique combination of software name, software version, and clade assignment database version when available. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.10. Type Assignment Software Name

7 distinct type assignment software configurations were reported for the FLU2 component. For this category, configurations were counted as unique combinations of software name, software version, and type assignment database version when available.

The type assignment plot is limited to configurations with evaluable concordance results, and therefore may contain fewer categories than the total number declared in metadata. The complete list of declared configurations and summary values is provided in Appendix Table 32.

<figure class="benchmark-figure landscape-benchmark-figure">
<img src="figures/FLU2/type_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by type assignment software configuration for FLU2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 40. Distribution of type concordance by declared type assignment software configuration for FLU2.** This boxplot summarises sample-level type concordance stratified by type assignment software configuration, where each configuration corresponds to a unique combination of software name, software version, and type assignment database version when available. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

#### 4.11. Subtype Assignment Software Name

8 distinct subtype assignment software configurations were reported for the FLU2 component. For this category, configurations were counted as unique combinations of software name, software version, and subtype assignment database version when available.

Subtype assignment configurations are plotted only when evaluable concordance data were available, so some declared categories may not be represented in the figure. Appendix Table 33 contains the full list of declarations and their summary values.

<figure class="benchmark-figure landscape-benchmark-figure">
<img src="figures/FLU2/subtype_assignment_metric_boxplots_by_pipeline.png" alt="Distribution of performance metrics by subtype assignment software configuration for FLU2." style="width: 96%; max-width: 96%;"/>
</figure>

**Figure 41. Distribution of subtype concordance by declared subtype assignment software configuration for FLU2.** This boxplot summarises sample-level subtype concordance stratified by subtype assignment software configuration, where each configuration corresponds to a unique combination of software name, software version, and subtype assignment database version when available. X-axis labels report the declared software configuration and the number of laboratories (`n`) contributing observations to each category. The central line indicates the median, boxes represent the interquartile range, whiskers denote the full observed range of sample-level observations across participating laboratories using each configuration, translucent points correspond to individual sample-level observations submitted by participating laboratories, and hollow circles beyond the whiskers indicate outliers.

## 5. Discussion

The benchmarking results show that observed workflow performance varied across analytical components and was associated with differences in software configuration, reference strategies, parameterisation, and reporting quality. The results do not point to a single universally best-performing pipeline. Instead, they highlight the importance of considering performance together with the consistency of the results across laboratories and the information available to describe the workflow.

The metadata confirm that a diverse range of analytical configurations is currently in use across the RELECOV network. This diversity is analytically valuable, as it allows different approaches to be compared within the same benchmarking exercise. At the same time, interpretation is limited by incomplete reporting. Only 60.1% of software-version fields were completed, while minimum coverage thresholds, variant calling parameters, and reference genome identifiers were reported for 45.6%, 38.9%, and 57.7% of submitted samples, respectively. This means that some plausible explanations for observed performance differences can only be considered as contributing factors rather than demonstrated causal effects.

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

### SARS1 (SARS-CoV-2, Illumina)

<h4 class="appendix-landscape-heading">Pipeline Benchmarking and Comparative Performance Supplementary Material</h4>

<h5 class="appendix-landscape-heading">De-hosting</h5>
**Appendix Table 1. Performance summary of declared de-hosting software for SARS1.**

| De-hosting software | Version | N labs | Median % host reads |
|---|---:|---:|---:|
| BWA-Mem | V2.2.1 | 1 | N/A |
| Kraken2 | 2.1.6 | 5 | 0 |
| Kraken2 | DRAGEN Microbial Amplicon ad hoc version | 1 | N/A |
| bowtie2 | 2.5.1 | 1 | 0.12 |

<h5 class="appendix-landscape-heading">Pre-processing</h5>
**Appendix Table 2. Performance summary of declared pre-processing software configurations for SARS1.**

| Pre-processing software | Version | N labs | Most common configuration | Median number of reads sequenced | Median reads passing filters |
|---|---:|---:|---:|---:|---:|
| Fastp | 0.20.1 | 1 | --detect_adapter_for_pe --cut_tail --cut_window_size 10 --cut_mean_quality 20 --length_required 35 --json json_file --html html_file --thread 34 | 3024406 | 2423968 |
| Fastp | 0.24.0 | 1 | N/A | 113415910 | N/A |
| Fastp | 1.0.1 | 6 | --detect_adapter_for_pe --cut_front --cut_tail --trim_poly_x --cut_mean_quality 30 --qualified_quality_phred 30 --unqualified_percent_limit 10 --length_required 50 | 3024406 | 2115268 |
| Trimmomatic | 0.39 | 2 | LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30 | 3024406 | 2950068 |
| Trimmomatic | 0.39-2 | 1 | N/A | N/A | N/A |
| Trimmomatic | 0.40 | 1 | N/A | 3024406 | 2631478 |

<h5 class="appendix-landscape-heading">Mapping</h5>
**Appendix Table 3. Performance summary of declared mapping software configurations for SARS1.**

| Mapping software | Version | N labs | Most common configuration | Most common depth of coverage threshold | Median % reads virus |
|---|---:|---:|---:|---:|---:|
| bowtie2 | 2.5.4 | 6 | --local --very-sensitive-local --seed 1 | 10x | 99.955 |
| bwa mem | 0.7.17-r1188 | 1 | -Y -M -t 34 | 30x | 99.64 |
| bwa mem | 0.7.19 | 1 | N/A | N/A | N/A |
| bwa mem | 0.7.19-r1273 | 1 | -p -Y -v 2 | 10x | N/A |
| bwa mem | V2.2.1 | 1 | -t 4 covid-ba2 | 10X | N/A |
| minimap2 | DRAGEN Microbial Amplicon ad-hoc version | 1 | N/A | N/A | N/A |

<h5 class="appendix-landscape-heading">Assembly</h5>
**Appendix Table 4. Performance summary of declared assembly software configurations for SARS1.**

| Assembly software | Version | N labs | Most common configuration | Median consensus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
| Megahit | DRAGEN Microbial Amplicon ad hoc version | 1 | N/A | N/A | 99.5697 | 2.5 |
| SPAdes | 4.0.0 | 1 | --rnaviral --threads 24 | 29808 | 99.6049 | 5 |
| SPAdes | 4.1.0 | 1 | N/A | 29801 | 99.6905 | 6 |

<h5 class="appendix-landscape-heading">Consensus software</h5>
**Appendix Table 5. Performance summary of declared consensus software configurations for SARS1.**

| Consensus software | Version | N labs | Most common configuration | Median consensus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
| DRAGEN consensus | DRAGEN Microbial Amplicon ad hoc version | 2 | N/A | N/A | 99.5697 | 2.5 |
| Ivar consensus |  | 1 | N/A | N/A | 99.7287 | 2.5 |
| Ivar consensus | 1.4.2 | 1 | -t 0.01 | 29903 | 99.6584 | 5 |
| Ivar consensus | 1.4.3 | 1 | -q 20 -t 0.8 -m 20 -n N | 29869 | 99.7388 | 3 |
| Ivar consensus | 1.4.4 | 2 | N/A | N/A | 99.6385 | 5 |
| bcftools consensus | 1.2.2 | 1 | -m | 29812 | 99.6938 | 6 |
| bcftools consensus | 1.21 | 1 | mosdepth -x -Q 1 \| awk '$4<=10' | 29.868 | 99.6367 | 5 |
| bcftools consensus | 1.22 | 4 | AF > 0.75 | 29885 | 99.6938 | 6 |

<h5 class="appendix-landscape-heading">Variant calling</h5>
**Appendix Table 6. Performance summary of declared variant calling software configurations for SARS1.**

| Variant calling software | Version | N labs | Most common configuration | Median high and low frequency (%) | Median high frequency only (%) | Median low frequency only (%) | Median variants (AF >=75%) | Median variants in VCF (AF >=75%) | Median variants with effect | Median variants with effect in VCF | Median metadata-VCF discrepancies | Median effect discrepancies | Median successful hits | Median total discrepancies |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| DRAGEN VariantCaller | DRAGEN Microbial Amplicon ad hoc version | 2 | N/A | 0 | 100 | 0 | N/A | 75 | N/A | 52 | N/A | N/A | 69.5 | 2.5 |
| Ivar | 1.22 | 1 | N/A | 100 | 0 | 0 | N/A | 61 | N/A | 45 | N/A | N/A | 65 | 1 |
| Ivar | 1.4.2 | 1 | -t 0.03 -m 10 | 100 | 0 | 0 | 65 | 51 | N/A | 37 | 3 | N/A | 62 | 9 |
| Ivar | 1.4.3 | 1 | -q 15 -t 0.01 -m 1 -r NC_045512.2.fasta -g NC_045512.2.gff3 | 100 | 0 | 0 | 61 | 76 | 47 | N/A | 4 | N/A | 65 | 56 |
| Ivar | 1.4.4 | 6 | -t 0.25 -q 20 -m 10 | 100 | 0 | 0 | 61 | 61 | 46.5 | 45 | 0 | 2 | 65 | 1 |
| Octopus | 0.7.4 | 1 | -P 1 | 0 | 100 | 0 | N/A | 65 | N/A | 46 | N/A | N/A | 64 | 4 |
| freebayes | 1.1.0 | 1 | N/A | 100 | 0 | 0 | N/A | 57 | N/A | 40 | N/A | N/A | 63.5 | 41.5 |
| lofreq | 2.1.5 | 1 | indelqual --dindel // call-parallel --pp-threads 64 --call-indels | 100 | 0 | 0 | N/A | 51 | N/A | N/A | N/A | N/A | N/A | N/A |

<h5 class="appendix-landscape-heading">Clade assignment</h5>
**Appendix Table 7. Performance summary of declared clade assignment software configurations for SARS1.**

| Clade assignment software | Version | N labs | Database version | % of clade match |
|---|---:|---:|---:|---:|
| Nextclade | 2.12.0 | 1 | 2024-10-17--16-48-48Z | 100 |
| Nextclade | 3.11.0 | 1 | N/A | 0 |
| Nextclade | 3.13.1 | 1 | 2026-01-06--14-59-32Z | 0 |
| Nextclade | 3.18.0 | 1 | 2026-01-14--19-24-43Z | 100 |
| Nextclade | 3.18.1 | 1 | N/A | 100 |
| Nextclade | 3.18.1 | 6 | 2026-01-06--14-59-32Z | 83.33 |
| Nextclade | 3.18.1 | 1 | 2026-01-14--19-24-43Z | 100 |
| Nextclade | 3.18.1 | 1 | web (latest) | 100 |
| Nextclade | 3.8.1 | 1 | 2026-01-06--14-59-32Z | 100 |
| Nextclade | 3.9.0 | 1 | 2026-01-06--14-59-32Z | 0 |

<h5 class="appendix-landscape-heading">Lineage assignment</h5>
**Appendix Table 8. Performance summary of declared lineage assignment software configurations for SARS1.**

| Lineage Assignment software | Version | N labs | Database version | % of lineage match |
|---|---:|---:|---:|---:|
| Pangolin | 4.3.1 | 1 | 1.33 | 100 |
| Pangolin | 4.3.1 | 3 | 1.37 | 100 |
| Pangolin | 4.3.1 | 2 | 4.3.1 | 100 |
| Pangolin | 4.3.4 | 3 | 1.37 | 58.33 |
| Pangolin | 4.3.4 | 1 | 1.37.1 | 100 |
| Pangolin | 4.3.4 | 4 | 4.3.4 | 100 |

### SARS2 (SARS-CoV-2, Oxford Nanopore Technologies)

<h4 class="appendix-landscape-heading">Pipeline Benchmarking and Comparative Performance Supplementary Material</h4>

<h5 class="appendix-landscape-heading">De-hosting</h5>
**Appendix Table 9. Performance summary of declared de-hosting software for SARS2.**

| De-hosting software | Version | N labs | Median % host reads |
|---|---:|---:|---:|
| Kraken2 | 2.1.6 | 4 | 0.01 |

<h5 class="appendix-landscape-heading">Pre-processing</h5>
**Appendix Table 10. Performance summary of declared pre-processing software configurations for SARS2.**

| Pre-processing software | Version | N labs | Most common configuration | Median number of reads sequenced | Median reads passing filters |
|---|---:|---:|---:|---:|---:|
| Fastp | 1.0.1 | 1 | N/A | N/A | N/A |
| NanoFilt | 2.6.0 | 1 | -q 8 -l 50 --headcrop 30 --tailcrop 30 --maxlength 50000 | N/A | N/A |
| Porechop | 0.3.2pre | 1 | --discard_middle | 71109 | 56514 |

<h5 class="appendix-landscape-heading">Mapping</h5>
**Appendix Table 11. Performance summary of declared mapping software configurations for SARS2.**

| Mapping software | Version | N labs | Most common configuration | Most common depth of coverage threshold | Median % reads virus |
|---|---:|---:|---:|---:|---:|
| minimap2 | 2.28-r1209 | 3 | -a -x map-ont -t 12 | 20x | 92.16 |
| minimap2 | 2.30 | 1 | N/A | N/A | N/A |
| minimap2 | 2.30-r1287 | 3 | -x map-ont | 10x | 47.6089 |

<h5 class="appendix-landscape-heading">Consensus software</h5>
**Appendix Table 12. Performance summary of declared consensus software configurations for SARS2.**

| Consensus software | Version | N labs | Most common configuration | Median consensus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
| Medaka consensus | 2.2.0 | 1 | --model r941_min_sup_g507 --chunk_len 800 --chunk_ovlp 400 | N/A | 99.4591 | 13 |
| Medaka consensus | 2.2.1 | 1 | -m r941_min_hac_variant_g507 -r N | N/A | 99.7823 | 31 |
| bcftools consensus | 1.17 | 1 | --min-depth 20 | 29856 | 99.5049 | 2 |
| bcftools consensus | 1.2.2 | 1 | -m | 29260 | 99.8896 | 3 |
| bcftools consensus | 1.22 | 3 | --min-depth 20 | 29898.5 | 99.4647 | 16 |
| bcftools consensus | 2.5 | 1 | N/A | N/A | 99.8794 | 11 |

<h5 class="appendix-landscape-heading">Variant calling</h5>
**Appendix Table 13. Performance summary of declared variant calling software configurations for SARS2.**

| Variant calling software | Version | N labs | Model used | Median high and low frequency (%) | Median high frequency only (%) | Median low frequency only (%) | Median variants (AF >=75%) | Median variants in VCF (AF >=75%) | Median variants with effect | Median variants with effect in VCF | Median metadata-VCF discrepancies | Median effect discrepancies | Median successful hits | Median total discrepancies |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Clair3 | 1.0.10 | 1 | r941_prom_hac_g360+g422 | 0 | 100 | 0 | 94 | 94 | 66 | 66 | 0 | 0 | 94 | 0 |
| Clair3 | 1.0.11 | 1 | r941_prom_hac_g360+g422 | 0 | 100 | 0 | 72 | 72 | 46 | 47 | 0 | 1 | 72 | 22 |
| Clair3 | 1.2.0 | 1 | r1041_e82_400bps_sup_g615 | 0 | 100 | 0 | N/A | 156 | N/A | 101 | N/A | N/A | 87 | 1 |
| Clair3 | 1.6.2 | 1 | N/A | 0 | 100 | 0 | N/A | 176 | N/A | 123 | N/A | N/A | 94 | 0 |
| Clair3 | artic_1.8.5 | 1 | r1041_e82_400bps_hac_v430 | 0 | 100 | 0 | 93 | 93 | 71 | 66 | 0 | 5 | 93 | 1 |
| Medaka | 1.2.1 | 1 | N/A | 100 | 0 | 0 | N/A | 75 | N/A | 52 | N/A | N/A | 94 | 1.5 |
| Medaka | 2.2.0 | 1 | r941_min_sup_variant_g507 | 100 | 0 | 0 | N/A | N/A | N/A | N/A | N/A | N/A | 91 | 23 |
| Medaka | 2.2.1 | 1 | r941_min_hac_variant_g507 | N/A | N/A | N/A | N/A | 0 | N/A | 0 | N/A | N/A | 91 | 29 |

<h5 class="appendix-landscape-heading">Clade assignment</h5>
**Appendix Table 14. Performance summary of declared clade assignment software configurations for SARS2.**

| Clade assignment software | Version | N labs | Database version | % of clade match |
|---|---:|---:|---:|---:|
| Nextclade | 3.11.0 | 1 | N/A | 0 |
| Nextclade | 3.18.0 | 1 | 2026-01-14--19-24-43Z | 100 |
| Nextclade | 3.18.1 | 5 | 2026-01-06--14-59-32Z | 100 |
| Nextclade | 3.18.1 | 1 | 3.0.0 | 100 |
| Nextclade | 3.18.1 | 1 | web (latest) | 100 |
| Nextclade | 3.18.2 | 1 | 3.0.0 | 100 |
| Nextclade | 3.18.3 | 1 | 3.0.0 | N/A |
| Nextclade | 3.8.1 | 1 | 2026-01-06--14-59-32Z | 100 |

<h5 class="appendix-landscape-heading">Lineage assignment</h5>
**Appendix Table 15. Performance summary of declared lineage assignment software configurations for SARS2.**

| Lineage Assignment software | Version | N labs | Database version | % of lineage match |
|---|---:|---:|---:|---:|
| INSaFLU | 2.2.2 | 1 | N/A | N/A |
| INSaFLU | 2.2.2 | 1 | 4.3.4 | 66.67 |
| Nextclade | 3.18.1 | 1 | 1.37 | 100 |
| Nextclade | 3.18.1 | 1 | 4.3.4 | 66.67 |
| Pangolin | 4.3.1 | 1 | 1.37 | 66.67 |
| Pangolin | 4.3.1 | 1 | 4.3.1 | 33.33 |
| Pangolin | 4.3.4 | 3 | 1.37 | 44.44 |
| Pangolin | 4.3.4 | 1 | 1.37.1 | 66.67 |

### FLU1 (Influenza virus, Illumina)

<h4 class="appendix-landscape-heading">Pipeline Benchmarking and Comparative Performance Supplementary Material</h4>

<h5 class="appendix-landscape-heading">De-hosting</h5>
**Appendix Table 16. Performance summary of declared de-hosting software for FLU1.**

| De-hosting software | Version | N labs | Median % host reads |
|---|---:|---:|---:|
| Kraken2 | 2.1.3 | 1 | 25.04 |
| Kraken2 | 2.1.6 | 2 | N/A |
| Kraken2 | DRAGEN Microbial Amplicon | 1 | N/A |
| NCBI Human Read Scrubber | 2.0.0 | 1 | 23.69 |
| bowtie2 | 2.5.1 | 1 | 10.24 |

<h5 class="appendix-landscape-heading">Pre-processing</h5>
**Appendix Table 17. Performance summary of declared pre-processing software configurations for FLU1.**

| Pre-processing software | Version | N labs | Most common configuration | Median number of reads sequenced | Median reads passing filters |
|---|---:|---:|---:|---:|---:|
| Custom preprocessing script |  | 1 | FaQCs --ascii 33 --min_L 50 --avg_q 30 | 1841284 | 1060538 |
| Fastp | 0.20.0 | 1 | --cut_front --cut_tail --cut_mean_quality 15 --qualified_quality_phred 15 --trim_poly_x --length_required 50 --detect_adapter_for_pe | 1841284 | 1841260 |
| Fastp | 0.24.0 | 3 | --detect_adapter_for_pe --correction --length_required 45 --qualified_quality_phred 25 --average_qual 30 --cut_front --cut_tail --trim_front1 5 --trim_tail1 5 --trim_front2 5 --trim_tail2 5 | 1841359 | 1295030 |
| Fastp | 1.0.1 | 2 | N/A | N/A | N/A |
| IRMA custom script |  | 1 | N/A | N/A | N/A |
| Trimmomatic | 0.39 | 1 | LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30 | 1841284 | 535168 |

<h5 class="appendix-landscape-heading">Mapping</h5>
**Appendix Table 18. Performance summary of declared mapping software configurations for FLU1.**

| Mapping software | Version | N labs | Most common configuration | Most common depth of coverage threshold | Median % reads virus |
|---|---:|---:|---:|---:|---:|
| BLAT | IRMA 1.3.0-rc1 | 1 | Default (IRMA module FLU_AD) | 1x | 25.97 |
| BLAT | v.36 | 1 | N/A | N/A | N/A |
| BLAT | v35 | 2 | N/A | 10x | 63.72 |
| bwa mem | 0.7.19-r1273 | 3 | -p -Y -v 2 | 10x | N/A |
| minimap2 | 2.26-r1175 | 1 | -ax sr | 20x | N/A |
| minimap2 | DRAGEN Microbial Amplicon ad-hoc version | 1 | N/A | N/A | N/A |

<h5 class="appendix-landscape-heading">Assembly</h5>
**Appendix Table 19. Performance summary of declared assembly software configurations for FLU1.**

| Assembly software | Version | N labs | Most common configuration | Median consensus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
| LABEL | IRMA 1.3.0-rc1 | 1 | Default (IRMA module FLU_AD) | 13139 | 95.9825 | 18.5 |
| LABEL | v0.6.5 | 1 | N/A | 13134.5 | 95.5875 | 50 |
| Megahit | DRAGEN Microbial Amplicon ad hoc version | 1 | N/A | N/A | 98.0694 | 16.5 |
| SPAdes |  | 1 | N/A | N/A | 95.9825 | 18 |

<h5 class="appendix-landscape-heading">Consensus software</h5>
**Appendix Table 20. Performance summary of declared consensus software configurations for FLU1.**

| Consensus software | Version | N labs | Most common configuration | Median consensus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
| IRMA consensus |  | 1 | N/A | N/A | 95.566 | 48.5 |
| IRMA consensus | 1.2.0 | 2 | MIN_AMBIG=0.25 MIN_CONS_SUPPORT=10 | 13134.5 | 95.8439 | 27 |
| IRMA consensus | 1.3.1 | 1 | MIN_AMBIG=0.75; MIN_CONS_SUPPORT=9 | 13139 | 96.014 | 21.5 |
| Ivar consensus |  | 1 | N/A | N/A | 92.9431 | 178 |
| Ivar consensus | 1.4.3 | 1 | N/A | N/A | 97.7187 | 23 |
| bcftools consensus | 1.17 | 1 | AF > 0.50 | 13136 | 96.3999 | 64 |
| bcftools consensus | 1.21 | 1 | mosdepth -x -Q 1 \| awk '$4<=10' | 13.1835 | 96.5581 | 18 |

<h5 class="appendix-landscape-heading">Variant calling</h5>
**Appendix Table 21. Performance summary of declared variant calling software configurations for FLU1.**

| Variant calling software | Version | N labs | Most common configuration | Median high and low frequency (%) | Median high frequency only (%) | Median low frequency only (%) | Median variants (AF >=75%) | Median variants in VCF (AF >=75%) | Median variants with effect | Median metadata-VCF discrepancies | Median total variants in VCF |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Custom variant calling script |  | 1 | -f 0.1 -d 10 -t 10 | 100 | 0 | 0 | 471 | 465.5 | 104.5 | 3.5 | 721 |
| IRMA custom VariantCaller |  | 3 | N/A | 0 | 0 | 100 | 987 | 0 | 165 | 987 | 166 |
| Ivar | 1.4.3 | 1 | N/A | 100 | 0 | 0 | N/A | 591 | N/A | N/A | 623 |
| Octopus | 0.7.4 | 1 | -P 1 | 0 | 100 | 0 | N/A | 233.5 | N/A | N/A | 233.5 |
| lofreq | 2.1.5 | 1 | indelqual --dindel // call-parallel --pp-threads 64 --call-indels | 100 | 0 | 0 | N/A | 233 | N/A | N/A | 528.5 |

<h5 class="appendix-landscape-heading">Clade assignment</h5>
**Appendix Table 22. Performance summary of declared clade assignment software configurations for FLU1.**

| Clade assignment software | Version | N labs | Database version | % of clade match |
|---|---:|---:|---:|---:|
| Nextclade | 3.13.1 | 1 | N/A | 0 |
| Nextclade | 3.13.1 | 1 | 2026-01-14--08-53-00Z | 100 |
| Nextclade | 3.13.1 | 1 | 2026-01-14--19-24-43Z | 100 |
| Nextclade | 3.18.1 | 2 | N/A | 75 |
| Nextclade | 3.18.1 | 4 | 2025-09-09--12-13-13Z | 100 |
| Nextclade | 3.18.1 | 5 | 2026-01-14--19-24-43Z | 100 |
| Nextclade | 3.18.1 | 1 | web (latest) | 100 |
| Nextclade | 3.9.1 | 1 | 2025-09-09--12-13-13Z | 100 |
| Nextclade | 3.9.1 | 1 | 2026-01-14--19-24-43Z | 100 |

<h5 class="appendix-landscape-heading">Type assignment</h5>
**Appendix Table 23. Performance summary of declared type assignment software configurations for FLU1.**

| Type Assignment software | Version | N labs | Database version | % of type match |
|---|---:|---:|---:|---:|
| Custom typing script |  | 4 | v2 | 93.75 |
| INSaFLU |  | 1 | 2.2.2 | 75 |
| IRMA typing |  | 2 | N/A | 100 |
| IRMA typing |  | 1 | 1 | 100 |
| IRMA typing |  | 1 | 1.2.0 | 100 |
| IRMA typing |  | 1 | influenza A & B, genome sets version 2 | 100 |

<h5 class="appendix-landscape-heading">Subtype assignment</h5>
**Appendix Table 24. Performance summary of declared subtype assignment software configurations for FLU1.**

| Subtype Assignment software | Version | N labs | Database version | % of subtype match |
|---|---:|---:|---:|---:|
| Custom subtyping script |  | 3 | v2 | 91.67 |
| INSaFLU | 2.2.2 | 1 | 2.2.2 | 75 |
| IRMA subtyping |  | 1 | influenza A & B, genome sets version 2 | 100 |
| IRMA subtyping | 1.2.0 | 1 | N/A | 100 |
| IRMA subtyping | 1.2.0 | 1 | 1 | 100 |
| IRMA subtyping | 1.3.0-rc1 | 1 | 1.3.0-rc1 | 100 |
| IRMA subtyping | 1.3.1 | 1 | N/A | 100 |
| MASH | 2.3 | 1 | 2.3 | 100 |

### FLU2 (Influenza virus, Oxford Nanopore Technologies)

<h4 class="appendix-landscape-heading">Pipeline Benchmarking and Comparative Performance Supplementary Material</h4>

<h5 class="appendix-landscape-heading">De-hosting</h5>
**Appendix Table 25. Performance summary of declared de-hosting software for FLU2.**

| De-hosting software | Version | N labs | Median % host reads |
|---|---:|---:|---:|
| Kraken2 | 2.1.2 | 1 | 5.8 |
| Kraken2 | 2.1.3 | 1 | 2.7289 |
| Kraken2 | 2.1.6 | 1 | N/A |

<h5 class="appendix-landscape-heading">Pre-processing</h5>
**Appendix Table 26. Performance summary of declared pre-processing software configurations for FLU2.**

| Pre-processing software | Version | N labs | Most common configuration | Median number of reads sequenced | Median reads passing filters |
|---|---:|---:|---:|---:|---:|
| FiltLong | 0.2.1 | 1 | --min_length 100 --keep_percent 98 | 20489 | 17217 |
| IRMA custom script |  | 2 | Default IRMA | 20489 | 17942 |
| NanoFilt | 2.6.0 | 1 | -q 8 -l 50 --headcrop 30 --tailcrop 30 --maxlength 50000 | N/A | N/A |
| Porechop | v0.2.4 | 2 | -q 7 -l 300 --maxlength 2500 | 20489 | 14129.5 |

<h5 class="appendix-landscape-heading">Mapping</h5>
**Appendix Table 27. Performance summary of declared mapping software configurations for FLU2.**

| Mapping software | Version | N labs | Most common configuration | Most common depth of coverage threshold | Median % reads virus |
|---|---:|---:|---:|---:|---:|
| BLAT |  | 1 | Default IRMA | Default IRMA | 77.6 |
| BLAT | IRMA 1.3.0-rc1 | 1 | Default (IRMA module FLU_AD) | 1x | 33.57 |
| BLAT | v35 | 2 | N/A | 50x | 78.9 |
| minimap2 | 2.30-r1287 | 2 | -ax map-ont | N/A | N/A |

<h5 class="appendix-landscape-heading">Assembly</h5>
**Appendix Table 28. Performance summary of declared assembly software configurations for FLU2.**

| Assembly software | Version | N labs | Most common configuration | Median consensus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
| Flye | 2.9.6-b1802 | 1 | --nano-raw --iterations 2 | N/A | 96.0597 | 18 |
| LABEL | IRMA 1.3.0-rc1 | 1 | Default (IRMA module FLU_AD) | 13132 | 95.8353 | 18 |
| LABEL | v0.6.5 | 1 | N/A | 13123 | 95.2088 | 262 |

<h5 class="appendix-landscape-heading">Consensus software</h5>
**Appendix Table 29. Performance summary of declared consensus software configurations for FLU2.**

| Consensus software | Version | N labs | Most common configuration | Median consensus genome length | Median genome identity | Median number of discrepancies per sample |
|---|---:|---:|---:|---:|---:|---:|
| IRMA consensus |  | 2 | Default IRMA | N/A | 95.5613 | 27 |
| IRMA consensus | 1.2.0 | 1 | MIN_AMBIG=0.25 MIN_CONS_QUALITY=10 MIN_CONS_SUPPORT=20 QUAL_THRESHOLD=10 | 13123 | 95.2088 | 262 |
| IRMA consensus | 1.3.1 | 2 | MIN_AMBIG=0.75; MIN_CONS_SUPPORT=9 | 13123 | 95.8173 | 18 |
| Medaka consensus | 2.2.1 | 1 | -m r941_min_hac_variant_g507 -r N | N/A | 96.9448 | 33.5 |
| bcftools consensus | 1.23 | 1 | N/A | N/A | 95.4243 | 524 |

<h5 class="appendix-landscape-heading">Variant calling</h5>
**Appendix Table 30. Performance summary of declared variant calling software configurations for FLU2.**

| Variant calling software | Version | N labs | Most common configuration | Median high and low frequency (%) | Median high frequency only (%) | Median low frequency only (%) | Median variants (AF >=75%) | Median variants in VCF (AF >=75%) | Median variants with effect | Median metadata-VCF discrepancies | Median total variants in VCF |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Custom variant calling script |  | 1 | -f 0.1 -d 10 -t 10 | 100 | 0 | 0 | 463.5 | 419 | 95 | 18 | 1431 |
| IRMA custom VariantCaller |  | 3 | N/A | 0 | 0 | 100 | 462 | 0 | 107 | 462 | 4840 |
| Medaka |  | 1 | default MEDAKA | N/A | N/A | N/A | 475 | N/A | 130 | N/A | 1135 |
| Medaka | 2.2.1 | 1 | -m r941_min_hac_variant_g507 | N/A | N/A | N/A | N/A | N/A | N/A | N/A | 564 |
| bcftools | 1.23 | 1 | N/A | 100 | 0 | 0 | N/A | 481 | N/A | N/A | 804.5 |

<h5 class="appendix-landscape-heading">Clade assignment</h5>
**Appendix Table 31. Performance summary of declared clade assignment software configurations for FLU2.**

| Clade assignment software | Version | N labs | Database version | % of clade match |
|---|---:|---:|---:|---:|
| Nextclade | 3.18.1 | 2 | N/A | 66.67 |
| Nextclade | 3.18.1 | 3 | 2025-09-09--12-13-13Z | 100 |
| Nextclade | 3.18.1 | 5 | 2026-01-14--19-24-43Z | 94.12 |
| Nextclade | 3.18.1 | 1 | web (latest) | 75 |
| Nextclade | 3.9.1 | 1 | 2025-09-09--12-13-13Z | 100 |
| Nextclade | 3.9.1 | 1 | 2026-01-14--19-24-43Z | 100 |

<h5 class="appendix-landscape-heading">Type assignment</h5>
**Appendix Table 32. Performance summary of declared type assignment software configurations for FLU2.**

| Type Assignment software | Version | N labs | Database version | % of type match |
|---|---:|---:|---:|---:|
| ABRicate |  | 1 | 2.2.2 | 100 |
| Custom typing script |  | 1 | N/A | 40 |
| INSaFLU |  | 1 | 2.2.2 | 80 |
| IRMA typing |  | 1 | N/A | 100 |
| IRMA typing |  | 1 | 1 | 100 |
| IRMA typing |  | 1 | 1.3.1 | 100 |
| IRMA typing |  | 1 | influenza A & B, genome sets version 2 | 100 |

<h5 class="appendix-landscape-heading">Subtype assignment</h5>
**Appendix Table 33. Performance summary of declared subtype assignment software configurations for FLU2.**

| Subtype Assignment software | Version | N labs | Database version | % of subtype match |
|---|---:|---:|---:|---:|
| ABRicate | 1.0.1 | 1 | 2.2.2 | 100 |
| Custom subtyping script |  | 1 | N/A | 40 |
| INSaFLU | 2.2.2 | 1 | 2.2.2 | 80 |
| IRMA subtyping |  | 1 | influenza A & B, genome sets version 2 | 100 |
| IRMA subtyping | 1.2.0 | 1 | 1 | 100 |
| IRMA subtyping | 1.3.0-rc1 | 1 | 1.3.0-rc1 | 100 |
| IRMA subtyping | 1.3.1 | 1 | N/A | 100 |
| IRMA subtyping | 1.3.1 | 1 | 1.3.1 | 100 |

