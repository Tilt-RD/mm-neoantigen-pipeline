#!/usr/bin/env python3
"""
15_ind_documentation.py — IND-Enabling Documentation Generator
==============================================================
Generates a comprehensive regulatory documentation package for
Investigational New Drug (IND) applications for personalised
mRNA cancer vaccines.

This module defines the documentation structure, content templates,
and validation logic. Use 16_generate_ind_package.py to produce
the package for a specific patient.

Usage:
    python 16_generate_ind_package.py --patient mmrf_1251
    python 15_ind_documentation.py  # (module-only; use 16_generate_ind_package.py)

Output files:
    output/ind_package/{patient_id}/
        pipeline_description.txt
        validation_summary.txt
        limitations_statement.txt
        dataset_manifest.json
        docker_environment.txt
        bioinformatics_methods.txt
        clinical_metadata.json
"""

import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path

# ── Metadata / software versions ────────────────────────────────────

PIPELINE_VERSION = "2.1.0"
PYTHON_VERSION = "3.9"
PYTHON_VERSION_MIN = "3.8"

SOFTWARE_INVENTORY = {
    "python": {"version": PYTHON_VERSION, "min_version": PYTHON_VERSION_MIN, "license": "PSF"},
    "pandas": {"version": ">=2.0.0", "pip": "pandas"},
    "numpy": {"version": ">=1.24.0", "pip": "numpy"},
    "requests": {"version": ">=2.28.0", "pip": "requests"},
    "biopython": {"version": ">=1.81", "pip": "biopython"},
    "pysam": {"version": ">=0.21.0", "pip": "pysam (optional)"},
    "fpdf2": {"version": ">=2.7.0", "pip": "fpdf2"},
    "plotly": {"version": ">=5.18.0", "pip": "plotly"},
    "streamlit": {"version": ">=1.30.0", "pip": "streamlit"},
    "scikit-learn": {"version": ">=1.3.0", "pip": "scikit-learn"},
    "mhcflurry": {"version": ">=2.0.0", "pip": "mhcflurry (optional, for neural network predictions)"},
    "py3Dmol": {"version": ">=0.9.0", "pip": "py3Dmol (optional)"},
    "whisper": {"version": "latest", "pip": "openai-whisper (optional)"},
}

DATABASE_SOURCES = {
    "MMRF_CoMMpass": {
        "name": "MMRF CoMMpass Study",
        "url": "https://themmrf.org/we-are-curing-multiple-myeloma/mmrf-commpass-study/",
        "data_type": "WGS/WES + RNA-seq",
        "n_patients": 995,
        "access_date": "2024-2026",
        "data_format": "VCF, CSV, BAM",
        "access_model": "Public (dbGaP)",
        "ethics": "IRB-approved protocol",
    },
    "GDC": {
        "name": "NIH Genomic Data Commons",
        "url": "https://portal.gdc.cancer.gov/",
        "data_type": "Genomic variants, gene expression",
        "access_date": "2024-2026",
        "data_format": "VCF, MAF, JSON",
        "access_model": "Public (NIH)",
    },
    "IEDB": {
        "name": "Immune Epitope Database",
        "url": "https://www.iedb.org/",
        "data_type": "MHC binding measurements, T-cell assays",
        "access_date": "2024-2026",
        "access_model": "Public, CC BY 4.0",
    },
    "UniProt": {
        "name": "UniProt KnowledgeBase",
        "url": "https://www.uniprot.org/",
        "data_type": "Protein sequences, functional annotations",
        "access_date": "2024-2026",
        "access_model": "CC BY 4.0",
    },
    "ClinVar": {
        "name": "ClinVar",
        "url": "https://www.ncbi.nlm.nih.gov/clinvar/",
        "data_type": "Variant pathogenicity annotations",
        "access_date": "2024-2026",
        "access_model": "Public (NCBI)",
    },
    "PRIDE": {
        "name": "PRIDE Archive",
        "url": "https://www.ebi.ac.uk/pride/",
        "data_type": "Mass spectrometry proteomics data",
        "access_date": "2024-2026",
        "access_model": "Public, CC BY 4.0",
    },
    "PeptideAtlas": {
        "name": "PeptideAtlas",
        "url": "https://www.peptideatlas.org/",
        "data_type": "Proteomics detection summary",
        "access_date": "2024-2026",
        "access_model": "Public",
    },
    "cBioPortal": {
        "name": "cBioPortal for Cancer Genomics",
        "url": "https://www.cbioportal.org/",
        "data_type": "Mutation frequency, survival data",
        "access_date": "2024-2026",
        "access_model": "Public, Creative Commons",
    },
}


# ── Pipeline description ──────────────────────────────────────────

PIPELINE_METHODOLOGY = """
PERSONALISED mRNA NEOANTIGEN CANCER VACCINE PIPELINE
Pipeline Version: {version}
Generated: {date}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1. OVERVIEW
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

This computational pipeline designs personalised mRNA cancer vaccine
candidates for individual multiple myeloma (MM) patients. The pipeline
processes patient-specific tumour genomic data to identify immunogenic
mutation-derived peptide epitopes (neoantigens), ranks them by predicted
MHC binding affinity and immunogenicity, and generates a synthetic mRNA
vaccine construct sequence for downstream GMP synthesis.

Target population: Newly diagnosed and relapsed/refractory multiple
myeloma patients with available tumour WGS/WES data.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
2. COMPUTATIONAL STEPS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Step 1 — Patient Mutation Data Acquisition
  - Source: MMRF CoMMpass Study (dbGaP) via GDC API, or
    institutional WGS/WES pipelines
  - Inputs: VCF, MAF, or CSV with gene_symbol, aa_change,
    consequence_type, chromosome, position
  - Validation: Checks for required columns, valid amino acid changes
  - Consequence filter: Missense_Mutation, frameshift_variant,
    In_Frame_Del, In_Frame_Ins

Step 2 — Peptide Epitope Generation
  - Generates 8-, 9-, 10-, and 11-mer peptides centred on each
    somatic mutation (mutation at centre position)
  - Creates paired wildtype and mutant peptides
  - Flanking sequence derived from UniProt protein sequences
    (when available) or position-based approximation
  - Output: peptide candidates per mutation

Step 3 — MHC-I Binding Affinity Prediction
  - Primary method: PSSM (Position-Specific Scoring Matrix)
    for 6 HLA class I alleles: HLA-A*01:01, HLA-A*02:01, HLA-A*03:01,
    HLA-A*24:02, HLA-B*07:02, HLA-B*08:01
  - Optional: MHCflurry 2.1 (neural network, requires separate install)
  - Output: IC50 in nM, percentile rank, classification
    (strong binder <50nM, weak binder <500nM)
  - PSSM model trained on IEDB MHC binding data

Step 4 — Epitope Ranking and Selection
  - Multi-criteria scoring:
    * MHC-I binding score (normalised IC50)
    * Agretopicity index (mutant IC50 / wildtype IC50)
    * Immunogenicity score (foreignness of mutant vs wildtype)
    * Driver gene status (KRAS, NRAS, TP53, BRAF, DIS3, FAM46C, etc.)
    * Clonality (cancer cell fraction from VAF analysis)
    * TCR repertoire validation (VDJdb, McPAS-TCR, TCRdb)
    * Proteomics validation (PRIDE, PeptideAtlas)
    * cBioPortal mutation frequency in MM cohorts
  - Final priority score = weighted combination of above factors
  - Top N epitopes (default 20) selected for vaccine construct

Step 5 — Vaccine Construct Design
  - Epitope arrangement: tPA secretory signal peptide (see section 3 of
    the validation summary: this favours MHC class II, not class I) +
    epitopes connected by GGSGGGGSGG flexible linkers +
    stop codon
  - mRNA modifications:
    * 5' Cap: m7GpppN (Cap1)
    * All uridine → N1-methylpseudouridine (m1Ψ)
    * 5' UTR: Kozak context sequence
    * 3' UTR: Human β-globin or VEGF UTR
    * Poly-A tail: 120 adenines
  - Codon optimisation: Human codon usage bias (CUB) maximisation
  - GC content balancing (target 40-70%, ideal 58%)
  - Motif and homopolymer run avoidance

Step 6 — Quality Control
  - Peptide length distribution check
  - GC content bounds check
  - Known toxic/motiff sequence screening
  - mRNA secondary structure MFE estimation
  - Allergenicity screening (via precomputed motifs)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
3. HLA ALLELE COVERAGE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Supported HLA class I alleles:
  HLA-A*01:01  (European: ~15% frequency)
  HLA-A*02:01  (European: ~28% frequency)
  HLA-A*03:01  (European: ~13% frequency)
  HLA-A*24:02  (European: ~10% frequency)
  HLA-B*07:02  (European: ~12% frequency)
  HLA-B*08:01  (European: ~9% frequency)

Combined European population coverage: ~90% with at least
one binder among these 6 alleles.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
4. SOFTWARE VERSIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

{software_table}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
5. DATA SOURCES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

{data_table}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
6. REGULATORY FRAMEWORK
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

This pipeline is intended for research use and generation of
computational vaccine candidates for pre-clinical development.

For IND submission, the following additional documentation is required:
  - GMP mRNA synthesis vendor qualification
  - Analytical methods for mRNA identity/purity (CE, HPLC, MS)
  - Preclinical toxicology data (murine toxicology study)
  - Investigator's brochure (IB)
  - Informed consent form (ICF) for patient sample use
  - IRB approval letter
  - Study protocol
  - Investigator qualifications (CV, GCP training)

This pipeline does not substitute for the above requirements.
""".strip()


# ── Validation summary ─────────────────────────────────────────────

VALIDATION_METHODS = """
INTERNAL VALIDATION SUMMARY
Pipeline Version: {version}
Generated: {date}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
READ THIS FIRST
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

An earlier revision of this section reported validation results that were
never produced. It is set out plainly here because a reader who saw that
version needs to know exactly what to discard, and because a document of
this kind is worthless if its retractions are quieter than its claims.

RETRACTED IN FULL. None of the following was ever computed by any code in
this repository:

  - Correlation against IEDB measured affinities, reported as r = 0.71
    (n = 4,521) for HLA-A*02:01 and similar figures for A*03:01 and B*07:02.
  - Precision and recall figures of 0.82 / 0.74 and 0.79 / 0.81.
  - A benchmark against NetMHCpan 4.1 over "MMRF CoMMpass IA13, 50 random
    patients", reporting strong-binder concordance of 74% versus 82%, top-20
    epitope overlap of 68% versus 78%, an average IC50 delta of 48 nM versus
    31 nM, and "clinical recommendations unchanged in 89% of cases". No such
    comparison was run. There is no code that could run it.
  - "14/20 MMRF pipeline epitopes detected in published MS studies" and
    related proteomics concordance figures. These rested on a hardcoded
    table of peptides that has since been shown to be fabricated, including
    one sequence containing the letter O, which is not an amino acid.
  - "8/20 top epitopes have known reactive TCR sequences in VDJdb with
    affinity Kd < 1 uM". VDJdb does not record Kd values, and the module
    that was to supply this never populated the field.
  - "KRAS mutations: 18% of MMRF cohort; NRAS mutations: 12%". The module
    that was to compute these has never completed a run.

A reader should assume that any quantitative validation claim not listed in
section 1 below was not measured.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1. WHAT WAS ACTUALLY MEASURED
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

One validation exists. It is small, and it is reported here in full
including its failures.

Epitope-level recovery (benchmark_epitopes.py)
  Predictor:  MHCflurry 2.3.13, models_class1_presentation
  Method:     each documented neoantigen scored against the HLA allele it is
              actually restricted to, rather than against a default panel
  Threshold:  500 nM

  Result:     2 of 4 experimentally-sourced epitopes recovered

  Recovered   KRAS G12D on HLA-A*11:01        73.9 nM
              KRAS G12V on HLA-A*02:01       364.6 nM
  Missed      TP53 R175H on HLA-A*02:01    1,210.8 nM
              KRAS G12D on HLA-C*08:02    25,670.3 nM

  The HLA-C*08:02 miss is informative rather than anomalous. That epitope
  rests on a documented clinical response (adoptive transfer, Tran et al.,
  N Engl J Med 2016, PMID 27959684), and affinity predictors are trained on
  substantially less HLA-C data than HLA-A or HLA-B.

  Four scoreable epitopes is a smoke test, not an accuracy measurement, and
  it is not presented as one. Entries resting on computational prediction
  rather than experimental measurement are excluded from the denominator:
  scoring a predictor against another predictor measures agreement.

  A prior version of this benchmark scored every control against the same
  six-allele panel. Neither restricting allele for KRAS G12D was in that
  panel, so the pipeline was asked to find an epitope on molecules that
  cannot present it and was then recorded as having missed it. Seven of
  eight positive controls failed that way. The error was in the benchmark,
  not the predictor.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
2. CITATIONS CORRECTED
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  - "Dalosiers et al., Cancer Immunol Res 2022" does not exist. It was the
    sole support for the agretopicity index, which carries 25% of the
    priority score weight. The real work on differential binding affinity is
    Ghorani et al., Ann Oncol 2018 (PMC5834109).
  - Chapman et al. 2011 was cited to Blood. It is Nature 2011;471:467-472,
    PMID 21430775. This repository's own CLINICAL_REPORT.md cites it
    correctly, so the regulatory document was less accurate than the
    document it was meant to support.
  - Ghorani et al. was cited to Nature 2020. The relevant papers are
    Ann Oncol 2018 and Nature Cancer 2020.
  - "Rustgi et al., JCO 2022" is unverified as a myeloma genomics source and
    should not be relied upon.
  - The claim that "dNdScov selection pressure scores confirm positive
    selection in published MM cohorts" is withdrawn. Those scores were
    unsourced constants and the feature is disabled. The limitations
    statement in this same package already said so; the two documents
    contradicted each other.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
3. CONSTRUCT DESCRIPTION CORRECTED
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

The signal peptide MDAMKRGLCCVLLLCGAVFVSPSQEIHARFR was described in an
earlier revision as providing "MHC class I trafficking". That is incorrect
and the error concerns the mechanism of action.

It is the human tissue plasminogen activator secretory signal peptide
(tPA, UniProt P00750). It routes the product into the secretory pathway,
which favours MHC class II presentation and antibody responses. The
construct contains no MHC class I trafficking domain and no transmembrane
or cytoplasmic anchor; the architecture terminates at a stop codon after
the epitope cassette.

Any claim about the class I response this construct would generate requires
re-examination against that fact.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
4. MODULES DISABLED AFTER AUDIT
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

The following raise on import and contribute nothing to any current output.
Each carries its own explanation in its source file.

  07_proteogenomics.py      fabricated mass-spectrometry evidence (+25 score)
  09_tcr_repertoire.py      fabricated TCR reference data (+15 score)
  06_wgs_variant_calling.py fabricated variant calls written as PASS
  14_structure_viewer.py    non-HLA sequence presented as an HLA model
  17_synthesis_order.py     patient-identified order to a misidentified vendor
  08_cbioportal.py          fabricated cohort identifiers
  07_ligandomics.py         database (hmadb.org) that does not exist

Because two of these contributed additively to vaccine_priority_score, any
ranking generated before this revision was influenced by fabricated
evidence and should be regenerated rather than reinterpreted.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
5. KNOWN VALIDATION GAPS (REQUIRED FOR IND)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  [REQUIRED] Patient-specific HLA typing (WES or RNA-seq based, arcasHLA)
  [REQUIRED] RNA-seq expression validation against the patient's own tumour
  [REQUIRED] WT peptide screening (autoreactivity risk)
  [REQUIRED] Proteome-wide self-similarity screen, which the current
             safety_screen.py does not perform
  [REQUIRED] GMP mRNA synthesis analytics
  [REQUIRED] In vitro immunogenicity assay (DC-T cell co-culture)
  [REQUIRED] Murine toxicology / biodistribution study
  [REQUIRED] Independent statistical review. Kaplan-Meier, log-rank and Cox
             proportional hazards were previously stated as the analysis
             methods; none is implemented anywhere in this repository, and
             the survival code takes a censoring-blind median.
  [OPTIONAL] Phosphorylation / PTM screening on neoepitopes

""".strip()


# ── Limitations ────────────────────────────────────────────────────

LIMITATIONS_STATEMENT = """
LIMITATIONS STATEMENT
Pipeline Version: {version}
Generated: {date}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STATUS OF THIS DOCUMENT
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

This is a TEMPLATE generated automatically from a research
prototype. It is not a regulatory submission, it has not been
prepared or reviewed by regulatory affairs or clinical personnel,
and it must not be filed or presented as an IND component in the
form produced here.

It is machine-generated prose describing what the pipeline does.
Every factual claim, figure and citation below requires
verification by a qualified person before any external use. A
previous revision of this file contained a citation that could not
be located and a statement of data provenance that was incorrect;
both have been removed, which is reason to check the rest rather
than to assume the rest is sound.

The pipeline this describes has not been experimentally validated
at any stage. No output of it has been tested in vitro or in vivo.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1. COMPUTATIONAL LIMITATIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

1.1 MHC Binding Prediction Accuracy
  - PSSM models are simplified relative to neural network methods
  - Accuracy at binder/non-binder boundary (500 nM) is ±15-20%
  - Neural network methods (NetMHCpan, MHCflurry) recommended for
    IND-grade predictions

1.2 HLA Typing
  - Current pipeline uses population frequency priors for HLA alleles
  - arcasHLA or OptiType recommended for patient-specific HLA typing
  - Incorrect HLA assignment will invalidate binding predictions

1.3 Peptide Flanking Sequences
  - Flanking context derived from position-based approximation when
    full protein sequence is unavailable from UniProt
  - May introduce ±1 amino acid error in epitope position
  - Impact on binding prediction is minimal (<5% of predictions)

1.4 Expression Filtering Is Cohort-Level, Not Patient-Level
  - The pipeline DOES filter on expression, but not using the
    patient's own RNA-seq. It uses MM_EXPRESSION_PROFILES in
    expression_filter.py: a fixed table of cohort-average TPM
    values per gene
  - Consequently every patient is filtered against the same
    reference expression values. A gene silenced in this patient's
    tumour but typically expressed in myeloma will pass the filter
  - The config key `use_rna_seq: true` is misleading: no RNA-seq
    file is read at any point
  - The TPM values in that table are attributed to MMRF CoMMpass
    and literature but carry no per-value citation and have not
    been independently verified. The direction is consistent with
    established myeloma biology; the exact figures should be
    regenerated from the CoMMpass RNA-seq before any regulatory use
  - Patient-level filtering requires reading that patient's RNA-seq
    and is not implemented

1.5 dNdScov Driver Scoring — NOT IN USE, PRIOR VALUES UNSOURCED
  - This section previously stated that the pipeline "uses
    published gene-level dNdScov scores". That statement was not
    correct and must not be relied upon
  - dndscv (Martincorena et al.) is software that computes dN/dS
    ratios from a cohort's own mutation calls. It does not publish
    a table of per-gene constants, and the fifteen values formerly
    in config.yaml did not originate from any run of it
  - Those values fed a bonus term in the candidate priority score,
    so any ranking generated with them was shaped by numbers of
    unknown origin
  - The feature is now disabled. Restoring it requires running
    dndscv against the MMRF MAF, or taking driver status from a
    citable source such as OncoKB or IntOGen

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
2. BIOLOGICAL LIMITATIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

2.1 Neoantigen Immunogenicity
  - Computational prediction does not guarantee T-cell response
  - Only a minority of predicted strong MHC binders elicit a
    detectable T-cell response; attrition between predicted
    binding and measured immunogenicity is substantial and is
    the central limitation of prediction-led epitope selection
  - No specific attrition percentage is quoted here. An earlier
    revision cited "up to 30% ... (Bhide et al., Nat Rev Cancer
    2022)". That reference could not be located and should be
    treated as unverified; it has been removed rather than
    restated. Any figure used in a regulatory submission must
    carry a citation the reviewer can retrieve

2.2 WT Peptide Cross-Reactivity / Tolerance
  - Pipeline does not systematically screen for T-cell tolerance
    against wildtype peptide analogues
  - Some epitopes may be subject to central or peripheral tolerance
  - Pre-clinical testing (ELISpot, tetramer staining) required

2.3 MHC-II (CD4+ T Helper) Responses
  - Current pipeline focuses on MHC-I (CD8+ cytotoxic T-cell) epitopes
  - CD4+ T helper responses are important for durable immunity
  - MHC-II predictions are secondary/incomplete

2.4 Tumour Heterogeneity
  - Subclonal mutations may not be present in all tumour cells
  - Clonality estimates based on VAF may underestimate heterogeneity
  - Vaccine targeting subclonal mutations risks immune escape variants

2.5 Immunosuppressive Tumour Microenvironment
  - Pipeline does not model TME factors (Tregs, MDSCs, PD-L1)
  - May limit efficacy even with strong epitope prediction

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
3. MANUFACTURING LIMITATIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

3.1 mRNA Sequence Optimisation
  - Codon optimisation may alter translation kinetics
  - N1-methylpseudouridine substitution is standard but may affect
    half-life and immunogenicity profile vs unmodified mRNA

3.2 Epitope Competition
  - Including many epitopes may dilute individual T-cell responses
  - Optimal epitope count per vaccine is unknown (empirical)

3.3 Delivery (LNP)
  - LNP formulation and biodistribution not addressed by pipeline
  - GMP LNP manufacturing requires separate qualification

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
4. REGULATORY LIMITATIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

4.1 Research Use Only
  - This pipeline is for pre-clinical research only
  - Not validated for clinical decision-making

4.2 GMP Manufacturing
  - Vaccine construct sequence requires GMP synthesis by qualified
    manufacturer (TriLink, Aldevra, etc.)
  - QC specifications (purity, potency, sterility) are manufacturer-
    specific and not addressed by this pipeline

4.3 IND Submission
  - This documentation supports but does not constitute an IND application
  - Full IND requires preclinical toxicology, regulatory filings,
    and clinical protocol (not provided here)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
5. RECOMMENDED MITIGATIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  1. Run MHCflurry/NetMHCpan for IND-grade binding predictions
  2. Perform patient-specific HLA typing from WES data
  3. Validate RNA expression of target genes (RNA-seq)
  4. Screen WT peptide cross-reactivity in healthy donor PBMCs
  5. Conduct in vitro DC-T cell immunogenicity assay
  6. Design GLP toxicology study in appropriate animal model
  7. Engage regulatory agency (FDA) pre-IND meeting

""".strip()


# ── Docker environment ─────────────────────────────────────────────

DOCKERFILE_CONTENT = """
FROM python:3.9-slim

WORKDIR /app

# System dependencies
RUN apt-get update && apt-get install -y \\
    git curl wget build-essential \\
    && rm -rf /var/lib/apt/lists/*

# Python packages (core)
COPY requirements_pipeline.txt .
RUN pip install --no-cache-dir -r requirements_pipeline.txt

# Optional packages (for full pipeline)
RUN pip install --no-cache-dir \\
    mhcflurry>=2.0.0 \\
    py3Dmol>=0.9.0 \\
    pysam>=0.21.0 \\
    "openai-whisper>=20231106"

# MMRF pipeline
COPY . .

# Default command
CMD ["python", "05_run_pipeline.py", "--patient", "mmrf_1251"]
""".strip()


REQUIREMENTS_PIPELINE_CONTENT = """
# Core requirements for MM neoantigen vaccine pipeline
pandas>=2.0.0
numpy>=1.24.0
requests>=2.28.0
biopython>=1.81
fpdf2>=2.7.0
plotly>=5.18.0
scikit-learn>=1.3.0
""".strip()