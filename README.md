# MPXV-Pooling-VSP2-FieldStudy

This repository supports the manuscript:

**“Evaluating High-Plex Sample Pooling for the Illumina Viral Surveillance Panel v2 in Mpox Virus Genomic Detection under Field and Resource-Constrained Conditions”**

It includes bioinformatics workflows, scripts, metadata, and processed results from an evaluation of sample multiplexing strategies during the 2024 Mpox outbreak in Burundi. The study compares sequencing performance across different pooling levels (3, 6, 9, and 15 samples per capture) using the Illumina VSP2 panel and iSeq 100 platform.

---

## 📂 Repository Contents

| Folder/File         | Description                                                                 |
|---------------------|-----------------------------------------------------------------------------|
| `data/`             | Sample metadata, Ct values, genome coverage, and variant tables             |
| `scripts/`          | Shell scripts and Python notebooks for preprocessing, mapping, and analysis |
| `results/`          | Processed outputs including variant metrics, assembly stats, and coverage   |
| `figures/`          | Publication-ready plots (PNG and PDF formats)                               |
| `MANUSCRIPT.md`     | Draft manuscript or link to final publication                               |
| `LICENSE`           | Licensing information (MIT by default)                                      |

---

## 🔬 Methods Overview

- **Preprocessing**: FastP v0.20.1 for quality filtering
- **Host Depletion**: Kraken2, SeqKit filtering
- **Alignment**: BWA-MEM v0.7.17
- **Variant Calling**: iVar, bcftools
- **Assembly**: SPAdes v3.15, QUAST
- **Taxonomic Classification**: KrakenUniq, Kaiju, Centrifuge

All steps are portable and designed for use in field or low-resource settings.

---

## Highlights

- MPXV enrichment remains robust in pools ≤6
- Larger pools reduce coverage, variant yield, and assembly completeness
- High Ct (low viral load) samples underperform in high-pooling conditions
- SPAdes assembly quality degrades significantly beyond 6-sample pools

---


