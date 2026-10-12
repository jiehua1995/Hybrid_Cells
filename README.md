# Constrained resolution of genomic incompatibility in hybrid *Drosophila* cell lines

This repository contains the R/Quarto notebooks used for data analysis and figure generation for a study of *Drosophila melanogaster* x *Drosophila simulans* hybrid cell lines.

## Repository structure

Every notebook is named after the figure(s) it produces, and its folder uses the same figure-number name. Each notebook's first line states which figure(s) it produces.

| Figure(s) | Notebook |
|---|---|
| Fig. 1 | `Fig1/Fig1.qmd` |
| Fig. 2a | `Fig2a/Fig2a.qmd` |
| Fig. 2b | `Fig2bc/Fig2b.qmd` |
| Fig. 2c | `Fig2bc/Fig2c.qmd` |
| Fig. 3b | `Fig3b/Fig3b.qmd` |
| Fig. 3c | `Fig3c/Fig3c.qmd` |
| Fig. 3d, Extended Data Fig. 7 | `Fig3d_ExtFig7/Fig3d_ExtFig7.qmd` |
| Fig. 4 | `Fig4/Fig4.qmd` |
| Fig. 5, Extended Data Fig. 9 | `Fig5_ExtFig9/Fig5_ExtFig9.qmd` |
| Extended Data Fig. 1 | `ExtFig1/ExtFig1.qmd` |
| Extended Data Fig. 2 | `ExtFig2/ExtFig2.qmd` |
| Extended Data Fig. 3 | `ExtFig3/ExtFig3.qmd` |
| Extended Data Fig. 4 | `ExtFig4/ExtFig4.qmd` |
| Extended Data Fig. 5 | `ExtFig5/ExtFig5.qmd` |
| Extended Data Fig. 6 | `ExtFig6/ExtFig6.qmd` |
| Extended Data Fig. 8 | `ExtFig8/ExtFig8.qmd` |

Support scripts:
- `Fig3c/functions.R` and `Fig3d_ExtFig7/functions.R` — helper functions sourced by the Fig. 3c / Fig. 3d notebooks.
- `Fig3d_ExtFig7/Fig3d_Org_DB.qmd` — builds the orthogroup-GO annotation package used by Fig. 3d / Extended Data Fig. 7.
- `ExtFig1/ExtFig1_Build_BSgenome.R` — builds the BSgenome data package for the Ras3 assembly.

## Repository status

The associated manuscript is currently in preparation or under review. The notebooks and repository structure may therefore be updated before publication.

A finalized and fully documented version of the code will be released upon acceptance of the manuscript.

## Requirements

The analyses were performed using R and Quarto. Required R packages are listed within the individual notebooks.

## Citation

Citation information will be added after publication of the associated manuscript.
