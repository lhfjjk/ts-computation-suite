# ChemRxiv submission package

**Title:** Beyond Candidate Generation: A Validation-First Protocol for AI-Assisted Transition-State Discovery

## Files

| File | Purpose |
|---|---|
| `manuscript.tex` / `manuscript.pdf` | Main manuscript (11 pages) |
| `references.bib` | Verified bibliography (Crossref-checked DOIs) |
| `SI.tex` / `SI.pdf` | Supporting Information (5 pages) |
| `cover_letter.tex` / `cover_letter.pdf` | Cover letter |
| `figures/workflow.pdf` / `workflow.png` | Figure 1 |
| `make_figure.py` | Regenerates Figure 1 |

## Build

```bash
pdflatex manuscript && bibtex manuscript && pdflatex manuscript && pdflatex manuscript
pdflatex SI && pdflatex SI
pdflatex cover_letter && pdflatex cover_letter
python make_figure.py   # only if Figure 1 changes
```

## Before submission: fill these

- [ ] **Corresponding email** in `cover_letter.tex` (`[corresponding email]`).
- [ ] **Affiliation** in `manuscript.tex` and `cover_letter.tex` (confirm the
      FJIRSM/UCAS wording, or replace with your current affiliation).
- [ ] **ORCID** in the ChemRxiv submission form (not required in the PDF).
- [x] **GitHub URL** is set to https://github.com/lhfjjk/ts-computation-suite
      in the manuscript and cover letter.
- [ ] **Case-study approval**: the numbers in Table 3 come from internal
      projects; confirm they may be released publicly.
- [ ] **Coauthor approval**: the manuscript currently lists a single author.
      Add coauthors and their contributions if needed.
- [ ] **License**: ChemRxiv asks you to choose one; CC BY 4.0 is the usual
      choice for a methods/preprint contribution.
- [ ] **Category**: Computational Chemistry (primary), Catalysis (secondary).

## Notes

- The manuscript is a methods/perspective contribution. It does not claim a new
  predictive model or new experimental data.
- The 11-page main text, 5-page SI, and one-page cover letter are within
  ChemRxiv's normal scope.
- All journal references were checked against Crossref; the DOI list is in
  `references.bib`. The Gaussian 16 entry is a software citation without a DOI.
- The visual check was done by inspecting the rendered pages; re-check the PDF
  in your own viewer before submitting, especially the figure and tables.
- The GitHub repository is public at https://github.com/lhfjjk/ts-computation-suite.
