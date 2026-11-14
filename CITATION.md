# Data citation and licence

This repository contains code only. The data is downloaded by
`python -m drbm.download` and is never committed.

## Sources

| Data | Source | File |
| --- | --- | --- |
| Drug response | Genomics of Drug Sensitivity in Cancer (GDSC), GDSC2, release 8.5 | `GDSC2_fitted_dose_response_27Oct23.xlsx` |
| Drug response, older screen (Experiment 7 only) | GDSC1, release 8.5 | `GDSC1_fitted_dose_response_27Oct23.xlsx` |
| Cell line details | Sanger Cell Model Passports | `model_list_latest.csv.gz` |
| Mutations | Sanger Cell Model Passports | `mutations_summary_20230202.zip` |
| Copy number | Sanger Cell Model Passports | `cnv_summary_20230303.zip` |
| Gene expression | Sanger Cell Model Passports (Sanger and Broad sequencing runs, merged) | `rnaseq_merged_20250117.zip` |
| Protein abundance | Sanger Cell Model Passports (ProCan-DepMapSanger) | `Proteomics_20221214.zip` |

Checksums of every downloaded file are written to `data/raw/MANIFEST.json`.
`model_list_latest.csv.gz` is a "latest" file and may change; the checksum
records which version was used.

## Papers to cite

- Yang, W., et al. (2013). Genomics of Drug Sensitivity in Cancer (GDSC): a
  resource for therapeutic biomarker discovery in cancer cells. *Nucleic Acids
  Research*, 41(D1), D955-D961. https://doi.org/10.1093/nar/gks1111
- Iorio, F., et al. (2016). A landscape of pharmacogenomic interactions in
  cancer. *Cell*, 166(3), 740-754. https://doi.org/10.1016/j.cell.2016.06.017
- van der Meer, D., et al. (2019). Cell Model Passports: a hub for clinical,
  genetic and functional datasets of preclinical cancer models. *Nucleic Acids
  Research*, 47(D1), D923-D929. https://doi.org/10.1093/nar/gky872
- Gonçalves, E., et al. (2022). Pan-cancer proteomic map of 949 human cell
  lines. *Cancer Cell*, 40(8), 835-849. https://doi.org/10.1016/j.ccell.2022.06.010

The clinical background for the choice of drug:

- Sweeney, C., et al. (2021). Ipatasertib plus abiraterone and prednisolone in
  metastatic castration-resistant prostate cancer (IPATential150): a
  multicentre, randomised, double-blind, phase 3 trial. *The Lancet*,
  398(10295), 131-142.

The journal details above were checked against library and publisher records
on 13 November 2026. The papers themselves were not read in full.

TODO (Paul): the gene expression file merges sequencing runs from the Sanger
and Broad institutes. Find and add the citations the Cell Model Passports
download page asks for (the Broad runs come from the Cancer Cell Line
Encyclopedia).

## Licence

GDSC, checked on 14 November 2026. The bulk download page grants a
non-exclusive, non-transferable right to use the data files for internal
research and educational purposes. It excludes resale of the data, alone or
combined with other data or products, and use of the data to provide
commercial services. This study is non-commercial research and redistributes
no data, so it falls inside those terms. Anyone reusing the code for
commercial work needs their own agreement with Genome Research Limited.

Cell Model Passports: not yet checked.

TODO (Paul): open the terms of use on the Cell Model Passports download page,
confirm that the data may be used for this purpose, and paste the required
acknowledgement here.
