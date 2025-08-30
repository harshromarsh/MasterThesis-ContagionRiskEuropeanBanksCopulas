# Master-Thesis---Contagion-Risk-in-European-Banks-using-Copulas
Code for the Master Thesis "Contagion Risk in European Banks using Copulas"

The aim of the thesis was to investigate the dependence structure and dynamics
of financial contagion in the European banking sector using copula-based models,
focusing on tail dependencies as well as systemic risk measures. 

This repository contains the R code and results and is structured in the following way:

```
├── Code/
│ ├── DataRetrieval.R # Retrieve and preprocess stock data (source: Yahoo Finance via Quantmod)
│ ├── DataRetrieval_CDS.R # Preprocess CDS data from a sep. Excel (source: Bloomberg)
│ ├── RVineCopula.R # Estimation of R-vine copulas (stocks)
│ ├── RVineCopulaCDS.R # Estimation of R-vine copulas (CDS)
│ ├── config.R # Configuration parameters (folders, colors, etc.)
│ └── functions.R # Helper functions (preprocessing, estimation, visualisation, etc.)
│
├── Plots/ # Figures generated
│ ├── CDS/ # From running RVineCopulaCDS.R
│ └── Stock/ # From running RVineCopula.R
│
├── Tables/ # LaTeX-format tables generated
│ ├── CDS/ # From running RVineCopulaCDS.R
│ └── Stock/ # From running RVineCopula.R
│
├── .gitignore
├── Master-Thesis---Contagion-Risk-in-European-Banks-using-Copulas.Rproj
└── README.md # Project description
```
