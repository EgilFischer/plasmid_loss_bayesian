# Plasmid loss and conjugation estimation with Bayesian methods

This project contains files to estimate the loss of plasmids from bacterial population is conjection with uptake by conjugation. Two main parameters can be estimated "\rho" = loss parameter and "\gamma" = conjugation parameter. 

## Installation
Install Rstan: https://cran.r-project.org/web/packages/rstan/index.html

More information about Stan: https://mc-stan.org

R-libraries:
- rstan
- rstanarm
- ggplot2
- deSolve
- tidyverse
- easyABC

### Dependencies

Link to [dependencies](/renv.lock)

## Project Structure

The project structure distinguishes three kinds of folders:
- read-only (RO): not edited by either code or researcher
- human-writeable (HW): edited by the researcher only.
- project-generated (PG): folders generated when running the code; these folders can be deleted or emptied and will be completely reconstituted as the project is run.


```
.
├── .gitignore
├── CITATION.cff
├── LICENSE
├── README.md
├── requirements.txt
├── data               <- All project data, ignored by git (RO)
|
├── docs               <- Documentation notebook for users (HW)
│   ├── manuscript     <- Manuscript source, e.g., LaTeX, Markdown, etc. (HW)
│   └── reports        <- Other project reports and notebooks (e.g. Jupyter, .Rmd) (HW)
├── results
│   ├── figures        <- Figures for the manuscript or reports (PG)
│   └── output         <- Other output for the manuscript or reports (PG)
└── src                  <- Source code for this project (HW)
    ├── R              <- R-codes (HW)
    └── Stan         <- Stan code (HW)

```

## Add a citation file
Citation for this repository  [CITATION](/CITATION.cff)

## License

This project is licensed under the terms of the [MIT License](/LICENSE).
