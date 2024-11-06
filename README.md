
# The ITSNdb app
A Streamlit web based application was developed to enhance accessibility and exploration of the [The Immunenogenic Tumor Specific Neoantigen database](https://github.com/elmerfer/ITSNdb) providing access to neoantignes peptides, wildtype counterpart and their predicted protein structures.
The developed web application consists of three main sections:
The Home section provides an overview of ITSNdb along with statistical visualisations displaying all the information related to the stored peptides such as the amount, classification types, HLA
information, frequency and tumor types.
DataBase exploration: By different filters individual peptides can be explored in detail, providing wild-type and mutated peptide sequence, mutation position and localization type, gene accession,
tumor type, publication source and HLA characteristics. Then by mean of an interactive view, 3D visualization tool allows for side-by-side comparisons of wildtype and mutated predicted 3D protein
structures along with the ir pLDDT confidence metrics, superimposition capabilities and comparative metrics.
The Software Validation tab allows users to upload their immunogenicity predictions based on three datasets ([ITSNdb](https://github.com/elmerfer/ITSNdb), validation and neopeptides from immunotherapy cohorts). They can evaluate
performance using metrics like ROC curves, F1 scores, and confusion matrices, supporting in-depth analysis to improve the accuracy of immunogenicity prediction models.

The developed web application effectively addresses the need for a comprehensive tool to explore neoantigens and assess immunogenicity predictions. By providing an intuitive interface and advanced
visualization capabilities, it enhances the understanding of neoantigen characteristics and supports personalized immunotherapy approaches.
The importance of this project is clear in its contribution to cancer immunotherapy, offering a platform to improve the accuracy and reliability of immunogenicity prediction models

*Access to the [ITSNdb AppWeb](https://itsndb.streamlit.app/)*

## Cite

Nibeyro et al. [Unraveling Tumor Specific Neoantigen immunogenicity prediction: a comprehensive analysis](https://doi.org/10.3389/fimmu.2023.1094236) Front. Immunol.Sec. Cancer Immunity and Immunotherapy Volume 14 - 2023 | doi: 10.3389/fimmu.2023.1094236 

Nibeyro et al. [MHC-I binding affinity derived metrics fail to predict tumor specific neoantigen immunogenicity](https://doi.org/10.1101/2022.03.14.484285) BioRxiv

## Authors

* **Elmer Andrés Fernández** - *Idea and Initial work* - [Profile](https://www.researchgate.net/profile/Elmer_Fernandez) - [CONICET](http://www.conicet.gov.ar) - [Fundación para el Progreso de la Medicina - FPM](https://fpmlab.org.ar/)
* **Guadalupe Nibeyro** -*Idea and initial work* - [CONICET](http://www.conicet.gov.ar) - [Fundación para el Progreso de la Medicina - FPM](https://fpmlab.org.ar/)
* **Rocio Flesia** - *Developer and maintener* Undergraduate student at [FECFyN](https://fcefyn.unc.edu.ar/)


