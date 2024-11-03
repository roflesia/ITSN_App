You can access the ITSNdb AppWeb by following the next link: https://itsndb.streamlit.app/ 

Neoantigens are unique peptides arising from tumor-specific genomic alterations, presented on the cell surface by Major Histocompatibility Complex Class I (MHC-I) molecules, where they can be recognized
by T cells to trigger an immune response against cancer cells. Cancer immunotherapies, such as immune checkpoint inhibitors and neoantigen vaccines, exploit these immunogenic neoantigens, making their
identification essential for precision medicine.
Given the costly and time-intensive nature of experimental neoantigen identification, computational preselection of peptides has emerged, with the primary challenge being the prediction of neoantigen
immunogenicity. To provide a comprehensive framework for assessing immunogenicity predictors and prioritizers, the Immunogenic Tumor-Specific Neoantigen Database (ITSNdb) was created. ITSNdb is a
curated database of experimentally validated neoantigens, collected from peer-reviewed literature and meeting specific criteria to enhance the accuracy and relevance of computational predictions. In this
work, a user-friendly web application was developed to facilitate ITSNdb’s usage by scientists without bioinformatics expertise.

# The ITSNdb app
A Streamlit web based application was developed to enhance accessibility and exploration of the neoantignes peptides as well as their wildtype counterpart and predicted protein structures.
The developed web application consists of three main sections:
The Home section provides an overview of ITSNdb along with statistical visualisations displaying all the information related to the stored peptides such as the amount, classification types, HLA
information, frequency and tumor types.
DataBase exploration: By different filters individual peptides can be explored in detail, providing wild-type and mutated peptide sequence, mutation position and localization type, gene accession,
tumor type, publication source and HLA characteristics. Then by mean of an interactive view, 3D visualization tool allows for side-by-side comparisons of wildtype and mutated predicted 3D protein
structures along with the ir pLDDT confidence metrics, superimposition capabilities and comparative metrics.
The Software Validation tab allows users to upload their immunogenicity predictions based on three datasets (ITSNdb, validation and neopeptides from immunotherapy cohorts). They can evaluate
performance using metrics like ROC curves, F1 scores, and confusion matrices, supporting in-depth analysis to improve the accuracy of immunogenicity prediction models.

The developed web application effectively addresses the need for a comprehensive tool to explore neoantigens and assess immunogenicity predictions. By providing an intuitive interface and advanced
visualization capabilities, it enhances the understanding of neoantigen characteristics and supports personalized immunotherapy approaches.
The importance of this project is clear in its contribution to cancer immunotherapy, offering a platform to improve the accuracy and reliability of immunogenicity prediction models
