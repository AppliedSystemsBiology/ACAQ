ACAQ: a Fiji and R toolkit targeted for automated confrontation assay quantification
====================================================================================

Zoltan Cseresnyes (1,+), Naim Al-Zaben (1,2,+), Mohamed I. Abdelwahab Hassan (2,3,4), Ruman Gerst (1,2), Kerstin Voigt (3), Marc Thilo Figge (1,2,\*)

(1)  Applied Systems Biology, Leibniz Institute for Natural Product Research and Infection Biology – Hans Knöll Institute, Jena, Germany

(2)  Faculty of Biological Sciences, Friedrich Schiller University Jena, Jena, Germany

(3)  Jena Microbial Resource Collection, Leibniz Institute for Natural Product Research and Infection Biology – Hans Knöll Institute, Jena, Germany

(4)  Pests and Plant Protection Department, National Research Centre, Giza, Egypt 

(+) These authors contributed equally to this work.

(\*) Correspondence should be addressed to MTF: <thilo.figge@leibniz-hki.de>.

Host-pathogen interactions are ideally studied by imaging confrontation assays, where interacting cell types are monitored during and/or after their interaction. The resulting images are quantified via automated segmentation and classification algorithms, which help to unravel biological functions from host-pathogen interactions. Here we describe a modular and easily extendable segmentation and analysis framework designed both for fluorescently labeled and unlabeled cells participating in such interactions. Acquiring images of unlabeled biological samples has multiple advantages: it minimizes the side effects that fluorescence labeling may have on biological functions, and it saves the time and expense spent on fluorescence labeling. Our new high- throughput confrontation assay-analysis and -visualisation open-source toolkit, algorithm for confrontation assay quantification (ACAQ), is a Fiji macro extended by R. It consists of two parts: ACAQ Analyzer and ACAQ Visualizer. The former is designed to analyze cells in two-dimensional microscopy images, using either fluorescence labeling or transmitted light bright-field microscopy. In the absence of labeling, Hessian filter-based segmentation is used, which is able to distinguish between various types of cells that differ in size. The ACAQ toolkit computes a collection of phagocytosis measures and morphological descriptors. The results are summarized and exported by ACAQ Visualizer, providing immediate feedback to the user.

# ACAQ Analyzer

You can find the code of ACAQ analyzer in the `ACAQ_Analyzer` folder. See the [supplementary information](https://asbdata.hki-jena.de/publidata/CseresnyesEtAl_ACAQ/ACAQ_ScientificReports_SupplementaryInformation.pdf) to find instructions on how to use ACAQ Analyzer. We recommend to download the [binary package](https://asbdata.hki-jena.de/publidata/CseresnyesEtAl_ACAQ/Source%20code%20and%20packages.zip) if you are using Windows, as it comes pre-packaged with all necessary binaries to start the application.

# ACAQ Visualizer

The source code of ACAQ visualizer can be found in the `ACAQ_Visualizer` folder.
You can find instructions on how to use the application in the [supplementary information](https://asbdata.hki-jena.de/publidata/CseresnyesEtAl_ACAQ/ACAQ_ScientificReports_SupplementaryInformation.pdf).

We recommend to download the [binary package](https://asbdata.hki-jena.de/publidata/CseresnyesEtAl_ACAQ/Source%20code%20and%20packages.zip) if you are using Windows, as it comes pre-packaged with all necessary binaries to start the application.
