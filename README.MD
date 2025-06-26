# Space-time data-driven modeling of wildfire initiation in the mountainous region of Trentino--South Tyrol, Italy
This is the R code to reproduce the analysis in "Space-time data-driven modeling of wildfire initiation in the mountainous region of Trentino--South Tyrol, Italy." by
Mateo Moreno <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0002-9530-3076)</sup>, 
Stefan Steger <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0003-0886-5191)</sup>, 
Laura Bozzoli <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)(https://orcid.org/0000-0002-2648-2007)</sup>, 
Stefano Terzi<sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)(https://orcid.org/0000-0003-0491-6772)</sup>,
Andrea Trucchia <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)(https://orcid.org/0000-0001-7294-9061)</sup>
Cees van Westen <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0002-2992-902X)</sup> and 
Luigi Lombardo <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0003-4348-7288)</sup>,
[![DOI](https://zenodo.org/badge/887347458.svg)](https://doi.org/10.5281/zenodo.15033256)

> Moreno, M., Steger, S., Bozzoli, L., Terzi, S., Trucchia, A., van Westen, C. J. & Lombardo, L. (2025). Space-time data-driven modeling of wildfire initiation in the mountainous region of Trentino--South Tyrol, Italy. https://eartharxiv.org/repository/view/9257/


## Abstract
Wildfires are complex hazards occurring worldwide, leading to substantial economic losses, fatalities, and carbon emissions. The interplay of climate change, land use alterations, and socioeconomic pressures is expected to further increase the frequency and intensity of wildfires. In this context, developing reliable, dynamic prediction tools is essential for risk mitigation. This work presents a spatiotemporal wildfire prediction model for the Trentino-South Tyrol region (13,600 km²) in the northeastern Italian Alps. Leveraging generalized additive models, we integrate multitemporal wildfire records (2000--2023) with static and dynamic environmental controls (e.g., topography, land cover, daily precipitation, and temperature). The resulting model predictions change dynamically over space and time in response to static features, seasonal trends, and evolving meteorological conditions. Model outputs were evaluated using established performance metrics, enabling the derivation of dynamic spatial wildfire probability thresholds. These thresholds are illustrated for varying amounts of precipitation, temperature, and different combinations of static factors. Validation through multiple perspectives yielded performance scores generally exceeding 0.8, confirming the model strong generalization and transferability. To demonstrate the practical application, the model was used to hindcast past wildfire initiation between 1--15 July 2022, a period marked by elevated wildfire activity. By integrating static and dynamic environmental controls, this research advances the spatiotemporal prediction of wildfires in complex alpine regions, supporting the development of early warning systems.


## Acknowledgements
The research that led to these results is related to the EO4MULTIHA project (https://eo4society.esa.int/projects/eo4multihazards/), which received funding from the European Space Agency. We thank Dr. Serkan Girgin for his support in using the CRIB platform. We thank Dr. Paolo Fiorucci and Dr. Marj Tonini for their valuable feedback. Finally, we thank the Office for Forest Planning of the Autonomous Province of Bolzano, especially Alessandro Andriolo, for providing and revising the wildfire data. 


## Further dissemination
Previous work was presented at the EGU General Assembly 2025. The abstract is available at:

> Moreno, M., Steger, S., Bozzoli, L., Terzi, S., Trucchia, A., van Westen, C., & Lombardo, L. (2025). Space-time data-driven modeling of wildfire initiation in the mountainous region of Trentino-Alto Adige, Italy. Abstract from EGU General Assembly 2025, Vienna, Austria. https://doi.org/10.5194/egusphere-egu25-17023

## Repo structure
The general structure is as follows:
- `dat`: data sets
- `dev`: development (scripts)
- `plt`: plots
