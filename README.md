#  KnowYourCity data for research

## Background
90% of population growth in the next 25 years will be in African and Asian cities, and most people will be added in deprived “slum” areas. Due to informality, dynamism, and insecurity, data on “slums” and their residents are incomplete and dispersed, at best. Many “slums” are simply missing from official (e.g., census, survey) and other widely used datasets (e.g., OSM). Satellite imagery is often used to fill this data gap with “slums” being identified as clusters of small disorganized buildings and other morphologies of informality. Who generates these data, however, matters because deprivation manifests in myriad of ways within and between cities. Whether training AI models, evaluating the accuracy of modelled data sets, identifying areas at risk of hazards, etcetera, data about “slums” by residents themselves reflect important contextual variations. We process the KnowYourCity (KYC) Campaign dataset by Slum Dwellers International (SDI) to encourage its use for research.

## About SDI
SDI is a global network of community-based organizations (CBOs) founded by the urban poor more than 25 years ago to ensure that “slums” are integrated into cities, and that their residents enjoy basic human rights. A key activity among SDI-affiliated CBOs is to map their own settlement boundaries and collect key infrastructure and population data using mapping, survey, and community engagement methods. In 2016, SDI launched the KYC Campaign website (https://sdinet.org) with support from the United Cities and Local Governments of Africa, Cities Alliance, and other partners. Since, more than 7700 “slum” communities in 224 cities across Africa, Asia, and Latin America have submitted community profiles with the following information: settlement boundaries; population count; status (e.g., illegal); number of water taps, toilets, and other key infrastructure; common diseases; health care access; and commercial facilities. The KYC Campaign website intends to strengthen the network among community based organizations, provide visibility to systematic issues that affect the urban poor, engage city governments, and showcase the work of community-based data collectors, most of whom volunteer their time and expertise. 

## Motivation
While KYC settlement profiles are open, they are not cleaned, documented, and compiled into an easily downloadable format, making them challenging to use in research. The motivation for this Github page is to document and share this rich, field-referenced source of data to enable its use by researchers, while also highlighting the work of community data collectors and experts. 

## Data collection
Data collection methods can vary from one city to the next. However, generally, settlement boundaries are mapped by collecting GPS coordinates around the settlement perimeter. Population estimates are generally derived by physically marking and counting all front doors in the settlement, sampling every nth household to estimate average household size, and then multiplying number of front doors by the average household size in the settlement; this estimated number is then discussed and agreed by consensus in an open community forum.

Error in settlement boundaries can occur when too few points were collected around the settlement perimeter, leading to sharp angles that omit sections of the settlement. Another challenge for field data collectors is that some settlements are located in difficult to navigate areas such as next to, or on top of, a water body or swamp. Error in population estimates can happen due to simple data entry error, a miscalculation when multiplying household size by number of front doors, or misalignment between the area mapped and the area surveyed. In a few settlements where informal households were located among formal residential buildings (eg apartments), community data collectors appear to have only enumerated informal households. 

## Data processing
We applied a web-scraping algorithm to gather settlement profiles from the KYC website (as of Dec 2021). In 12 cities that had more than 12 profiled settlements, we manually reviewed and adjusted settlement boundaries to ensure they only included built-up areas (e.g., not water bodies), roughly followed morphological features such as roads (based on satellite imagery collected in the same year), and that settlements did not overlap. We then calculated population per square meter and visually inspected settlements in Google Earth over historical imagery, and dropped any settlement with an unrealistic population count/density. Visual comparison of nearby mapped settlements with similar morphology was a helpful gauge of (in)consistency, and thus accuracy, of population estimates.

## Provided data and files
**kyc_settlement_population_extract:** The python script used to scrape data from the KYC Campaign website.

**ori_...:** Original downloaded settlement boundaries and population counts formatted as a shapefile by city.

**cln_...:** Cleaned settlement boundaries and population counts formatted as a shapefile by city.

## Recommended citation
Slum Dwellers International Profiling Teams. 2022. KnowYourCity data for research. Data processed by Dana R. Thomson and Hazem Mahmoud. Availalbe at: https://github.com/hazemmahmoud88/KnowYourCity-data-for-research.
