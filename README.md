## CoastalWaves
MATLAB scripts for wave processing from RBR pressure loggers (.rsk) and Spotter buoys (.csv).
The code corrects atmospheric pressure, computes water levels, and estimates wave parameters (Hs, Tp, Tm01, Tm02, infragravity and sea-swell components).

## Requirements
- MATLAB R2013b or later
- [RSKtools (RBR Global)](https://rbr-global.com/support/matlab-tools/) for reading `.rsk` files

## Data availability

- All datasets used in this study—RBR raw files (.rsk), Spotter buoy files (.csv), and meteorological data (.csv)—are archived on Zenodo (DOI: 10.5281/zenodo.XXXXXXX
).

The original meteorological observations were obtained from the official  
[Environment and Climate Change Canada portal](https://climat.meteo.gc.ca/historical_data/search_historic_data_f.html).  

## Usage
This repository provides MATLAB classes and functions to process raw wave data from RBR pressure loggers (.rsk) and Spotter buoys (.csv). Inline comments describe the processing steps.

## Author
**Faten Zouaghi** *(University of Quebec at Rimouski)*  
ORCID:[0009-0003-0722-4947]((https://orcid.org/0009-0003-0722-4947)  
Contact: [faten_zouaghi@uqar.ca](mailto:faten_zouaghi@uqar.ca)

## License

This repository is licensed under Apache License 2.0.
