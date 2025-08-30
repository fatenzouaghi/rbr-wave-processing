## rbr-wave-processing
MATLAB scripts for RBR pressure data processing (.rsk files).
The code corrects atmospheric pressure, computes water levels, and estimates wave parameters (Hs, Tp, Tm01, Tm02, infragravity and sea-swell components).

## Requirements
- MATLAB R2013b or later
- [RSKtools (RBR Global)](https://rbr-global.com/support/matlab-tools/) for reading `.rsk` files

## Data availability

- **Raw RBR file (.rsk)** used in this study is archived on Zenodo: DOI: 10.5281/zenodo.XXXXXXX  
- **Meteorological data (.csv)** for this example are also archived on Zenodo to allow full reproducibility.  

The original meteorological observations were obtained from the official  
[Environment and Climate Change Canada portal](https://climat.meteo.gc.ca/historical_data/search_historic_data_f.html).  

## Usage
This repository contains MATLAB code to process RBR pressure (.rsk) into wave parameters (HS, Hs_IG, Hs_SW, Tp, Tm01, Tm02) and water levels (ABSOLUTE_WL, WL_CGVD2013). Further instructions on using the functions are provided via inline comments and function help blocks throughout the code.
## Author
**Faten Zouaghi** *(University of Quebec at Rimouski)*  
ORCID: [0009-0003-0722-4947]((https://orcid.org/0009-0003-0722-4947)  
Contact: [faten_zouaghi@uqar.ca](mailto:faten_zouaghi@uqar.ca)

## License

This project is released under the MIT License.
