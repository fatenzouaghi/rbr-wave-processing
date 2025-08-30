## rbr-wave-processing
MATLAB scripts for RBR pressure data processing.
The code corrects atmospheric pressure, computes water levels, and estimates wave parameters (Hs, Tp, Tm01, Tm02, infragravity and sea-swell components).
## Repository structure
- `README.md` — project description (this file)
- `LICENSE` — license (MIT)
- `CITATION.cff` — citation metadata (GitHub/Zenodo)
- `config.m` — configuration file (paths, constants)
- `src/` — MATLAB functions 
- `examples/` — example script(s)
- `data/` — raw RBR file (`.rsk`) and `meteo.mat`
- `docs/code-availability.md` — statement for the manuscript

## Requirements
- MATLAB R2013b or later
- [RSKtools (RBR Global)](https://rbr-global.com/support/matlab-tools/) for reading `.rsk` files

## Data availability

- **Raw RBR file (.rsk)** used in this study is archived on Zenodo: DOI: 10.5281/zenodo.XXXXXXX  
- **Meteorological data** for this example are  archived on Zenodo to allow full reproducibility.  

The original meteorological observations were obtained from the official  
[Environment and Climate Change Canada portal](https://climat.meteo.gc.ca/historical_data/search_historic_data_f.html).  
The `meteo.mat` file contains:  
- `Time_UTC` – time vector in UTC  
- `Press` – atmospheric pressure in kPa (converted to Pa in the script)  
- `Temperature` – air temperature in °C  

In this example, the station altitude used for barometric correction is **44.50 m** above sea level (specific to the Grise Fiord station). 


## Code availability

The custom MATLAB scripts are mirrored on GitHub:
GitHub repository: https://github.com/fatenzouaghi/rbr-wave-processing


## Citation

If you use this code, please cite:

@software{zouaghi2025rbr,
  author       = {Zouaghi, Faten},
  title        = {RBR Wave Processing (MATLAB)},
  year         = {2025},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.YYYYYYY},
  url          = {https://github.com/fatenzouaghi/rbr-wave-processing}
}


## License

This project is released under the MIT License.
