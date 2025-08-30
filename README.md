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

## Quick start

Clone the repository:

```bash
git clone https://github.com/fatenzouaghi/rbr-wave-processing.git
cd rbr-wave-processing
```
Run with the provided .rsk file:
```
addpath(genpath("src"));

% Site parameters used in this example:
zmembrane = -1.191;      % m (example)
zbottom = -1.597;      % m
HAB       = zmembrane - bottom;
alti      = 44.5;        % m (met station altitude for this example)

opts = struct('alti',alti,'HAB',HAB,'fs',4,'nfft',1024, ...
              'use_attenuation',true,'minFreq',0.0083,'igCutoff',0.05, ...
              'maxFreq',0.5,'transferStop',0.40);
out  = process_rbr_pressure("data/206599_20220715_2119.rsk","data/meteo.mat",opts);
[out.Hs, out.HsIG, out.HsSW, out.Tp, out.Tm01, out.Tm02]
```
## Data availability

- **Raw RBR file (.rsk)** used in this study is provided in the `data/` folder and archived on Zenodo: DOI: 10.5281/zenodo.XXXXXXX  
- **Meteorological data** for this example are provided as `data/meteo.mat` (also archived on Zenodo) to allow full reproducibility.  

The original meteorological observations were obtained from the official  
[Environment and Climate Change Canada portal](https://climat.meteo.gc.ca/historical_data/search_historic_data_f.html).  
The `meteo.mat` file contains:  
- `Time_UTC` – time vector in UTC  
- `Press` – atmospheric pressure in kPa (converted to Pa in the script)  
- `Temperature` – air temperature in °C  

In this example, the station altitude used for barometric correction is **44.50 m** above sea level (specific to the Grise Fiord station). 


## Code availability

The custom MATLAB scripts are archived on Zenodo and mirrored on GitHub:

Code DOI: 10.5281/zenodo.YYYYYYY

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
