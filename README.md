# fmri_preproc
----------------

[![Documentation Status](https://readthedocs.org/projects/fmri-preproc/badge/?version=latest)](https://fmri-preproc.readthedocs.io/en/latest/?badge=latest)         


A resting-state fMRI preprocessing pipeline for neonatal neuroimaging data.

# Quick start
----------------

## External dependencies        

| Dependency  | Environmental variable (if applicable)  |
|---|---|
| [`FSL`](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/)  | `FSLDIR`  |
| [`ANTs`](http://stnava.github.io/ANTs/)  | `ANTSPATH`  |
| [`C3D`](http://www.itksnap.org/pmwiki/pmwiki.php?n=Convert3D.Documentation)  | ***Must be in system path**  |
| [`FIX`](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX)  | ***Must be in system path**  |
| [The dHCP structural pipeline](https://github.com/BioMedIA/dhcp-structural-pipeline)  | N/A  |

> If at CCHMC:    
>            
> ```bash
> # Load modules
> module load anaconda3/1.0.0
> module load fsl/6.0.4
> module load ants/2.3.1
> module load c3d/1.0.0
> module load R/4.0.2
> module load octave/6.2.0
> module load fix/1.6.15
> module load cuda/9.1
> ```        
>
>             
> `FIX` uses [`GNU Octave`](https://www.gnu.org/software/octave/index). Make sure the executable can be called from the command line prior to running the pipeline.

## Install python dependencies

If not done so already, install the required dependencies by typing:

```bash
git clone https://github.com/AdebayoBraimah/fmri_preproc.git

pkg_path=$(realpath fmri_preproc)

cd ${pkg_path}

pip install -r requirements.txt
```

## Setting environmental variables

**NOTE**: Only append to the `LD_LIBRARY_PATH` **IF** python issues arise when calling FSL executables.        

```bash
git clone https://github.com/AdebayoBraimah/fmri_preproc.git
pkg_path=$(realpath fmri_preproc)

# Set path variables
export ANTSPATH=$(dirname $(which antsRegistrationSyN.sh))/ # If not set already

# Append to Linker and Python paths
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FSLDIR}/fslpython/envs/fslpython/lib
export PYTHONPATH=${PYTHONPATH}:$(which python3):${pkg_path}

# Source FIX settings file
source $(dirname $(which fix))/settings.sh

# Append package executable to PATH
export PATH=${PATH}:${pkg_path}
```

# Usage
----------------

Once the environmental variables are set, simply typing `fmri_preproc.py` or `fmri_preproc.py -h` should pull up the main help menu. 

The pipeline is split into four main sections (plus the `run-all` option) - each accessible via a subcommand (the first positional argument) with their own separate help and usage menus.

The main help menu is shown below:        

```
usage: fmri_preproc.py [-h] [-s <str>] [-r <str>] [-o <str>] [-e <str>] [-version]
                       [--write-config <file>] [-v] [--log-level <str>] [--config <file>]
                       {import-data,pre-mcdc,mcdc,post-mcdc,run-all} ...

Performs preprocessing of neonatal resting-state functional MR neuroimages.

optional arguments:
  -h, --help                                     show this help message and exit

subcommands:
  Data preprocessing pipeline sections for 'fmri_preproc'.

  {import-data,pre-mcdc,mcdc,post-mcdc,run-all}  Type the 'subcommand' name followed by '-h' or '--
                                                 help' for more information.
    import-data                                  Imports the subject's imaging data into a working
                                                 directory for image preprocessing.
    pre-mcdc                                     Pre-Motion Correction and Distortion Correction
                                                 (MCDC).
    mcdc                                         Motion Correction and Distortion Correction (MCDC).
    post-mcdc                                    Post-MCDC preprocessing, which includes: func ->
                                                 template registration, ICA-based denoising, and QC
                                                 (web page).
    run-all                                      Perform/Run all preprocessing steps of the
                                                 pipeline.

Required Arguments:
  -s <str>, --sub <str>                          REQUIRED: Subject ID.
  -r <str>, --run <str>                          REQUIRED: Run number.
  -o <str>, --outdir <str>                       REQUIRED: Output parent directory.

Optional Arguments:
  -e <str>, --ses <str>                          OPTIONAL: Session ID.
  -version, --version                            Prints version, then exits.
  --write-config <file>                          Writes a template configuration (JSON) to file.
  -v, --verbose                                  Enables verbose output to the command and log
                                                 files.
  --log-level <str>                              Log level to be written to output log files.
                                                 Acceptable inputs include: 'info', 'debug',
                                                 'critical', 'warning', and 'error'
                                                 [default='info'].
  --config <file>, --settings <file>             Configuration (JSON) file that contains additional
                                                 input arguments not specified at the command line.
```

# How to cite
----------------

This preprocessing pipeline is based off of the [`dhcp-fmri preprocessing pipeline`](https://git.fmrib.ox.ac.uk/seanf/dhcp-neonatal-fmri-pipeline).

When using this pipeline, please be sure to cite the appropriate reference shown below:

**Primary citation**:

Fitzgibbon, SP, Harrison, SJ, Jenkinson, M, Baxter, L, Robinson, EC, Bastiani, M, Bozek, J, Karolis, V, Cordero Grande, L, Price, AN, Hughes, E, Makropoulos, A, Passerat-Palmbach, J, Schuh, A, Gao, J, Farahibozorg, S, O'Muircheartaigh, J, Ciarrusta, J, O'Keeffe, C, Brandon, J, Arichi, T, Rueckert, D, Hajnal, JV, Edwards, AD, Smith, SM, *Duff, E, *Andersson, J  **The developing Human Connectome Project automated functional processing framework for neonates**., *NeuroImage (2020), 223: 117303*, 2020. doi: [https://doi.org/10.1016/j.neuroimage.2020.117303](https://doi.org/10.1016/j.neuroimage.2020.117303) *Authors contributed equally.

```bibtex
@article {Fitzgibbon766030,
	author = {Fitzgibbon, Sean P. and Harrison, Samuel J. and Jenkinson, Mark and Baxter, Luke and Robinson, Emma C. and Bastiani, Matteo and Bozek, Jelena and Karolis, Vyacheslav and Grande, Lucilio Cordero and Price, Anthony N. and Hughes, Emer and Makropoulos, Antonios and Passerat-Palmbach, Jonathan and Schuh, Andreas and Gao, Jianliang and Farahibozorg, Seyedeh-Rezvan and O{\textquoteright}Muircheartaigh, Jonathan and Ciarrusta, Judit and O{\textquoteright}Keeffe, Camilla and Brandon, Jakki and Arichi, Tomoki and Rueckert, Daniel and Hajnal, Joseph V. and Edwards, A. David and Smith, Stephen M. and Duff, Eugene and Andersson, Jesper},
	title = {The developing Human Connectome Project (dHCP) automated resting-state functional processing framework for newborn infants},
	elocation-id = {117303},
	year = {2020},
	doi = {10.1016/j.neuroimage.2020.117303},
	publisher = {Elsevier},
	URL = {https://doi.org/10.1016/j.neuroimage.2020.117303},
	eprint = {https://www.sciencedirect.com/science/article/pii/S1053811920307898/pdfft?md5=18806cf190a26f783de4bef456fe28b6&pid=1-s2.0-S1053811920307898-main.pdf},
	journal = {NeuroImage}
}
```
