<h1 align="center">Compute SNR on Corpus Callosum</h1>

<p align="center">
  The signal-to-noise ratio of a DWI in the corpus callosum, for the b0 and every gradient direction, with DIPY, plus the masks used to measure it.
</p>

<p align="center">
  <a href="https://doi.org/10.25663/brainlife.app.120"><img src="https://img.shields.io/badge/brainlife.io-run%20this%20App-2b6cb0" alt="Run on brainlife.io"></a>
  <a href="https://doi.org/10.25663/brainlife.app.120"><img src="https://img.shields.io/badge/DOI-10.25663%2Fbrainlife.app.120-blue" alt="DOI"></a>
  <a href="https://github.com/brainlife/abcd-spec"><img src="https://img.shields.io/badge/ABCD%20spec-v1.1-green" alt="ABCD spec v1.1"></a>
  <a href="https://dipy.org/"><img src="https://img.shields.io/badge/DIPY-0.16.0-orange" alt="DIPY 0.16.0"></a>
  <img src="https://img.shields.io/badge/license-MIT-lightgrey" alt="MIT license">
</p>

## 🧠 At a glance

<!-- brainlife-app:begin (generated from the App manifest; edit the manifest, not this block) -->

> **Computes the signal-to-noise ratio of DWI in the corpus callosum for the b0 and each gradient direction with DIPY, and outputs SNR statistics and the corpus callosum and noise masks.**

Fits a tensor, segments the corpus callosum from the colour-FA map within a central bounding box (DIPY segment_from_cfa), estimates noise from a background mask (median_otsu) and reports SNR for b0 and for each gradient direction, including the directions closest to x, y and z. Use it as a QC measure of DWI data quality. Follows the DIPY snr_in_cc example.

| | |
|---|---|
| 📥 **Inputs** | `neuro/dwi` |
| 📤 **Outputs** | `raw` · `neuro/mask` ×2 · `neuro/snr-stats` |
| 🐳 **Container** | `brainlife/dipy:0.16.0` |
| ⚙️ **Resources** | 8 CPUs |
| 🔗 **Run it** | [doi.org/10.25663/brainlife.app.120](https://doi.org/10.25663/brainlife.app.120) |

#### 📥 Inputs

| Input | Datatype | Required | Description |
|---|---|---|---|
| `dwi` | `neuro/dwi` | yes | Diffusion-weighted image with bvals and bvecs. |

#### 📤 Outputs

| Output | Datatype | Description |
|---|---|---|
| `output` | `raw` | Raw SNR results (JSON/CSV) and intermediate files. |
| `cc_mask` | `neuro/mask` | Corpus callosum mask in DWI space used to measure signal. |
| `noise_mask` | `neuro/mask` | Background noise mask in DWI space used to estimate noise. |
| `snr-stats` | `neuro/snr-stats` | SNR statistics for b0 and each gradient direction. |

<!-- brainlife-app:end -->

## 🔬 How it works

The App follows the DIPY [SNR in the corpus callosum](https://dipy.org/documentation/1.1.0./examples_built/snr_in_cc/) example:

1. A tensor model is fitted and the corpus callosum is segmented from the colour-FA map inside a central bounding box (`segment_from_cfa`).
2. The noise is estimated from a background mask (`median_otsu`).
3. The SNR is reported for the b0 and for each gradient direction, and for the directions closest to the x, y and z axes.

Diffusion is freer along white matter fibres, so the diffusion-weighted signal is attenuated most in the direction the fibres run. The fibres of the corpus callosum run mainly left-right, so the x-most directions should generally show a lower SNR than the y- and z-most directions. Use the App as a quality-control measure of DWI data.

## 🚀 Running the App

### ☁️ On brainlife.io

Open the App at [doi.org/10.25663/brainlife.app.120](https://doi.org/10.25663/brainlife.app.120), choose a `neuro/dwi` dataset from your project and press **Execute**. The App writes the SNR statistics (`neuro/snr-stats`), the corpus callosum and noise masks (`neuro/mask`, tagged `CC` and `noise`) and the raw results (`raw`, tagged `snr-cc`), and shows a bar graph of the SNR in each direction.

### 💻 Locally

You need [Singularity/Apptainer](https://apptainer.org/) and [`jq`](https://jqlang.github.io/jq/). `main` runs the SNR computation in the `brainlife/dipy:0.16.0` container and the csv conversion in `brainlife/freesurfer-stats:1.2`.

1. Clone this repository and check out the App's branch: `git checkout plot`.
2. Create `config.json` with the paths to your DWI files:

   ```json
   {
       "dwi": "./input/dwi/dwi.nii.gz",
       "bvals": "./input/dwi/dwi.bvals",
       "bvecs": "./input/dwi/dwi.bvecs"
   }
   ```

3. Run the App:

   ```bash
   ./main
   ```

### 🧪 Sample data

Download a sample DWI from brainlife.io with the [brainlife CLI](https://github.com/brainlife/cli):

```bash
npm install -g brainlife
bl login
mkdir input
bl dataset download 5bfee5acc203920043d43ddd && mv 5bfee5acc203920043d43ddd input/dwi
```

## 📦 Outputs

All output files are written under the current working directory:

```text
output/
└── snr.json            SNR in b0 and the x-, y- and z-most directions, SNR in every
                        direction, and the unit vector of each direction
snr-stats/
└── snr.csv             the same SNR values as a table
cc_mask/
└── mask.nii.gz         corpus callosum mask used to measure the signal
noise_mask/
└── mask.nii.gz         background mask used to estimate the noise
product.json            bar graph of the SNR in each direction, shown on brainlife.io
```

## 📚 Citing

> [!TIP]
> If you use this App in your work, please cite the App, the software it builds on and brainlife.io (below). GitHub's **Cite this repository** button gives the App citation in APA and BibTeX.

- Descoteaux M, Deriche R, Le Bihan D, Mangin J-F, Poupon C. Multiple q-shell diffusion propagator imaging. *Medical Image Analysis* 15(4):603 (2011). [doi:10.1016/j.media.2010.07.001](https://doi.org/10.1016/j.media.2010.07.001)
- Jones DK, Knösche TR, Turner R. White matter integrity, fiber count, and other fallacies: the do's and don'ts of diffusion MRI. *NeuroImage* 73:239 (2013). [doi:10.1016/j.neuroimage.2012.06.081](https://doi.org/10.1016/j.neuroimage.2012.06.081)
- Compute SNR on Corpus Callosum. brainlife.io App. [doi:10.25663/brainlife.app.120](https://doi.org/10.25663/brainlife.app.120)
- Garyfallidis et al. Dipy, a library for the analysis of diffusion MRI data. *Frontiers in Neuroinformatics* 8 (2014). [doi:10.3389/fninf.2014.00008](https://doi.org/10.3389/fninf.2014.00008)
- Hayashi S, Caron BA, Heinsfeld AS, et al. brainlife.io: a decentralized and open-source cloud platform to support neuroscience research. *Nature Methods* 21:809–813 (2024). [doi:10.1038/s41592-024-02237-2](https://doi.org/10.1038/s41592-024-02237-2)
- Avesani P, McPherson B, Hayashi S, et al. The open diffusion data derivatives, brain data upcycling via integrated publishing of derivatives and reproducible open cloud services. *Scientific Data* 6:69 (2019). [doi:10.1038/s41597-019-0073-y](https://doi.org/10.1038/s41597-019-0073-y)

## 👥 Authors

- Franco Pestilli (project director)
- David Hunt
- Soichi Hayashi

## 💰 Funding

[![Wellcome-226486/Z/22/Z](https://img.shields.io/badge/Wellcome-226486%2FZ%2F22%2FZ-purple.svg)](https://wellcome.org/)
[![NIH-NINDS-U24NS140384](https://img.shields.io/badge/NIH_NINDS-U24NS140384-green.svg)](https://reporter.nih.gov/project-details/11261037)
[![NIH-NINDS-UM1NS132207](https://img.shields.io/badge/NIH_NINDS-UM1NS132207-green.svg)](https://reporter.nih.gov/project-details/11411328)
[![NIH-NIBIB-R01EB029272](https://img.shields.io/badge/NIH_NIBIB-R01EB029272-green.svg)](https://reporter.nih.gov/project-details/10262925)
[![NIH-NIBIB-R01EB030896](https://img.shields.io/badge/NIH_NIBIB-R01EB030896-green.svg)](https://reporter.nih.gov/project-details/10428625)
[![NIH-NIMH-R01MH126699](https://img.shields.io/badge/NIH_NIMH-R01MH126699-green.svg)](https://reporter.nih.gov/project-details/10460628)
[![NIH-NIDCD-R01DC020088](https://img.shields.io/badge/NIH_NIDCD-R01DC020088-green.svg)](https://reporter.nih.gov/project-details/11140284)
[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-2203524](https://img.shields.io/badge/NSF_BCS-2203524-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2203524)
[![NSF-IIS-1636893](https://img.shields.io/badge/NSF_IIS-1636893-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1636893)
[![NSF-OAC-2148729](https://img.shields.io/badge/NSF_OAC-2148729-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2148729)
[![NSF-IIS-1912270](https://img.shields.io/badge/NSF_IIS-1912270-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1912270)
[![NSF-IIS-2148700](https://img.shields.io/badge/NSF_IIS-2148700-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2148700)
[![NSF-OAC-1916518](https://img.shields.io/badge/NSF_OAC-1916518-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1916518)

<sub>brainlife.io is publicly funded, and acknowledging that support helps keep it running. This App was developed with support from Wellcome, the U.S. National Institutes of Health (including the BRAIN Initiative BRAIN CONNECTS program: The Axonal Projectome EXchange, U24NS140384, and the Center for Mesoscale Connectomics, UM1NS132207) and the U.S. National Science Foundation.</sub>

## ⚖️ License

MIT, Copyright (c) 2020 brainlife.io, The University of Texas at Austin and Indiana University. DIPY is distributed under its own license.
