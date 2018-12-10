[![Abcdspec-compliant](https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg)](https://github.com/brain-life/abcd-spec)
[![Run on Brainlife.io](https://img.shields.io/badge/Brainlife-bl.app.1-blue.svg)](https://doi.org/10.25663/bl.app.1)

# app-snr_in_cc
Brainlife.io app that computes the signal-to-noise ratio in the corpus callosum in all b-vec directions.

This app receives DWI datatype as input and outputs the signal-to-noise ratio in the corpus callosum in each b-vec direction (and b0) as well as the SNR in the X-most, Y-most, and Z-most directions.  Since the tracts of the corpus callosum run primarily left-to-right/right-to-left, the X-most directions should generally see lower SNRs than the Y- and Z-most directions.

In the white matter, examining diffusion-weighted images (DWI) of the brain can be a useful measure of the local orientation of white matter fascicles, as diffusion in the WM is more free in the directions the fibers run.  The readout signal of a DWI image will be attenuated (quieted) in the directions of greatest water molecule diffusion, and so we can approximate the local orientation of WM fascicles in each voxel of the brain by examining the signal-to-noise ratios (SNR) in different directions.

### Authors
- David Hunt (davhunt@indiana.edu)
- Soichi Hayashi (hayashis@iu.edu)

### Project director
- Franco Pestilli (franpest@indiana.edu)

### Funding 
[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-1636893](https://img.shields.io/badge/NSF_BCS-1636893-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1636893)

## Running the App 

### On Brainlife.io

You can submit this App online at [https://doi.org/10.25663/brainlife.app.120](https://doi.org/10.25663/brainlife.app.120) via the "Execute" tab.

### Running Locally (on your machine)

1. git clone this repo.
2. Inside the cloned directory, create `config.json` with something like the following content with paths to your input files.

```json
{
	"dwi": "./input/dtiinit/dwi_aligned_trilin_noMEC.nii.gz",
	"bvals": "./input/dtiinit/dwi_aligned_trilin_noMEC.nii.bvals",
	"bvecs": "./input/dtiinit/dwi_aligned_trilin_noMEC.nii.bvecs"
}
```

3. Launch the App by executing `main`

```bash
./main
```

### Sample Datasets

If you don't have your own input file, you can download sample datasets from Brainlife.io, or you can use [Brainlife CLI](https://github.com/brain-life/cli).

```
npm install -g brainlife
bl login
mkdir input
bl dataset download 5bfee5acc203920043d43ddd && mv 5bfee5acc203920043d43ddd input/dwi
```

## Output

All output files will be generated under the current working directory (pwd). The main outputs of this App are the subject's corpus callosum Nifti image called `cc.nii.gz` and the subject's brain mask Nifti image called `mask_noise.nii.gz`.

#### Product.json

This app will also produce a file called `product.json` that contains all the SNR data in each b-vec direction, as well as the SNR data in the X-most, Y-most, and Z-most directions. This file allows web interfaces, DB and API calls on the results of the processing.

### Dependencies

This App requires the following libraries to run locally.


  - NIBABEL: https://github.com/nipy/nibabel
  - NUMPY: https://github.com/numpy/numpy
  - jsonlab: https://github.com/fangq/jsonlab.git

