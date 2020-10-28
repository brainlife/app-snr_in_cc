[![Abcdspec-compliant](https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg)](https://github.com/brain-life/abcd-spec)
[![Run on Brainlife.io](https://img.shields.io/badge/Brainlife-bl.app.120-blue.svg)](https://doi.org/10.25663/bl.app.120)

# app-snr_in_cc
Brainlife.io app that computes the signal-to-noise ratio in the corpus callosum in all b-vec directions.

This app receives DWI datatype as input and outputs the signal-to-noise ratio in the corpus callosum in each b-vec direction (and b0) as well as the SNR in the X-most, Y-most, and Z-most directions.  Since the tracts of the corpus callosum run primarily left-to-right/right-to-left, the X-most directions should generally see lower SNRs than the Y- and Z-most directions.

In the white matter, examining diffusion-weighted images (DWI) of the brain can be a useful measure of the local orientation of white matter fascicles, as diffusion in the WM is more free in the directions the fibers run.  The readout signal of a DWI image will be attenuated (quieted) in the directions of greatest water molecule diffusion, and so we can approximate the local orientation of WM fascicles in each voxel of the brain by examining the signal-to-noise ratios (SNR) in different directions.

See [here](https://dipy.org/documentation/1.1.0./examples_built/snr_in_cc/) for more information

### Authors
- David Hunt (davhunt@indiana.edu)
- Soichi Hayashi (hayashis@iu.edu)

### Project director
- Franco Pestilli (franpest@indiana.edu)

### Funding Acknowledgement
brainlife.io is publicly funded and for the sustainability of the project it is helpful to Acknowledge the use of the platform. We kindly ask that you acknowledge the funding below in your publications and code reusing this code.

[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-1636893](https://img.shields.io/badge/NSF_BCS-1636893-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1636893)
[![NSF-ACI-1916518](https://img.shields.io/badge/NSF_ACI-1916518-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1916518)
[![NSF-IIS-1912270](https://img.shields.io/badge/NSF_IIS-1912270-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1912270)
[![NIH-NIBIB-R01EB029272](https://img.shields.io/badge/NIH_NIBIB-R01EB029272-green.svg)](https://grantome.com/grant/NIH/R01-EB029272-01)

### Citations
We kindly ask that you cite the following articles when publishing papers and code using this code.

1. Descoteaux, M., Deriche, R., Le Bihan, D., Mangin, J.-F., and Poupon, C. Multiple q-shell diffusion propagator imaging. Medical Image Analysis, 15(4), 603, 2011. [https://doi.org/10.1016/j.media.2010.07.001](https://doi.org/10.1016/j.media.2010.07.001)

2. Jones, D. K., Knosche, T. R., & Turner, R. White Matter Integrity, Fiber Count, and Other Fallacies: The Dos and Don'ts of Diffusion MRI. NeuroImage, 73, 239, 2013. [https://doi.org/10.1016/j.neuroimage.2012.06.081](https://doi.org/10.1016/j.neuroimage.2012.06.081)

3. Avesani, P., McPherson, B., Hayashi, S. et al. The open diffusion data derivatives, brain data upcycling via integrated publishing of derivatives and reproducible open cloud services. Sci Data 6, 69 (2019). [https://doi.org/10.1038/s41597-019-0073-y](https://doi.org/10.1038/s41597-019-0073-y)

#### MIT Copyright (c) 2020 brainlife.io The University of Texas at Austin and Indiana University


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

All output files will be generated under the current working directory (pwd). The main outputs of this App are the subject's corpus callosum mask `cc.nii.gz` and the subject's brain mask Nifti image called `mask_noise.nii.gz`, as well as a bar graph of the SNR in each direction.

#### SNR.json

This app will also produce a file called `SNR.json` in `output` directory that contains all the SNR data in the X-most, Y-most, and Z-most directions ("SNR in b0, X, Y, Z) as well as in each b-vec direction ("SNR in all directions"). The unit vectors that correspond to the direction #s are found in "direction vectors".

### Dependencies

This App requires the following libraries to run locally.


  - NIBABEL: https://github.com/nipy/nibabel
  - NUMPY: https://github.com/numpy/numpy
  - jsonlab: https://github.com/fangq/jsonlab.git

