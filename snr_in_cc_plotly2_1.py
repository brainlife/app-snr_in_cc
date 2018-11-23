#!/usr/bin/env python

"""

=============================================
SNR estimation for Diffusion-Weighted Images
=============================================

Computing the Signal-to-Noise-Ratio (SNR) of DW images is still an open
question, as SNR depends on the white matter structure of interest as well as
the gradient direction corresponding to each DWI.

In classical MRI, SNR can be defined as the ratio of the mean of the signal
divided by the standard deviation of the underlying Gaussian noise, that is
$SNR = mean(signal) / std(noise)$. The noise standard deviation can be computed
from the background in any of the DW images. How do we compute the mean of the
signal, and what signal?

The strategy here is to compute a 'worst-case' SNR for DWI. Several white
matter structures such as the corpus callosum (CC), corticospinal tract (CST),
or the superior longitudinal fasciculus (SLF) can be easily identified from the
colored-FA (CFA) map. In this example, we will use voxels from the CC, which
have the characteristic of being highly red in the CFA map since they are
mainly oriented in the left-right direction. We know that the DW image closest
to the X-direction will be the one with the most attenuated diffusion signal.
This is the strategy adopted in several recent papers (see [Descoteaux2011]_
and [Jones2013]_). It gives a good indication of the quality of the DWI data.

First, we compute the tensor model in a brain mask (see the :ref:`reconst_dti`
example for further explanations).

"""

#from __future__ import division, print_function
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.morphology import binary_dilation

#from dipy.data import fetch_stanford_hardi, read_stanford_hardi
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.segment.mask import median_otsu
from dipy.reconst.dti import TensorModel

from dipy.segment.mask import segment_from_cfa
from dipy.segment.mask import bounding_box

import sys
import json

## load data
#fetch_stanford_hardi()
#img, gtab = read_stanford_hardi()

img = nib.load(sys.argv[1])
bvals, bvecs = read_bvals_bvecs(sys.argv[2], sys.argv[3])
gtab = gradient_table(bvals, bvecs)

data = img.get_data()
affine = img.affine

print('Computing brain mask...')
b0_mask, mask = median_otsu(data)

print('Computing tensors...')
tenmodel = TensorModel(gtab)
tensorfit = tenmodel.fit(data, mask=mask)

"""Next, we set our red-green-blue thresholds to (0.6, 1) in the x axis
and (0, 0.1) in the y and z axes respectively.
These values work well in practice to isolate the very RED voxels of the cfa map.

Then, as assurance, we want just RED voxels in the CC (there could be
noisy red voxels around the brain mask and we don't want those). Unless the brain
acquisition was badly aligned, the CC is always close to the mid-sagittal slice.

The following lines perform these two operations and then saves the computed mask.
"""

print('Computing worst-case/best-case SNR using the corpus callosum...')

threshold = (0.6, 1, 0, 0.1, 0, 0.1)
CC_box = np.zeros_like(data[..., 0])

mins, maxs = bounding_box(mask)
mins = np.array(mins)
maxs = np.array(maxs)
diff = (maxs - mins) // 4
bounds_min = mins + diff
bounds_max = maxs - diff

CC_box[bounds_min[0]:bounds_max[0],
       bounds_min[1]:bounds_max[1],
       bounds_min[2]:bounds_max[2]] = 1

mask_cc_part, cfa = segment_from_cfa(tensorfit, CC_box, threshold,
    return_cfa=True)

cfa_img = nib.Nifti1Image((cfa*255).astype(np.uint8), affine)
mask_cc_part_img = nib.Nifti1Image(mask_cc_part.astype(np.uint8), affine)
nib.save(mask_cc_part_img, 'cc.nii.gz')

#region = 40
# fig = plt.figure('Corpus callosum segmentation')
# plt.subplot(1, 2, 1)
# plt.title("Corpus callosum (CC)")
# plt.axis('off')
# red = cfa[..., 0]
# plt.imshow(np.rot90(red[region, ...]))

# plt.subplot(1, 2, 2)
# plt.title("CC mask used for SNR computation")
# plt.axis('off')
# plt.imshow(np.rot90(mask_cc_part[region, ...]))
#fig.savefig("CC_segmentation.png", bbox_inches='tight')

"""Now that we are happy with our crude CC mask that selected voxels in the x-direction,
we can use all the voxels to estimate the mean signal in this region.

"""

mean_signal = np.mean(data[mask_cc_part], axis=0)

"""Now, we need a good background estimation. We will re-use the brain mask
computed before and invert it to catch the outside of the brain. This could
also be determined manually with a ROI in the background.
[Warning: Certain MR manufacturers mask out the outside of the brain with 0's.
One thus has to be careful how the noise ROI is defined].
"""

mask_noise = binary_dilation(mask, iterations=10)
mask_noise[..., :mask_noise.shape[-1]//2] = 1
mask_noise = ~mask_noise
mask_noise_img = nib.Nifti1Image(mask_noise.astype(np.uint8), affine)
nib.save(mask_noise_img, 'mask_noise.nii.gz')

noise_std = np.std(data[mask_noise, :])
print('Noise standard deviation sigma= ', noise_std)

"""We can now compute the SNR for each DWI. For example, report SNR
for DW images with gradient direction that lies the closest to
the X, Y and Z axes.
"""

# Exclude null bvecs from the search
idx = np.sum(gtab.bvecs, axis=-1) == 0
gtab.bvecs[idx] = np.inf
#axis_X = np.argmin(np.sum((gtab.bvecs-np.array([1, 0, 0]))**2, axis=-1))
#axis_Y = np.argmin(np.sum((gtab.bvecs-np.array([0, 1, 0]))**2, axis=-1))
#axis_Z = np.argmin(np.sum((gtab.bvecs-np.array([0, 0, 1]))**2, axis=-1))
dist_x = []
dist_y = []
dist_z = []

for i in range(0, len(gtab.bvecs)):
	if gtab.bvecs[i][0] >= 0.0:
		dist_x.append(np.sqrt((gtab.bvecs[i][0]-1)**2 + gtab.bvecs[i][1]**2 + gtab.bvecs[i][2]**2))
	elif gtab.bvecs[i][0] < 0.0:
		dist_x.append(np.sqrt((gtab.bvecs[i][0]+1)**2 + gtab.bvecs[i][1]**2 + gtab.bvecs[i][2]**2))

for i in range(0, len(gtab.bvecs)):
	if gtab.bvecs[i][1] >= 0.0:
		dist_y.append(np.sqrt(gtab.bvecs[i][0]**2 + (gtab.bvecs[i][1]-1)**2 + gtab.bvecs[i][2]**2))
	elif gtab.bvecs[i][1] < 0.0:
		dist_y.append(np.sqrt(gtab.bvecs[i][0]**2 + (gtab.bvecs[i][1]+1)**2 + gtab.bvecs[i][2]**2))

for i in range(0, len(gtab.bvecs)):
	if gtab.bvecs[i][2] >= 0.0:
		dist_z.append(np.sqrt(gtab.bvecs[i][0]**2 + gtab.bvecs[i][1]**2 + (gtab.bvecs[i][2]-1)**2))
	elif gtab.bvecs[i][2] < 0.0:
		dist_z.append(np.sqrt(gtab.bvecs[i][0]**2 + gtab.bvecs[i][1]**2 + (gtab.bvecs[i][2]+1)**2))

x_list, y_list, z_list = [], [], []
for i in range(0, len(gtab.bvecs)):
	if dist_x[i] < dist_y[i] and dist_x[i] < dist_z[i]:
		x_list.append(i)
	elif dist_y[i] < dist_x[i] and dist_y[i] < dist_z[i]:
		y_list.append(i)
	elif dist_z[i] < dist_x[i] and dist_z[i] < dist_y[i]:
		z_list.append(i)

x_sorted, y_sorted, z_sorted = [], [], []
for i in x_list:
	x_sorted.append((i, dist_x[i]))
x_sorted = sorted(x_sorted, key=lambda tup: tup[1])

for i in y_list:
	y_sorted.append((i, dist_y[i]))
y_sorted = sorted(y_sorted, key=lambda tup: tup[1])

for i in z_list:
	z_sorted.append((i, dist_z[i]))
z_sorted = sorted(z_sorted, key=lambda tup: tup[1])

bvecs_sorted = []
bvecs_sorted.append(x_sorted)
bvecs_sorted.append(y_sorted)
bvecs_sorted.append(z_sorted)

colors = []
direction = []
for i in range(0, len(bvecs_sorted[0])):
	colors.append("#FF0000")
	direction.append("X, ")
for i in range(0, len(bvecs_sorted[1])):
	colors.append("#0000FF")
	direction.append("Y, ")
for i in range(0, len(bvecs_sorted[2])):
	colors.append("#00FF00")
	direction.append("Z, ")

SNR_output = []
directions = []

for j in range(0, len(bvecs_sorted)):
	for i in range(0, len(bvecs_sorted[j])):
		SNR = mean_signal[bvecs_sorted[j][i][0]]/noise_std
		SNR_output.append(SNR)
		directions.append(gtab.bvecs[bvecs_sorted[j][i][0]])

dirxs = []
for i in range(0, len(directions)):
	dirxs.append(str(directions[i]))
drrxs = []
for i in range(0, len(dirxs)):
	drrxs.append(direction[i] + dirxs[i])





#SNR_output = []
#directions = []
#b0 = True
#for direction in range(0, len(gtab.bvecs)):
    #SNR = mean_signal[direction]/noise_std
    #if gtab.bvecs[direction][0] == np.inf and b0 == True:
        #print("SNR for the b=0 image is :", SNR)
        #SNR_output.append(SNR)
        #directions.append('b0')
        #b0 = False
    #elif gtab.bvecs[direction][0] != np.inf:
        #print("SNR for direction", direction, " ", gtab.bvecs[direction], "is :", SNR)
        #SNR_output.append(SNR)
        #directions.append(str(direction) + ', ' + str(gtab.bvecs[direction]))





#data = []
results = {"brainlife": []}

results['brainlife'].append({
	"type": "plotly",
	"layout": {
		"yaxis": {
			"title": "SNR"
		},
		"xaxis": {
			"type": "category"
		},
		"title": "SNRs in each directions"
	},
	"name": "SNRs in different directions",
	"data": [
	{
		"opacity": 0.6,
		"text": [
			drrxs
		],
		"marker": {
			"color": colors,
			"line": {
				"color": "rgb(8,48,107)",
				"width": 1.5
			}
		},
		"y": SNR_output,
		"x": directions,
		"type": "bar"
	}
	]
})

#data.append({
            #"SNR data": SNR_output,
            #"directions": directions
            #})

with open("product.json", "w") as out_file:
    json.dump(results, out_file)

"""SNR for the b=0 image is : ''42.0695455758''"""
"""SNR for direction 58  [ 0.98875  0.1177  -0.09229] is : ''5.46995373635''"""
"""SNR for direction 57  [-0.05039  0.99871  0.0054406] is : ''23.9329492871''"""
"""SNR for direction 126 [-0.11825  -0.039925  0.99218 ] is : ''23.9965694823''"""

"""

Since the CC is aligned with the X axis, the lowest SNR is for that gradient
direction. In comparison, the DW images in the perpendical Y and Z axes have a
high SNR. The b0 still exhibits the highest SNR, since there is no signal
attenuation.

Hence, we can say the Stanford diffusion data has a 'worst-case' SNR of
approximately 5, a 'best-case' SNR of approximately 24, and a SNR of 42 on the
b0 image.

"""

"""
References
----------

.. [Descoteaux2011] Descoteaux, M., Deriche, R., Le Bihan, D., Mangin, J.-F.,
   and Poupon, C. Multiple q-shell diffusion propagator imaging. Medical Image
   Analysis, 15(4), 603, 2011.

.. [Jones2013] Jones, D. K., Knosche, T. R., & Turner, R. White Matter
   Integrity, Fiber Count, and Other Fallacies: The Dos and Don'ts of Diffusion
   MRI. NeuroImage, 73, 239, 2013.

"""
