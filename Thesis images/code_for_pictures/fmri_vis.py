from nilearn import plotting
from nilearn.image.image import mean_img

img_file = '/Users/matthew/nilearn_data/cobre/fmri_0040012.nii.gz'
img_mean = mean_img(img_file)

display = plotting.plot_epi(img_mean, colorbar=True)

display.savefig('fmri_visualisation.png')   
display.close() 

img_file = '/Users/matthew/Desktop/PhD/PhD/014_misc_data/0040012_mprage.nii.gz'

display = plotting.plot_anat(img_mean)

display.savefig('mri_visualisation.png')   
display.close() 