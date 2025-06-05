# PCNN3D
The 3D Pulse-coupled Neural Networks (PCNN)  method allows automatic brain extraction of the rodent (rats or mice) brain. It is optimized for T2-weighted MRI (eg, 2D multislice or 3D Fast Spin-Echo or FLASH). The default setting is for mouse brain. It can be adapted to rat brain by changing the threshold of brain volume to account for the larger brain. This software is free for academic and non-commercial use.

Download

Matlab code (ver 1.2) and an easy code (i.e., PCNN3D_run_v1_3.m) for running PCNN3D on Nifti data.

* This requires a Nifti toolbox (eg, https://au.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image).

Usage


First load the image data as a 3D matrix into Matlab. For Analyze/NIFTI format, one can use the PCNN3D_run_v1_3.m in the package.

Then use the following command to generate the brain mask:

[I_border, G_I, optG] = PCNN3D(I, p, voxeldim, BrSize, maxiter, OptMask)

I: 3D matrix of the brain image

p: radius of the structural element (in pixel)

voxeldim: [x, y, z] voxel dimensions of the data set (in mm)

BrSize: [min max] brain volume sizes (in mm3) for estimating optimal iteration; default values [100 550] is for mouse brain. For adult rat brain, one may try [1200 4400]

maxiter: maximum iteration. Default is 200

OptMask: output only the mask based on the optimal guess. Note: the guess may not always give the best result.

The outputs:

I_border: A cell array that stores the binary mask at each iteration.

G_I: the plot of volume of mask (in mm3) vs. iterations, like the following.

optG: the suggested optimal iteration. 

For example, if the optimal iteration is #37, I_border{37} will be the binary mask which can be multiplied with the input image, I, to produce the extracted brain.

For more information of each command, type help command in Matlab command line. For example, type help PCNN3D will show all the adjustable options.

Tips
1. Set proper brain volume range is most critical for getting good result. PCNN3D will fail to converge if the volume range is not right. If your data has only partial brain coverage, reduce the max volume.

2. If the data has voxel size magnified by 10x, the volume range should be magnified by 1000x. If you use the easy code to run, then just set the ZoomFactor to 10 instead of changing the volume range.

3. For optimal results, coil B1 field inhomogeneity should be corrected. Tools like N3 or N4 typically works well for rodent brain.

4. Effect of "Radius of the structural element" can be seen below. Note: larger radius, slower the processing.

5. To run Matlab version in command-line or inside a shell script, use another easy code and run command like the following:

matlab -nodisplay -nodesktop -r "datpath='/my_data_path/epi_mean.nii.gz'; BrSize=[600,1000]; run /pcnn_path/PCNN3D_run_v1_3.m; exit";

6. The max brain size is a key factor to prevent the algorithm expanding beyond the size you expected. You can determine the optimal range by looking at the brain mask size of good extraction results:

fslstats brain_mask.nii.gz -V

Version History

[15 Apr, 2016] v.1.2: bug fix for larger brain volume like rat. Add 'OptMask' input option. Add 'optG' output.

[Aug, 2017] easy code for processing Nifti format is provided.

[May, 2019] Update easy codes for Matlab compatibility and running in command-line.

Citation

N Chou, J Wu, J Bai, A Qiu, KH Chuang. “Robust automatic rodent brain extraction using 3D Pulse-coupled NeuralNetworks (PCNN).” IEEE Trans Imag Proc 20(9):2554-64, 2011.

