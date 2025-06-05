%% PCNN3D auto brain extraction
% save as *_mask.nii.gz
% requires nifti toolbox
%
% 2017/07/27 ver1.0
% 2017/08/29 ver1.1 bug fix; add ZoomFactor; use mean image for 4D data
% 2017/09/12 ver1.2 save as 8-bit mask
% 2018/01/23 ver1.3 use with command line bash script
% Kai-Hsiang Chuang, QBI/UQ

%% init setup
warning('off','all');
addpath('/Users/uqkchua4/Documents/MATLAB/NIFTI')
addpath('/Users/uqkchua4/Documents/MATLAB/PCNN')
%datpath='/Users/uqkchua4/Documents/MyData/glym/IgG_DPZ/DTI/20230829_MouseDTI_b16/b0.nii.gz';
%datpath='/Users/uqkchua4/Documents/MyData/CaFMRI/EPI_n4.nii.gz'
%BrSize=[300,450]; % brain size range for MOUSE (mm3).
%BrSize=[1500,3000]; % brain size range for RAT (mm3)
if ~exist('StrucRadius')
	StrucRadius=3; % use =3 for low resolution, use 5 or 7 for highres data
end
if ~exist('ZoomFactor')
	ZoomFactor=10; % resolution magnification factor
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run PCNN
[nii] = load_untouch_nii(datpath);
mtx=size(nii.img);
if length(mtx)==4
    disp('Data is 4D, use the average image to generate mask')
    nii.img=mean(nii.img,4);
end
voxdim=nii.hdr.dime.pixdim(2:4);
try
[I_border, G_I, optG] = PCNN3D(single(nii.img), StrucRadius, voxdim, BrSize*ZoomFactor^3);
catch
beep
disp('Error! Cannot converge. Please adjust the BrainSize option')
disp(['No mask created for ',datpath])
exit
end
V=zeros(mtx);
for n=1:mtx(3)
    V(:,:,n)=I_border{optG}{n};
end

%% save data
disp(['Saving mask at ',datpath(1:end-7),'_mask.nii.gz....'])
nii.img=V;
nii.hdr.dime.dim(1)=3; nii.hdr.dime.dim(5)=1;
nii.hdr.dime.datatype=2; nii.hdr.dime.bitpix=8; % save as unsigned char
p=strfind(datpath,'.nii');
save_untouch_nii(nii,[datpath(1:p-1),'_mask.nii.gz'])

disp('Done')
