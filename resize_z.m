function out = resize_z(in,factor)

%=====================================================================
% Resizes (downsamples) a 3D image in the z-direction (through-plane)
%
% usage : out = resize_z(in,factor);
%
%         out    : resized image
%         in     : 3D input image
%         factor : downsampling factor
%
% e.g. I_resz2 = resize_z(I,2);
% 
% Nigel Chou, Feb 2010
%=====================================================================

for k=1:floor(size(in,3)/factor)
    temp = in(:,:,(factor*(k-1)+1):(factor*k));
    out(:,:,k)=mean(temp,3);
end


