function out = resize_slice(in,factor)

%=====================================================================
% resizes (downsamples) a 3D image along x and y directions 
% (i.e. each slice) leaving z-direction unchanged
%
% usage : out = resize_slice(in,factor);
%
%         out    : resized image
%         in     : 3D input image
%         factor : downsampling factor
%
% e.g. I_res2 = resize_slice(I,2);
% 
% Nigel Chou, Feb 2010
%=====================================================================

for k=1:floor(size(in,1)/factor)
    temp = in((factor*(k-1)+1):(factor*k),:,:);
    in2(k,:,:)=mean(temp,1);
end

for k=1:floor(size(in2,2)/factor)
    temp = in2(:,(factor*(k-1)+1):(factor*k),:);
    out(:,k,:)=mean(temp,2);
end