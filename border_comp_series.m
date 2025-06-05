function out = border_comp_series(in,slices,crop,varargin)

% =========================================================================
%
% Displays borders of masks from various automatic methods
% on selected coronal slices. Shows all slices in one figure.
%
% usage: Out = border_comp_series(in,slices,crop,varargin);
% 
%   in     : original image volume (must be 3D input image)
%   slices : list slice numbers to be compared
%   crop   : no. of pixels on each side to be cropped ( [left right top
%            bottom]) e.g. [10 10 20 20]
%   col    : color in RGB e.g. [1 0 0] for red
%
% e.g.
% I_comparison = border_comp_series(I,[2 30 48],[10 10 20 20],I_PCNNmasks{1,54},[1 0 0]);
%
% Nigel Chou. 20 June 2010
%==========================================================================

% ----- check no. of arguments -------
vararg_l = length(varargin);
if mod(vararg_l,2)~=0
    disp('bw image arguments must be in sets of 2');
end

% ----- Adjust for contrast -------
in = in(crop(1)+1:end-crop(2),crop(3)+1:end-crop(4),:);
in = cont_adj(in);
in = flipdim(in,1);

% ----- add space between figures -------
space = 4;

comp_no = vararg_l/2; % number of comparisons

% ----- initialize RGB output matrix -------
[ydim,xdim,zdim] = size(in);
out = ones(ydim*comp_no +space*(comp_no-1),xdim*length(slices)+space*(length(slices)-1),3);

% ----- read binary mask volumes -------
binmask = cell(1,comp_no); colour = binmask;
for l = 1:comp_no
    binmask{1,l} = varargin{2*l-1};
    if iscell(binmask{1,l})
        binmask{1,l} = readsparse3d(binmask{1,l});
    end
    binmask{1,l}=binmask{1,l}(crop(1)+1:end-crop(2),crop(3)+1:end-crop(4),:);
    binmask{1,l}=flipdim(binmask{1,l},1);
    colour{1,l} = varargin{2*l};
end

for k = 1:length(slices)
    for n = 1:comp_no
        sliceR = in(:,:,slices(k)); sliceG = sliceR; sliceB = sliceR;
        border = bwperim(binmask{1,n}(:,:,slices(k)));
        sliceR(border) = colour{1,n}(1);
        sliceG(border) = colour{1,n}(2);
        sliceB(border) = colour{1,n}(3);
        startpix_x = (k-1)*xdim+1+(k>1)*(k-1)*space;
        startpix_y = (n-1)*ydim+1+(n>1)*(n-1)*space;
        out(startpix_y:startpix_y+ydim-1,startpix_x:startpix_x+xdim-1,:) = cat(3,sliceR,sliceG,sliceB);
    end
end


imshow(out);

