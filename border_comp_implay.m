function out = border_comp_implay(in,varargin)

% =========================================================================
%
% Displays borders of mask on individual coronal slices. up to 3 borders
% can be drawn and viewed using image-processing toolbox's 'implay'
%
%
% usage: Out = border_comp_implay(OrigVol,bwvol,col,bwvol2,col2,...)
% 
%   in     : original image volume (must be 3D input image)
%   binmask: binary mask volume
%   col    : color in RGB e.g. [1 0 0] for red
%
% e.g. I_dispborder = border_comp_implay(I,I_PCNNmasks{1,54},[1 0 0]);
%
% Nigel Chou. 15 June 2010
%==========================================================================

% ----- check no. of arguments -------
vararg_l = length(varargin);
if mod(vararg_l,2)~=0
    disp('bw image arguments must be in sets of 2');
end

% ----- Adjust for contrast -------
in = cont_adj(in);

% ----- initialize RGB output matrix -------
[xdim,ydim,zdim] = size(in);
out = zeros(xdim,ydim,3,zdim);

% ----- read binary mask volumes -------
binmask = cell(1,vararg_l/2); col = binmask;
for l = 1:vararg_l/2
    binmask{1,l} = varargin{2*l-1};
    if iscell(binmask{1,l})
        binmask{1,l} = readsparse3d(binmask{1,l});
    end
    col{1,l} = varargin{2*l};
end

for k = 1:zdim
    sliceR = in(:,:,k); sliceG = sliceR; sliceB = sliceR;
    for n = 1:vararg_l/2
        border = bwperim(binmask{1,n}(:,:,k));
        sliceR(border) = col{1,n}(1);
        sliceG(border) = col{1,n}(2);
        sliceB(border) = col{1,n}(3);   
    end
    out(:,:,:,k) = cat(3,sliceR,sliceG,sliceB);
end


implay(out);

