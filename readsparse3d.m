function out = readsparse3d(in)

%=============================================================
%
% Reads 3D binary image sets saved as slicewise sparse matrices.
% 
% usage: out = readsparse3d(in);
%   
%       in : input image in cell format with each cell as one
%            (coronal) slice in sparse format
%       out: output as double 3D matrix
%
% e.g. bw = readsparse3d(bw_cell);
% 
% Nigel Chou, March 2010
%=============================================================

zl = size(in,2);
[xl,yl] = size(in{1,1});

out = zeros(xl,yl,zl);
for k=1:zl
    out(:,:,k) = in{1,k};
end