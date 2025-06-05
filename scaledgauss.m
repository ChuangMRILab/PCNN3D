function h = scaledgauss( r,voxdim,sigma)

%=====================================================================
% creates a scaled gaussian weighting matrix for 3D PCNN
% intended to correct for anisotropy in voxel dimensions
%
% usage : h = scaledgauss( r,voxdim,sigma)
%
%         h     : scaled gaussina weighting matrix of size
%                 [2r+1 2r+2 2r+1]
%         r     : radius
%         voxdim: voxel dimensions
%         sigma : spread of gaussian function
%
% e.g. h = scaledgauss( 3,[.1 .1 .3],.65);
% 
% Nigel Chou, July 2009
%=====================================================================

% Scale pixel dimensions so that smallest dimension is 1 
voxdim = voxdim./min(voxdim);

if length(voxdim)==2
    xl=voxdim(1); yl=xl;
    zl=voxdim(2);
elseif length(voxdim)==3
     xl=voxdim(1); 
     yl=voxdim(2);
    zl=voxdim(3);
elseif length(voxdim)<2 || length(voxdim)>3
    disp('please enter a 3-vector for pixel dimensions')
end

% ---- initialize matrix -------
R2 = zeros( 2*r+1, 2*r+1, 2*r+1 );

% ---- find distance from center voxel --------
for x=1:2*r+1
    for y=1:2*r+1
        for z=1:2*r+1
            R2(x,y,z) =  (xl.*(r+1-x)).^2 + (yl.*(r+1-y)).^2 + (zl.*(r+1-z)).^2 ;
        end
    end
end

% ------ gaussian function ---------
h = exp(-R2./(2.*sigma^2));
h(h<eps*max(h(:))) = 0; % remove very small values

% ------ scale gaussian such that sum of all elements =1 ---------
sumh = sum(h(:));
if sumh ~= 0,
    h  = h/sumh;
end;