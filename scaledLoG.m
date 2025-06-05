function h = scaledLoG( r,voxdim,sigma)

%=====================================================================
% creates a scaled laplacian of gaussian weighting matrix for 3D PCNN
% (not used in paper)
%
% usage : h = scaledLoG( r,voxdim,sigma)
%
%         h     : scaled weighting matrix of size
%                 [2r+1 2r+2 2r+1]
%         r     : radius
%         voxdim: voxel dimensions
%         sigma : spread of function
%
% e.g.  h = scaledLoG( 2,[.1 .1 .3],3);
% 
% Nigel Chou, Feb 2009
%=====================================================================

% Scale pixel dimensions so that smallest dimension is 1 
pixdim = pixdim./min(pixdim);

if length(pixdim)==2
    xl=pixdim(1); yl=xl;
    zl=pixdim(2);
elseif length(pixdim)==3
     xl=pixdim(1); 
     yl=pixdim(2);
    zl=pixdim(3);
elseif length(pixdim)<2 || length(pixdim)>3
    disp('please enter a 3-vector for pixel dimensions')
end

R2 = zeros( 2*r+1, 2*r+1, 2*r+1 );

for x=1:2*r+1
    for y=1:2*r+1
        for z=1:2*r+1
            R2(x,y,z) =  (xl.*(r+1-x)).^2 + (yl.*(r+1-y)).^2 + (zl.*(r+1-z)).^2 ;
        end
    end
end

h = (1-R2./(2.*sigma^2)).*exp(-R2./(2.*sigma^2));
%h(h<eps*max(h(:))) = 0;

sumh = sum(h(:));
if sumh ~= 0,
    h  = h/sumh;
end;

if h(r+1,r+1,r+1)<0
    h=-h;
end