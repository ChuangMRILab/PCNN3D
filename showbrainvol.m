function data_scaled = showbrainvol(data,vox_xy,vox_z,viewangle)

%==========================================================================
% 
% Displays 3D surface of a brain mask with anisotropic voxel resolution
% 
% usage: showbrainvol(data,vox_xy,vox_z,viewangle);
%
%   data      : binary input volume 
%   vox_xy    : pixel dimensions in x and y (assumed to be
%               the same)
%   vox_z     : pixel dimension in z
%   viewangle : arguments for view command, default [74,16]
%
% e.g. showbrainvol(Out{1,47},0.09765625,0.30000925);
%
% Nigel Chou; June 23, 2009
%==========================================================================

if ~exist('viewangle'),
    viewangle = [74,16];
end

[X,Y,Z] = size(data);
data_scaled = zeros(X,Y,ceil(Z*vox_z/vox_xy));

for k=1:X
    for l=1:Y
        temp=data(k,l,:);
        data_vec = temp(:);
        temp2=interp1( 0:vox_z:(Z*vox_z-0.5*vox_z), data_vec, 0:vox_xy:Z*vox_z-.5*vox_xy,'nearest');
        data_scaled(k,l,:)=temp2;
    end
end

figure;
data_scaled = smooth3(data_scaled,'box',5);
% ------- Choose colour -----------
p1 = patch(isosurface(data_scaled,.5), 'FaceColor',[253/256 224/256 187/256],'EdgeColor','none');
p2 = patch(isocaps(data_scaled,.5), 'FaceColor','interp','EdgeColor','none');
isonormals(data_scaled,p1);
view(0,45); axis equal vis3d tight
camlight; lighting phong

% ------- Choose final viewing angle ----------
view(viewangle(1),viewangle(2));
