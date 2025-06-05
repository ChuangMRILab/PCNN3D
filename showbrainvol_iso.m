function showbrainvol_iso(data,viewangle)

%==========================================================================
%
% Displays 3D surface of a brain mask with isotropic voxel resolution
%
% usage: showbrainvol_iso(data,viewangle)
%   
%   data      : binary input volume 
%   viewangle : arguments for view command, default [74,16]
%
% e.g. showbrainvol(Out{1,47},[74,16])
%
% Nigel Chou, June  2010
%==========================================================================

if iscell(data)
    data_z = zeros([ size(data{1,1}) length(data)]);
    for k=1:length(data)
        data_z(:,:,k) = data{1,k};
    end
    data = data_z;
end

if ~exist('viewangle'),
    viewangle = [74,16];
end

data = flipdim(data,1);

figure;
data = smooth3(data,'box',3);

% ------- Choose colour -----------
%p1 = patch(isosurface(data,.5), 'FaceColor',[253/256 224/256 187/256],'EdgeColor','none');
p1 = patch(isosurface(data,.5), 'FaceColor',[200/256 170/256 160/256],'EdgeColor','none');
p2 = patch(isocaps(data,.5), 'FaceColor','interp','EdgeColor','none');
isonormals(data,p1);
view(0,45); axis equal vis3d tight
camlight; lighting phong

% ------- Choose final viewing angle ----------
view(viewangle(1),viewangle(2));


%>> showbrainvol2(CT1_28dec09_2_3D_manmask_nc);xlim([0 140]);ylim([5 88]);zlim([0 207]);
