function drawborderslice(in,bw,col,lw)

% =========================================================================
%
% Displays border of a binary mask on a 2D image.
%
% usage: Out = drawborderslice(in,bw,col,lw);
% 
%   in     : original image
%   bw     : binary mask
%   col    : color in RGB e.g. [1 0 0] for red
%   lw     : linewidth of border
%
% e.g. drawborderslice(I_slice,I_slice_mask,[1 0 0],2);
%
% Nigel Chou. June 2009
%==========================================================================

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
imshow(in,[0 max(in(:))/2],'InitialMagnification','fit');

hold on
drawborder(bw,col,lw);
hold off
