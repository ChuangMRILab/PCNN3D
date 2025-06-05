function drawborder(bw,col,lw)

% =========================================================================
%
% Draws border of a binary mask. Can be drawn onto a 2D image which is 
% already displayed in window using command 'hold on'.
%
% usage: drawborder(bw,col,lw);
% 
%   bw     : binary mask
%   col    : color in RGB e.g. [1 0 0] for red
%   lw     : linewidth of border
%
% e.g. drawborder(I_slice_mask,[1 0 0],2);
%
% Nigel Chou. June 2009
%==========================================================================

%hold on;

if issparse(bw)
    bw=full(bw);
end

b = bwboundaries(bw);
for k=1:numel(b);
    plot(b{k}(:,2),b{k}(:,1),'color',col,'Linewidth',lw);
end;

%hold off;