function Out = border_comp(OrigVol,lw,varargin)

% =========================================================================
%
% Displays borders of PCNN mask on individual coronal slices as a line (color
% can be specified). Multiple borders can be drawn on the same image. Displays
% each coronal slice as seperate figure.
%
% usage: Out = drawbordervol_comp(OrigVol,lw,bwvol,col,iterno,bwvol2,col2,iterno2,...)
% 
%   OrigVol: original image volume
%   lw     : linewidth of border
%   bwvol  : series of binary mask volumes (as matlab cell format)
%   col    : color in RGB e.g. [1 0 0] for red
%   iterno : iteration no. for PCNN (if only one mask, then set iterno=0)
%
% e.g.  drawbordervol_comp(I,2,I_border,[1 0 0],54);
%
% Nigel Chou, Oct 2009
%==========================================================================

scrsz = get(0,'ScreenSize');

vararg_l = length(varargin);

if mod(vararg_l,3)~=0
    disp('bw image arguments must be in sets of 3');
end

for k=1:size(OrigVol,3)
    if max(max(OrigVol(:,:,k)))~=0
        figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
        imshow(OrigVol(:,:,k),[min(min(OrigVol(:,:,k))) max(max(OrigVol(:,:,k)))],'InitialMagnification','fit');
        hold on;
        for n = 1:vararg_l/3
            bwVol = varargin{3*n-2};
            col = varargin{3*n-1};
            iterno = varargin{3*n};
            if iterno==0
                drawborder(bwVol(:,:,k),col,lw);
            else
                drawborder(bwVol{1,iterno}{1,k},col,lw);
                %drawborder(bwVol{1,iterno}(:,:,k),col,lw);
            end
        end
        %Out(k)=getframe;
        hold off;
    end
end
