function its = findit(G,range)

%==========================================================================
%
% finds 'best' iteration as well as range with best iterations for PCNN
%
% usage: its=findit(G,range);
%
%   its    : vector of 3 values: itmin: beginning of plateau region
%                                itmid: middle of plateau
%                                       (most likely the optimal iteration)
%                                itmax: end of plateau
%   G      : plot of volume of mask (in mm^3) vs. iterations
%   range  : range of brain sizes expected e.g. [100 550] for C57BL6 mouse
%            for lower bound preferable choose a smaller number
%
% e.g. its = findit(I_PCNN_G,[100 550])
%
% Nigel Chou, last updated 1 July 2010
%==========================================================================

its =zeros(1,3);

range = G>range(1) & G<range(2);
if isempty(range)
    display('warning: no pleateau region found');
end
Gshort=G(range);
index=find(range);
Gshort = [Gshort G(index(end)+1)];

%---- find increase for each iteration ----
dGshort=Gshort(2:end)-Gshort(1:end-1);
%---- find 'rate' of increase -------
ddGshort=dGshort(2:end)-dGshort(1:end-1);

%---- max rate of increase ----
maxgradinc=find(ddGshort==max(ddGshort));
its(3) = index(maxgradinc+1);

%---- min rate of increase ----
mingradinc=find(ddGshort==min(ddGshort));
its(1) = index(mingradinc+1);

%---- middle of plateau ---- 
its(2) = its(1) + round((its(3)-its(1))/2);

