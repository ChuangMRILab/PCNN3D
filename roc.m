function out = roc(in,series, manmask,type,plotstr,ss,es)

%=====================================================================
% Plots reciever-operator charactersitic (ROC) curve
%
% usage : out = roc(in,series, manmask,type,plotstr,ss,es);
%
%         out     : ROC plot
%         in      : 3D input image
%         series  : series of binary masks
%         manmask : manual/reference mask
%         type    : use 'abs' or 'corr' FPR (default = 'corr')
%         plotstr : string to define line e.g. ('k-' for black line)
%         ss      : starting slice
%         es      : ending slice (subtracting from end), (default = 0)
%
% e.g. ROC = roc(I,I_PCNNborders, I_manmask,'corr','kx',3);
% 
% Nigel Chou, March 2010
%=====================================================================

if ~exist('ss','var')
    ss=1;
end

if ~exist('es','var')
    es=0;
end

if ~exist('type','var')
    type='corr';
end

if strcmp(type,'abs')
    for k=1:length(series)
        temp=jaccard(series{1,k},manmask,in,ss,es);
        out(1,k)=temp(2); % true positive
        out(2,k)=temp(3); % false negative
    end
elseif strcmp(type,'corr')
    for k=1:length(series)
        temp=jaccard(series{1,k},manmask,in,ss,es);
        out(1,k)=temp(2); % true positive
        out(2,k)=temp(4); % false negative corrected
    end
end

%---- plot ROC curve ---------
plot(out(2,:),out(1,:),plotstr,'LineWidth',2)
xlabel(['FPR (' type ')']);
ylabel('TPR');


