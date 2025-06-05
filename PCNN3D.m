function [Out , G, optG]= PCNN3D(In,p,voxdim,BrSize,maxiter,OptMask)

%==========================================================================
% PCNN3D 
%
% Performs skull-stripping  of MRI volumes using 3D PCNN on the entire 
% volume
%
% usage: [Out , G]= PCNN3D(In,p,voxdim,BrSize,maxiter,OptMask);
%
%   Out    : series of binary mask volumes (MATLAB cell format)
%   G      : plot of volume of mask (in mm^3) vs. iterations
%   In     : 3D input image
%   p      : radius of structural element
%   voxdim : vector containing the x, y, and z voxel dimensions
%   BrSize : (optional) vector of min and max of brain size (in mm^3). Default=[100 550] for mouse brain. This will be used to estimate optimal iteration.
%   maxiter: (optional) max iteration. Default=200.
%   OptMask: (optional) =1: output only the mask of the optimal guess. Default=0.
%
% e.g. [I_border, G_I]= PCNN3D(I,4,[.1 .1 .25],[100 550]);
%
% Nigel Chou, last updated 1 July 2010
% Kai-Hsiang Chuang, 20 Apr 2012, Ver 1.1
%                    16 Mar 2016, Ver 1.2 bug fixed for rat brain. Add OptMask & optG.
%==========================================================================

%close all;
if ~exist('maxiter','var')
    maxiter = 200;
end
if ~exist('BrSize','var')
	BrSize=[100 550];
end
if length(BrSize) ~= 2
	error('BrSize should be in the form of [min max]')
end
if ~exist('OptMask','var')
	OptMask=0;
end

Out = cell(1,maxiter);

% ----------------------------------
%  Set parameters
% ----------------------------------
beta = 0.2;
tauF=0.3;
tauL=1; %.3
tauT=10; %10, 30
VF=0.01;
VL=.2;
VT=20;
r=3;

% ----------------------------------
%  Initialize
% ----------------------------------
minI=min(In(:));
maxI=max(In(:));
In = (In-minI)./(maxI-minI);

[N1,N2,N3] = size(In);
tot_vol = N1*N2*N3;

pad = 2;
In = cat(3, zeros(N1,N2,pad), In, zeros(N1,N2,pad));
N3 = N3 + 2*pad;

F = zeros(N1,N2,N3);
L = F; U = F; Y = F; 
A = F; Abreak = A;
T = 1*ones(N1,N2,N3);


% ----------------------------------
%  Set-up structural elements
% ----------------------------------

% Set Radius of structural element
%p = 3;

% ==== Cube =====
strel_box = ones(2*p+1,2*p+1,2*p+1);


% ==== Diamond =====
seed=cat(3, [0 0 0; 0 1 0; 0 0 0],...
    [0 1 0; 1 1 1; 0 1 0],...
    [0 0 0; 0 1 0; 0 0 0]);
strel_diam=seed;
if p>1
    for k=1:p-1
        strel_diam=imdilate(strel_diam,seed,'full');
    end
    strel_diam(isinf(strel_diam))=0;
end

% ==== Sphere ======
strel_sph = zeros(2*p+1,2*p+1,2*p+1);
for x=1:2*p+1
    for y=1:2*p+1
        for z=1:2*p+1
            if sqrt( (x-(p+1)).^2 + (y-(p+1)).^2 + (z-(p+1)).^2 ) <= p
                strel_sph(x,y,z)=1;
            end
        end
    end
end

%voxdim = [0.09765625 0.09765625 0.30000925];
mindim = min(voxdim);
xl=voxdim(1)/mindim ; yl=voxdim(2)/mindim ; zl=voxdim(3)/mindim ;

% ==== Sphere Scaled ======
strel_sphsc = zeros(2*p+1,2*p+1,2*p+1);
for x=1:2*p+1
    for y=1:2*p+1
        for z=1:2*p+1
            if sqrt( (xl.*(x-(p+1))).^2 + (yl.*(y-(p+1))).^2 + (zl.*(z-(p+1))).^2 ) <= p
                strel_sphsc(x,y,z)=1;
            end
        end
    end
end


% ----------------------------------
%  pre-calculate exponential factors
% ----------------------------------

eF=exp(-log(2)./tauF);
eL=exp(-log(2)./tauL);
eT=exp(-log(2)./tauT);

% ----------------------------------
%  do niter iterations
% ----------------------------------
n = 0; G = zeros(1,maxiter); brainsize=0;
fraction = 0; tic;

%h = waitbar(0,'PCNN (3D)');

gauss = scaledgauss(r,voxdim,.5);

while (( n<maxiter ) && ( brainsize<BrSize(2) ) )
        
    %convMY = smooth3(Y,'gaussian',2*r+1 ,.6);
    convMY = convn(Y,gauss,'same');
    F = eF.*F + In + VF.*convMY;
    %convNY = smooth3(Y,'gaussian',2*r+1 ,.6);
    convNY = convn(Y,gauss,'same');
    L = eL.*L + VL.*convNY;
    U = F.*(1 + beta*L);
    T = eT.*T + VT.*Y;
    Y = U > T;
    
    A = A + Y;
    %A = double(A>0);
    %Y = double(Y);
    %G(n)=sum(sum(sum(Y)));
    
%     for k=p+1:N1-p
%         for l=p+1:N2-p
%             for m=p+1:N3-p
%                 if sum(A(k,l-p:l+p,m))>=(2*p+1) && sum(A(k-p:k+p,l,m))>=(2*p+1) && sum(A(k,l,m-p:m+p))>=(2*p+1)
%                     Abreak(k,l,m)=1;
%                 else
%                     Abreak(k,l,m)=0;
%                 end
%             end
%         end
%     end

     %A = imfill(A,6,'holes'); % fill before breaking
    
     type = 'strel_sphsc';
     
     
    % Break connections --------
    Abreak = eval(sprintf('imerode(A,%s)',type));
    
    [L,ncomp] = bwlabeln(Abreak,6);
    
    if ncomp==0
        Afinal = Abreak;
    end
    
    s = regionprops(L,'Area');
    area_values = [s.Area];
    maxcomp_ind = find(area_values==max(area_values));
    
    if length(maxcomp_ind)==1
        Afinal = ismember(L,maxcomp_ind);
    else
        Afinal = Abreak;
    end
    
    % Dilate Mask --------
    Afinal = eval(sprintf('imdilate(Afinal,%s)',type));
    
    % Fill Holes ---------
    Afinal = imfill(Afinal,6,'holes');
    
    %     else
    %         numpix=zeros(1,ncomp);
    %         for q=1:ncomp
    %             numpix(q) = sum(sum(sum(L==q)));
    %         end
    %         numpix
    %         biggestconarea = max(numpix)
    %
    %         maxcompind = find(numpix == max(numpix));
    %
    %         if length(maxcompind)==1
    %             Afinal = L==maxcompind;
    %         else
    %             Afinal=Abreak;
    %         end
    %     end
    
    for k=1:N3
        Afinal(:,:,k) = imfill(Afinal(:,:,k),'holes');
    end
    
    Afinal = Afinal(:,:,pad+1:end-pad);
    pixsize = sum(Afinal(:));
    voxvol = voxdim(1).*voxdim(2).*voxdim(3);
    brainsize = pixsize.*voxvol;
    %== Convert to sparse =====
    Out{1,n+1} = cell(1,size(Afinal,3));
    for k=1:size(Afinal,3)
        Out{1,n+1}{1,k} = sparse(Afinal(:,:,k));
    end
    
    
    fraction = pixsize/tot_vol;
    G(n+1)=brainsize;
    
    n = n + 1;
    
    fprintf('iteration %2d, number of regions:%6d  , fraction of tot vol is %3.2f \n',n, ncomp, fraction);
    
    %waitbar(n/maxiter,h,sprintf('PCNN (3D) iteration %d',n));
    
end

%===============================
% Estimate Optimal Iteration
%===============================

its = findit(G,BrSize);
optG=its(2);
disp('=============================================');
fprintf('\nGuess for best iteration is %d. \n', its(2) );
fprintf('Best iterations are in the range [%d %d]. \n \n', its(1), its(3) );
disp('=============================================');
    
timetaken = toc;

fprintf('\nTotal time taken: %d min %0.0f s \n\n', floor(timetaken/60), mod(timetaken,60) )

% A16b=zeros(size(A16));
% for k=1:45;
%     temp=imfill(A16(:,:,k),'holes');
%     A16b(:,:,k)
% [L,ncomp] = bwlabel(temp,4);
%
%     if ncomp==0
%         A16b(:,:,k)=temp;

%     else
%
%         numpix=zeros(1,ncomp);
%         for q=1:ncomp
%             numpix(q) = sum(sum(sum(L==q)));
%         end
%
%         maxcompind = find(numpix == max(numpix))
%         if length(maxcompind)==1
%         A16b(:,:,k) = L==maxcompind;
%         else
%         A16b(:,:,k)=temp;
%         end
%     end
% end
if OptMask == 1
    Out = Out(1,its(2));
else
    Out = Out(1,1:n-1);
end
%figure; plot(1:n-1,G(1:n-1),'k-');
%xlim([1 max(G)]);

