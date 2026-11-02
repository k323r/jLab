function[varargout]=isomax(varargin)
% ISOMAX  Returns those transform maxima that are isolated from others.
%  
%   This function is part of 'element analysis' described in Lilly (2017), 
%   "Element analysis: a wavelet-based method for analyzing time-localized
%   events in noisy time series", available at www.jmlilly.net.
%  
%   BOOL=ISOMAX(SIZ,INDEX,WW,FF,GAMMA,BETA,MU,LAMBDA) returns an boolean
%   array is true for transform maxima that are isolated from other maxima.
% 
%   SIZ is the original size of the wavelet transform array, while the 
%   index into maxima locations INDEX, transform values WW, and transform
%   frequencies FF, are all as computed by TRANSMAX. 
%
%   ISOMAX returns a boolean array BOOL that is true for those maxima that 
%   are isolated at level LAMBDA.  The maxima are interpreted as being from  
%   the transform *of* a (MU,GAMMA) wavelet *with* a (BETA,GAMMA) wavelet.
%
%   Specifically, LAMBDA defines a contour around a transform maxima at 
%   which the transform modulus is expected to have decayed to a fraction
%   LAMBDA of its peak value.  A maxima is 'isolated' if it is larger in 
%   magnitude that all other points encompassed that region of influence.
%
%   Note that if the third or fourth element of SIZ is greater than one, 
%   this indicates a set of original time series organized as a 2D or 3D 
%   array respectively, according to the WAVETRANS convention.  In this 
%   case the isolation criterion is applied only to events from the same 
%   time series; the extra dimensions are simply looped over.  
%
%   [BOOL,Z]=ISOMAX(...) also returns a matrix Z with LENGTH(INDEX) columns
%   containing the regions of influence T+1i*F surrounding all of the input
%   maxima.  One may then use 'plot(z)' to plot these regions.   
%
%   ISOMAX(...,'parallel') optionally uses a PARFOR loop when computing the
%   regions of influence, which can be useful for large datasets.
%
%   See also MAXPROPS, TRANSMAX, TRANSMAXDIST, MAX2EDDY.
%
%   Usage: bool=isomax(siz,index,ww,ff,gamma,beta,mu,lambda);
%          [bool,z]=isomax(siz,index,ww,ff,gamma,beta,mu,lambda);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2017 J.M. Lilly --- type 'help jlab_license' for details        


if strcmpi(varargin{1},'--f')
%   This doesn't work right now.  
%   'isomax --f' generates a sample figure.
%    isomax_fig;
    return
end
if strcmpi(varargin{1},'--t')
    isomax_test,return
end

parstr='serial';
str='band';
for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'par')||strcmpi(varargin{end}(1:3),'ser')
            parstr=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'ban')||strcmpi(varargin{end}(1:3),'ene')
            str=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end
        
siz=varargin{1};
index=varargin{2};
ww=varargin{3};
ff=varargin{4};
ga=varargin{5};
be=varargin{6};
mu=varargin{7};
lambda=varargin{8};

[ii,jj,nn,kk]=ind2sub(siz,index);
   
bool=false(size(ii));

if nargout == 2
    z=zeros(1000,length(ii));
end

if length(siz)==2
    siz(3)=1;
    siz(4)=1;
elseif length(siz)==3
    siz(4)=1;
end

[tz,fz]=morseregion(lambda,ga,be,mu,1,str);
%figure,plot(tz,fz)

for n=1:siz(3)
    %n
    booln=(nn==n);
    for k=1:siz(4)
        index=find(booln&(kk==k));
        %length(find(index))
        if ~isempty(index)
            if nargout ==1
                bool(index)=isomax1(ii(index),ww(index),ff(index),ga,be,mu,lambda,tz,fz,parstr,str);
            elseif nargout ==2
                [bool(index),z(:,index)]=isomax1(ii(index),ww(index),ff(index),ga,be,mu,lambda,tz,fz,parstr,str);
            end
        end
    end
end


varargout{1}=bool;
if nargout==2
    varargout{2}=z;
end

function[bool,z]=isomax1(ii,ww,ff,ga,be,mu,lambda,tz,fz,parstr,str)

originalposition=[1:length(ww)]';
%Sort by magnitude
if ~isempty(ww)
    [temp,sorter]=sort(-abs(ww));
    %sorter
    vindex(ww,ii,ff,originalposition,sorter,1);
end

bool=[];
[c,rho,frho]=maxprops(1,1,ga,be,mu,str);
fact=1./frho;
%fact=frac(morsefreq(ga,be),morsefreq(ga,mu)).*frac(mu+1,be).^(1./ga);
%aresame(fact,1./frho,1e-10)

if isempty(ga)||isempty(ii)
    z=[];
else
    z=vrep(tz+1i*fz,length(ff),2);
    fmat=vrep(ff(:)',length(tz),1);
    z=real(z)./fmat.*fact+1i*imag(z).*fmat./fact;
%     if strcmpi(parstr(1:3),'par')
%         parfor i=1:length(ff)
%             z(:,i)=tz./ff(i).*fact+sqrt(-1)*fz.*ff(i)./fact;
%         end
%     else
%         for i=1:length(ff)
%             z(:,i)=tz./ff(i).*fact+sqrt(-1)*fz.*ff(i)./fact;
%         end
%     end
%     
%     
    %z2=morseregion(fs(jj),lambda,ga,be,mu);
    %aresame(z,z2,1e-6)
    %plot(z,'k'),hold on,plot(z2,'r--')
    %z=real(z)+sqrt(-1)*imag(z); % figure,plot(z)
    %tmat=vrep(t(ii)',size(z,1),1);
    tmat=vrep(ii',size(z,1),1);
    %vsize(tmat,t,ii,z)
    z=z+tmat;
    bool=true(size(ii));

    if strcmpi(parstr(1:3),'par')
        parfor i=2:length(ww)
            ti=ii([1:i-1]);
            fi=ff([1:i-1]);
            bool(i)=~anyany(inpolygon(ti,fi,real(z(:,i)),imag(z(:,i))));
        end
    else
        for i=2:length(ww)
            ti=ii([1:i-1]);
            fi=ff([1:i-1]);
            bool(i)=~anyany(inpolygon(ti,fi,real(z(:,i)),imag(z(:,i))));
        end
    end

    bool(originalposition)=bool;
    z(:,originalposition)=z;
    %vindex(ii,jj,ff,ww,bool,1);
    %vindex(z,bool,2);
    %if ~isempty(ii)
    %    [ii,sorter]=sort(ii);
    %    vindex(jj,ff,ww,sorter,1);
    %    vindex(z,sorter,2);
    %end
end

%if isempty(ii)
%    disp('ISOMAX found no isolated maxima with the specified properties.')
%end

function[]=isomax_test

ga=3;be=2;mu=0;
fs=2*pi./logspace(1,3,20)';
siz=[1000 20];
ii=[300 305 800]';
jj=[10 10 10]';
ww=[1 0.5 0.5]';
ff=fs(jj);

[bool,z]=isomax(siz,sub2ind(siz,ii,jj),ww,ff,ga,be,mu,1/2);
b1=aresame(bool(:),[true false true]');
b2=inpolygon(ii(2),ff(2),real(z(:,1)),imag(z(:,1)));
reporttest('ISOMAX smaller maximum within region of influence is not isolated',b1&&b2);

siz=[1000 20 2];
bool=isomax(siz,sub2ind(siz,ii,jj,[1 2 1]'),ww,ff,ga,be,mu,1/2);
reporttest('ISOMAX maxima from different time series are treated separately',allall(bool));
