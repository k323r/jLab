function[v,psi,zeta,psif]=simpleddy(varargin)
%SIMPLEDDY  Streamfunction, velocity, and vorticity for various eddy profiles. 
%
%   [V,PSI,ZETA]=SIMPLEDDY(R,RO,VO) returns the streamfunction PSI,
%   azimuthal velocity V, and relative vorticity ZETA at radial locations R
%   for the Rankine vortex.
%
%   The eddy core radius is RO, in km, and at this radius the eddy 
%   currents achieve their maximum value of |VO| cm/s.  The sign of VO
%   gives the rotation sense of the eddy, positive for counterclockwise.
%
%   PSI has units of km^2/s, V has units of cm/s, and ZETA of 1/s. 
%
%   RO and VO should either be scalars or an arrays of the same size as R. 
%   The output arrays will have size LENGTH(R) x N.
%   _______________________________________________________________________
%
%   Eddy models
%
%   SIMPLEDDY can also return profiles associated with various other types
%   of model eddies.  The following forms are available:
%
%         SIMPLEDDY(R,RO,VO,'rankine')        -- Rankine vortex   [default]
%         SIMPLEDDY(R,RO,VO,'gaussian')       -- Gaussian vortex
%
%   RO and VO may then each be either scalars or arrays of the same length. 
%   _______________________________________________________________________
% 
%   Eddy 'slices'
%
%   SIMPLEDDY can also be used to return an eddy 'slice', that is, the
%   eddy streamfunction, currents, and vorticity observed along a curve.
%  
%   CV=SIMPLEDDY(CX,RO,VO,...) where CX=X+iY is an array of locations on 
%   the complex plane, with units of kilometers, returns the complex-valued
%   currents CV=U+iV due to an eddy located at the origin.
%
%   [CV,PSI,ZETA]=SIMPLEDDY(CX,RO,VO,...) also returns the streamfunction
%   and vorticity observed along this curve.  
%   _______________________________________________________________________
%
%   'simpleddy --t' runs a test.
%   'simpleddy --f' makes some sample figures.
%
%   Usage: [v,psi,zeta]=simpleddy(r,RO,VO);
%          [v,psi,zeta]=simpleddy(r,RO,VO,'gaussian');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015--2020 J.M. Lilly --- type 'help jlab_license' for details
 

%   See also EDDYCURVE, EDDYFIT, EDDYGUESS.



%   XXX what if X has more than one column.  


if strcmpi(varargin{1}, '--f')
    rankineddy_fig
    gaussian_fig
    return
end

if strcmpi(varargin{1}, '--t')
    simpleddy_test;
    return
end

x=varargin{1};
if size(x,1)==1
    disp('Transposing the first argument to SIMPLEDDY to make it a column vector.')
    x=conj(x');
end

ro=varargin{2};
vo=varargin{3};

eddystr='gaussian';
alpha=[];
        
if nargin>3
    if ~ischar(varargin{4})
        alpha=varargin{4};
        if alpha==0
             eddystr='rankine';
             alpha=1;
        elseif alpha==inf
            eddystr='gaussian';
        else
            error(['Sorry, the value ALPHA=' num2str(alpha) ' is not supported.'])
        end
    else
        if ischar(varargin{end})
            eddystr=varargin{end};
            varargin=varargin(1:end-1);
        elseif ischar(varargin{end-1})
            eddystr=varargin{end-1};
            alpha=varargin{end};
            varargin=varargin(1:end-2);
        end
    end
end

vo=vo/100/1000;
if strcmpi(eddystr(1:3),'gau')
    [v,psi,zeta]=gaussianeddy(ro,vo,abs(x));
elseif strcmpi(eddystr(1:3),'ran')
    [v,psi,zeta]=rankineddy(alpha,ro,vo,abs(x));
else
    error(['Sorry, the eddy model ''' eddystr ''' is not supported.'])
end


if size(x,2)==1
    x=vrep(x,size(v,2),2);
end
if ~isreal(x)||anyany(x<0)
    v=sqrt(-1)*frac(x,abs(x)).*v;
end

v(~isfinite(v))=0;
psi(~isfinite(psi))=0;
zeta(~isfinite(zeta))=0;



function[v,psi,zeta]=rankineddy(alpha,ro,vo,r)
%RANKINEDDY Velocity and streamfunction for a generalized Rankine vortex.
%
%   RANKINEDDY computes the velocity and streamfunction observed at
%   the origin due to a Rankine vortex located at a specified point,
%   as described in Lilly and Rhines (2002).
%
%   [v,psi]=rankineddy(eta,ro,vo);
%
%     Input
%	eta: complex-valued position of eddy (km)
%	ro:  size of eddy (km)
%	vo:  velocity at edge, negative for anticyclone (cm/s) 
%
%     Output
%	v,psi:	azimuthal velocity and streamfunction at origin
%
%   'rankineddy --f' generates a sample figure
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2013 J.M. Lilly --- type 'help jlab_license' for details    

if isempty(alpha)||~isfinite(alpha)
    alpha=1;
end

iO=find(r>ro);
v=vo.*r./ro;
vO=vo.*(ro.^alpha)./(r.^alpha);
if ~isempty(iO)
 %size(r),size(v), size(ro),size(iO), size(vo)
    v(iO)=vO(iO);
end
v=v*100*1000;

psi=(1/2).*(r./ro).^2;
if alpha==1
    psiO=log(r./ro)+1/2;
else
    psiO=1/2+frac(1,alpha-1)-frac(1,alpha-1).*frac(ro.^(alpha-1),r.^(alpha-1));
end

if ~isempty(iO)
    psi(iO)=psiO(iO);
end
psi=psi.*vo.*ro;

zeta=2*vo./ro+0*r;
zetaO=frac(vo.*ro,r.^(alpha+1)).*(1-alpha);
if ~isempty(iO)
    zeta(iO)=zetaO(iO);
end

function[v,psi,zeta]=gaussianeddy(ro,vo,r)
%GAUSSIANEDDY Velocity and streamfunction for a Gaussian vortex.
%
%   GAUSSIANEDDY computes the velocity and streamfunction observed at
%   the origin due to an eddy having a Gaussian streamfunction profile
%   located at a specified point.
%
%   [v,psi]=gaussianeddy(eta,ro,vo);
%
%     Input
%	eta: complex-valued position of eddy (km)
%	ro:  size of eddy (km)
%	vo:  velocity at edge, negative for anticyclone (cm/s) 
%
%     Output
%	v,psi:	azimuthal velocity and streamfunction at origin
%
%   'gaussianeddy --f' generates a sample figure
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2012 J.M. Lilly --- type 'help jlab_license' for details    

psi=-1.*vo.*ro.*exp(1/2.*(1-r.^2./ro.^2));
v=100*1000*r.*(vo./ro).*exp(1/2.*(1-r.^2./ro.^2));
zeta=(vo./ro).*(2-r.^2./ro.^2).*exp(1/2.*(1-r.^2./ro.^2));

function[]=simpleddy_test

x=[-100:.1:100]';
y=[0:2:16];
clear cx
for i=1:length(y)
    cx(:,i)=x+sqrt(-1)*y(i);
end

Ro=10;
Vo=15;

bool(1)=aresame(simpleddy(cx,Ro,Vo,'gaussian'),simpleddy(cx,Ro,Vo,inf));
bool(2)=aresame(simpleddy(cx,Ro,Vo,'rankine'),simpleddy(cx,Ro,Vo,0));
reporttest('SIMPLEDDY alternate input format matches string format',all(bool))

bool=false(1,2);
try
    simpleddy(cx,Ro,Vo,'matern',.75);
catch
    bool(1)=true;
end
try
    simpleddy(cx,Ro,Vo,.75);
catch
    bool(2)=true;
end
reporttest('SIMPLEDDY unsupported Matern eddy gives an error',all(bool))

function[]=eddyslice_figures

x=[-100:.1:100]';
y=-[0:2:16];
clear cx
for i=1:length(y)
    cx(:,i)=x+sqrt(-1)*y(i);
end

Ro=10;
Vo=15;


[cv,psi,zeta]=simpleddy(cx,Ro,Vo,'gaussian');
%[cv,psi,zeta]=simpleddy(cx,Ro,Vo,inf);
%[cv,psi,zeta]=simpleddy(cx,Ro,Vo);

figure,
subplot(2,2,1),plot(cx),hold on,axis([-1 1 -1 1]*Ro*2),axis equal
h=plot(Ro*rot(0:.01:2*pi+.1));
linestyle -h h 2G
subplot(2,2,2),hodograph(cv)
subplot(2,2,3),plot(x,real(cv)),vlines([-Ro Ro],'k:');hlines(0,'k:');
title('Along-Stream Currents')
subplot(2,2,4),plot(x,imag(cv)),vlines([-Ro Ro],'k:');hlines(0,'k:');
title('Cross-Stream Currents')


function[]=gaussian_fig
ro=10;
vo=-10;
eta=-sqrt(-1)*5+(-50:.1:50)';

%[v,psi]=simpleddy(eta,ro,vo,inf);
[v,psi]=simpleddy(eta,ro,vo,'gaussian');

figure,
subplot(121),
uvplot(real(eta),v),
title('Anticyclone sliced north of center, halfway to edge')
subplot(122)
polar(angle(v),abs(v))
title('Hodograph of Gaussian eddy')

function[]=rankineddy_fig
ro=10;
vo=-10;
eta=-sqrt(-1)*5+(-50:.1:50)';
[v,psi]=simpleddy(eta,ro,vo,'rankine',2);

figure,
subplot(121),
uvplot(real(eta),v),
title('Anticyclone sliced north of center, halfway to edge')
subplot(122)
polar(angle(v),abs(v))
title('Hodograph of Rankine eddy')
