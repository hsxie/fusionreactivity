% Hua-sheng XIE, huashengxie@gmail.com, 2023-02-12 07:37
% calculate the <sigma*v> using Monte-Carlo, for drift tri-Maxwellian
% Nath13 model
% 09:43 test OK, imethod=1,3 yield same results
% however, to improve the imethod=2 of the weight
% 23-02-17 07:53 update
function [sgmv,stdsgmv]=fsgmvmcdtm(vt1x,vt1y,vt1z,vd1x,vd1y,vd1z,...
    vt2x,vt2y,vt2z,vd2x,vd2y,vd2z,N,imethod)

% constants
% kB=1.3807e-23; % J/K
qe=1.6022e-19; % C
% me=9.1094e-31; % kg
mp=1.6726e-27; % kg
% epsilon0=8.8542e-12; % F/m
% mu0=4e-7*pi; % H/m
% c=2.99792458e8; % m/s

md=2*mp;
mt=3*mp;
% mhe=3*mp;
% mb=11*mp;

m1=md; m2=mt;
mr=m1*m2/(m1+m2);

nrpt=3; % repeat to obtain standard error
sgmvv=zeros(nrpt,1);
for jrpt=1:nrpt % 23-02-17 08:11
if(imethod==1 || imethod==2)
    v1x=vt1x*randn(N,1)+vd1x; % randn -> exp(-x^2/2)
    v1y=vt1y*randn(N,1)+vd1y;
    v1z=vt1z*randn(N,1)+vd1z;
    v2x=vt2x*randn(N,1)+vd2x;
    v2y=vt2y*randn(N,1)+vd2y;
    v2z=vt2z*randn(N,1)+vd2z;
    wgt=1+0.*v1x;
elseif(imethod==3)
    
    vt=sqrt(vt1x^2/6+vt1y^2/6+vt1z^2/3+vt2x^2/6+vt2y^2/6+vt2z^2/3+...
        vd1x^2+vd1y^2+vd1z^2+vd2x^2+vd2y^2+vd2z^2); % choose a vt
    v1x=vt*randn(N,1); % randn -> exp(-x^2/2)
    v1y=vt*randn(N,1);
    v1z=vt*randn(N,1);
    v2x=vt*randn(N,1);
    v2y=vt*randn(N,1);
    v2z=vt*randn(N,1);
    
    % weight of Monte-Carlo integral
    wgt=vt^6/(vt1x*vt1y*vt1z*vt2x*vt2y*vt2z...
        )*exp(-(v1x-vd1x).^2/(2*vt1x^2)+v1x.^2/(2*vt^2)...
        -(v1y-vd1y).^2/(2*vt1y^2)+v1y.^2/(2*vt^2)...
        -(v1z-vd1z).^2/(2*vt1z^2)+v1z.^2/(2*vt^2)...
        -(v2x-vd2x).^2/(2*vt2x^2)+v2x.^2/(2*vt^2)...
        -(v2y-vd2y).^2/(2*vt2y^2)+v2y.^2/(2*vt^2)...
        -(v2z-vd2z).^2/(2*vt2z^2)+v2z.^2/(2*vt^2));
end

if(imethod==1 || imethod==3)
    v=sqrt((v1x-v2x).^2+(v1y-v2y).^2+(v1z-v2z).^2);
elseif(imethod==2)
    v=zeros(N*N,1);
    for i=1:N
        for j=1:N
            v((i-1)*N+j)=sqrt((v1x(i)-v2x(j)).^2+(v1y(i)-v2y(j)).^2+(v1z(i)-v2z(j)).^2);
        end
    end
    wgt=1+0.*v;
end

E=0.5*mr*v.^2;
EkeV=E/(qe*1e3);

icase=1;
if(icase==1) % DT
    sgm=fsgmdt(EkeV);
end

sgmv0=mean(sgm.*v.*wgt,'omitnan');

sgmvv(jrpt)=sgmv0;

end
sgmv=mean(sgmvv);
stdsgmv=std(sgmvv);

end
