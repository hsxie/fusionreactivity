% Hua-sheng XIE, huashengxie@gmail.com, 2023-02-12 07:37
% calculate the <sigma*v> using Monte-Carlo, for drift ring Maxwellian
% 13:43 test OK, imethod=1,3 yield same results
% however, to improve the imethod=3 of the weight
% 23-02-13 19:53 update imethod=1, to avoid interp1 error
% 23-02-19 10:33 update
function [sgmv,stdsgmv]=fsgmvmcdrm(vt1x,vt1z,vd1x,vd1y,vd1z,vd1r,...
    vt2x,vt2z,vd2x,vd2y,vd2z,vd2r,N,imethod)

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
    
    As1=exp(-vd1r^2/(2*vt1x^2))+sqrt(pi/2)*vd1r/vt1x*erfc(-vd1r/(sqrt(2)*vt1x));
    As2=exp(-vd2r^2/(2*vt2x^2))+sqrt(pi/2)*vd2r/vt2x*erfc(-vd2r/(sqrt(2)*vt2x));

    ua1=rand(N,1);
    ua2=rand(N,1);
    if(vd1r==0)
        vpa1=sqrt(2)*vt1x*sqrt(-log(1-ua1));
    else % to improve
        xx=(0:0.001:5)*(vt1x+vd1r);
        yy=1/As1*(exp(-vd1r^2/(2*vt1x^2))+...
            sqrt(pi)*vd1r/(sqrt(2)*vt1x)*erf(vd1r/(sqrt(2)*vt1x)))-...
            1/As1.*(exp(-(xx-vd1r).^2./(2*vt1x^2))+...
            sqrt(pi)*vd1r/(sqrt(2)*vt1x)*erf((vd1r-xx)./(sqrt(2)*vt1x)));
        ind1=find(yy<1e-10); % 23-02-13 19:29 to avoid interp1 error
        ind2=find((1-yy)<1e-10);
        yy=yy(ind1(end):ind2(1));
        xx=xx(ind1(end):ind2(1));
        vpa1=interp1(yy,xx,ua1); % use interp to calculate the inverse function
    end
    if(vd2r==0)
        vpa2=sqrt(2)*vt2x*sqrt(-log(1-ua2));
    else % to improve
        xx=(0:0.001:5)*(vt2x+vd2r);
        yy=1/As2*(exp(-vd2r^2/(2*vt2x^2))+...
            sqrt(pi)*vd2r/(sqrt(2)*vt2x)*erf(vd2r/(sqrt(2)*vt2x)))-...
            1/As2.*(exp(-(xx-vd2r).^2./(2*vt2x^2))+...
            sqrt(pi)*vd2r/(sqrt(2)*vt2x)*erf((vd2r-xx)./(sqrt(2)*vt2x)));
        ind1=find(yy<1e-10); % 23-02-13 19:29 to avoid interp1 error
        ind2=find((1-yy)<1e-10);
        yy=yy(ind1(end):ind2(1));
        xx=xx(ind1(end):ind2(1));
        vpa2=interp1(yy,xx,ua2); % use interp to calculate the inverse function
    end

    phia1=rand(N,1)*2*pi;
    v1x=vpa1.*sin(phia1)+vd1x;
    v1y=vpa1.*cos(phia1)+vd1y;
    phia2=rand(N,1)*2*pi;
    v2x=vpa2.*sin(phia2)+vd2x;
    v2y=vpa2.*cos(phia2)+vd2y;
    
    % Gaussian rand para velocity with drift
    v1z=vt1z*randn(N,1)+vd1z; % randn -> exp(-x^2/2)
    v2z=vt2z*randn(N,1)+vd2z;
    
    wgt=1+0.*v1x;
elseif(imethod==3)
    
    vt=sqrt(vt1x^2/3+vt1z^2/6+vt2x^2/3+vt2z^2/6+vd1r^2+vd2r^2+...
        vd1x^2+vd1y^2+vd1z^2+vd2x^2+vd2y^2+vd2z^2); % choose a vt
    v1x=vt*randn(N,1); % randn -> exp(-x^2/2)
    v1y=vt*randn(N,1);
    v1z=vt*randn(N,1);
    v2x=vt*randn(N,1);
    v2y=vt*randn(N,1);
    v2z=vt*randn(N,1);
    
    % weight of Monte-Carlo integral
    As1=exp(-vd1r^2/(2*vt1x^2))+sqrt(pi/2)*vd1r/vt1x*erfc(-vd1r/(sqrt(2)*vt1x));
    As2=exp(-vd2r^2/(2*vt2x^2))+sqrt(pi/2)*vd2r/vt2x*erfc(-vd2r/(sqrt(2)*vt2x));
    wgt=vt^6/(vt1x^2*vt1z*vt2x^2*vt2z*As1*As2...
        )*exp(-(sqrt((v1x-vd1x).^2+(v1y-vd1y).^2)-vd1r).^2/(2*vt1x^2)+...
        (v1x.^2+v1y.^2)/(2*vt^2)-(v1z-vd1z).^2/(2*vt1z^2)+v1z.^2/(2*vt^2)...
        -(sqrt((v2x-vd2x).^2+(v2y-vd2y).^2)-vd2r).^2/(2*vt2x^2)+...
        (v2x.^2+v2y.^2)/(2*vt^2)-(v2z-vd2z).^2/(2*vt2z^2)+v2z.^2/(2*vt^2));
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
