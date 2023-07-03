% Hua-sheng XIE, huashengxie@gmail.com, 2023-02-12 21:00
% calculate the <sigma*v> using Monte-Carlo, for slowing down
% 23-02-19 10:33 update
% method=3 vs method=1, difference ~2%, to check the reason
function [sgmv,stdsgmv]=fsgmvmcsd(vb1,vc1,vb2,vc2,N,imethod)

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
if(imethod==1|| imethod==2) % 23-02-16 21:39
    
    v1=rand(N,1);
    u1=vc1*(exp(v1*log(1+vb1^3/vc1^3))-1).^(1/3);
    
%     ffu1=@(u)3/log(1+vb1^3/vc1^3)*u.^2./(u.^3+vc1^3).*heaviside(vb1-u);
%     maxfv1=ffu1(min(2^(1/3)*vc1,0.999*vb1));
%     u1=vb1*rand(N,1);
%     fu1=maxfv1*rand(N,1);
%     u1=u1(fu1<=ffu1(u1));
%     N1=size(u1,1);
    
    v2=rand(N,1);
    u2=vc2*(exp(v2*log(1+vb2^3/vc2^3))-1).^(1/3);

    phi1=rand(N,1)*2*pi;
    theta1=rand(N,1)*pi;
    v1x=u1.*sin(theta1).*cos(phi1);
    v1y=u1.*sin(theta1).*sin(phi1);
    v1z=u1.*cos(theta1);
    
    phi2=rand(N,1)*2*pi;
    theta2=rand(N,1)*pi;
    v2x=u2.*sin(theta2).*cos(phi2);
    v2y=u2.*sin(theta2).*sin(phi2);
    v2z=u2.*cos(theta2);
    
    wgt=1+0.*v1x;
elseif(imethod==3)
    
    vt=1/2*sqrt(vc1^2/3+vb1^2/3+vc2^2/3+vb2^2/3); % choose a vt
    v1x=vt*randn(N,1); % randn -> exp(-x^2/2)
    v1y=vt*randn(N,1);
    v1z=vt*randn(N,1);
    v2x=vt*randn(N,1);
    v2y=vt*randn(N,1);
    v2z=vt*randn(N,1);
    
    % weight of Monte-Carlo integral
    As1=3/(4*pi*log(1+vb1^3/vc1^3));
    As2=3/(4*pi*log(1+vb2^3/vc2^3));
    wgt=((2*pi)^3*vt^6)*(As1*As2...
        ).*heaviside(vb1-sqrt(v1x.^2+v1y.^2+v1z.^2))./(...
        (sqrt(v1x.^2+v1y.^2+v1z.^2)).^3+vc1^3).*heaviside(...
        vb2-sqrt(v2x.^2+v2y.^2+v2z.^2))./(...
        (sqrt(v2x.^2+v2y.^2+v2z.^2)).^3+vc2^3).*exp(...
        (v1x.^2+v1y.^2+v1z.^2+v2x.^2+v2y.^2+v2z.^2)/(2*vt^2));
elseif(imethod==4) % use uniform rand number, 23-02-19 16:28
    
    v1x=vb1*(2*rand(N,1)-1);
    v1y=vb1*(2*rand(N,1)-1);
    v1z=vb1*(2*rand(N,1)-1);
    v2x=vb2*(2*rand(N,1)-1);
    v2y=vb2*(2*rand(N,1)-1);
    v2z=vb2*(2*rand(N,1)-1);
    
    % weight of Monte-Carlo integral
    As1=3/(4*pi*log(1+vb1^3/vc1^3));
    As2=3/(4*pi*log(1+vb2^3/vc2^3));
    wgt=(8*vb1^3*8*vb2^3)*(As1*As2...
        ).*heaviside(vb1-sqrt(v1x.^2+v1y.^2+v1z.^2))./(...
        (sqrt(v1x.^2+v1y.^2+v1z.^2)).^3+vc1^3).*heaviside(...
        vb2-sqrt(v2x.^2+v2y.^2+v2z.^2))./(...
        (sqrt(v2x.^2+v2y.^2+v2z.^2)).^3+vc2^3);
end

if(imethod==1 || imethod==3 || imethod==4)
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
