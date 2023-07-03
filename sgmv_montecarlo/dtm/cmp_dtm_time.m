% 2023-02-12 09:06 Hua-sheng XIE, huashengxie@gmail.com, ENN
% tst drift tri-Maxwellian sigmv*v
% 23-02-19 08:40 update
% 23-02-19 17:01 cpuime, error vs N

close all;clear;clc;

load('../../bench/benchdat.mat');
dat=dat.';
Tr1d=dat(1,:);
sgmv1d_all=dat(2:6,:);
pars=[0,1.02;
    0,2;
    20,2;
    0,0.5;
    20,0.5];
jE=5;

% constants
kB=1.3807e-23; % J/K
qe=1.6022e-19; % C
% me=9.1094e-31; % kg
mp=1.6726e-27; % kg
% epsilon0=8.8542e-12; % F/m
% mu0=4e-7*pi; % H/m
% c=2.99792458e8; % m/s

md=2*mp;
mt=3*mp;

m1=md; m2=mt;
mr=m1*m2/(m1+m2);

TT=5:5:100;

% Rt=0.3; % Tper/Tpar
% EdkeV=20; % keV

Rt=pars(jE,2); % Tper/Tpar
EdkeV=pars(jE,1); % keV

vd=sqrt(2*EdkeV*(qe*1e3)/mr);

sgmvv1=0.*TT; sgmvv2=0.*TT; sgmvv3=0.*TT;
stdsgmvv1=0.*TT; stdsgmvv2=0.*TT; stdsgmvv3=0.*TT;
runtime1=0; runtime2=0; runtime3=0;
for j=1:length(TT)
    TrkeV=TT(j); % keV
    
    % sgmvdt=fsgmv(TrkeV,1);
    
    Tr=TrkeV*qe*1e3/kB; % keV -> K
    Trper=3*Rt*Tr/(2*Rt+1);
    Trpar=3*Tr/(2*Rt+1);
    T1x=Trper; T1y=1.0*Trper; T1z=Trpar; T2x=Trper; T2y=1.0*Trper; T2z=Trpar;
    vt1x=sqrt(kB*T1x/m1); vt1y=sqrt(kB*T1y/m1);vt1z=sqrt(kB*T1z/m1);
    vt2x=sqrt(kB*T2x/m2); vt2y=sqrt(kB*T2y/m2); vt2z=sqrt(kB*T2z/m2);
    
    vd1x=0*vd; vd1y=0; vd1z=1*vd;
    vd2x=0; vd2y=0.0*vd; vd2z=0;
    
    N=1000; N2=10*N; N3=100*N;
    tmp1=cputime;
    [sgmv1,stdsgmv1]=fsgmvmcdtm(vt1x,vt1y,vt1z,vd1x,vd1y,vd1z,vt2x,vt2y,vt2z,vd2x,vd2y,vd2z,N,1);
    tmp1=cputime-tmp1; runtime1=runtime1+tmp1;
    tmp2=cputime;
    [sgmv2,stdsgmv2]=fsgmvmcdtm(vt1x,vt1y,vt1z,vd1x,vd1y,vd1z,vt2x,vt2y,vt2z,vd2x,vd2y,vd2z,N2,1);
    tmp2=cputime-tmp2; runtime2=runtime2+tmp2;
    tmp3=cputime;
    [sgmv3,stdsgmv3]=fsgmvmcdtm(vt1x,vt1y,vt1z,vd1x,vd1y,vd1z,vt2x,vt2y,vt2z,vd2x,vd2y,vd2z,N3,1);
    tmp3=cputime-tmp3; runtime3=runtime3+tmp3;
    
    sgmvv1(j)=sgmv1; sgmvv2(j)=sgmv2; sgmvv3(j)=sgmv3;
    stdsgmvv1(j)=stdsgmv1; stdsgmvv2(j)=stdsgmv2; stdsgmvv3(j)=stdsgmv3;
end

%%
close all;
if(1==0)
figure('unit','normalized','DefaultAxesFontSize',16,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,...
    'position',[0.01,0.05,0.45,0.55]);

% plot(TT,sgmvv1,'--',TT,sgmvv2,'x','linewidth',3); hold on;

plot(Tr1d,sgmv1d_all(jE,:),':','linewidth',3); hold on;
errorbar(TT,sgmvv1,stdsgmvv1,'-','linewidth',3); hold on;
errorbar(TT,sgmvv2,stdsgmvv2,'--','linewidth',3); hold on;
errorbar(TT,sgmvv3,stdsgmvv3,'x','linewidth',3); hold on;
xlabel('T_r [keV]');
ylabel('<\sigma{}v> (m^3/s)');
legend('1D analytical','method=1, using f_{1,2}, O(N)','method=2, using f_{1,2}, O(N^2)',...
    'method=3, using Gaussian with weight, O(N)',...
    'location','best');
legend('boxoff');

title(['R_t=',num2str(Rt),', E_d=',num2str(EdkeV),'keV',', N_1=',num2str(N),...
    ', N_2^2=',num2str(N2^2),', N_3=',num2str(N3),10,'v_{d1}/v_d=[',num2str(vd1x/vd),',',...
    num2str(vd1y/vd),',',num2str(vd1z/vd),'], v_{d2}/v_d=[',...
    num2str(vd2x/vd),',',num2str(vd2y/vd),',',num2str(vd2z/vd),']']);
else

figure('unit','normalized','DefaultAxesFontSize',13,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,...
    'position',[0.01,0.05,0.75,0.4]);

% plot(TT,sgmvv1,'--',TT,sgmvv2,'x','linewidth',3); hold on;

subplot(131);
% plot(Tr1d,sgmv1d_all(jE,:),':','linewidth',3); hold on;
errorbar(TT,sgmvv1,stdsgmvv1,':x','linewidth',2); hold on;
ylim([0,10e-22]);
xlabel('T_r [keV]'); ylabel('<\sigma{}v> (m^3/s)');
hleg1=legend(['method=1, using f_{1,2}, O(N),',...
    10,'N_1=',num2str(N),', runtime=',num2str(runtime1,2),...
    's,',10,'error=',num2str(100*mean(stdsgmvv1./sgmvv1),2),'%'],'location','best');
legend('boxoff');
set(hleg1,'Fontsize',10);
text(2,0.95e-21,'(a)','Fontsize',10, 'FontWeight','bold');

text(25,0.5e-21,['R_t=',num2str(Rt),', E_d=',num2str(EdkeV),'keV,',...
    10,'v_{d1}/v_d=[',num2str(vd1x/vd),',',...
    num2str(vd1y/vd),',',num2str(vd1z/vd),'], v_{d2}/v_d=[',...
    num2str(vd2x/vd),',',num2str(vd2y/vd),',',num2str(vd2z/vd),']'],...
    'Fontsize',10, 'FontWeight','bold');

subplot(132);
% plot(Tr1d,sgmv1d_all(jE,:),':','linewidth',3); hold on;
errorbar(TT,sgmvv2,stdsgmvv2,':x','linewidth',2); hold on;
ylim([0,10e-22]);
xlabel('T_r [keV]'); %ylabel('<\sigma{}v> (m^3/s)');
hleg2=legend(['method=1, using f_{1,2}, O(N),',...
    10,'N_1=',num2str(N2),', runtime=',num2str(runtime2,2),...
    's,',10,'error=',num2str(100*mean(stdsgmvv2./sgmvv2),2),'%'],...
    'location','best');
legend('boxoff');
set(hleg2,'Fontsize',10);
text(2,0.95e-21,'(b)','Fontsize',10, 'FontWeight','bold');

subplot(133);
% plot(Tr1d,sgmv1d_all(jE,:),':','linewidth',3); hold on;
errorbar(TT,sgmvv3,stdsgmvv3,':x','linewidth',2); hold on;
ylim([0,10e-22]);
xlabel('T_r [keV]'); %ylabel('<\sigma{}v> (m^3/s)');
hleg3=legend(['method=1, using f_{1,2}, O(N),',...
    10, 'N_1=',num2str(N3),', runtime=',num2str(runtime3,2),...
    's,',10,'error=',num2str(100*mean(stdsgmvv3./sgmvv3),2),'%'],...
    'location','best');
legend('boxoff');
set(hleg3,'Fontsize',10);
text(2,0.95e-21,'(c)','Fontsize',10, 'FontWeight','bold');

end

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
% set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
%   'PaperSize',[screenposition(3:4)]);
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[screenposition(3:4)]);

% print(gcf,'-dpdf','-painters','tst.pdf');
% print(gcf,'-dpdf',['smgv_mc_cmp_dtm_Rt=',num2str(Rt),'_Ed=',num2str(EdkeV),'_N1=',num2str(N),...
%     '_N2=',num2str(N2),'_N3=',num2str(N3),'.pdf']);
print(gcf,'-dpng',['smgv_mc_cmp_dtm_Rt=',num2str(Rt),'_Ed=',num2str(EdkeV),'_N1=',num2str(N),...
    '_N2=',num2str(N2),'_N3=',num2str(N3),'.png']);