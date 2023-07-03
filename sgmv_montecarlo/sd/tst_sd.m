% 2023-02-12 21:13 Hua-sheng XIE, huashengxie@gmail.com, ENN
% tst slowing down sigmv*v
% 23-02-19 11:10 update

close all;clear;clc;

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

TT=5:5:150;

EdkeV=100; % keV
vd=sqrt(2*EdkeV*(qe*1e3)/mr);

sgmvv1=0.*TT; sgmvv2=0.*TT; sgmvv3=0.*TT; sgmvv4=0.*TT;
stdsgmvv1=0.*TT; stdsgmvv2=0.*TT; stdsgmvv3=0.*TT; stdsgmvv4=0.*TT;
runtime1=0; runtime2=0; runtime3=0; runtime4=0;
for j=1:length(TT)
    TrkeV=TT(j); % keV
    
    % sgmvdt=fsgmv(TrkeV,1);
    
    Tr=TrkeV*qe*1e3/kB; % keV -> K
%     T1r=Tr; T2r=2.0*Tr;
    vt=sqrt(kB*Tr/m1);
    vc1=vt; vc2=sqrt(2)*vt;
    vb1=1*vd; vb2=2.0*vd;
    
    N=40000; N3=10*N; N2=5*floor(sqrt(N)); N4=10*N;
    tmp1=cputime;
    [sgmv1,stdsgmv1]=fsgmvmcsd(vb1,vc1,vb2,vc2,N,1);
    tmp1=cputime-tmp1; runtime1=runtime1+tmp1;
    tmp2=cputime;
    [sgmv2,stdsgmv2]=fsgmvmcsd(vb1,vc1,vb2,vc2,N2,2);
    tmp2=cputime-tmp2; runtime2=runtime2+tmp2;
    tmp3=cputime;
    [sgmv3,stdsgmv3]=fsgmvmcsd(vb1,vc1,vb2,vc2,N3,3);
    tmp3=cputime-tmp3; runtime3=runtime3+tmp3;
    tmp4=cputime;
    [sgmv4,stdsgmv4]=fsgmvmcsd(vb1,vc1,vb2,vc2,N4,4);
    tmp4=cputime-tmp4; runtime4=runtime4+tmp4;
    
    sgmvv1(j)=sgmv1; sgmvv2(j)=sgmv2; sgmvv3(j)=sgmv3; sgmvv4(j)=sgmv4;
    stdsgmvv1(j)=stdsgmv1; stdsgmvv2(j)=stdsgmv2; stdsgmvv3(j)=stdsgmv3;
    stdsgmvv4(j)=stdsgmv4;
end

%%
close all;
if(1==0)
figure('unit','normalized','DefaultAxesFontSize',16,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,...
    'position',[0.01,0.05,0.45,0.55]);

% plot(TT,sgmvv1,'--',TT,sgmvv2,'x','linewidth',3); hold on;
% plot(TT,sgmvv1,'--','linewidth',3); hold on;
errorbar(TT,sgmvv1,stdsgmvv1,'--','linewidth',3); hold on;
errorbar(TT,sgmvv2,stdsgmvv2,'x','linewidth',3); hold on;
errorbar(TT,sgmvv3,stdsgmvv3,'x','linewidth',3); hold on;
xlabel('T_r [keV]');
ylabel('<\sigma{}v> (m^3/s)');
legend('method=1, using f_{1,2}','method=2, using Gaussian with weight',...
    'method=3, N\times N','location','best');
legend('boxoff');
ylim([0,1e-21]);

title(['E_d=',num2str(EdkeV),'keV',', N_1=',num2str(N),...
    ', N_2=',num2str(N2),', N_3^2=',num2str(N3^2),10,'v_{1}=[',num2str(vc1/vd),',',...
    num2str(vb1/vd),'], v_{2}=[',...
    num2str(vc2/vd),',',num2str(vb2/vd),']']);
else
    
figure('unit','normalized','DefaultAxesFontSize',13,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,...
    'position',[0.01,0.05,0.75,0.4]);

% plot(TT,sgmvv1,'--',TT,sgmvv2,'x','linewidth',3); hold on;

subplot(131);
errorbar(TT,sgmvv1,stdsgmvv1,':x','linewidth',2); hold on;
% errorbar(TT,sgmvv4,stdsgmvv4,'.','linewidth',2); hold on;
xlabel('T_r [keV]'); ylabel('<\sigma{}v> (m^3/s)');
ylim([0,10e-22]);
hleg1=legend(['method=1, using f_{1,2}, O(N),',...
    10,'N_1=',num2str(N),', runtime=',num2str(runtime1,2),...
    's,',10,'error=',num2str(100*mean(stdsgmvv1./sgmvv1),2),'%'],'location','best');
legend('boxoff');
set(hleg1,'Fontsize',10);
text(2,9.5e-22,'(a)','Fontsize',10, 'FontWeight','bold');

text(20,0.5e-21,['E_d=',num2str(EdkeV),'keV,',10,...
    '[v_{c1}/v_t,v_{b1}/v_d]=[',num2str(vc1/vt,3),',',...
    num2str(vb1/vd,3),'],',10,'[v_{c2}/v_t,v_{b2}/v_d]=[',...
    num2str(vc2/vt,3),',',num2str(vb2/vd,3),']']);

% text(10,0.5e-21,['E_d=',num2str(EdkeV),'keV,',...
%     10,'v_{d1}/v_d=[',num2str(vd1x/vd),',',...
%     num2str(vd1y/vd),',',num2str(vd1z/vd),',',num2str(vd1r/vd),'],',10,'v_{d2}/v_d=[',...
%     num2str(vd2x/vd),',',num2str(vd2y/vd),',',num2str(vd2z/vd),',',num2str(vd2r/vd),']'],...
%     'Fontsize',10, 'FontWeight','bold');

subplot(132);
errorbar(TT,sgmvv2,stdsgmvv2,':x','linewidth',2); hold on;
xlabel('T_r [keV]'); %ylabel('<\sigma{}v> (m^3/s)');
ylim([0,10e-22]);
hleg2=legend(['method=2, using f_{1,2}, O(N^2),',...
    10,'N_2=',num2str(N2),', runtime=',num2str(runtime2,2),...
    's,',10,'error=',num2str(100*mean(stdsgmvv2./sgmvv2),2),'%'],...
    'location','best');
legend('boxoff');
set(hleg2,'Fontsize',10);
text(2,0.95e-21,'(b)','Fontsize',10, 'FontWeight','bold');

subplot(133);
errorbar(TT,sgmvv3,stdsgmvv3,':x','linewidth',2); hold on;
errorbar(TT,sgmvv4,stdsgmvv4,'.','linewidth',2); hold on;
xlabel('T_r [keV]'); %ylabel('<\sigma{}v> (m^3/s)');
ylim([0,10e-22]);
hleg3=legend(['method=3, using Gaussian',10,...
    'with weight, O(N),',10, 'N_3=',num2str(N3),', runtime=',num2str(runtime3,2),...
    's,',10,'error=',num2str(100*mean(stdsgmvv3./sgmvv3),2),'%'],...
    ['method=3, using uniform rand',10,...
    'with weight, O(N),',10, 'N_3=',num2str(N3),', runtime=',num2str(runtime4,2),...
    's,',10,'error=',num2str(100*mean(stdsgmvv4./sgmvv4),2),'%'],...
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
% print(gcf,'-dpdf',['smgv_mc_sd_Ed=',num2str(EdkeV),'_N1=',num2str(N),...
%     '_N2=',num2str(N2),'_N3=',num2str(N3),'.pdf']);
print(gcf,'-dpng',['smgv_mc_sd_Ed=',num2str(EdkeV),'_N1=',num2str(N),...
    '_N2=',num2str(N2),'_N3=',num2str(N3),'.png']);