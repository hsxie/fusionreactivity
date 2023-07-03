% Hua-sheng XIE, huashengxie@gmail.com, 2018-11-22 13:37
% Ref: [1] Bosch & Hale, 1992, fusion cross-sections
% [2] Nevins & Swain, 2000, p- 11 B
% 2020-04-15 23:42 rewrite to function
% 22-03-15 12:35 add range 0.55-4.7MeV
% 22-11-16 08:18 fixed a bug when only E>550keV

function sigmadt=fsgmdt(E)

% E=10.^(0:0.002:3.0); % keV

% 1. D + T -> n + 4He + 17.4MeV
BGdt  =  34.3827;

% 22-11-16 08:20 wrong
% ind1=find(abs(E-550)==min(abs(E-550)));
% % ind1=find(abs(E-5500)==min(abs(E-5500)));
% E1=E(1:ind1); E2=E((ind1+1):end);
ind1=find(E<550); ind2=find(E>=550 & E<=4700);
E1=E(ind1); E2=E(ind2);

% 0.5-550 keV, S_error<1.9%
Adt1   =  [6.927e4,   7.454e8,    2.050e6,    5.2002e4,   0.0   ];
Bdt1   =  [6.38e1,    -9.95e-1,   6.981e-5,   1.728e-4    ];
Sdt1=(Adt1(1)+Adt1(2)*E1+Adt1(3)*E1.^2+Adt1(4)*E1.^3+Adt1(5)*E1.^4)./(1.0+...
    Bdt1(1)*E1+Bdt1(2)*E1.^2+Bdt1(3)*E1.^3+Bdt1(4)*E1.^4);

% 550-4700 keV, S_error<2.5%
Adt2   =  [-1.4714e6,   0.0   ];
Bdt2   =  [-8.4127e-3,    4.7983e-6,   -1.0748e-9,   8.85184e-14    ];
Sdt2=(Adt2(1))./(1.0+...
    Bdt2(1)*E2+Bdt2(2)*E2.^2+Bdt2(3)*E2.^3+Bdt2(4)*E2.^4);

% Sdt=[Sdt1,Sdt2];
Sdt=0.*E;  % 22-11-16 08:27 to check use 0 or NaN
Sdt(ind1)=Sdt1;
Sdt(ind2)=Sdt2;
sigmadt=Sdt./(E.*exp(BGdt./sqrt(E)))*1e-31; % unit: mb (1e-31 m^2) -> m^2


% sigmadt(E>4700 | E<0.5)=NaN;
