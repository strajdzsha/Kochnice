% J1=scatteredInterpolant(TezistaTrouglova(1,:)',TezistaTrouglova(2,:)',B0(2,:)');
% pdeplot(model,'XYData',J1(p1(1,:)',p1(2,:)'))
% R=(ksil^2+psil^2+h^2)^(1/2);
% B0x=@(x,y)(((x^2+y^2+h^2)^(1/2))^(-5))*(3*h*x);
SilaDrag=[-4.9453*10^(-4) -4.8377*10^(-4) -4.7420*10^(-4) -4.6542*10^(-4) -4.5647*10^(-4) -4.4510*10^(-4) -4.3711*10^(-4) -4.2632*10^(-4) -4.1125*10^(-4) -3.8991*10^(-4) -3.5959*10^(-4) -3.1728*10^(-4) -2.5998*10^(-4) -1.8619*10^(-4) -1.6952*10^(-4) -1.5227*10^(-4) -1.3449*10^(-4) -1.1623*10^(-4) -9.7531*10^(-5) -7.8469*10^(-5) -5.9112*10^(-5) -3.9527*10^(-5) -1.9797*10^(-5) -1.9804*10^(-6)];
SilaLift=[3.3603*10^(-4) 3.2696*10^(-4) 3.1535*10^(-4) 3.0018*10^(-4) 2.7991*10^(-4) 2.5229*10^(-4) 2.3479*10^(-4) 2.1422*10^(-4) 1.9015*10^(-4) 1.6230*10^(-4) 1.3076*10^(-4) 9.643*10^(-5) 6.1599*10^(-5) 3.0329*10^(-5) 2.4974*10^(-5) 2.0038*10^(-5) 1.5545*10^(-5) 1.1545*10^(-5) 8.0931*10^(-6) 5.2167*10^(-6) 2.9477*10^(-6) 1.3109*10^(-6) 3.2459*10^(-7) 2.5255*10^(-9)];
SilaTrans=[1.6384*10^(-6) 1.5783*10^(-6) 1.5101*10^(-6) 1.4218*10^(-6) 1.3039*10^(-6) 1.1372*10^(-6) 1.0373*10^(-6) 9.1842*10^(-7) 7.8162*10^(-7) 6.3087*10^(-7) 4.7096*10^(-7) 3.0885*10^(-7) 1.7483*10^(-7) 7.1626*10^(-8) 5.6696*10^(-8) 4.3248*10^(-8) 3.1940*10^(-8) 2.1959*10^(-8) 1.4636*10^(-8) 8.4311*10^(-9) 4.0995*10^(-9) 1.1074*10^(-9) -1.6900*10^(-10) -7.8352*10^(-11)];
Brzine=[100 90 80 70 60 50 45 40 35 30 25 20 15 10 9 8 7 6 5 4 3 2 1 0.1];
h=[0.13 0.11 0.1 0.09 0.08 0.075 0.07 0.065 0.06 0.055 0.05];
% C=((mi0)^2)*(57*10^(6))/(128*pi*h^3);
SilaDragH=[-1.0928*10^(-4) -2.1513*10^(-4) -3.1728*10^(-4) -4.8675*10^(-4) -7.8425*10^(-4) -0.0010 -0.0013 -0.0018 -0.0025 -0.0035 -0.0051];
SilaLiftH=[3.2806*10^(-5) 6.4760*10^(-5) 9.6430*10^(-5) 1.4854*10^(-4) 2.3899*10^(-4) 3.0926-10^(-4) 4.0668*10^(-4) 5.4473*10^(-4) 7.4547*10^(-4) 0.0010 0.0015];
SilaTransH=[8.6073*10^(-9) 1.8184*10^(-7) 3.0885*10^(-7) 5.6846*10^(-7) 1.1104*10^(-6) 1.6233*10^(-6) 2.3877*10^(-6) 3.6586*10^(-6) 5.6765*10^(-6) 9.1275*10^(-6) 1.4949*10^(-5) ];
plot(Brzine,-SilaDrag,'LineWidth',0.9,'Marker','o')
hold on
plot(Brzine,SilaLift,'LineWidth',0.9,'Marker','*')
hold off
plot(Brzine,SilaTrans)
plot(h,-SilaDragH,'LineWidth',0.9,'Marker','o')
SilaDragH1=(-SilaDragH).^(-1/3);
plot(h,SilaDragH1,'LineStyle','none','Marker','o')
hold on
xq=linspace(min(h),max(h),20);
P = polyfit(h,SilaDragH1,1);
yfit = P(1)*xq+P(2);
plot(xq,yfit);
hold off