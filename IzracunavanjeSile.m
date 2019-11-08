v=20;
h=0.13;
clear J;
clear B0;
clear Biz;
clear f;
clear F;
J=zeros(2,m);
B0=zeros(3,m);
Biz=zeros(1,m);
IndukovanoPoljepom=(IndukovanoPolje(:,iternum))';
fdrag=zeros(1,m);
flift=zeros(1,m);
ftrans=zeros(1,m);
a=10^(-3);
sigma=57*10^6;
Fdrag=0;
Flift=0;
Ftrans=0;
mi0=4*pi*10^(-7);

% dB0zy=@(x,y)-(mi0/(4*pi))*v*(75*x^3+(75*y^2-3)*x)/(25*(x^2+y^2+0.01)^(7/2));%izvod B0z po y
% E=scatteredInterpolant(p1(1,:)',p1(2,:)',resultsX(iternum,:)');

for l=1:m
    nodes=T(:,l);
    hil=TezistaTrouglova(1,l);%x i y koordinate tezista trouglova
    psil=TezistaTrouglova(2,l);
    R=(hil^2+psil^2+h^2)^(1/2);

    Biz(1,l)=(IndukovanoPoljepom(nodes(1,1))+IndukovanoPoljepom(nodes(2,1))+IndukovanoPoljepom(nodes(3,1)))/3;%posto od ranije imamo vrednosti polja u temenima, sad se racuna polje u tezistima
    B0(1,l)=(R^(-5))*(3*h*hil)*(mi0/(4*pi)); 
    B0(2,l)=(R^(-5))*(3*h*psil)*(mi0/(4*pi));
    B0(3,l)=(R^(-3))*(3*h^2*(R^(-2))-1)*(mi0/(4*pi));
    
    GradientX=resultsY(iternum,:);%nisam zamenio X i Y sve je dobro 
    GradientY=resultsX(iternum,:);
    J(1,l)=-sigma*(GradientX(nodes(1,1))+GradientX(nodes(2,1))+GradientX(nodes(3,1)))/3;%x i y komponente vektora struje
    J(2,l)=-sigma*((GradientY(nodes(1,1))+GradientY(nodes(2,1))+GradientY(nodes(3,1)))/3+v*(Biz(1,l)+B0(3,l)));

    fdrag(1,l)=J(2,l)*(B0(3,l)+Biz(1,l));%da li je + ili -
    flift(1,l)=J(1,l)*B0(2,l)-J(2,l)*B0(1,l);
    ftrans(1,l)=-(J(1,l)*(B0(3,l)+Biz(1,l)));
end
f1=scatteredInterpolant(TezistaTrouglova(1,:)',TezistaTrouglova(2,:)',fdrag');
Fdrag = a*integral2(@(x,y) f1(x,y), -r/3,r/3,-r/3,r/3);%trostruki integral je sveden na dvostruki jer nema promene f po z osi
f2=scatteredInterpolant(TezistaTrouglova(1,:)',TezistaTrouglova(2,:)',flift');
Flift=a*integral2(@(x,y) f2(x,y), -r/3,r/3,-r/3,r/3);
f3=scatteredInterpolant(TezistaTrouglova(1,:)',TezistaTrouglova(2,:)',ftrans');
Ftrans=a*integral2(@(x,y) f3(x,y), -r/3,r/3,-r/3,r/3);
