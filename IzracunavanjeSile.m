clear J;
clear B0;
clear Bez;
clear f;
clear F;
clear h;
p=30;
J=zeros(2,m);
B0=zeros(3,m);
Bez=zeros(1,m);
IndukovanoPoljepom=(IndukovanoPolje(:,iternum))';
f=zeros(1,m);
a=10^(-3);
sigma=57*10^6;
F=0;
h=0.1;
v=5;
for l=1:m
    nodes=T(:,l);
    ksil=TezistaTrouglova(1,l);%x i y koordinate tezista trouglova
    psil=TezistaTrouglova(2,l);
    R=(ksil^2+psil^2+h^2)^(1/2);

    Bez(1,l)=(IndukovanoPoljepom(nodes(1,1))+IndukovanoPoljepom(nodes(2,1))+IndukovanoPoljepom(nodes(3,1)))/3;%posto od ranije imamo vrednosti polja u temenima, sad se racuna polje u tezistima
    B0(1,l)=(R^(-5))*(3*h*ksil); 
    B0(2,l)=(R^(-5))*(3*h*psil);
    B0(3,l)=(R^(-3))*(3*h^2*R^(-2)-1);
    
    GradientX=-resultsY(iternum,:);%nisam zamenio X i Y sve je dobro 
    GradientY=-resultsX(iternum,:);
    J(1,l)=-sigma*(GradientX(nodes(1,1))+GradientX(nodes(2,1))+GradientX(nodes(3,1)))/3;%x i y komponente vektora struje
    J(2,l)=-sigma*((GradientY(nodes(1,1))+GradientY(nodes(2,1))+GradientY(nodes(3,1)))/3+v*(Bez(1,l)+B0(3,l)));
    %f(1,l)=J(1,l)*B0(2,l)-J(2,l)*B0(1,l);
    f(1,l)=J(2,l)*B0(3,l);
end
f1=scatteredInterpolant(TezistaTrouglova(1,:)',TezistaTrouglova(2,:)',f');
F = a*integral2(@(x,y) f1(x,y), -r,r,-r,r);%trostruki integral je sveden na dvostruki jer nema promene f po z osi