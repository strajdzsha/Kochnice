J=zeros(2,m);
B0=zeros(1,m);
Bez=zeros(1,m);
IndukovanoPoljepom=(IndukovanoPolje(:,iternum))';
f=zeros(1,m);
a=10^(-3);
d=20;
F=0;
v=20;
for l=1:m
    nodes=T(:,l);
    ksil=TezistaTrouglova(1,l);
    psil=TezistaTrouglova(2,l);
    R=(ksil^2+psil^2+h^2)^(1/2);

    Bez(1,l)=(IndukovanoPoljepom(nodes(1,1))+IndukovanoPoljepom(nodes(2,1))+IndukovanoPoljepom(nodes(3,1)))/3;
    B0(1,l)=(R^(-5))*(3*h*ksil);
    B0(2,l)=(R^(-5))*(3*h*psil);
    B0(3,l)=(R^(-3))*(3*h^2*R^(-2)-1);
    
    GradientX=-resultsY(iternum,:);
    GradientY=-resultsX(iternum,:);
    J(1,l)=-(GradientX(nodes(1,1))+GradientX(nodes(2,1))+GradientX(nodes(3,1)))/3;%fale komponente od B0
    J(2,l)=-(GradientY(nodes(1,1))+GradientY(nodes(2,1))+GradientY(nodes(3,1)))/3+v*(Bez(1,l)+B0(3,l));
    f(1,l)=J(1,l)*B0(2,l)-J(2,l)*B0(1,l);
    F=f(1,l)*PovrsineTrouglova(1,l)*a;
end
f1=scatteredInterpolant(TezistaTrouglova(1,:)',TezistaTrouglova(2,:)',f');
fMatrix=f1(I,Q);