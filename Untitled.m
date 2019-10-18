
J1=scatteredInterpolant(TezistaTrouglova(1,:)',TezistaTrouglova(2,:)',B0(2,:)');
pdeplot(model,'XYData',J1(p1(1,:)',p1(2,:)'))
R=(ksil^2+psil^2+h^2)^(1/2);
B0x=@(x,y)(((x^2+y^2+h^2)^(1/2))^(-5))*(3*h*x);
for q=1:p
    Max(1,q)=max(f2{1,q}(p1(1,:)',p1(2,:)'));
end
f2=[6.79 34 68.2 103 137 172]*10^8;
f3=[8.0506 3.54 2.68 2.01 1.34 0.671]*10^10;
