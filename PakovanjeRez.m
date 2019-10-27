d=500;
xq=linspace(-r,r,d);
yq=linspace(-r,r,d);
[I,J] = meshgrid(xq,yq);

AproxSol=scatteredInterpolant(p1(1,:)',p1(2,:)',results.XGradients);
AproxSolMatrix=AproxSol(I,J);
PropperSol=scatteredInterpolant(p1(1,:)',p1(2,:)',resultsX(iternum,:)');
PropperSolMatrix=PropperSol(I,J);
Bi=scatteredInterpolant(p(1,:)',p(2,:)',IndukovanoPolje(:,iternum));
BiMatrix=Bi(I,J);
CurrentDensityMatrix=PropperSolMatrix+(BiMatrix)*v;
plot(xq,AproxSolMatrix(d/2,:))
grid on
hold on
plot(xq,CurrentDensityMatrix(d/2,:))
xlabel('x koordinata za y=0')
ylabel('Gustina struje')
title('Grafik zavisnosti y komponente gustine struje od x koordinate za y=0')
pdeplot(model,'XYData',Bi(p1(1,:)',p1(2,:)')*v+PropperSol(p1(1,:)',p1(2,:)'))
pdeplot(model,'XYData',results1.NodalSolution)
contour(CurrentDensityMatrix,20)

IterDiff=zeros(1,iternum-1);
for i=1:iternum-1
    IterDiff(1,i)=min(resultsX(i+1,:)-resultsX(i,:));
end
niz=linspace(1,iternum-1,iternum-1);
plot(niz,IterDiff)

%3d plot indukovanog polja%
surfBi=surf(I,J,BiMatrix);
set(surfBi,'LineStyle','none')
colormap(parula)
hold on
for i=1:d/10
    plot3(xq,ones(size(xq))*yq(i*10),BiMatrix(i*10,:) ,'LineWidth',0.1,'Color','k')
    plot3(ones(size(xq))*xq(i*10),yq,BiMatrix(:,i*10) ,'LineWidth',0.1,'Color','k')
end
hold off

%3d plot gustine struje (samo y komponenta)%
surfCurrent=surf(I,J,CurrentDensityMatrix);
set(surfCurrent,'LineStyle','none');
colormap(parula)
hold on
for i=1:d/10
    plot3(xq,ones(size(xq))*yq(i*10),CurrentDensityMatrix(i*10,:) ,'LineWidth',0.1,'Color','k')
    plot3(ones(size(xq))*xq(i*10),yq,CurrentDensityMatrix(:,i*10) ,'LineWidth',0.1,'Color','k')
end
hold off

%3d plot elektricnog polja (samo y komponenta)%
surfEField=surf(I,J,PropperSolMatrix);
set(surfEField,'LineStyle','none');
colormap(parula)
hold on
for i=1:d/10
    plot3(xq,ones(size(xq))*yq(i*10),PropperSolMatrix(i*10,:) ,'LineWidth',0.1,'Color','k')
    plot3(ones(size(xq))*xq(i*10),yq,PropperSolMatrix(:,i*10) ,'LineWidth',0.1,'Color','k')
end
hold off

%3d plot elektricnog potencijala%
ElektricniPotencijal=scatteredInterpolant(p1(1,:)',p1(2,:)',results1n(iternum,:)');
ElektricniPotencijalMatrix=ElektricniPotencijal(I,J);
surfEPField=surf(I,J,ElektricniPotencijalMatrix);
set(surfEPField,'LineStyle','none');
colormap(parula)
hold on
for i=1:d/10
    plot3(xq,ones(size(xq))*yq(i*10),ElektricniPotencijalMatrix(i*10,:) ,'LineWidth',0.1,'Color','k')
    plot3(ones(size(xq))*xq(i*10),yq,ElektricniPotencijalMatrix(:,i*10) ,'LineWidth',0.1,'Color','k')
end
hold off
