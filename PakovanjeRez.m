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
CurrentDensity=PropperSolMatrix+(BiMatrix)*v;
imagesc(CurrentDensity)
plot(xq,CurrentDensity(d/2,:))
grid on
hold on
plot(xq,AproxSolMatrix(d/2,:))
pdeplot(model,'XYData',Bi(p1(1,:)',p1(2,:)')*v+PropperSol(p1(1,:)',p1(2,:)'))
pdeplot(model,'XYData',results1.NodalSolution)
contour(CurrentDensity)
X=zeros(2,n1);
X(1,1)=results1;
% hold on
% plot(xq,v*BiMatrix(d/2,:));
