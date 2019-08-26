d=500;
xq=linspace(-r,r,d);
yq=linspace(-r,r,d);
[I,J] = meshgrid(xq,yq);
AproxSol=scatteredInterpolant(p1(1,:)',p1(2,:)',results.XGradients);
AproxSolMatrix=AproxSol(I,J);
PropperSol=scatteredInterpolant(p1(1,:)',p1(2,:)',resultsX(iternum,:)');
PropperSolMatrix=PropperSol(I,J);
Bi=scatteredInterpolant(p(1,:)',p(2,:)',X(:,iternum));
BiMatrix=Bi(I,J);
CurrentDensity=PropperSolMatrix+(BiMatrix)*v;
imagesc(CurrentDensity)
plot(xq,CurrentDensity(d/2,:))
grid on
% hold on
% plot(xq,V2(d/2,:));
