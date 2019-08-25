d=500;
xq=linspace(-r,r,d);
yq=linspace(-r,r,d);
[I,J] = meshgrid(xq,yq);
% uintrp = interpolateSolution(results1n(7,:),I,J);
U=scatteredInterpolant(p(1,:)',p(2,:)',resultsY(iternum,:)');
V=U(I,J);
Bi=scatteredInterpolant(p(1,:)',p(2,:)',X(:,iternum));
V1=Bi(I,J);
K=V+(V1)*v;
imagesc(K)
plot(xq,K(:,d/2))
grid on
