%----------GEOMETRIJA-----------------------------------------------------------------------------

r=0.5;
circ1=[1
    0
    0
    r
    0
    0
    0
    0
    0
    0];
gd = [circ1];
ns = char('circ1');
ns = ns';
sf = 'circ1';
[dl,bt] = decsg(gd,sf,ns);
pdegplot(dl,'FaceLabels','on','EdgeLabels','on')
axis equal
g=decsg(gd,sf,ns);

%REŠAVANJE PDE (APROKSIMATIVNI SLU?AJ)

model=createpde();
geometryFromEdges(model,g);
pdegplot(model,'EdgeLabels','on')
hmax=0.02;
mesh=generateMesh(model,'Hmax',hmax);
[p,e,t]=meshToPet(mesh);
pdeplot(model); 
applyBoundaryCondition(model,'dirichlet','Edge',[1,2,3,4],'u',0);
figure; 
pdeplot(model); 
axis equal
v=20;
mi0=4*pi*10^(-7);
f=@(location,state)-mi0*v*(75*location.y.^3+(75*location.x.^2-3).*location.y)./(25*(location.y.^2+location.x.^2+0.01).^(7/2));
specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',f,'face',1);
results = solvepde(model);
pdeplot(model,'XYData',results.NodalSolution)

%PAKOVANJE REZULTATA I RA?UNANJE J

w=150;
d=2*w+1;
xq=linspace(-r,r,d);
yq=linspace(-r,r,d);
[X,Y] = meshgrid(xq,yq);
U=scatteredInterpolant(p(1,:)',p(2,:)',results.YGradients);
V=U(X,Y);
plot(xq,V(:,w));



