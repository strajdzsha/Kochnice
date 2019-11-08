%----------GEOMETRIJA-----------------------------------------------------------------------------
v=20;
h=0.13;
r=0.6;
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
pdegplot(model,'EdgeLabels','on','FaceLabels','on')
hmax=0.02;
mesh=generateMesh(model,'Hmax',hmax);
[p,e,t]=meshToPet(mesh);
pdeplot(model); 
applyBoundaryCondition(model,'dirichlet','Edge',[1,2,3,4],'u',0);
figure; 
pdeplot(model); 
axis equal
mi0=4*pi*10^(-7);
f=@(location,state)-(mi0*v/(4*pi)).*(3*location.x.*(location.x.^4+(2*location.y.^2-3*h^2).*location.x.^2+location.y.^4-3*(h*location.y).^2-4*h^4))./((location.x.^2+location.y.^2+h^2).^(9/2));
specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',f,'face',1);
results = solvepde(model);
pdeplot(model,'XYData',results.NodalSolution)

%PAKOVANJE REZULTATA I RACUNANJE J/rezultati se pakuju uz pomoc
%scatteredInterpolant koji interpolira niz na ceo prostor zadat
%koordinatama tacaka u kojima su vrednosti niza racunate/ovo se jasnije
%vidi u skripti PakovanjeRez

w=150;
d=2*w+1;
xq=linspace(-r,r,d);
yq=linspace(-r,r,d);
[X,Y] = meshgrid(xq,yq);
U=scatteredInterpolant(p(1,:)',p(2,:)',results.XGradients);
V=U(X,Y);
plot(xq,V(w,:));



