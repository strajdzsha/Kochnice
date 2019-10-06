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
pdegplot(model,'EdgeLabels','on','FaceLabels','on')
hmax=0.4;
mesh=generateMesh(model,'Hmax',hmax);
[p1,e1,t]=meshToPet(mesh);
pdeplot(model,'NodeLabels','on','ElementLabels','on'); 
applyBoundaryCondition(model,'dirichlet','Edge',[1,2,3,4],'u',0);

M=size(t(1,:));
m=M(1,2);%broj trouglova
T=t([1,2,3],:);%prve tri linije matrice t1 predstavljaju nodove u temenima trouglova/matlab nodama naziva i sredista ivica-nama irelevantni podaci
n1=length(p1);
n=max(T(:));%broj nodova (koji su temena trougla, n1 je ukupan broj nodova-pogledati komentar na liniji 37)
p=zeros(2,n);%nodovi koji su temena trouglova
for i=1:n
    p(:,i)=p1(:,i);
end

PovrsineTrouglova=zeros(1,m);
TezistaTrouglova=zeros(2,m);
for i=1:m
    xy=p(:,T(:,i));
    x=xy(1,:);%koordinate nodova
    y=xy(2,:);%-||-
    TezistaTrouglova(1,i)=(x(1,1)+x(1,2)+x(1,3))/3;
    TezistaTrouglova(2,i)=(y(1,1)+y(1,2)+y(1,3))/3;
    S1=[x',y',[1,1,1]'];
    PovrsineTrouglova(1,i)=det(S1)/2;
end

NodoviPoTrouglovima=cell(n,1);%lista nodova po trouglovima kojima pripadaju
for l=1:m
    for h=1:3
        node=T(h,l);
        NodoviPoTrouglovima{node,1}=[NodoviPoTrouglovima{node,1} l];
    end
end

v=20;
mi0=4*pi*10^(-7);
asigma=57*10^3;

f=@(location,state)-mi0*v*(75*location.x.^3+(75*location.y.^2-3).*location.x)./(25*(location.x.^2+location.y.^2+0.01).^(7/2));
specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',f,'face',1);
results = solvepde(model);
pdeplot(model,'XYData',results.NodalSolution)
%provera scatteredInterpolant
F=scatteredInterpolant(p1(1,:)',p1(2,:)',results.NodalSolution);
xq=linspace(-r,r,20);
yq=linspace(-r,r,20);
[X,Y] = meshgrid(xq,yq);
F1=F(X,Y);