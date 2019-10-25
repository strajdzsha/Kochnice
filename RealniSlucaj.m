%GEOMETRIJA%
v=100;
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

%MESH%

model=createpde();
geometryFromEdges(model,g);
hmax=0.015;%ne sitniti!
mesh=generateMesh(model,'Hmax',hmax);


%DRNDANJE SA TROUGLOVIMA%

[p1,e1,t]=meshToPet(mesh);%p1-matrica koordinata svih nodova;e1-nebitno;t-matrica 7*(broj trouglova),sadrzi koordinate nodova koje pripadaju tom trouglu 
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

%REŠAVANJE PDE%

applyBoundaryCondition(model,'dirichlet','Edge',[1,2,3,4],'u',0);
mi0=4*pi*10^(-7);
asigma=57*10^3;

B0z=@(x,y)(3*0.01./(0.01+x.^2+y.^2).^(5/2)-1./(0.01+x.^2+y.^2).^(3/2))*(mi0/(4*pi));%z-komponenta polja dipola 
iternum=10;%broj iteracija
IndukovanoPolje=zeros(n,iternum);%resenja za Bez

resultsX=zeros(iternum,n1);
resultsY=zeros(iternum,n1);
resultsX(1,:)=results.XGradients;
resuktsY(1,:)=results.YGradients;

DesnaStranaJednacine=Desno(m,resultsY(1,:),resultsX(1,:),TezistaTrouglova,n,v,B0z,asigma,T,p,PovrsineTrouglova,mi0);
Alpha=SklapanjeMatrice(n,asigma,mi0,NodoviPoTrouglovima,p,TezistaTrouglova,PovrsineTrouglova,v);
IndukovanoPolje(:,1)=Alpha\DesnaStranaJednacine;

dBezy=zeros(1,m);%izvod Bez po y
results1n=zeros(iternum,n1);

for i=1:iternum
    
    for j=1:m%racunanje BYO na osnovu predhodnog resenja/racuna se u tezistima trouglova
        nodes=T(:,j);
        IndukovanoPoljepom=IndukovanoPolje(:,i);
        bez=IndukovanoPoljepom(nodes');
        x=p(1,nodes);
        y=p(2,nodes);
        XY=[ones(1,3)',x',y'];
        ABC=XY\bez;%ova sintaksa sledi iz FEM fromulacije/vrednost neke funkcije se aproksimira kao U=a+b*x+c*y
        dBezy(1,j)=ABC(3,1);%koeficijent C je izvod po y
    end
    
    F=scatteredInterpolant(TezistaTrouglova(1,:)',TezistaTrouglova(2,:)',dBezy');%interpoliranje dBezy po celom mesh-u
    f=@(location,state)-F(location.x,location.y)-(mi0/(4*pi))*v*(75*location.x.^3+(75*location.y.^2-3).*location.x)./(25*(location.x.^2+location.y.^2+0.01).^(7/2));
    
    specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',f,'face',1);
    
    results1 = solvepde(model);
    results1n(i,:)=results1.NodalSolution;
    pdeplot(model,'XYData',results1.NodalSolution)
    resultsX(i+1,:)=results1.XGradients;
    resultsY(i+1,:)=results1.YGradients;
    
    if i<iternum
        DesnaStranaJednacine=Desno(m,resultsY(i+1,:),resultsX(i+1,:),TezistaTrouglova,n,v,B0z,asigma,T,p,PovrsineTrouglova,mi0);
%         Alpha=SklapanjeMatrice(n,asigma,mi0,NodoviPoTrouglovima,p,TezistaTrouglova,PovrsineTrouglova,v);
        IndukovanoPolje(:,i+1)=Alpha\DesnaStranaJednacine;
    end
    
end


function Beta=BioSavart(k,K,m,p,TezistaTrouglova,PovrsineTrouglova,mi0)%k-broj node za koje se poziva beta;K-vektor gustine struje dimenzija 2xm
Bezk=0;
for l=1:m
    xk=p(1,k);
    yk=p(2,k);
    ksil=TezistaTrouglova(1,l);
    psil=TezistaTrouglova(2,l);
    dx=xk-ksil;
    dy=yk-psil;
    Sl=PovrsineTrouglova(1,l);
    cross=K(1,l).*dy-K(2,l).*dx;
    Bezk=Bezk+cross*Sl/(dx^2+dy^2)^(3/2);
end
Beta=Bezk*mi0/(4*pi);
end

function DesnaStranaJednacine=Desno(m,FX,FY,TezistaTrouglova,n,v,B0z,asigma,T,p,PovrsineTrouglova,mi0)%racunanje desne strane (12)

 K2=zeros(2,m);%ovo je gradijent fi
 DesnaStranaJednacine=zeros(n,1);
 
 for l=1:m
     nodes=T(:,l);
     K2(1,l)=(FX(nodes(1,1))+FX(nodes(2,1))+FX(nodes(3,1)))/3;
     K2(2,l)=(FY(nodes(1,1))+FY(nodes(2,1))+FY(nodes(3,1)))/3;
 end
 
 K3=zeros(2,m);%ovo je B0z
 K3(1,:)=0;
 
 for l=1:m
     K3(2,l)=B0z(TezistaTrouglova(1,l),TezistaTrouglova(2,l))*v;
 end
 
 for k=1:n
     DesnaStranaJednacine(k,1)=asigma*(BioSavart(k,-K2-K3,m,p,TezistaTrouglova,PovrsineTrouglova,mi0));
 end
 
end

function Alpha=SklapanjeMatrice(n,asigma,mi0,NodoviPoTrouglovima,p,TezistaTrouglova,PovrsineTrouglova,v)

Alpha=zeros(n,n);
C=-mi0*v*asigma/(12*pi);

for k=1:n
    for j=1:n
        B=NodoviPoTrouglovima{j};%trazimo kojim sve trouglovima pripada noda j
        b=length(B);
        for q=1:b
            xk=p(1,k);
            yk=p(2,k);
            ksil=TezistaTrouglova(1,B(q));
            psil=TezistaTrouglova(2,B(q));
            dx=xk-ksil;
            dy=yk-psil;
            dl=(dx.^2+dy.^2).^(3/2);
            Sl=PovrsineTrouglova(1,B(q));
            Alpha(k,j)=Alpha(k,j)+dx*Sl/dl;
        end
    end
end
Alpha=Alpha.*C;

for i=1:n%dodavanje jedinica po dijagonali
    Alpha(i,i)=Alpha(i,i)+1;
end

end
