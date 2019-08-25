%GEOMETRIJA%

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
hmax=0.02;%ne sitniti!
mesh=generateMesh(model,'Hmax',hmax);


%DRNDANJE SA TROUGLOVIMA%

[p,e,t]=meshToPet(mesh);
M=size(t(1,:));
N=size(p(1,:));
n=N(1,2);%broj nodova
m=M(1,2);%broj trouglova
O=zeros(2,m);
T=t([1,2,3],:);
S=zeros(1,m);

for i=1:m
    xy=p(:,T(:,i));
    x=xy(1,:);%koordinate nodova
    y=xy(2,:);%-||-
    O(1,i)=(x(1,1)+x(1,2)+x(1,3))/3;
    O(2,i)=(y(1,1)+y(1,2)+y(1,3))/3;
    S1=[x',y',[1,1,1]'];
    S(1,i)=det(S1);
end

A=cell(n);%lista nodova po trouglovima kojima pripadaju
for l=1:m
    for h=1:3
        node=T(h,l);
        A{node}=[A{node} [l]];
    end
end

%REŠAVANJE PDE%

applyBoundaryCondition(model,'dirichlet','Edge',[1,2,3,4],'u',0);
v=20;
mi0=4*pi*10^(-7);
asigma=57*10^3;

B0z=@(x,y)(3*0.01./(0.01+x.^2+y.^2).^(5/2)-1./(0.01+x.^2+y.^2).^(3/2))*(mi0/4*pi);
iternum=10;%broj iteracija
X=zeros(n,iternum);%resenja za Bez

resultsX=zeros(iternum,n);
resultsY=zeros(iternum,n);
resultsX(1,:)=results.XGradients;
resuktsY(1,:)=results.YGradients;

Znamo=Desno(m,resultsX(1,:),resuktsY(1,:),O,n,v,B0z,asigma,T,p,S,mi0);
Alpha=SklapanjeMatrice(n,asigma,mi0,A,p,O,S,v);
X(:,1)=Alpha\Znamo;

BYO=zeros(1,m);%izvod Bez po y
results1n=zeros(iternum,n);

for i=1:iternum
    
    for j=1:m%racunanje BYO na osnovu predhodnog resenja/racuna se u tezistima trouglova
        nodes=T(:,j);
        Xpom=X(:,i);
        bez=Xpom(nodes');
        x=p(1,nodes);
        y=p(2,nodes);
        XY=[ones(1,3)',x',y'];
        ABC=XY\bez;
        BYO(1,j)=ABC(3,1);%koeficijent C je izvod po y
    end
    
    F=scatteredInterpolant(O(1,:)',O(2,:)',BYO');%interpoliranje BYO po celom mesh-u
    f=@(location,state)-F(location.x,location.y)-mi0*v*(75*location.y.^3+(75*location.x.^2-3).*location.y)./(25*(location.y.^2+location.x.^2+0.01).^(7/2));
    
    specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',f,'face',1);
    
    results1 = solvepde(model);
    results1n(i,:)=results1.NodalSolution;
    pdeplot(model,'XYData',results1.NodalSolution)
    resultsX(i+1,:)=results1.XGradients;
    resultsY(i+1,:)=results1.YGradients;
    
    if i<iternum
        Znamo=Desno(m,resultsX(i+1,:),resultsY(i+1,:),O,n,v,B0z,asigma,T,p,S,mi0);
        Alpha=SklapanjeMatrice(n,asigma,mi0,A,p,O,S,v);
        X(:,i+1)=Alpha\Znamo;
    end
    
end


function beta=BioSavart(k,K,m,p,O,S,mi0)%k-broj node za koje se poziva beta;K-vektor gustine struje dimenzija 2xm
Bezk=0;
for l=1:m
    xk=p(1,k);
    yk=p(2,k);
    ksil=O(1,l);
    psil=O(2,l);
    dx=xk-ksil;
    dy=yk-psil;
    Sl=S(1,l);
    cross=K(1,l).*dy-K(2,l).*dx;
    Bezk=Bezk+cross*Sl/(dx^2+dy^2)^(3/2);
end
beta=Bezk*mi0/(4*pi);
end

function Znamo=Desno(m,FX,FY,O,n,v,B0z,asigma,T,p,S,mi0)%racunanje desne strane (12)

 K2=zeros(2,m);%ovo je gradijent fi
 Znamo=zeros(n,1);
 
 for l=1:m
     nodes=T(:,l);
     K2(1,l)=(FX(nodes(1,1))+FX(nodes(2,1))+FX(nodes(3,1)))/3;
     K2(2,l)=(FY(nodes(1,1))+FY(nodes(2,1))+FY(nodes(3,1)))/3;
 end
 
 K3=zeros(2,m);%ovo je B0z
 K3(1,:)=0;
 
 for l=1:m
     K3(2,l)=B0z(O(1,l),O(2,l))*v;
 end
 
 for k=1:n
     Znamo(k,1)=asigma*(BioSavart(k,-K2,m,p,O,S,mi0)+BioSavart(k,-K3,m,p,O,S,mi0));
 end
 
end

function Alpha=SklapanjeMatrice(n,asigma,mi0,A,p,O,S,v)

Alpha=zeros(n,n);
C=-mi0*v/(12*pi*asigma);

for k=1:n
    for j=1:n
        B=A{j};
        b=length(B);
        for q=1:b
            xk=p(1,j);
            yk=p(2,j);
            ksil=O(1,B(q));
            psil=O(2,B(q));
            dx=xk-ksil;
            dy=yk-psil;
            dl=(dx.^2+dy.^2).^(3/2);
            Sl=S(1,B(q));
            Alpha(k,j)=Alpha(k,j)+dx*Sl/dl;
        end
    end
end
Alpha=Alpha.*C;

for i=1:n%dodavanje jedinica po dijagonali
    Alpha(i,i)=Alpha(i,i)+1;
end

end
