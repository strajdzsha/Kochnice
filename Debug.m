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
hmax=0.015;
mesh=generateMesh(model,'Hmax',hmax);
[p1,e1,t]=meshToPet(mesh);
pdeplot(model,'NodeLabels','on','ElementLabels','on'); 
applyBoundaryCondition(model,'dirichlet','Edge',[1,2,3,4],'u',0);

m=length(t(1,:));%broj trouglova
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
h=0.1;

iternum=3;
Rezultati=zeros(n1,iternum);
RezultatiX=zeros(n1,iternum);
RezultatiY=zeros(n1,iternum);

f=@(location,state)(mi0*v/(4*pi)).*(3*location.y.*(location.y.^4+(2*location.x.^2-3*h^2).*location.y.^2+location.x.^4-3*(h*location.x).^2-4*h^4))./((location.y.^2+location.x.^2+h^2).^(9/2));
specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',f,'face',1);
results = solvepde(model);
pdeplot(model,'XYData',results.YGradients)
Rezultati(:,1)=results.NodalSolution;
RezultatiX(:,1)=results.XGradients;
RezultatiY(:,1)=results.YGradients;
B0z = @(x,y) (mi0/(4*pi))*((x^2+y^2+h^2)^(-3/2))*((3*(h^2)/(x^2+y^2+h^2))-1);

for i=1:iternum
    Alpha = sklapanjeMatrice(n, mi0 ,v, p, asigma, TezistaTrouglova, PovrsineTrouglova, NodoviPoTrouglovima);
    D=PopunjavanjeMatrice(mi0,n, m, p, v, T, B0z, asigma, TezistaTrouglova, PovrsineTrouglova, RezultatiX(:,i), RezultatiY(:,i));
    Biz=Alpha\D;
    dBizy=zeros(1,m);
    for l=1:m
        node=T(:,l);
        x=p(1,node);
        y=p(2,node);
        XY=[x',y',[1 1 1]'];
        Bizpom=Biz(node);
        ABC = XY\Bizpom;
        dBizy(1,l)=ABC(2,1);
    end
    dBizyInterpolirano=scatteredInterpolant(TezistaTrouglova(1,:)',TezistaTrouglova(2,:)',dBizy');
    f=@(location,state)dBizyInterpolirano(location.x,location.y)+(mi0*v/(4*pi)).*(3*location.y.*(location.y.^4+(2*location.x.^2-3*h^2).*location.y.^2+location.x.^4-3*(h*location.x).^2-4*h^4))./((location.y.^2+location.x.^2+h^2).^(9/2));
    specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',f,'face',1);
    results = solvepde(model);
    Rezultati(:,i+1)=results.NodalSolution;
    RezultatiX(:,i+1)=results.XGradients;
    RezultatiY(:,i+1)=results.YGradients; 
end

BizInterpolirano=scatteredInterpolant(p(1,:)',p(2,:)',Biz);
pdeplot(model,'XYData',BizInterpolirano(p1(1,:)',p1(2,:)'))

function Alpha = sklapanjeMatrice(n, mi0 ,v, p, asigma, TezistaTrouglova, PovrsineTrouglova, NodoviPoTrouglovima)
Alpha = zeros(n,n);
for k = 1:n
    for j = 1:n
        Trouglovi = NodoviPoTrouglovima{j};
        BrojTrouglova = length(Trouglovi);
        for q = 1:BrojTrouglova
            hil=TezistaTrouglova(1,Trouglovi(q));
            psil=TezistaTrouglova(2,Trouglovi(q));
            Sl=PovrsineTrouglova (1, Trouglovi(q));
            xk=p(1,k);
            yk=p(2,k);
            dxk=xk-hil;
            dyk=yk-psil;
            dlk=(dxk.^2+dyk.^2).^3/2;
            Alpha(k,j)=Alpha(k,j)+dxk*Sl/dlk;
        end
    end
end
Alpha=Alpha.*asigma.*(-mi0*v/(12*pi));
for i = 1:n
    Alpha(i,i)=Alpha(i,i)+1;
end
end 

function D = PopunjavanjeMatrice (mi0,n, m, p, v, T, B0z, asigma, TezistaTrouglova,PovrsineTrouglova, XGradients, YGradients)
    K = zeros(2,m);
    D = zeros(n,1);
    for l = 1:m
        node = T(:,l);
        hil = TezistaTrouglova(1,l);
        psil = TezistaTrouglova(2,l);
        K (1,l) = -asigma*(XGradients(node(1,1))+XGradients(node(2,1))+XGradients(node(3,1)))/3;
        K (2,l) =-asigma*((YGradients(node(1,1))+YGradients(node(2,1))+YGradients(node(3,1)))/3-v*B0z(hil,psil));
    end
    for k=1:n
        D(k,1)=BioSavart(k,K,TezistaTrouglova,p,PovrsineTrouglova,mi0,m);
    end
end

function beta = BioSavart(k,K,TezistaTrouglova,p,PovrsineTrouglova,mi0,m)
beta=0;
for l = 1:m
    hil=TezistaTrouglova(1,l);
    psil=TezistaTrouglova(2,l);
    Sl=PovrsineTrouglova (1,l);
    xk=p(1,k);
    yk=p(2,k);
    dxk=xk-hil;
    dyk=yk-psil;
    dlk=(dxk.^2+dyk.^2).^(3/2);
    beta=beta+(K(1,l).*dyk-K(2,l).*dxk)*Sl./dlk;
end
beta=beta*mi0/(4*pi);
end