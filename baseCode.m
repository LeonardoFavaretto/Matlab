Anexo 1 – Código em malab.: versão 2017.a

clearvars %===============================================================
%A primeira etapa é gerar uma matriz livre de escala para determinar as===
%conexões entre os agentes. Para tal utilizaremos um algoritmo B-A.=======
%=========================================================================
%modelo(100,100,10,5,0.1,0.05,0.005,1.1)           funcao salva para teste
%0.061506462    0.034349226 0.252879718 0.520964498 25.39376737 0.609493363 parametros calibrados
n=100;  %número de agentes
t=100;      %número de rodadas de tempo
kc=10;      %knowlegde categories
ln=5;      %location    number
p1= 0.061506462;%0.00300584244790664;     %probability of activate migration function    
p2=0.034349226;      %chance conexão na mesma localidade
p3=0.252879718;%chance conexão aleatória
p4=0.720964498;    %0.520964498;;   %redução efeito popução
p5 = 25.39376737;%24.9137680673842;   % importância diversidade na migração
p6 = 0.609493363; %time coefficient reduction
       %diversity multiplier
%r = (b-a).*rand() + a;     %random knowledge share between 1 and 5% a=0.01; b=0.05;
%=========================================================================
 
Nodes=n;
mlinks=1;
seed=1;
seed = full(seed);
pos = length(seed);
rand('state',sum(100*clock));
Net = zeros(Nodes, Nodes, 'single');
Net(1:pos,1:pos) = seed;
sumlinks = sum(sum(Net));
while pos < Nodes
    pos = pos + 1;
    linkage = 0;
    while linkage ~= mlinks
        rnode = ceil(rand * pos);
        deg = sum(Net(:,rnode)) * 2;
        rlink = rand * 1;
        if rlink < deg / sumlinks && Net(pos,rnode) ~= 1 && Net(rnode,pos) ~= 1
            Net(pos,rnode) = 1;
            Net(rnode,pos) = 1;
            linkage = linkage + 1;
            sumlinks = sumlinks + 2;
        end
    end
end
clear Nodes deg linkage pos rlink rnode sumlinks mlinks
SFNet = Net;
 
%=========================================================================
%Gráfico da geometria de rede=============================================
%=========================================================================
%format compact
%format long e
%theta = linspace(0,2*pi,length(Net)+1);
%xy = zeros(length(Net)+1,2);
%x = cos(theta);
%y = sin(theta);
%xy(1:length(Net)+1,1) = x(1:length(Net)+1);
%xy(1:length(Net)+1,2) = y(1:length(Net)+1);
%figure, gplot(Net,xy,'.-');
%set(gcf, 'Color', [1 1 1]);
%axis('equal');
%xlim([-1.1 1.1]);
%ylim([-1.1 1.1]);
%axis off;
 
%========================================================================
%========================================================================
%Gera números aleatórios com distribuição log-normal com média 0 e dp 1.
%Location = randi(ln,n,1);
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Looping@@@@@@@@@@@@@
K = lognrnd(0,1,n,kc);      %K stands for knowledge
L=[1:n; randi(ln,1,n)];   %localização por agente
%============================================
 
 
 
 
%================================ Criação variável diversidade
%================Distribuição normal, 3 variáveis
Div = normrnd(0,1,n,3);         % sexualidade, orientação política, grupo étnico
%Cria variável de 1 3 5, dependendo dos quadrantes.
% entre 1 e -1, assume 1. Entre 2 e 1 e -2 e -1 assume 3. maior que 2,e
% menor que -2 assume 5.
 
 
for i=1:n
    for j=1:3
if Div(i,j)<= 1 && Div(i,j) >= -1
    Div1(i,j) = 1;
elseif Div(i,j) > 1 && Div(i,j) <=2
    Div1(i,j) = 2;
elseif Div(i,j) < -1 && Div(i,j) >=-2
 Div1(i,j)= 2;
else
Div1(i,j)=3;
end
    end
end
 
Div2= sum(Div1,2);
 
 
 
 
for leo=1:t
 
%=============== Knowledge exchange======================================
 
for i=1:n
    for j=1:n
if Net(i,j)==1 && rand()<=0.05 + (L(2,i)-L(2,j))/100          %Care here
K1(i,:)= K(i,:)+ 0.05*K(j,:);    %###############################
else
    K1(i,:)= K(i,:);
 
end
    end
end
 
 
%========================================================================
% Migration Padron ======================================================
%========================================================================
 
%Kcoef=sum(K,2)/100;              %  Knowledge coeficient
 
for i=1:n
for j=1:n
T(i,j)=sum(abs(K(i,:)-K(j,:)));  %afinity matrix, the smaller the bigger 
  
end
end
 
 
[sim, index]=sort(T,2);
for i=1:n
    for j=i:n
index1(i,j)=L(2,index(i,j));                     %converte matrix proximidade em localidade
    end
end
index2=index1+index1'- eye(size(index1,1)).*diag(index1);
Classif=index2(1:n,2:n);
 
%==========================================NOVO=============================
%===========================================================================
%===========================================================================
 
for i=1:n
    
        for h=1:ln
if L(2,i)==h
numj(i,h)=1;
else 
    numj(i,h)=0;
end
        end
end
 
poploc0 = sum(numj);
 
 
%============================================
% gera indice diversidade agregada por localidade
for i=1:n
    
        for h=1:ln
if L(2,i)==h
num_teste(i,h)=Div2(i,1);
else 
    num_teste(i,h)=0;
end
        end
end
 
Div_per_loc0 = sum(num_teste,1);
% normalizamos de 1 a 100 a diversidade
Div_per_loc = (p5*((Div_per_loc0 - min(Div_per_loc0))/(max(Div_per_loc0)- min(Div_per_loc0))))/100;
 
%=================================================== :)
%=================================================== 
%=================================================== 
 
 
sorteio0=[repelem(1,round((p4*poploc0(1,1)*(1+Div_per_loc(1,1))/(p6*t)))) repelem(2,round((p4*poploc0(1,2)*(1+Div_per_loc(1,2))/(p6*t)))) repelem(3,round((p4*poploc0(1,3)*(1+Div_per_loc(1,3))/(p6*t)))) repelem(4,round((p4*poploc0(1,4)*(1+Div_per_loc(1,4))/(p6*t)))) repelem(5,round((p4*poploc0(1,5)*(1+Div_per_loc(1,5))/(p6*t))))];
sorteio1=repmat(sorteio0,n,1);
 
 
sorteio=[repmat(Classif(:,1),1,ln+round(1)), repmat(Classif(:,2),1,round(ln/2)), repmat(Classif(:,3),1,round(ln/3)), repmat(Classif(:,4),1,round(ln/4)),repmat(Classif(:,5),1,round(ln/5)),randi(ln,n,ln) sorteio1];
 
tet=sum(K,2);
%pond=round(1./sim(1:n,2:n)).*tet; %ponderação da migração por afinidade e conhecimento
 
 
 
for i=1:n
    L1(1,i)=L(1,i);
if rand()>=(1-p1)-(tet(i,1)/100)
    L1(2,i)=sorteio(i,randi(size(sorteio,2)));
else
    L1(2,i)=L(2,i);
end
 
end
%========================================================================
%Dynamic Network=========================================================
%========================================================================
for i=1:n
    for j=1:n
if rand()<=p2  &&  L1(2,i)==L1(2,j);       %4chance conexão dentro da mesma localidade ###############################
    Net1(i,j)=1;
elseif rand()<=p3                 %chance fazer uma conexão aleatória
    Net1(i,randi(n))=1;
elseif rand>=0.99
    Net(i,randi(n))=0;      % CHANCE DESFAZER CONEXÃO
else
    Net1(i,j)=Net(i,j);
end
    end
end
%====================================================================
%=====================PARAMETROS#####################################
%====================================================================
%potencial criativo
for i=1:n
    for j=1:kc
if K1(i,j)>=4
    KO(i,j)=1;
else
    KO(i,j)=0;
end
    end
end
 
for i=1:n
if sum(KO(i,:))>=2
    PC(i,1)=1;
else
    PC(i,1)=0;
end
end
 
 
 
% Análise pMigração
for i=1:n
    
        for h=1:ln
if L1(2,i)==h
numi(i,h)=1;
else
    numi(i,h)=0;
end
        end
end
 
poploc(leo,:)= sum(numi);
%@@@ @@@@@@@@@@@@@@@@@@
criloc(1,:)=sum(PC,2);
criloc1=[criloc;L1(2,:)];
criloc2=criloc1(1,:).*criloc1(2,:);
for h=1:ln
    for i=1:n
if criloc2(1,i)==h
    criloc3(h,i)=1;
else
    criloc3(h,i)=0;
end
    end
end
 POTCRI(leo,1)=sum(PC);           %potencial criativo total
criativelocation(leo,:)=sum(criloc3,2)'; % number of criative per location
percriloc=criativelocation./poploc;    %percentage of creative people per location
div_locat(leo,:) = Div_per_loc0;
Net=Net1;
L=L1;
K=K1;
 
end
 
 
%Gráfico final da rede
%format compact
%format long e
%theta = linspace(0,2*pi,length(Net)+1);
%xy = zeros(length(Net)+1,2);
%x = cos(theta);
%y = sin(theta);
%xy(1:length(Net)+1,1) = x(1:length(Net)+1);
%xy(1:length(Net)+1,2) = y(1:length(Net)+1);
%figure, gplot(Net,xy,'.-');
%set(gcf, 'Color', [1 1 1]);
%axis('equal');
%xlim([-1.1 1.1]);
%ylim([-1.1 1.1]);
%axis off;

