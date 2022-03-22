//------------------------------------------------------------------------------
// DINAMICA ATOMICA E MOLECULAR COM TEMPERATURA - DINATEMP - VERSAO 153 (17-06-2014)
// PROGRAMA PARA SIMULAÇÃO DE SISTEMAS ATOMICOS VIA DINAMICA MOLECULAR, EM AMBIENTE
// COM TEMPERATURA CONTROLADA (OPCIONAL) E DENTRO DE CAIXA 3D (OPCIONAL)
// RODANDO EM PROCESSAMENTO PARALELO. ENSEMBLES NVE E NVT 
// COPYRIGHT: PROF. GILMAR G. FERREIRA / GRAD. ROCHELLE DINIZ - 2014
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// CARREGANDO SUBROTINAS E DEFININDO VARIAVEIS GLOBAIS
//------------------------------------------------------------------------------
clear; clc; stacksize('max') ;mode(0); ieee(1); funcprot(0); tic; 
global('np','m','kb','Q','Temp_equil','Temp_equil_n','L1','L2','eps','Time','Temperatura')

exec epot.sce; exec depot.sce; exec epp.sce; exec depp.sce; 
exec cometa.sce; exec histograma.sce;
exec ODE_V150_ROCHELE.sce;
//exec ODE_V153_PARALLEL4.sce; exec DINA3.sce; exec DVDR1.sce;
disp('INICIANDO DINATEMP (VERSAO Nº 155_ROCHELE)');

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
disp('CARREGANDO PARAMETROS DA SIMULAÇÃO...')
//------------------------------------------------------------------------------
exec PARAMETROS.sce;
np=2;           // Número de Partículas
m=1;              // MassParanaíbaa (em UNidades Reduzidas)
Q=1;            // Fator de Transferencia de energia do/para reservatorio. Cuidado!
L1=0; L2=10; // Localização das paredes da caixa, em cada eixo (x,y,z)
Temp_equil =200;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
disp('GERANDO CONDIÇÕES INICIAIS...')
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// POSICOES [   x   y   z   ]
//MP = rand(np,3)*L2;
//13
  
MP =[2.7831211    2.33857      1.302439   
    1.7182633    2.7451285    1.4133933];
//    2.5961411    3.18523      2.0387368  
//    3.2600231    2.3180241    2.3419101  
 //   1.9106687    1.6229539    1.3166708  
//    2.8579645    1.3674849    1.8872663  
 //   1.6032056    1.9754951    3.0920316  
//   1.5552503    2.959385     2.5274605  
 //   1.1509337    2.0066901    2.079874   
 //   2.5040957    2.7017412    3.1060571  
 //   2.687329     1.5721567    2.9892673  
//    1.8361666    1.1436161    2.3281971  
//    2.2086547    2.1631869    2.2040399
 
//------------------------------------------------------------------------------
// MOMENTOS [   px   py   pz   ]
MM = zeros(np,3);
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// VETOR DE ENTRADA
//------------------------------------------------------------------------------
//    x        y        z       px       py       pz
u=[MP(:,1); MP(:,2); MP(:,3); MM(:,1); MM(:,2); MM(:,3)];
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
disp('INTEGRANDO EQUAÇÕES DE MOVIMENTO...')
//------------------------------------------------------------------------------
PASSO=0.0001; TFINAL=250;         //  unidades reduzidas
t=[0:PASSO:TFINAL];
u=ode(u,0,t,ODE_V150_ROCHELE);
//u=ode(u,0,t,ODE_V153_PARALLEL4); 
//proximos trabalhos: 
// Rodar com Caixa; Rodar em paralelo...
//------------------------------------------------------------------------------

toc; Tempo_analise_min=ans/60

//------------------------------------------------------------------------------
disp('GERANDO CONDIÇÕES FINAIS...')
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// CALCULANDO "TEMPO REAL" DE SIMULAÇÃO (Ps)
//------------------------------------------------------------------------------
TS =  PASSO*TFINAL; // "Tempo de Simulação" do Argonio em Unid. Red. Tempo (1 URT = 2.17X10-12 s) 
TR = TS*Time;                   // "Tempo Real" em pico-segundos (x10-12 s) - ARGONIO
TR = 0:TR/(length(t)-1):TR;     // Criando vetor de "Tempo Real"
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// POSICOES [   x   y   z   ]
//------------------------------------------------------------------------------
x = u(1:np,:); y = u(np+1:2*np,:);z = u(2*np+1:3*np,:);
//save('COORDENADAS.dat',np,u,x,y,z,L1,L2,TR,Tempo_analise_min);
//--------------------------t----------------------------------------------------

//------------------------------------------------------------------------------
// MOMENTOS [  px  py  pz   ]
//------------------------------------------------------------------------------
px = u(3*np+1:4*np,:);py = u(4*np+1:5*np,:);pz = u(5*np+1:6*np,:);
MM = [px; py; pz];
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// ENERGIA CINETICA  (em "ϵ" )
//------------------------------------------------------------------------------
ECIN = (sum(MM.^2,1))./(2*m);
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// TEMPERATURA (ATENÇÃO: VERIFICAR SE ESTA USANDO UNIDADES REDUZIDAS OU ATOMICAS)
//------------------------------------------------------------------------------
// EM UNIDADES REDUZIDAS
dim = 3; // n de dimensões
temp = ((ECIN.*2)./(dim.*np.*kb)).*eps;    // vindo em "ϵ" saindo em "K"  (1 eps = 1.65e−21 J para AR)
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// ENERGIA POTENCIAL PARA N PARTICULAS (EPN)
//------------------------------------------------------------------------------
EPN=0;
for i=1:np
    for j=i+1:np
        r=sqrt((x(i,:)-x(j,:)).^2 + (y(i,:)-y(j,:)).^2 + (z(i,:)-z(j,:)).^2);
        EPN = EPN + epot(r);
    end
end
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// ENERGIA POTENCIAL COM PAREDES REPULSIVAS (EPPR)
//------------------------------------------------------------------------------
epx = epp(abs(x-L1)) + epp(abs(x-L2));
epy = epp(abs(y-L1)) + epp(abs(y-L2));
epz = epp(abs(z-L1)) + epp(abs(z-L2));
EPPR = sum(epx,1) + sum(epy,1) + sum(epz,1);
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// ENERGIA POTENCIAL TOTAL  (em "ϵ" )
//------------------------------------------------------------------------------
EPT = EPN + EPPR;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// ENERGIA ETOTAL (CINETICA + POTENCIAL)
//------------------------------------------------------------------------------
ETOTAL = ECIN + EPT ; // (vindo em "ϵ" )
ETOTAL = ETOTAL*eps*6.241506363e+18; // vindo em "ϵ"saindo em "ev"
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// CALCULANDO FUNCAO DISTRIBUIÇÃO RADIAL (Gr)
//------------------------------------------------------------------------------
//for k=1:length(u(1,:)) // tempo de simulação
   // for i=1:np
       // for j=1:np
          //  if i~=j // evitando duplo calculo para cada par de partículas
            //    r(j)=sqrt((x(i,k)-x(j,k))^2+(y(i,k)-y(j,k))^2+(z(i,k)-z(j,k))^2);
            //end
       // end
       // [outx,outy]=histograma([0:0.01:10],r); // histograma para particula "i"
       // ACUMULANDO(i,:)=outy';// registrando matriz x de histogramas das particlulas

    //end
//MEDIA(k,:) = mean(ACUMULANDO,1); // calculando media da matriz x em tempo "k"
//end
//close;
//Gr = mean(MEDIA,1); // calculando media do histograma no tempo total da simluação

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
disp('PLOTANDO RESULTADOS DA SIMULAÇÃO...')
//-----------------------------------------------------------------------------
// energia total (em "AU" ou "ϵ" ou "eV" )
f2 = scf(2);
plot(TR,ETOTAL,'r')
//plot(ETOTAL,'r')
title(['VARIAÇÃO ENERGIA TOTAL']);
xlabel("Tempo (Ps)", "fontsize", 2);
ylabel("Eexec PARAMETROS.sce;
np=2;           // Número de Partículas
m=1;              // MassParanaíbaa (em UNidades Reduzidas)
Q=1;            // Fator de Transferencia de energia do/para reservatorio. Cuidado!
L1=0; L2=10; // Localização das paredes da caixa, em cada eixo (x,y,z)
Temp_equil =50;nergia (eV)", "fontsize", 2);


// temperatura (em K)
f1 = scf(1);
plot(TR,temp,'r')
title(['VARIAÇÃO TEMPERATURA']);
xlabel("Tempo (Ps)", "fontsize", 2);
ylabel("Temperatura (K)", "fontsize", 2);


// função radial
//f3 = scf(3);
//title(['DISTRIBUIÇÃO RADIAL SISTEMA AR13']);
//plot(outx,Gr','g');
//plot(outx,MEDIA(1,:),'b');
//plot(outx,MEDIA(length(MEDIA(:,1)),:),'K');
//xlabel("Distância (x 3.4 Å)", "fontsize", 2);
//ylabel("g(r)", "fontsize", 2);

// filmando trajetórias
load('COORDENADAS.dat')
f4 = scf(4);
cometa(x(:,1:length(x(1,:)))',y(:,1:length(y(1,:)))',z(:,1:length(z(1,:)))');
a=gca(); a.isoview="on"; // isoview mode
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// REGISTRANDO ULTIMAS COORDENADAS PARA POSTERIOR SIMULAÇÃO DE CLUSTERS...
//------------------------------------------------------------------------------
ultimax= x(:,length(x(1,:)));ultimay= y(:,length(y(1,:)));ultimaz= z(:,length(z(1,:)));
MU = [ultimax ultimay ultimaz];
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
disp('SIMULAÇÃO CONCLUIDA COM SUCESSO.'); 
//------------------------------------------------------------------------------
