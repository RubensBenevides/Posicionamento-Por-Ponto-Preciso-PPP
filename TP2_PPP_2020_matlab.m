clc
clear
D = 0 ; H = 11; M = 13; S = 45;     % Dia, hora, minuto, segundo da época GPS
GM_T = 3.986008E14;                 % cte. geocêntrica TERRA
GM_S = 1.3271244004193938E20;       % cte. geocêntrica SOL
GM_L = 4.902800066E12;              % cte. geocêntrica LUA
lat  =  (22+(7/60)+(36.7/3600));    % lat.  da  ppte1971: 22° 07' 36.7'' S 
long = -(51+(24/60)+(30.7/3600));   % long. da  ppte1971: 51° 24' 30.7'' W
H_elip = 431.049;                   % altitude elipsoidal em m do pornto
N_geo  = -5.03;                     % ondulalção geoidal em m do ponto
Tr = 18+273;                        % temperatura padrão em Kelvin
Pr = 1013.25;                       % pressão padrão em mbar (760 mmHg)
Hur= 50;                            % umidade relativa padrão (50%)
omega = 7.292115147*(10^-5);        % velocidade angular da Terra em rad/s
c  = 299792458;                     % velocidade da luz em m/s
m1 = 2.5457;                        % coef. ion-free PD e fase na L1
m2 = 1.5457;                        % coef. ion-free fase na L1
m2f= 1.9837;                        % coef. ion-free fase na L2
F  = -4.442807633*10^-10;           % parâmetro do erro relativístico (dtR)
timezone = -3;                      % fuso horário do ponto
timehack = [2012 07 15 11 13 45];   % hora sideral em Greenwich (0 TU)
UA = 149597870700;                  % unidade astronômica em metros
sats = [18;21;22;25;29;30;31];      % PRNs observados na época
SATs = length(sats);
%********************************************************
% VETOR DE PARÂMETROS APROXIMADOS [Xo] E OBSERVAÇÕES [Lb]
%********************************************************
Xo=[3606986.0630;         % Coordenada X0
    -4345293.2446;        % Coordenada Y0
    -2956654.2117;        % Coordenada Z0
    0;                    % Valor inicial do dtr
    0;                    % Ambiguidade do sat. 1 
    0;                    % Ambiguidade do sat. 2
    0;                    % Ambiguidade do sat. 3
    0;                    % Ambiguidade do sat. 4
    0;                    % Ambiguidade do sat. 5
    0;                    % Ambiguidade do sat. 6
    0];                   % Ambiguidade do sat. 7
% Pseudodistâncias da L1
C1=[20673352.930;
    20462608.391;
    23710479.281;
    21504467.148;    
    21432292.453;
    22830292.500;
    23714758.398];
% Pseudodistâncias da L2
P2=[20673358.508;20462613.207;
    23710484.977;21504476.008;    
    21432298.734;22830299.895;
    23714764.621];
% Fases da L1
L1=[108639303.4188;107531967.0188;
    124599350.4087;113006744.4607;
    112627436.5598;119973979.7467;
    124622029.1557];
% Fases da L2
L2=[84653910.02547;83790961.72746;
    97090462.70744;88057229.08846;
    87761627.83746;93486219.20545;
    97107947.98744];
% Combinação ion-free
PD_L0   = (m1*C1)-(m2*P2);
Fase_L0 = (m1*L1)-(m2f*L2);
Lb      = [PD_L0;Fase_L0];
RMSEx = zeros(100,1);
% Pesos das observáveis
P_PD    = eye(7)*inv(3)^2;   % Inverso do quadrado da precisão da PD
P_fase  = eye(7)*inv(0.1)^2; % Inverso do quadrado da precisão da Fase
Peso    = [P_PD    zeros(7); % Peso = matriz diagonal 14x14 com pesos
           zeros(7) P_fase];
           
% ATRASO NOS HARDWARES (TGD) (tabelado no ppte1971.12n)
TGD = [-1.07102096081E-08;
       -1.11758708954E-08;
       -1.67638063431E-08;
        5.58793544769E-09;
       -9.31322574616E-09;
       -4.19095158577E-09;
       -1.30385160446E-08];

% INSTANTE DE RECEPÇÃO DA ONDA (tr)
tr=(D*24+H)*3600+M*60+S;

% TEMPO DE TRANSMISSÃO DE CADA SATÉLITE (tt_GPS)
tt_GPS = tr-(PD_L0./c);

%*********************************************************
% INTERPOLAÇÃO DE XsYsZs E dts EM CADA tt_GPS
%*********************************************************        
XYZsp3 = dlmread('XsYzZs_sp3.csv'); % Ler coordenadas XsYsZs do arquivo excel
epocas_sp3 = [(36900:900:43200)]; % Define 8 epocas variando 900s (4 antes, 4 depois)
dts_clk = dlmread('ts_clk05s.csv'); % Lê tempos do relógio do arquivo excel
epocas_clk_5 = [(40410:5:40445)]; % Define 8 epocas variando 5s (4 antes, 4 depois)
for n=1:SATs
    % Organizar os dados lidos para usar na interpolação
    coorX=XYZsp3(n:SATs:8*SATs,1)*1000; 
    coorY=XYZsp3(n:SATs:8*SATs,2)*1000;
    coorZ=XYZsp3(n:SATs:8*SATs,3)*1000;
    % Interpolação de XsYsZs (IGS.sp3) em cada tt_GPS para cada satélite
    Xs(n,1)=spline(epocas_sp3,coorX,tt_GPS(n));
    Ys(n,1)=spline(epocas_sp3,coorY,tt_GPS(n));
    Zs(n,1)=spline(epocas_sp3,coorZ,tt_GPS(n));
    % Interpolação do dts em cada ttGPS para cada satélite
    dts_inter=dts_clk(n:SATs:8*SATs);
    %dts_GPS(n,1)=spline(epocas_clk_5,dts_inter,tt_GPS(n)); 
end
dts_GPS = dlmread('CLK_exato_11_13_45.csv');
dts_GPS = dts_GPS(1:7)
% *******************************************************
% TROPOSFERA - MODELO DE SAASTAMOINEN
% *******************************************************
H_ort = (H_elip-N_geo)/1000;           % Alt. ortométrica em m
D  = 0.0026*cosd(2*lat)+0.00028*H_ort; % Coeficiente do modelo
P  = Pr*(1-0.0000226*(H_ort))^5.225;   % Pressão calculada em H_ort
T  = Tr-0.0065*H_ort;                  % Temperatura calculada em H_ort
Hu = Hur*exp(-0.0006396*H_ort);        % Umidade relativa calculada em H_ort
e  = 6.108*exp((17.15*T-4684)/(T-38.45))*(Hu/100); % Pressão de vapor
% Interpolação do coeficiente B na H_ort;
vetor_H=[0:0.5:3,4,5];
vetor_B=[1.156;1.079;1.006;0.938;0.874;0.813;0.757;0.654;0.563];
B=spline(vetor_H,vetor_B,H_ort);       % B é um fator de correção
% Elevação de cada PRN
r_est = [3687624.3674;-4620818.6827;-2386880.3805];
R_est = norm(r_est);
for n=1:SATs
     r_sat = [Xs(n);Ys(n);Zs(n)];
     l2    = norm(r_sat);
     r_rel = r_sat - r_est ;
     l3    = norm(r_rel);
     alfa  = acos(-(l2^2-R_est^2-l3^2)/(2*R_est*l3));
     zen(n,1)  = pi-alfa;
     elev(n,1) = rad2deg((pi/2)-zen(n)); 
end
% Modelo de Saastamoinen sem dR (como Bernise 5.0)
Trop = 0.002277*inv(1+D)*sec(zen).*(P+((1255/T)+0.05)*e-B*tan(zen).^2);

%**********************************************************
% PESOS DAS OBSERVÁVEIS PROPORCIONAIS AO SEN² DA ELEVAÇÃO
%*******************************************************
%P_PD   = inv(diag(sind(elev).^2)).*inv(3^2);
%P_fase = inv(diag(sind(elev).^2)).*inv(0.1^2); 
%Peso   = [P_PD,zeros(7);zeros(7),P_fase];

%*******************************************************
% MODELAGEM DA MARÉ TERRESTRE
%*******************************************************
% Efemérides da Lua (Ascensão Reta e Declinação)
RA_DEC = dlmread('horizons_luna_matlab.txt'," ",[67,6,74,11]);
RA  = RA_DEC(:,1) + RA_DEC(:,2)/60 + RA_DEC(:,3)/3600 ;
DEC = RA_DEC(:,4) + RA_DEC(:,5)/60 + RA_DEC(:,6)/3600; 
RA_lua  = spline(8:15,RA,tr/3600)*15; DEC_lua = spline(8:15,DEC,tr/3600);
% Efemérides do Sol (Ascensão Reta e Declinação)
dlmread('horizons_sun_matlab.txt'," ",[66,6,73,11]);
RA  = RA_DEC(:,1) + RA_DEC(:,2)/60 + RA_DEC(:,3)/3600; 
DEC = RA_DEC(:,4) + RA_DEC(:,5)/60 + RA_DEC(:,6)/3600;
RA_sol = spline(8:15,RA,tr/3600)*15; DEC_sol = spline(8:15,DEC,tr/3600);
% Hora sideral em Greenwich (0 TU) na época 15/07/2012 11:13:45
GAST = sideraltime(timehack,timezone); % Saída em horas decimais
GAST = GAST*15; Hsol = GAST-RA_sol; Hlua = GAST-RA_lua;
% SOL: coordenadas, vetor, versores e magnitude
Xsol = UA*cosd(DEC_sol)*cosd(Hsol);
Ysol = UA*cosd(DEC_sol)*sind(Hsol);
Zsol = UA*sind(DEC_sol); Rsol = UA; rsol = [Xsol,Ysol,Zsol]/Rsol;
% Calcular distância média terra-lua
a = 384400; e = 0.05490; b = a*sqrt(1-e^2); Rlua = sqrt(a*b)*1000;
% LUA: coordenadas, vetor, versores e magnitude
Xlua = Rlua*cosd(DEC_lua)*cosd(Hlua);
Ylua = Rlua*cosd(DEC_lua)*sind(Hlua);
Zlua = Rlua*sind(DEC_lua); rlua = [Xlua,Ylua,Zlua]/Rlua;
% Coordenadas da estação
r_est = [3687624.3674;-4620818.6827;-2386880.3805];
R_est = norm(r_est); rest = r_est/R_est;
% Números de Love
h0 = 0.6078 ; l0 = 0.0847; h2 = -0.0006; l2 = 0.0002;
h2 = h0+h2*(3*(sind(lat)^2)-1)/2; l2 = l0+l2*(3*(sind(lat)^2)-1)/2;
% Correção SOL
for n=1:3
    AA = (GM_S/GM_T)*(R_est^4/Rsol^3);
    BB = (3*l2*rsol(n)*rest(n)*rsol(n));
    CC = (3*(h2/2-l2))*(rsol(n)*rest(n))^2;
    DD = h2/2;
    EE = -0.025*sind(lat)*cosd(lat)*sind(GAST+long)*rest(n);
    MareSol(n,1) = AA*(BB+(CC-DD)*rest(n))+EE;
end
% Correção LUA
for n=1:3
    AA = (GM_L/GM_T)*(R_est^4/Rlua^3);
    BB = (3*l2*rlua(n)*rest(n)*rlua(n));
    CC = (3*(h2/2-l2))*(rsol(n)*rest(n))^2;
    DD = h2/2;
    EE = -0.025*sind(lat)*cosd(lat)*sind(GAST+long)*rest(n);
    MareLua(n,1) = AA*(BB+(CC-DD)*rest(n))+EE;
end
mare=MareLua+MareSol;

%*******************************************************
% ERRO RELATIVÍSTICO
%*******************************************************
% 1 - Carregar parâmetros orbitais de cada PRN
a_e_toe_vn_m0 = dlmread('a_e_toe_vn_m0.csv');
a_e_toe_vn_m0;
sqr_a = a_e_toe_vn_m0(:,1); % a^1/2
excen = a_e_toe_vn_m0(:,2); % exentricidade
toe   = a_e_toe_vn_m0(:,3); % toe = 360000
vn    = a_e_toe_vn_m0(:,4); % variação do mov. ang. médio
M0    = a_e_toe_vn_m0(:,5); % anomalia média inicial
% 2 - Propagar os parâmetros para tt_GPS
for n=1:SATs 
    nk(n,1)=(((GM_T/(sqr_a(n))^6))^(1/2))+vn(n); % mov. ang. no toe
    tk(n,1)=tt_GPS(n)-toe(n);                
    Mk(n,1)=M0(n)+(nk(n)*tk(n));              
    % Anomalia média de cada PRN pelo método de Newton
    Ek0(n,1)=Mk(n,1)+excen(n,1)*Mk(n,1);
    Ek1(n,1)=Mk(n,1)+excen(n,1)*Ek0(n,1);
    while abs(Ek1(n,1)-Ek0(n,1))>1e-15; 
        Ek0(n,1)=Mk(n,1)+excen(n,1)*Ek1(n,1);
        Ek1(n,1)=Mk(n,1)+excen(n,1)*Ek0(n,1);
    end
    % 3 - Salvar correção
    dtR(n,1)=F*excen(n)*sqr_a(n)*sin(Ek1(n));
end
% CORREÇÃO DE RELATIVIDADE E ATRASO DE HARDWARE NOS dts (dtsC)
dtsC_L1 = dts_GPS+dtR-TGD;
dtsC_L2 = dts_GPS+dtR-(TGD.*1.647); 
% Combinação do dts em L0
dtsC_L0 = (dtsC_L1.*m1)+(dtsC_L2.*m2);

% ************************************************************
% AJUSTAMENTO + INTERPOLAÇÕES
% ************************************************************
dtr=0; itera = 0; fim = false;
while (~fim)
    itera = itera+1;
    % tempo de voo da onda (tal_L0) combinação na PD
    tal_L0 = (PD_L0./c)-dtr+dtsC_L1;
    % tempo de transmissão da cada GPS refinado
    ttGPS_L0_ref = tr-dtr+dtsC_L1-tal_L0;
    % Interpolar XsYsZs e dts novamente em ttGPS_L0_ref
    for n=1:SATs
        coorX=XYZsp3(n:SATs:8*SATs,1)'*1000;
        coorY=XYZsp3(n:SATs:8*SATs,2)'*1000; 
        coorZ=XYZsp3(n:SATs:8*SATs,3)'*1000;
        Xs_ref(n,1)=spline(epocas_sp3,coorX,ttGPS_L0_ref(n));
        Ys_ref(n,1)=spline(epocas_sp3,coorY,ttGPS_L0_ref(n));
        Zs_ref(n,1)=spline(epocas_sp3,coorZ,ttGPS_L0_ref(n));
        % Interpolação do dts
        dts_inter=dts_clk(n:SATs:8*SATs,1)';
        %dtsC_L0(n,1) = spline(epocas_clk_5,dts_inter,ttGPS_L0_ref(n));      
    end
    % CORREÇÂO SAGNAC
    alfa = tal_L0.*omega;
    Xs_rot = cos(alfa).*Xs_ref + sin(alfa).*Ys_ref;
    Ys_rot = cos(alfa).*Ys_ref - sin(alfa).*Xs_ref;
    Zs_rot = Zs_ref; % Zs não sofre correção Sagnac
    % Distâncias (parte geométrica)
    psr = sqrt((Xs_rot-Xo(1)+mare(1)).^2+(Ys_rot-Xo(2)+mare(2)).^2+(Zs_rot-Xo(3)+mare(3)).^2);
    % Calcular L0 = f(X0)
    PD =   psr + Xo(4) - dtsC_L1*c + Trop;
    FASE = psr + Xo(4) - dtsC_L1*c + Xo(5:end)*0.19 + Trop;
    L0 = [PD;FASE];
    L  = L0 - Lb;
    % Componentes da matriz A
    A = [(Xo(1)-Xs_rot(1))/psr(1) (Xo(2)-Ys_rot(1))/psr(1) (Xo(3)-Zs_rot(1))/psr(1) 1 0 0 0 0 0 0 0; 
         (Xo(1)-Xs_rot(2))/psr(2) (Xo(2)-Ys_rot(2))/psr(2) (Xo(3)-Zs_rot(2))/psr(2) 1 0 0 0 0 0 0 0;
         (Xo(1)-Xs_rot(3))/psr(3) (Xo(2)-Ys_rot(3))/psr(3) (Xo(3)-Zs_rot(3))/psr(3) 1 0 0 0 0 0 0 0;
         (Xo(1)-Xs_rot(4))/psr(4) (Xo(2)-Ys_rot(4))/psr(4) (Xo(3)-Zs_rot(4))/psr(4) 1 0 0 0 0 0 0 0;
         (Xo(1)-Xs_rot(5))/psr(5) (Xo(2)-Ys_rot(5))/psr(5) (Xo(3)-Zs_rot(5))/psr(5) 1 0 0 0 0 0 0 0;
         (Xo(1)-Xs_rot(6))/psr(6) (Xo(2)-Ys_rot(6))/psr(6) (Xo(3)-Zs_rot(6))/psr(6) 1 0 0 0 0 0 0 0;
         (Xo(1)-Xs_rot(7))/psr(7) (Xo(2)-Ys_rot(7))/psr(7) (Xo(3)-Zs_rot(7))/psr(7) 1 0 0 0 0 0 0 0;
         (Xo(1)-Xs_rot(1))/psr(1) (Xo(2)-Ys_rot(1))/psr(1) (Xo(3)-Zs_rot(1))/psr(1) 1 0.19 0 0 0 0 0 0;
         (Xo(1)-Xs_rot(2))/psr(2) (Xo(2)-Ys_rot(2))/psr(2) (Xo(3)-Zs_rot(2))/psr(2) 1 0 0.19 0 0 0 0 0;
         (Xo(1)-Xs_rot(3))/psr(3) (Xo(2)-Ys_rot(3))/psr(3) (Xo(3)-Zs_rot(3))/psr(3) 1 0 0 0.19 0 0 0 0;
         (Xo(1)-Xs_rot(4))/psr(4) (Xo(2)-Ys_rot(4))/psr(4) (Xo(3)-Zs_rot(4))/psr(4) 1 0 0 0 0.19 0 0 0;
         (Xo(1)-Xs_rot(5))/psr(5) (Xo(2)-Ys_rot(5))/psr(5) (Xo(3)-Zs_rot(5))/psr(5) 1 0 0 0 0 0.19 0 0;
         (Xo(1)-Xs_rot(6))/psr(6) (Xo(2)-Ys_rot(6))/psr(6) (Xo(3)-Zs_rot(6))/psr(6) 1 0 0 0 0 0 0.19 0;
         (Xo(1)-Xs_rot(7))/psr(7) (Xo(2)-Ys_rot(7))/psr(7) (Xo(3)-Zs_rot(7))/psr(7) 1 0 0 0 0 0 0 0.19];
    N   = A'*Peso*A;
    U   = A'*Peso*L;
    Xc  = -inv(N)*U; % Correção dos parâmetros
    Xa  = Xo+Xc;
    Xo  = Xa;
    dtr = Xo(4)/c;
    erro_max = max(abs(Xc));
    fim = (max(abs(Xc))<1.0E-4) || (itera > 10);
end
itera
gl=14-11; % graus de liberdadde
V=A*Xc+L  % resíduos das observações (Pseudodistâncias e Fases)
La=Lb+V;  % Observações ajustadas    
var_posteriori=V'*Peso*V/gl  % Variância a posteriori
mvcXa=var_posteriori*inv(N); % MVC dos parâmetros (X,Y,Z,dt,fase1,...,faseN)

% Teste Qui-quadrado
qui_quadrado = gl*var_posteriori;
qui_quadrado_tabelado_95 = 7.815;
if qui_quadrado < qui_quadrado_tabelado_95
   fprintf('Passou no teste qui-quadrado (valor calculado = %.2f < 7.82 valor tabelado)\n', qui_quadrado);
else
   fprintf('Falhou no teste qui-quadrado (valor calculado = %.2f > 7.82 valor tabelado)\n', qui_quadrado);
end
% Métricas de Erro
Erro_x = Xa(1)-(3687624.3674)
Erro_y = Xa(2)-(-4620818.6827)
Erro_z = Xa(3)-(-2386880.3805)
Erro_xyz = norm([Erro_x,Erro_y,Erro_z])
EQM = sqrt((Erro_xyz^2+var_posteriori^2)/2) % média quadrática
