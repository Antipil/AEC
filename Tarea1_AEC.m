%% Primera Tarea de AEC - 22.05.23
%  Nataly Antipil %
%% Detalles a considerar antes de:
% Se eligieron las variables temperatura del aire a 2 m de altura desde la superficie, y la tasa de
% precipitación. Son datos mensuales, entre enero de 1948 y marzo de 2023.

% Por alguna razón el archivo de la tasa de precipitación, contiene valores negativos (pequeños).
% Dejarlos en cero antes de realizar cualquier análisis.

%No cometer el error de promediar o sumar los meses 1, 2 y 12. En realidad es 12, 13 y 14, etc.
%Perderán el invierno boreal de 1948 y ganarán el invierno boreal de 2023.

%% Ploteo de todos los datos

lon = ncread('air.2m.mon.mean.nc','lon');
lat = ncread('air.2m.mon.mean.nc','lat');
YY=length(lat);
XX=length(lon);
aux=ncread('air.2m.mon.mean.nc','air');
contourf(lon,lat,reshape(squeeze(aux(:,:,1)),XX,YY)') % Solo enero
%% Cargar archivos y los entiendo mejor
 ncdisp('air.2m.mon.mean.nc') % variables:lon(192),lat(94),air(192,94,903),Unidad:K

% Extraer la región del continente euro-asiático: 70°N a 0°; 11°E a 140°E
 %clear all
lon = ncread('air.2m.mon.mean.nc','lon',[7],[69],[1]);
lat = ncread('air.2m.mon.mean.nc','lat',[11],[47],[1]);
YY=length(lat);
XX=length(lon);
air=ncread('air.2m.mon.mean.nc','air',[7 11 1],[69 47 903],[1 1 1]);
%contourf(lon,lat,reshape(squeeze(air(:,:,1)),XX,YY)') 

promair=nanmean(air(:,:,1:1:903),3); % Promedio de la temp del aire en todo ese tiempo
% me explayo: tiene dimencion (lon x lat) porque promedia la tercera
% dimencion, es decir, el tiempo; es como si promediara la tercera
% dimencion de cada grilla y la dejara en un solo mapa de lon por lat.

%   Mapa promedio de temperatura del aire
[lati,long] = meshgrid(double(lat),double(lon)); %grilla para graficar
figure()
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,promair,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('Mapa promedio de temperatura del aire')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off

%% 1. En el caso del continente euro-asiático, comparar las distribuciones de probabilidad de 
% cada grilla, entre el invierno boreal (DEF) y el verano boreal (JJA).

clear all
clc
% EXTRAIGO LA REGIÓN del continente euro-asiático: 70°N a 0°; 11°E a 140°E
lon = ncread('air.2m.mon.mean.nc','lon',[7],[69],[1]);
lat = ncread('air.2m.mon.mean.nc','lat',[11],[47],[1]);
YY=length(lat);
XX=length(lon);
air=ncread('air.2m.mon.mean.nc','air',[7 11 1],[69 47 903],[1 1 1]); 
%promair=nanmean(air(:,:,1:1:903),3); % Promedio de la temp del aire en todo ese tiempo

puntero=find(isnan(reshape(squeeze(air(:,:,1)),XX*YY,1))==0); % guarda la pocición
        % de los puntos en el mapa que tienen datos (es decir ve si hay o no hay nan).

% Extrae datos de 'air' sin los continentes
% transformar la matriz de 3 dimensiones a 2 dimensiones
for i=1:903 % Para todo el tiempo
    auy=reshape(squeeze(air(:,:,i)),XX*YY,1); % 
    aaa(:,i)=auy(puntero); % Es la iteracion de todo los datos
end
aaa=aaa-273.15; % a celcius

% SELECCIONO SOLO LOS MESES DESEADOS
% Invierno HN
a1=aaa(:,12:12:end-1); % D % o diciembre es del dato 3?
a2=aaa(:,13:12:end);   % E
a3=aaa(:,14:12:end);   % F
tsairdef=(a1+a2+a3)/3;
[M,N]=size(tsairdef);
% su verano (JJA)
a4=aaa(:,6:12:end);   % J =9(?)
a5=aaa(:,7:12:end);    % J
a6=aaa(:,8:12:end);    % A
tsairjja=(a4+a5+a6)/3;
[M,N]=size(tsairjja);

[lati,long] = meshgrid(double(lat),double(lon)); %grilla para graficar

%%%%%%%%%%%%%%               Promedio          %%%%%%%%%%%%%%%%%%%%
promdef=reshape(nanmean(tsairdef(:,:),2),length(lon),length(lat));
promjja=reshape(nanmean(tsairjja(:,:),2),length(lon),length(lat));
% ojo que nanmean(tsairdef(:,:),2) es el promedio de tsairdef(3243x75)
% temporalmente tal que queda (3243x1) que luego se reorganiza (reshape) de
% una columna a un mapa lonxlat.


%   Mapa promedio de temperatura del aire
figure()
subplot(1,2,1)
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,promdef,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('Promedio de la temperatura del aire en DEF de 1948-2022')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off
subplot(1,2,2)
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,promjja,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('Promedio de la temperatura del aire en JJA de 1948-2022')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off

%%%%%%%%%%%%%%%%%%%5 Promedio significativo
% pormedios
promair_def=nanmean(tsairdef(:,:),2) % Promedio de la temp aire def - tx

% p(1) es beta y p(2) es alfa
p=polyfit((1:length(promair_def))',promair_def',1)
Promdef_lin=polyval(p,(1:length(promair_def))'); % componente lineal

% se extrae tendencia lineal
a_promdef=promair_def-Promdef_lin;

% tendencia lineal
Bobs=p(1)*10; % C/decada

% significancia
for m=1:1000
    ser=remuestreo(promair_def);
    p=polyfit((1:length(promair_def))',ser',1);
    Bmon(m)=p(1)*10;
end
[Bobs prctile(Bmon,2.5) prctile(Bmon,97.5)] % con alfa=5%
[Bobs prctile(Bmon,5) prctile(Bmon,95)] % con alfa=10%



% significancia estadistica
for j=1:100
    for i=1:M
        remuest_prom(i,j)=(nanmean(remuestreo(a_promdef(i,:)),2));  
    end
end

for i=1:length(puntero)
    if promair_def(i)<=prctile(remuest_prom(i,:),2.5) | promair_def(i)>=prctile(remuest_prom(i,:),97.5)
        def_sig(i)=promair_def(i);
    else
        def_sig(i)=NaN;
    end
end

%figura
[lati,long] = meshgrid(double(lat),double(lon)); %grilla para graficar
def_sigMap=reshape(def_sig',length(lon),length(lat));

figure()
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,def_sigMap,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('correlacion: temp. aire - Precipitacion, MAM, 1948-2022')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off


%%%%%%%%%%%             Desviacion estandar        %%%%%%%%%%%%%%%%%%
stddef=reshape(nanstd(tsairdef(:,:),0,2),length(lon),length(lat));
stdjja=reshape(nanstd(tsairjja(:,:),0,2),length(lon),length(lat));
figure()
subplot(1,2,1)
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,stddef,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('Desviación estándar de la temperatura en DEF de 1948-2022')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off
subplot(1,2,2)
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,stdjja,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('Desviación estándar de la temperatura en JJA de 1948-2022')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off

%%%%%%%%%%%%%%               Skewness          %%%%%%%%%%%%%%%%%%%%
skwdef=reshape(skewness(tsairdef(:,:),0,2),length(lon),length(lat));
skwjja=reshape(skewness(tsairjja(:,:),0,2),length(lon),length(lat));
figure()
subplot(1,2,1)
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,skwdef,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('Skewness de la temperatura del aire en DEF de 1948-2022')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off
subplot(1,2,2)
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,skwjja,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('Skewness de la temperatura del aire en JJA de 1948-2022')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off


%% 1.2 Obtener significancia estadística del estadístico que decidan utilizar para la comparación.


for j=1:100
    for i=1:M
        remCorr_AirPrt(i,j)=corr(a_air(i,:)',remuestreo(a_prate(i,:))');  
    end
end


for i=1:length(puntero)
    if r_PrateAir(i)<=prctile(remCorr_AirPrt(i,:),2.5) | r_PrateAir(i)>=prctile(remCorr_AirPrt(i,:),97.5)
        r_sig(i)=r_PrateAir(i);
    else
        r_sig(i)=NaN;
    end
end

%figura
[lati,long] = meshgrid(double(lat),double(lon)); %grilla para graficar
r_sigMap=reshape(r_sig',length(lon),length(lat));
figure()
subplot(1,2,1)
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,skwdef,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('Skewness Significativo de temperatura del aire en DEF')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off
subplot(1,2,2)
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,skwjja,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('Skewness Significativo de temperatura del aire en JJA')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off

%% 2. En el caso del continente euro-asiático, determinar la tendencia lineal en el periodo de 
% análisis, para el invierno boreal (DEF) y el verano boreal (JJA). 

[M,N,T]=size(air); % dimensiones del tensor.

% Reshape del tensor en una matriz bidimensional (M x N, T)
%reshaped_tensor = reshape(tensor, [], size(tensor, 3));

% Crear un vector de tiempo o índices de tiempo (puede ser meses, años, etc.)
tiempo = 1:size(air, 3);

% Aplicar la regresión lineal para cada punto en la matriz
tendencia = zeros(size(aaa, 1), 2);  % Vector para almacenar las pendientes y las intersecciones
for i = 1:size(aaa, 1)
    % Obtener los coeficientes de la regresión lineal
    coeficientes = polyfit(tiempo, aaa(i, :), 1);
    
    % Almacenar la pendiente y la intersección en el vector 'tendencia'
    tendencia(i, :) = coeficientes;
end

% Reshape de nuevo los resultados a la forma original del tensor
tendencia_reshaped = reshape(tendencia(:, 1), M, N);

figure()
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,tendencia_reshaped,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('Tendencia lineal anual de la temperatura del aire')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off

%% tendencia lineal en el periodo de análisis, para el invierno boreal (DEF) y el verano boreal (JJA). 
[M,N,T]=size(air); % dimensiones del tensor.

%%%%%% DEF
% Crear un vector de tiempo o índices de tiempo (puede ser meses, años, etc.)
tiempo = 1:size(tsairdef, 2);
% Aplicar la regresión lineal para cada punto en la matriz
pdef = zeros(size(tsairdef, 1), 2);  % Vector para almacenar las pendientes y las intersecciones
for i = 1:size(tsairdef, 1)
    % Obtener los coeficientes de la regresión lineal
    coef_def = polyfit(tiempo, tsairdef(i, :), 1);
    
    % Almacenar la pendiente y la intersección en el vector 'tendencia'
    pdef(i, :) = coef_def;
end

% Reshape de nuevo los resultados a la forma original del tensor
p_def = reshape(pdef(:, 1), M, N);

%%%%% JJA
tiempo = 1:size(tsairjja, 2);
% Aplicar la regresión lineal para cada punto en la matriz
pjja = zeros(size(tsairjja, 1), 2);  % Vector para almacenar las pendientes y las intersecciones
for i = 1:size(tsairjja, 1)
    % Obtener los coeficientes de la regresión lineal
    coef_jja = polyfit(tiempo, tsairjja(i, :), 1);
    
    % Almacenar la pendiente y la intersección en el vector 'tendencia'
    pjja(i, :) = coef_jja;
end
% Reshape de nuevo los resultados a la forma original del tensor
p_jja = reshape(pjja(:, 1), M, N);

figure()
subplot(1,2,1)
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,p_def,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('Tendencia lineal de la temperatura del aire en DEF')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off
subplot(1,2,2)
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,p_jja,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('Tendencia lineal de la temperatura del aire en JJA')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off

%% 2.2 Obtener significancia 
% estadística del estadístico que decidan utilizar para obtener la tendencia lineal temporal.

for m=1:100
 for i=1:length(tsairdef)
    ser=remuestreo(tsairdef(i,:));
    p=polyfit((1:length(tsairdef))',ser',1);
    Bmon(m)=p(1)*10;
 end
end

[Bobs prctile(Bmon,2.5) prctile(Bmon,97.5)] % con alfa=5%




for j=1:1000
    for i=1:M
        remCorr_AirPrt(i,j)=corr(a_air(i,:)',remuestreo(a_prate(i,:))');  
    end
end
%%%%%

tiempo = 1:size(tsairdef, 2);
p_remdef= zeros(size(tsairdef, 1), 2);
for j=1:100;
   for i = 1:size(tsairdef, 1);
      % ser(i,:)=remuestreo(tsairdef(i,:))
      coef_rem= polyfit(tiempo',remuestreo(tsairdef(i,:))', 1); % Obtener los coeficientes de la regresión lineal
      p_remdef(i, j) =coef_rem; % Almacenar la pendiente y la intersección en el vector 'tendencia'
   end
end



% Reshape de nuevo los resultados a la forma original del tensor
p_def = reshape(pdef(:, 1), M, N);

for i=1:length(puntero)
    if pdef(i, 1)<=prctile(p_remdef(i, :),2.5) | pdef(i, 1)>=prctile(p_remdef(i, :),97.5) % alfa 0.5 %??????? y estos prctile?
        p_sig(i)=r_PrateAir(i);
    else
        p_sig(i)=NaN;
    end
end


%% 3. En el caso de Sudamérica tropical-subtropical, determinar la correlación punto a punto entre 
% la temperatura a 2 m y la precipitación, durante el otoño austral (MAM). 
% Obtener significancia estadística.

% next
 ncdisp('prate.sfc.mon.mean.nc') % variables:lon(192),lat(94),prate(192,94,903),Unidad:kg m-2 s-1
% Extraer la región tropical-subtropical de Sudamérica: 14°N - -30°S; 279° (81°W) - 330° (30°W).

%%%%%%%
clear all
clc
% EXTRAIGO LA REGIÓN de la región tropical-subtropical de Sudamérica: 14°N - -30°S; 279° (81°W) - 330° (30°W)
lon = ncread('air.2m.mon.mean.nc','lon',[149],[28],[1]);
lat = ncread('air.2m.mon.mean.nc','lat',[40],[44],[1]);
YY=length(lat);
XX=length(lon);
air=ncread('air.2m.mon.mean.nc','air',[149 40 1],[28 44 903],[1 1 1]);
prate=ncread('prate.sfc.mon.mean.nc','prate',[149 40 1],[28 44 903],[1 1 1]);
prate(prate < 0) = 0; % Igualar todos los valores negativos a cero
%contourf(lon,lat,reshape(squeeze(air(:,:,1)),XX,YY)')
%promair=nanmean(air(:,:,1:1:903),3); % Promedio de la temp del aire en todo ese tiempo

puntero=find(isnan(reshape(squeeze(air(:,:,1)),XX*YY,1))==0); % donde no hay nan.

%%%%%% Extrae datos de 'air' sin los continentes
for i=1:903 % Para todo el tiempo
    auy=reshape(squeeze(air(:,:,i)),XX*YY,1); % 
    aaa(:,i)=auy(puntero); % Es la iteracion de todo los datos
end
aaa=aaa-273.15; % a celcius

% seleciono solo Nuestro Otoño
a1=aaa(:,3:12:end-1);   % M
a2=aaa(:,4:12:end);   % A
a3=aaa(:,5:12:end);   % M
airmam=(a1+a2+a3)/3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[M,N]=size(airmam);

% anomalias
for i=1:N % iteramos en cada columna año de los 75
    a_air(:,i)=(airmam(:,i)-mean(airmam')');
end


%%%%% Extrae datos de 'prate' sin los continentes
for i=1:903
    auy=reshape(squeeze(prate(:,:,i)),XX*YY,1); % 
    aa(:,i)=auy(puntero); % Es la iteracion de todo los datos
end
% Nuestro Otoño
a4=aa(:,3:12:end-1);   % M
a5=aa(:,4:12:end);   % A
a6=aa(:,5:12:end);   % M
pratemam=(a4+a5+a6)/3;

% anomalias
for i=1:N
    a_prate(:,i)=(pratemam(:,i)-mean(pratemam')');
end


%%%%%%%%%% point-to-point correlation (correlacion grilla-a-grilla) %%%%%%%%%

for i=1:M % Para todas las pocisiones M(lon*lat), se toma la primera pocion y correlacionandola dato da dato entre la aslp y atsm
    r_PrateAir(i)=corr(a_prate(i,:)',a_air(i,:)');
end

%figura
[lati,long] = meshgrid(double(lat),double(lon)); %grilla para graficar
r_PrateAirMap=reshape(r_PrateAir',length(lon),length(lat));

figure()
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,r_PrateAirMap,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('Correlación entre temp. aire y precipitación,MAM de 1948-2022')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off

%%%%%%%%%% Significancia estadistica

for j=1:100
    for i=1:M
        remCorr_AirPrt(i,j)=corr(a_air(i,:)',remuestreo(a_prate(i,:))');  
    end
end


for i=1:length(puntero)
    if r_PrateAir(i)<=prctile(remCorr_AirPrt(i,:),2.5) | r_PrateAir(i)>=prctile(remCorr_AirPrt(i,:),97.5) % alfa 0.5 %??????? y estos prctile?
        r_sig(i)=r_PrateAir(i);
    else
        r_sig(i)=NaN;
    end
end

%figura
[lati,long] = meshgrid(double(lat),double(lon)); %grilla para graficar
r_sigMap=reshape(r_sig',length(lon),length(lat));

figure()
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,r_sigMap,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
colormap('jet')
colorbar
hold on
title('Correlación significativa entre temp. aire y precipitación,MAM de 1948-2022')
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',2)
set(gca,'FontSize',15)
hold off


% %%%%%%%%     Test pa ver los graficos de valores Promedio     %%%%%%%%%%
% promair_mam=reshape(nanmean(airmam(:,:),2),length(lon),length(lat));
% promprate_mam=reshape(nanmean(pratemam(:,:),2),length(lon),length(lat));
% % ojo que nanmean(tsairdef(:,:),2) es el promedio de tsairdef(3243x75)
% % temporalmente tal que queda (3243x1) que luego se reorganiza (reshape) de
% % una columna a un mapa lonxlat.
% 
% %   Mapa promedio de temperatura del aire
% figure()
% subplot(1,2,1)
% m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
% m_contourf(long,lati,promair_mam,...
%  'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
% colormap('jet')
% colorbar
% hold on
% title('Mapa promedio de temperatura del aire en otoño')
% m_grid('Box','Fancy','LineStyle','none','FontSize',14);
% m_gshhs_c('color','k','linewidth',2)
% set(gca,'FontSize',15)
% hold off
% subplot(1,2,2)
% m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
% m_contourf(long,lati,promprate_mam,...
%  'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
% colormap('jet')
% colorbar
% hold on
% title('Mapa promedio de la precipitacion en otoño')
% m_grid('Box','Fancy','LineStyle','none','FontSize',14);
% m_gshhs_c('color','k','linewidth',2)
% set(gca,'FontSize',15)
% hold off


%%  PREGUNTAS

% 1.1 Estan bien los mpas stm y skw?

% Enero comienza del mes 1 en aaa? => D=12 E=13...

% .2 Cómo saco la significancia estadística?

% 2.1 Hago un mapa de tendencia lineal, de ser así que significaría?

% 