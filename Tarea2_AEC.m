%% Tarea2 / AEC / 8.7.23 / EOF Compuesta
%  Nataly Antipil %
clear all
% close all
clc
close all

%% Detalles a considerar antes de:

% A. Campos globales con información mensual desde enero de 1871 a diciembre de 2012

% -----
% i) Datos de viento
% ncdisp('uwnd.sig995.mon.mean.nc');
% ncdisp('vwnd.sig995.mon.mean.nc');
% veo todos los datos de latitud y longitud para despues cortarlos:
% lon = ncread('uwnd.sig995.mon.mean.nc','lon');
% lat = ncread('uwnd.sig995.mon.mean.nc','lat');

% Con la componente zonal (uwnd) y meridional (vwnd) del viento cerca de la superficie 
% (nivel sigma 0,995; 1 es la superficie), en m/s, CALCULAR LA MAGNITUD DEL VIENTO.

% CORTAR ENTRE enero de 1980 a diciembre de 2012. Extraer la región: 0°-30°S; 180°-280° (80°W).
lon = ncread('uwnd.sig995.mon.mean.nc','lon',[91],[50],[1]);
lat = ncread('uwnd.sig995.mon.mean.nc','lat',[46],[15],[1]);
YY=length(lat);
XX=length(lon);
uwnd=ncread('uwnd.sig995.mon.mean.nc','uwnd',[91 46 1308],[50 15 384],[1 1 1]);
vwnd=ncread('vwnd.sig995.mon.mean.nc','vwnd',[91 46 1308],[50 15 384],[1 1 1]);
% contourf(lon,lat,reshape(squeeze(uwnd(:,:,1)),XX,YY)') % para comprobar que es la zona grafico enero1980

mag_wnd=(uwnd.^2+vwnd.^2).^(1/2); %Magnitud del viento

% ------
% ii) Flujo de calor latente
% Campo reanalizado del flujo de calor latente en superficie, en W/m². 
% Un valor positivo (negativo) indica un flujo hacia arriba (abajo).
% ncdisp('lhtfl.mon.mean.nc') % observo detalles de dimenciones, tiempo, etc
% 
% CORTAR ENTRE enero de 1980 a diciembre de 2012. Extraer la región: 
% 0,9524°S a 29,5234°S; 180° a 279,3750° (80,625°W)
% lon2 = ncread('lhtfl.mon.mean.nc','lon'); % veo lat y lon total para
% lat2 = ncread('lhtfl.mon.mean.nc','lat'); % seleccionar la zona deseada

lon2 = ncread('lhtfl.mon.mean.nc','lon',[97],[53],[1]);
lat2 = ncread('lhtfl.mon.mean.nc','lat',[47],[16],[1]);
YY2=length(lat2);
XX2=length(lon2);
LH=ncread('lhtfl.mon.mean.nc','lhtfl',[97 47 1308],[53 16 384],[1 1 1]);
% contourf(lon2,lat2,reshape(squeeze(LH(:,:,1)),XX2,YY2)')

% B. Índices E y C de El Niño-Oscilación del Sur (ENOS). Desde enero de 1888 a febrero de 2022.
% Mayor información en: http://met.igp.gob.pe/elnino/indices.html
% 
% El índice E muestra eventos El Niño extraordinarios, con un gran desarrollo 
% en la zona costera occidental de Sudamérica tropical-subtropical. 
% 
% El índice C es para eventos El Niño/La Niña que tienen un mayor desarrollo 
% en el Pacífico ecuatorial central.

load EC.txt

%% 1.-
% Obtener el promedio climatológico (o promedio de largo plazo) de la
% magnitud del viento (wind speed) y del flujo de calor latente (latent heat). 
% El promedio climatológico corresponde a la media de los 396 meses (33 años). 
% No confundir con el ciclo anual, que corresponde al promedio climatológico 
% de cada meses, donde hay 33 elementos para calcular el promedio de cada mes.
         % ¿Qué tipo de relación entre ambas variables, %%%%%%%%%%%%%%%%
         % sugieren los patrones espaciales (visualmente)? %%%%%%%%%%%%%

prom_wnd=nanmean(mag_wnd(:,:,:),3);              % Solo este promedio simple ?????????????????????                       
prom_LH=nanmean(LH(:,:,:),3);

figure()
%barra=[]
[lati,long] = meshgrid(double(lat),double(lon)); %grilla para graficar
subplot 211
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,prom_wnd,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
hold on
m_contour(long, lati, prom_wnd, [0 0], 'LineColor', 'k', 'LineWidth', 2);
hold off
colormap(jet(30))
colorbar
hold on
title({'Magnitud del viento [m s^{-1}]'})
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',1)
%xlabel('Longitud','color','k');
%ylabel('Latitud','color','k');
set(gca,'FontSize',15)
set(gcf,'color','w');
subplot 212
[lati2,long2] = meshgrid(double(lat2),double(lon2));
m_proj('mercator','lon2',[min(long2(:)) max(long2(:))],'lat2',[min(lati2(:)) max(lati2(:))])
m_contourf(long2,lati2,prom_LH,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
hold on
m_contour(long2, lati2, prom_LH, [0 0], 'LineColor', 'k', 'LineWidth', 2);
hold off
colormap(jet(30))
colorbar
hold on
title({'Flujo de calor latente [W m^{-2}]'})
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',1)
%xlabel('Longitud','color','k');
%ylabel('Latitud','color','k');
set(gca,'FontSize',15)
set(gcf,'color','w');
sgtitle('Promedio climatológico (1980-2012)')

%% 2.-
% A partir de los campos de anomalías (extraer ciclo anual), obtener los modos 
% de covarianza acoplados o combinados de los campos de rapidez del viento y 
% del flujo turbulento de calor latente, que no están mezclados, de acuerdo al 
% criterio de North et al. 1982. (EOF combinada)

c = 0;                                                             % FECHAS
for i = 1981:2012
    for j = 1:12
        c = c+1;
        fecha(c) = datenum(i,j,1);
    end
end
                                                 % DIMENCIONES DE LOS CAMPOS
[X1,Y1,T1]=size(mag_wnd);
[X2,Y2,T2]=size(LH);

%Notas:
% Estamos trabajando con anomalias (en  EOF)
% Esta pesca los datos y ve como varian en general.
% Estandarizamos las variables para dar a cada serie de tiempo el mismo peso.

% tabajando directamente con los datos
% an_wnd(:,:,:) = (mag_wnd(:,:,:) - nanmean(mag_wnd(:,:,:),3))./nanstd(mag_wnd(:,:,:),0,3);
% an_LH(:,:,:) = (LH(:,:,:) - nanmean(LH(:,:,:),3))./nanstd(LH(:,:,:),0,3);

% Extrayendo el ciclo anual                          % ANOMALIAS A TRABAJAR
for i = 1:12
    an_wnd(:,:,i:12:size(mag_wnd,3)) = (mag_wnd(:,:,i:12:end) - nanmean(mag_wnd(:,:,i:12:end),3))./nanstd(mag_wnd(:,:,i:12:end),0,3);
end
for i = 1:12
    an_LH(:,:,i:12:size(LH,3)) = (LH(:,:,i:12:end) - nanmean(LH(:,:,i:12:end),3))./nanstd(LH(:,:,i:12:end),0,3);
end

F1 = reshape(permute(an_wnd,[3 1 2]),T1,X1*Y1);
F2 = reshape(permute(an_LH,[3 1 2]),T2,X2*Y2);
F=[F1 F2]; % (384x1598)                   % VECTOR F, EN ESTE CASO COBINADO

[L,E,A,err] = EOF(F);                                             % EOF
% % L: varianza de cada modo (384x1) (=length de fecha)
% % E: mapa = correlación con el modo (1598x1598)
% % A: modo % ¿Qué es A? ¿por que tiene dimencion (384x1598) siendo varianza 384x1? ?????????????????????????????????'
% % err: (1x1598)

varianzas1=(L(:)/sum(L))*100;

N=length(F(1,:)) % mapas o elementos
M=length(F(:,1)) % grillas o series de tiempo
NN=N;

figure()                                       % CRITERIO North et al. 1982.
plot(L(1:30),'ob','LineWidth',1) % varianza de cada modo
hold on
grid on
% varianza de cada modo +- la aproximación del intervalo de confianza:
plot(L(1:30)+L(1:30)*sqrt(2/NN),'+r','LineWidth',1)
plot(L(1:30)-L(1:30)*sqrt(2/NN),'+r','LineWidth',1) 
title('Criterio de North et al.(1982)','FontSize',15)
legend('Varianza de cada modo','Intervalo de confianza','Location','northeast','FontSize',15)
xlabel('Modo','color','k','FontSize',15); 
ylabel('Valores propios','color','k','FontSize',15);
set(gcf,'color','w');
axis tight;
grid on

% modos de covarianza acoplados o combinados de los campos de rapidez 
% del viento y del flujo turbulento de calor latente (?)             

for k=1:3; % poner aqí la cantidad de modos estadisticamente significatívos.
                                                                     % CORR
for x = 1:X1
    for y = 1:Y1
        wndcorr(x,y) = corr(squeeze(an_wnd(x,y,:)),A(:,k)/std(A(:,k)));
    end
end
for x = 1:X2
    for y = 1:Y2
        LHcorr(x,y) = corr(squeeze(an_LH(x,y,:)),A(:,k)/std(A(:,k)));
    end
end

modo=A(:,k)/std(A(:,k)); % componente principal del modo

figure(k) %                                                               Esta bien este plot?
subplot(2,2,[1,2])
 plot(fecha,-modo), axis tight
 datetick('x','yyyy','keepticks')
 hold on
 title(['Modo ' num2str(k) ': ' num2str((L(k)/sum(L))*100) '% de varianza total'])
 set(gca,'FontSize',23)
 grid on
   xlabel('Fecha')
ylabel('\sigma')
axis tight

subplot(2,2,3)
m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_contourf(long,lati,-wndcorr,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
hold on
m_contour(long, lati, wndcorr, [0 0], 'LineColor', 'k', 'LineWidth', 2);
hold off
colormap(jet(30))
colorbar
hold on
title({"Mapa de correlación modo " + k + " "," y magnitud del viento"})
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',1)
set(gca,'FontSize',15)
set(gcf,'color','w');

subplot(2,2,4)
[lati2,long2] = meshgrid(double(lat2),double(lon2));
m_proj('mercator','lon2',[min(long2(:)) max(long2(:))],'lat2',[min(lati2(:)) max(lati2(:))])
m_contourf(long2,lati2,-LHcorr,...
 'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
hold on
m_contour(long2, lati2, LHcorr, [0 0], 'LineColor', 'k', 'LineWidth', 2);
hold off
colormap(jet(30))
colorbar
hold on
title({"Mapa de correlación modo " + k + " "," y flujo de calor latente"})
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_c('color','k','linewidth',1)
set(gca,'FontSize',15)
set(gcf,'color','w');
end
%sgtitle('Promedio climatológico (1980-2012)')


%% 3. Para cada modo, evaluar:

% a. ¿Se observa tendencia lineal de la componente principal del modo? 
%    Calcular la significancia estadística.

k=1
modo=A(:,k)/std(A(:,k)); % componente principal del modo

% plot(fecha,-modo) % Componente principal.

xx=(1:length(modo));
% p(1) es beta y p(2) es alfa
p=polyfit(xx',-modo',1)
ttx=polyval(p,xx'); % componente lineal

figure()
plot(fecha,-modo,'LineWidth',1)
hold on
plot(fecha,ttx,'r','LineWidth',1)
 datetick('x','yyyy','keepticks')
 hold on
 title(['Modo ' num2str(k) ': ' num2str((L(k)/sum(L))*100) '% de varianza total'])
 legend('Serie de tiempo','Tendencia','Location','best','FontSize',12) %,'Orientation','horizontal'
 set(gca,'FontSize',23)
 grid on
xlabel('Fecha')
ylabel('\sigma')
axis tight
set(gcf,'color','w')

% se extrae tendencia lineal
txa=-modo-ttx';

% tendencia lineal
Bobs=p(1)*10; % C/decada

% significancia
for m=1:1000
    ser=remuestreo(-modo);
    p=polyfit(xx',ser',1);
    Bmon(m)=p(1)*10;
end
[Bobs prctile(Bmon,2.5) prctile(Bmon,97.5)] % con alfa=5%
[Bobs prctile(Bmon,5) prctile(Bmon,95)] % con alfa=10%

% [M,N,T]=size(wndcorr); % dimensiones del tensor.                        ????????????????
% 
% % Reshape del tensor en una matriz bidimensional (M x N, T)
% %reshaped_tensor = reshape(tensor, [], size(tensor, 3));
% 
% % Crear un vector de tiempo o índices de tiempo (puede ser meses, años, etc.)
% tiempo = 1:size(wndcorr, 3);
% 
% for i=1:384 % Para todo el tiempo
%     aaa(:,i)=reshape(squeeze(wndcorr(:,:,i)),XX*YY,1); % cambio las dimenciones a (lon*lat,tiempo)
% end
% 
% % Aplicar la regresión lineal para cada punto en la matriz
% tendencia = zeros(size(aaa, 1), 2);  % Vector para almacenar las pendientes y las intersecciones
% for i = 1:size(aaa, 1)
%     % Obtener los coeficientes de la regresión lineal
%     coeficientes = polyfit(tiempo, aaa(i, :), 1);
% 
%     % Almacenar la pendiente y la intersección en el vector 'tendencia'
%     tendencia(i, :) = coeficientes;
% end
% 
% % Reshape de nuevo los resultados a la forma original del tensor
% tendencia_reshaped = reshape(tendencia(:, 1), M, N);
% 
% figure()
% m_proj('mercator','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
% m_contourf(long,lati,tendencia_reshaped,...
%  'LabelSpacing',600,'ShowText','on','LineWidth',0.1, 'color',[0 0 0])
% colormap('jet')
% colorbar
% hold on
% title('Tendencia lineal anual de la temperatura del aire')
% m_grid('Box','Fancy','LineStyle','none','FontSize',14);
% m_gshhs_c('color','k','linewidth',2)
% set(gca,'FontSize',15)
% hold off

% b. Analizar el modo, poniendo atención en el tipo de relación que se observa:
%    si la rapidez del viento aumenta o disminuye qué ocurre con el flujo de 
%    calor latente (o al revés)

% Puedo correlacionar los dos mapas (subplot 2,23 y 2,2,4)               ?????????????
   % Aquí solo hay que tramuyar como al correlacionar el modo X con la
   % magnitud del viento o con el flujo de calor latente ambos dan valores
   % similares existendo mucha correlación entre los dos mapas, es
   % decir,cuando mag_wnd tiene valores positucos LH también y cuando esta
   % negativo el otro igual disminuye en escalas incluso similares ---v
   % 
   % (la escala de la correlación implica que la escala de ambos aumentos es
   % similar o es solo porque esta normalizada, y si esta normalizada, que
   % implica que ambas tengan la misma escakla? que la proporcion de vaiación respecto
   % al total de cada variable es igual?)                         ????????????????????????????????



% c. Cómo se relaciona linealmente (o qué tipo de relación lineal tiene) la 
%    componente principal del modo con los índices E y C del ENSO. ----> Relación lineal = correlació
%    Calcular significancia estadística

% Los indices E y C resumen la variabilidad asociada a El Niño y La Niña, representando el calentamiento superficial anómalo en el Pacífico este y centro, respectivamente. Debido a la forma en que fueron calculados (usando componentes principales) la correlación lineal entre ellos es baja, por lo que permite distinguir mejor la variabilidad propia de cada una de estas regiones.
E_ec=EC(1201:1584,3);
C_ec=EC(1201:1584,4);
fecha = datenum(EC(1201:1584,1),EC(1201:1584,2), 1)

figure()
subplot 211
plot(fecha,E_ec,'LineWidth',1.5)
datetick('x', 'yyyy');
title(['E index'])
set(gca,'FontSize',16)
grid on
xlabel('Fecha','FontSize',12); %ylabel('E');
axis tight
set(gcf,'color','w')
subplot 212
plot(fecha,C_ec,'LineWidth',1.5)
datetick('x', 'yyyy');
title(['C index'])
set(gca,'FontSize',16)
grid on
xlabel('Fecha','FontSize',12); %ylabel('E');
axis tight
set(gcf,'color','w')

for k=1:3
modo=A(:,k)/std(A(:,k))
% Calcular el coeficiente de correlación lineal
[E_corr(k),p_E(k)] = corr(-modo,E_ec)
[C_corr(k),p_C(k)] = corr(-modo,C_ec)
end


