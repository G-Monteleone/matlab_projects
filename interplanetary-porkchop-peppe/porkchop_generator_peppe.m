clc; clear; close all;

tic;

% % Window Earth-Asteroids
% departure_1 = '31-Mar-2025'; %Data di partenza minima 
% departure_2 = '30-Apr-2026'; %Data di partenza massima (Differenza di 395 giorni in questo caso)
% arrival_1 = '13-Aug-2026'; %Data di arrivo minima
% arrival_2 = '12-Sep-2027'; %Data di arrivo massima (Differenza di 395 giorni anche in questo caso)

% % Window Asteroid-Asteroids Prima Iterazione
% departure_1 = '31-Mar-2025'; %Data di partenza minima 
% departure_2 = '30-Apr-2026'; %Data di partenza massima (Differenza di 395 giorni in questo caso)
% arrival_1 = '13-Aug-2026'; %Data di arrivo minima
% arrival_2 = '17-Oct-2028'; %Data di arrivo massima (Differenza di 395 giorni anche in questo caso)

% Window Terra-Asteroide Con date per iterazione
departure_1 = '3-Jan-2029'; %Data di partenza minima 
departure_2 = '3-Jan-2030'; %Data di partenza massima (Differenza di 395 giorni in questo caso)
arrival_1 = '1-Mar-2030'; %Data di arrivo minima
arrival_2 = '28-Feb-2031'; %Data di arrivo massima (Differenza di 395 giorni anche in questo caso)

% %% M2020: depart 30-Jul-2020 arrive 18-Feb-2021
% departure_m2020 = '25-Jul-2020';
% arrival_m2020 = '15-Feb-2021';
% 
% departure_m2020 = '30-Jul-2020';
% arrival_m2020 = '18-Feb-2021';

% departure_m2020 = '20-May-2018';
% arrival_m2020 = '6-Dec-2018';

jd_departure_1 = juliandate(departure_1,'dd-mm-yyyy');
jd_departure_2 = juliandate(departure_2,'dd-mm-yyyy');
jd_arrival_1 = juliandate(arrival_1,'dd-mm-yyyy');
jd_arrival_2 = juliandate(arrival_2,'dd-mm-yyyy');

% jd_departure_m2020 = juliandate(departure_m2020,'dd-mm-yyyy');
% jd_arrival_m2020 = juliandate(arrival_m2020,'dd-mm-yyyy');

%%
% %N_Relativi Fissi
% N_departure_dates = 395; %Prima avevamo 500
% N_arrival_dates = 395; %Prima avevamo 500
% departure_dates = linspace(jd_departure_1,jd_departure_2,N_departure_dates); 
% arrival_dates = linspace(jd_arrival_1,jd_arrival_2,N_arrival_dates);

% %N_Relativi Variabili con le date scelte con step di 1 giorno 
N_departure_dates = jd_departure_2-jd_departure_1; 
N_arrival_dates = jd_arrival_2-jd_arrival_1; 
departure_dates = linspace(jd_departure_1,jd_departure_2,N_departure_dates); 
arrival_dates = linspace(jd_arrival_1,jd_arrival_2,N_arrival_dates);

%N_Relativi Variabili con le date scelte con step di 5 giorni 
% N_departure_dates = (jd_departure_2-jd_departure_1)/5; 
% N_arrival_dates = (jd_arrival_2-jd_arrival_1)/5; 
% departure_dates = linspace(jd_departure_1,jd_departure_2,N_departure_dates); 
% arrival_dates = linspace(jd_arrival_1,jd_arrival_2,N_arrival_dates);

% Initialize arrays dove salveró DV di Lambert
DV_depart_short = zeros(N_arrival_dates,N_departure_dates); 
DV_depart_long = zeros(N_arrival_dates,N_departure_dates);
DV_arrive_short = zeros(N_arrival_dates,N_departure_dates);
DV_arrive_long = zeros(N_arrival_dates,N_departure_dates);

% departure planet
planet1 = '2014EK24'; %Cambiato da Earth con identificativo asteroide

% arrival planet
planet2 = '2012FC71'; %Cambiato da Mars con identificativo asteroide

% initialize data matrices
DV_depart = zeros(N_arrival_dates,N_departure_dates,2); %DeltaV totali dei due casi short e long
DV_arrive = zeros(N_arrival_dates,N_departure_dates,2);
DT = zeros(N_arrival_dates,N_departure_dates);

for j = 1:N_departure_dates
    for i = 1:N_arrival_dates %Considera ogni giorno di arrivo per un dato giorno di partenza (come il codice giá fatto)
        
        % calculate transfer at (i,j) transfer
        [DV_depart(i,j,1),DV_arrive(i,j,1)] = ... %Prograde
            interplanetary_transfer(departure_dates(j),...
            arrival_dates(i),planet1,planet2,1);
        [DV_depart(i,j,2),DV_arrive(i,j,2)] = ... %Retrograde
            interplanetary_transfer(departure_dates(j),...
            arrival_dates(i),planet1,planet2,-1);
        
        % populate transfer time matrix
        DT(i,j) = arrival_dates(i) - departure_dates(j);
        
        % progress updates
        num_complete = N_arrival_dates*(j-1) + i;
        num_total = N_arrival_dates * N_departure_dates;
        clc;
        fprintf('Progress: %.1f%% \n',100*num_complete/num_total);
    end
end

disp('CALCULATION COMPLETE');


%% REMOVE NEGATIVE TIME TRANSFERS
for j = 1:N_departure_dates %Controllo per evitare date di partenza successive all'arrivo
    for i = 1:N_arrival_dates
        
        if departure_dates(j) > arrival_dates(i)
            DV_depart(i,j,:) = 10000 * ones(1,1,2);
            DV_arrive(i,j,:) = 10000 * ones(1,1,2);
        end
    end
end


%% PLOTTING
% plotdefaults(16,14,2);

X = datenum(datetime(departure_dates,'convertfrom','juliandate'));
Y = datenum(datetime(arrival_dates,'convertfrom','juliandate'));
Z = min(DV_depart + DV_arrive,[],3);
V = linspace(0,15,100);

% plot total Delta-V porkchop plot
figure(1)
contour(X,Y,Z,V);
hold on
% plot(datenum(datetime(jd_departure_m2020,'convertfrom','juliandate')),...
%     datenum(datetime(jd_arrival_m2020,'convertfrom','juliandate')),'k+')
colormap(jet);
hold off

% colorbar
cb = colorbar;
cb.Title.String ='Total $\Delta V$ [km/s]';
cb.Title.Interpreter = 'latex';
set(cb,'ticklabelinterpreter','latex');

% date ticks
set(gca,'ytick',linspace(min(Y),max(Y),7))
set(gca,'xtick',linspace(min(X),max(X),5))
datetick('x',2,'keeplimits','keepticks')
datetick('y',2,'keeplimits','keepticks')

% labels
xlabel('Departure Date')
ylabel('Arrival Date')
% legend('Total $\Delta V$ Contour','NASA Mars 2020 Mission','location','northwest')
title('2014EK24 to 2012FC71 2025-2030 Launch Window') %Cambiare in funzione di date e corpi celesti

% save figure
% saveas(gca,'porkchop_plot.pdf')
% exportgraphics(gcf,'porkchop.jpg','resolution','500')







% %% Plot DV_Prograde
% % plotdefaults(16,14,2);
% 
% X = datenum(datetime(departure_dates,'convertfrom','juliandate'));
% Y = datenum(datetime(arrival_dates,'convertfrom','juliandate'));
% Z_1 = DV_depart(:,:,1) + DV_arrive(:,:,1);
% V = linspace(0,15,100);
% 
% % plot total Delta-V porkchop plot
% figure(2)
% contour(X,Y,Z,V);
% hold on
% % plot(datenum(datetime(jd_departure_m2020,'convertfrom','juliandate')),...
% %     datenum(datetime(jd_arrival_m2020,'convertfrom','juliandate')),'k+')
% colormap(jet);
% hold off
% 
% % colorbar
% cb = colorbar;
% cb.Title.String ='Total $\Delta V$ [km/s]';
% cb.Title.Interpreter = 'latex';
% set(cb,'ticklabelinterpreter','latex');
% 
% % date ticks
% set(gca,'ytick',linspace(min(Y),max(Y),7))
% set(gca,'xtick',linspace(min(X),max(X),5))
% datetick('x',2,'keeplimits','keepticks')
% datetick('y',2,'keeplimits','keepticks')
% 
% % labels
% xlabel('Departure Date')
% ylabel('Arrival Date')
% % legend('Total $\Delta V$ Contour','NASA Mars 2020 Mission','location','northwest')
% title('Earth to 2012FC71 2025-2026 Launch Window') %Cambiato 2040 con 2025-2026
% 
% % save figure
% % saveas(gca,'porkchop_plot.pdf')
% % exportgraphics(gcf,'porkchop.jpg','resolution','500')






% %% Plot DV_Retrograde
% % plotdefaults(16,14,2);
% 
% X = datenum(datetime(departure_dates,'convertfrom','juliandate'));
% Y = datenum(datetime(arrival_dates,'convertfrom','juliandate'));
% Z_2 = DV_depart(:,:,2) + DV_arrive(:,:,2);
% V = linspace(0,15,100);
% 
% % plot total Delta-V porkchop plot
% figure(3)
% contour(X,Y,Z,V);
% hold on
% % plot(datenum(datetime(jd_departure_m2020,'convertfrom','juliandate')),...
% %     datenum(datetime(jd_arrival_m2020,'convertfrom','juliandate')),'k+')
% colormap(jet);
% hold off
% 
% % colorbar
% cb = colorbar;
% cb.Title.String ='Total $\Delta V$ [km/s]';
% cb.Title.Interpreter = 'latex';
% set(cb,'ticklabelinterpreter','latex');
% 
% % date ticks
% set(gca,'ytick',linspace(min(Y),max(Y),7))
% set(gca,'xtick',linspace(min(X),max(X),5))
% datetick('x',2,'keeplimits','keepticks')
% datetick('y',2,'keeplimits','keepticks')
% 
% % labels
% xlabel('Departure Date')
% ylabel('Arrival Date')
% % legend('Total $\Delta V$ Contour','NASA Mars 2020 Mission','location','northwest')
% title('Earth to 2012FC71 2025-2026 Launch Window') %Cambiato 2040 con 2025-2026
% 
% % save figure
% % saveas(gca,'porkchop_plot.pdf')
% % exportgraphics(gcf,'porkchop.jpg','resolution','500')

%% plot separate contour plots
Z1 = min(DV_depart,[],3);
Z2 = min(DV_arrive,[],3);
V1 = linspace(0,10,8);
V2 = linspace(0,10,8);
V3 = linspace(min(DT,[],'all'),max(DT,[],'all'),10);
% plotdefaults(16,14,3);


figure(4) %prima era 2
[C1,h1] = contour(X,Y,Z1,V1,'r');
hold on
[C2,h2] = contour(X,Y,Z2,V2,'b');
% [C3,h3] = contour(X,Y,DT,V3,'k');
hold off

clabel(C1,h1,'interpreter','latex','fontsize',8)
clabel(C2,h2,'interpreter','latex','fontsize',8)
% clabel(C3,h3,'interpreter','latex','fontsize',14)


% date ticks
set(gca,'ytick',linspace(min(Y),max(Y),7))
set(gca,'xtick',linspace(min(X),max(X),5))
datetick('x',2,'keeplimits','keepticks')
datetick('y',2,'keeplimits','keepticks')

% labels
xlabel('Departure Date')
ylabel('Arrival Date')
legend('Departure $\Delta V$ [km/s]','Arrival $\Delta V$ [km/s]','location','northwest')
title('2014EK24 to 2012FC71 2025-2030 Launch Window') %cambiato 2040 con 2025-2026

exportgraphics(gca,'porkchop_separate.jpg','resolution','500') %Cambiato 500 con 395


elapsedTime = toc;
fprintf('La simulazione ha impiegato %.2f secondi.\n', elapsedTime);