close all
clear all
clc



%% DATI
tic
partenza0 = 'Earth';
partenza1 = '2003YN107';
partenza2 = '2006JY26';
partenza3 = '2006RH120';
partenza4 = '2014EK24';
partenza5 = '2012FC71';
partenza6 = '2013BS45';

bodies = {partenza0, partenza1, partenza2, partenza3, partenza4, partenza5, partenza6};
N_asteroids = length(bodies) - 1; % 6 asteroidi più la Terra

% Tutte le permutazioni possibili dei 6 asteroidi
perm_combos = perms(2:length(bodies)); % Permutazioni dei numeri da 2 a 7 (indice dei corpi celesti esclusa la Terra)

% Struttura per memorizzare il miglior risultato tra tutte le combinazioni
best_global = struct('DVTOTBest', inf, 'PBest', [], 'SeqBest', []);

% Parametri PSO
Npar = 5; % Numero delle particelle
Nmax = 30; % Numero massimo di iterazioni
CI = 1.2;
CC = 1.4995;
CS = 1.4995;

% Definizione del bound (finestre temporali)
bound = [juliandate('15-Jun-2025','dd-mm-yyyy'),juliandate('30-Jun-2025','dd-mm-yyyy'); ... % Partenza da Terra
         juliandate('15-Jul-2026','dd-mm-yyyy'),juliandate('30-Jul-2026','dd-mm-yyyy');...
         juliandate('15-Jul-2027','dd-mm-yyyy'),juliandate('30-Jul-2027','dd-mm-yyyy');...
         juliandate('29-Nov-2028','dd-mm-yyyy'),juliandate('15-Dec-2028','dd-mm-yyyy');...
         juliandate('29-Jun-2029','dd-mm-yyyy'),juliandate('05-Jul-2029','dd-mm-yyyy');...
         juliandate('30-Aug-2030','dd-mm-yyyy'),juliandate('31-Aug-2030','dd-mm-yyyy');];
Vlim = [-5, +5; -5, +5; -5, +5; -5,+5; -5, +5; -5,+5];

% Loop su tutte le combinazioni possibili
for combo_idx = 1:size(perm_combos, 1)
    Seq = [1, perm_combos(combo_idx, :)]; % Inserisci la Terra come punto di partenza

    % PSO per la combinazione corrente
    particle = struct('Seq', [], 'Part0', [], 'Arr1', [], 'Part1', [], 'Arr2', [], 'Part2', [], 'Arr3', [], 'Part3', [], 'Arr4', [], 'Part4', [], 'Arr5', [], 'DV0', [], 'DV1', [], 'DV2', [], 'DV3', [], 'DV4', [], 'DVTOT', inf, 'DVTOTBest', inf, 'PBest', [], 'SeqBest', []);
    gbest = struct('DVTOTBest', inf, 'PBest', [], 'SeqBest', []);
    
    % Inizializzazione
    for part = 1:Npar
        particle(part).Seq = Seq;
        particle(part).Part0 = round(bound(1, 2) + rand * (bound(1, 1) - bound(1, 2)));
        particle(part).Part0stringa = datestr(datetime(particle(part).Part0, 'convertfrom', 'juliandate'), 'dd-mmm-yyyy');
        particle(part).Arr1 = ceil(bound(2, 2) + rand * (bound(2, 1) - bound(2, 2)));
        particle(part).Arr1stringa = datestr(datetime(particle(part).Arr1, 'convertfrom', 'juliandate'), 'dd-mmm-yyyy');
        particle(part).Part1 = particle(part).Arr1 + 10;
        particle(part).Part1stringa = datestr(datetime(particle(part).Part1, 'convertfrom', 'juliandate'), 'dd-mmm-yyyy');
        particle(part).Arr2 = ceil(bound(3, 2) + rand * (bound(3, 1) - bound(3, 2)));
        particle(part).Arr2stringa = datestr(datetime(particle(part).Arr2, 'convertfrom', 'juliandate'), 'dd-mmm-yyyy');
        particle(part).Part2 = particle(part).Arr2 + 10;
        particle(part).Part2stringa = datestr(datetime(particle(part).Part2, 'convertfrom', 'juliandate'), 'dd-mmm-yyyy');
        particle(part).Arr3 = ceil(bound(4, 2) + rand * (bound(4, 1) - bound(4, 2)));
        particle(part).Arr3stringa = datestr(datetime(particle(part).Arr3, 'convertfrom', 'juliandate'), 'dd-mmm-yyyy');
        particle(part).Part3 = particle(part).Arr3 + 10;
        particle(part).Part3stringa = datestr(datetime(particle(part).Part3, 'convertfrom', 'juliandate'), 'dd-mmm-yyyy');
        particle(part).Arr4 = ceil(bound(5, 2) + rand * (bound(5, 1) - bound(5, 2)));
        particle(part).Arr4stringa = datestr(datetime(particle(part).Arr4, 'convertfrom', 'juliandate'), 'dd-mmm-yyyy');
        particle(part).Part4 = particle(part).Arr4 + 10;
        particle(part).Part4stringa = datestr(datetime(particle(part).Part4, 'convertfrom', 'juliandate'), 'dd-mmm-yyyy');
        particle(part).Arr5 = ceil(bound(6, 2) + rand * (bound(6, 1) - bound(6, 2)));
        particle(part).Arr5stringa = datestr(datetime(particle(part).Arr5, 'convertfrom', 'juliandate'), 'dd-mmm-yyyy');

        particle(part).velocity = zeros(1, 10);

        for i = 1:N_asteroids-1 % Cambia N_asteroids a N_asteroids-1 per evitare di superare Part4
            [dv, dp, ad] = Pork(particle(part).(sprintf('Part%d', i-1)), particle(part).(sprintf('Part%d', i-1))+10, particle(part).(sprintf('Arr%d', i)), particle(part).(sprintf('Arr%d', i))+10, bodies{Seq(i)}, bodies{Seq(i+1)});
            [colMin, rowIdx] = min(dv);
            [dv, colIdx] = min(colMin);
            rowIdx = rowIdx(colIdx);
            EffPart = dp(colIdx);
            EffArr = ad(rowIdx);
            particle(part).(sprintf('DV%d', i-1)) = dv;
            particle(part).(sprintf('Part%d', i-1)) = EffPart;
            particle(part).(sprintf('Part%dstringa', i-1)) = datestr(datetime(EffPart, 'convertfrom', 'juliandate'), 'dd-mmm-yyyy');
            particle(part).(sprintf('Arr%d', i)) = EffArr;
            particle(part).(sprintf('Arr%dstringa', i)) = datestr(datetime(EffArr, 'convertfrom', 'juliandate'), 'dd-mmm-yyyy');
        end
        
        particle(part).DVTOT = sum([particle(part).DV0, particle(part).DV1, particle(part).DV2, particle(part).DV3, particle(part).DV4]);
        particle(part).DVTOTBest = particle(part).DVTOT;
        particle(part).PBest = [particle(part).Part0, particle(part).Arr1, particle(part).Part1, particle(part).Arr2, particle(part).Part2, particle(part).Arr3, particle(part).Part3, particle(part).Arr4, particle(part).Part4, particle(part).Arr5];
        particle(part).SeqBest = particle(part).Seq;
        if particle(part).DVTOTBest < gbest.DVTOTBest
            gbest.DVTOTBest = particle(part).DVTOTBest;
            gbest.PBest = particle(part).PBest;
            gbest.SeqBest = particle(part).SeqBest;
        end
    end

    iterations = 1;
    while iterations <= Nmax
        for part = 1:Npar
            r1 = rand(1, 10);
            r2 = rand(1, 10);
            particle(part).velocity = ceil(CI * particle(part).velocity + CC * r1 .* (particle(part).PBest - [particle(part).Part0, particle(part).Arr1, particle(part).Part1, particle(part).Arr2, particle(part).Part2, particle(part).Arr3, particle(part).Part3, particle(part).Arr4, particle(part).Part4, particle(part).Arr5]) + CS * r2 .* (gbest.PBest - [particle(part).Part0, particle(part).Arr1, particle(part).Part1, particle(part).Arr2, particle(part).Part2, particle(part).Arr3, particle(part).Part3, particle(part).Arr4, particle(part).Part4, particle(part).Arr5]));
            particle(part).Part0 = particle(part).Part0 + particle(part).velocity(1);
            particle(part).Arr1 = particle(part).Arr1 + particle(part).velocity(2);
            particle(part).Part1 = particle(part).Part1 + particle(part).velocity(3);
            particle(part).Arr2 = particle(part).Arr2 + particle(part).velocity(4);
            particle(part).Part2 = particle(part).Part2 + particle(part).velocity(5);
            particle(part).Arr3 = particle(part).Arr3 + particle(part).velocity(6);
            particle(part).Part3 = particle(part).Part3 + particle(part).velocity(7);
            particle(part).Arr4 = particle(part).Arr4 + particle(part).velocity(8);
            particle(part).Part4 = particle(part).Part4 + particle(part).velocity(9);
            particle(part).Arr5 = particle(part).Arr5 + particle(part).velocity(10);

            particle(part).Part0 = ceil(max(bound(1,1), min(particle(part).Part0, bound(1,2))));
            particle(part).Arr1 = ceil(max(bound(2,1), min(particle(part).Arr1, bound(2,2))));
            particle(part).Arr2 = ceil(max(bound(3,1), min(particle(part).Arr2, bound(3,2))));
            particle(part).Arr3 = ceil(max(bound(4,1), min(particle(part).Arr3, bound(4,2))));
            particle(part).Arr4 = ceil(max(bound(5,1), min(particle(part).Arr4, bound(5,2))));
            particle(part).Arr5 = ceil(max(bound(6,1), min(particle(part).Arr5, bound(6,2))));

            for i = 1:N_asteroids-1 % Cambia N_asteroids a N_asteroids-1 per evitare di superare Part4
                [dv, dp, ad] = Pork(particle(part).(sprintf('Part%d', i-1)), particle(part).(sprintf('Part%d', i-1))+10, particle(part).(sprintf('Arr%d', i)), particle(part).(sprintf('Arr%d', i))+10, bodies{Seq(i)}, bodies{Seq(i+1)});
                [colMin, rowIdx] = min(dv);
                [dv, colIdx] = min(colMin);
                rowIdx = rowIdx(colIdx);
                EffPart = dp(colIdx);
                EffArr = ad(rowIdx);
                particle(part).(sprintf('DV%d', i-1)) = dv;
                particle(part).(sprintf('Part%d', i-1)) = EffPart;
                particle(part).(sprintf('Part%dstringa', i-1)) = datestr(datetime(EffPart, 'convertfrom', 'juliandate'), 'dd-mmm-yyyy');
                particle(part).(sprintf('Arr%d', i)) = EffArr;
                particle(part).(sprintf('Arr%dstringa', i)) = datestr(datetime(EffArr, 'convertfrom', 'juliandate'), 'dd-mmm-yyyy');
            end

            particle(part).DVTOT = sum([particle(part).DV0, particle(part).DV1, particle(part).DV2, particle(part).DV3, particle(part).DV4]);
            if particle(part).DVTOT < particle(part).DVTOTBest
                particle(part).DVTOTBest = particle(part).DVTOT;
                particle(part).PBest = [particle(part).Part0, particle(part).Arr1, particle(part).Part1, particle(part).Arr2, particle(part).Part2, particle(part).Arr3, particle(part).Part3, particle(part).Arr4, particle(part).Part4, particle(part).Arr5];
                particle(part).SeqBest = particle(part).Seq;
            end
            if particle(part).DVTOTBest < gbest.DVTOTBest
                gbest.DVTOTBest = particle(part).DVTOTBest;
                gbest.PBest = particle(part).PBest;
                gbest.SeqBest = particle(part).SeqBest;
            end
        end
        iterations = iterations + 1;
    end

    % Aggiorna il miglior risultato globale se necessario
    if gbest.DVTOTBest < best_global.DVTOTBest
        best_global.DVTOTBest = gbest.DVTOTBest;
        best_global.PBest = gbest.PBest;
        best_global.SeqBest = gbest.SeqBest;
    end
end

% Risultato finale
fprintf('Miglior \x{0394}V totale: %f\n', best_global.DVTOTBest);
fprintf('Sequenza migliore: %s\n', strjoin(bodies(best_global.SeqBest), ' -> '));
toc
