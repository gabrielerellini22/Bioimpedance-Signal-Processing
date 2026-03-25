clear all % Clears the workspace
clc % Clears the command window

% === SETTINGS ===
% Path to the folder containing the data
cartella_cellule = 'C:\Users\Gabriele Rellini\OneDrive\Desktop\Università\2. Magistrale\Tesi\Sperimentazione\4 - Misure Perugia 26-11-25\Misure\Esperimento';
% Saving path
salvaPath = 'C:\Users\Gabriele Rellini\OneDrive\Desktop\Università\2. Magistrale\Tesi\Sperimentazione\4 - Misure Perugia 26-11-25\Elaborazione';

% Find all subfolders in the main directory
folders = dir(cartella_cellule); 
folders = folders([folders.isdir]); 
folders = folders(~ismember({folders.name},{'.','..'})); 
dati_per_cartella = struct; 

% Loop through all folders to process the data
for k = 1:length(folders)
    cartella = fullfile(cartella_cellule, folders(k).name); 
    nome_campo = matlab.lang.makeValidName(folders(k).name); 
    
    files_csv = dir(fullfile(cartella,'*.csv')); 
    files_a = files_csv(endsWith({files_csv.name},'a.csv')); 
    files_b = files_csv(endsWith({files_csv.name},'b.csv')); 
    
    if numel(files_a) ~= numel(files_b)
        warning('Number of a.csv and b.csv files does not match in %s', cartella); 
        continue 
    end
    
    [~, idx_a] = sort_nat({files_a.name}); 
    [~, idx_b] = sort_nat({files_b.name}); 
    files_a = files_a(idx_a); 
    files_b = files_b(idx_b); 
    
    dati_per_cartella.(nome_campo) = {}; 
    
    for j = 1:numel(files_a)
        file_a_path = fullfile(files_a(j).folder, files_a(j).name); 
        file_b_path = fullfile(files_b(j).folder, files_b(j).name); 
        disp(['Processing pair: ', files_a(j).name, ' + ', files_b(j).name]); 
        
        % === READING FILE B (calibration) ===
        opts_b = detectImportOptions(file_b_path); opts_b.VariableNamingRule = 'preserve'; 
        fid_b = fopen(file_b_path,'r'); prima_riga = fgetl(fid_b); fclose(fid_b);
        pattern = ',(?=\D)|(?<=\D),'; 
        Calibrazione = regexp(prima_riga, pattern, 'split'); 
        Calibrazione = strrep(Calibrazione, ',', '.'); 
        
        % === READING FILE A (measurement) ===
        opts_a = detectImportOptions(file_a_path); opts_a.VariableNamingRule = 'preserve'; 
        C = readtable(file_a_path, opts_a); 
        
        fid_a = fopen(file_a_path,'r'); for r=1:5, riga=fgetl(fid_a); end; fclose(fid_a);
        Start_time = regexp(riga, pattern, 'split'); Start_time = strrep(Start_time, ',', '.'); 
        
        % === BIOIMPEDANCE CALCULATION ===
        Iraw = table2array(C(:,3)); 
        Qraw = table2array(C(:,4)); 
        time_vector = ((table2array(C(:,1))) - str2double(Start_time{2})) / 1000; 
        time_vector = time_vector(~isnan(time_vector)); 
        
        I_coef = str2double(Calibrazione{2}); Q_coef = str2double(Calibrazione{4});
        I_phase_coef = str2double(Calibrazione{6}); Q_phase_coef = str2double(Calibrazione{8});
        I_offset = str2double(Calibrazione{10}); Q_offset = str2double(Calibrazione{12});
        
        I_cal_real = @(i)(((i - I_offset) / I_coef) * cosd(I_phase_coef));
        I_cal_imag = @(i)(((i - I_offset) / I_coef) * sind(I_phase_coef));
        Q_cal_real = @(q)(((q - Q_offset) / Q_coef) * sind(Q_phase_coef));
        Q_cal_imag = @(q)(((q - Q_offset) / Q_coef) * cosd(Q_phase_coef));
        
        I_cal_r = arrayfun(I_cal_real, Iraw);
        I_cal_i = arrayfun(I_cal_imag, Iraw);
        Q_cal_r = arrayfun(Q_cal_real, Qraw);
        Q_cal_i = arrayfun(Q_cal_imag, Qraw);
        
        Load_real = (I_cal_r - Q_cal_r) * str2double(Calibrazione{14});
        Load_imag = (I_cal_i + Q_cal_i) * str2double(Calibrazione{14});
        
        Load_real = Load_real(~isnan(Load_real));
        Load_imag = Load_imag(~isnan(Load_imag));
        
        Load_mag = sqrt(Load_real.^2 + Load_imag.^2); 
        Load_angle = rad2deg(atan2(Load_imag, Load_real)); 
        time_vector = time_vector(1:length(Load_real)); 
        
        % === TABLE CREATION ===
        % ADDING REAL/IMAG COLUMNS (for future compatibility)
        Bioimpedance = table(time_vector, Load_real, Load_imag, Load_mag, Load_angle, ...
            'VariableNames', {'Time_s','Real','Imag','Magnitude','Phase'}); 
        
        dati_per_cartella.(nome_campo){end+1} = Bioimpedance;
    end
end
disp('Processing completed'); 

% Save the structure containing processed data
if ~exist(salvaPath, 'dir'), mkdir(salvaPath); end
save(fullfile(salvaPath,'dati_Perugia_26_11.mat'),'dati_per_cartella');

% === Plots for 1-minute raw data ===
folder_save = 'C:\Users\Gabriele Rellini\OneDrive\Desktop\Università\2. Magistrale\Tesi\Sperimentazione\4 - Misure Perugia 26-11-25\Elaborazione\Grafici dati grezzi';
if ~exist(folder_save, 'dir'), mkdir(folder_save); end
cartelle = fieldnames(dati_per_cartella);
cat_16_9 = {}; cat_4_5 = {}; cat_ctrl = {};

for k = 1:numel(cartelle)
    nome_cartella = cartelle{k};
    signals = dati_per_cartella.(nome_cartella);
    if numel(signals) == 3 
        segnali_1min = {signals{1}, signals{3}};
    elseif numel(signals) == 2 
        segnali_1min = signals;
    else
        warning('Folder with unexpected number of signals: %s', nome_cartella); continue
    end
    
    if contains(nome_cartella,'16_9UApp','IgnoreCase',true)
        cat_16_9{end+1} = segnali_1min;
    elseif contains(nome_cartella,'4_5UApp','IgnoreCase',true)
        cat_4_5{end+1} = segnali_1min;
    elseif contains(nome_cartella,'Controllo','IgnoreCase',true)
        cat_ctrl{end+1} = segnali_1min;
    else
        warning('Unclassified folder: %s', nome_cartella);
    end
end

plot_1min(cat_16_9, folder_save, 'Stimulation at 16 Hz and 16.9 uApp');
plot_1min(cat_4_5, folder_save, 'Stimulation at 16 Hz and 4.5 uApp');
plot_1min(cat_ctrl, folder_save, 'Controls');

% === 20-minute Plots ===
cat_16_9_20min = {}; cat_4_5_20min = {};
for k = 1:numel(cartelle)
    nome_cartella = cartelle{k};
    signals = dati_per_cartella.(nome_cartella);
    if numel(signals) == 3
        if contains(nome_cartella,'16_9UApp','IgnoreCase',true)
            cat_16_9_20min{end+1} = signals{2};
        elseif contains(nome_cartella,'4_5UApp','IgnoreCase',true)
            cat_4_5_20min{end+1} = signals{2};
        end
    end
end

plot_20min(cat_16_9_20min, folder_save, 'Stimulation at 16 Hz and 16.9 uApp - 20min');
plot_20min(cat_4_5_20min, folder_save, 'Stimulation at 16 Hz and 4.5 uApp - 20min');

function plot_1min(cat_cell, folder_save, nome_png)
    n = numel(cat_cell);  
    nrows = ceil(sqrt(n)); ncols = ceil(n/nrows);  
    % Magnitude Subplot
    f = figure('Visible', 'off');
    for i = 1:n
        segnali = cat_cell{i};
        subplot(nrows, ncols, i); hold on; legend_labels = {};
        for j = 1:numel(segnali)
            if j == 1, col='b'; lbl='Pre'; else, col='r'; lbl='Post'; end
            plot(segnali{j}.Time_s, segnali{j}.Magnitude, col, 'LineWidth', 1);
            legend_labels{end+1} = lbl;
        end
        legend(legend_labels); xlabel('Time [s]'); ylabel('Mag [\Omega]'); title(sprintf('Folder: %s', nome_png));
    end
    saveas(f, fullfile(folder_save, [nome_png, '_subplot_mag.png'])); close(f);
    
    % Phase Subplot
    f = figure('Visible', 'off');
    for i = 1:n
        segnali = cat_cell{i};
        subplot(nrows, ncols, i); hold on; legend_labels = {};
        for j = 1:numel(segnali)
            if j == 1, col='b'; lbl='Pre'; else, col='r'; lbl='Post'; end
            plot(segnali{j}.Time_s, segnali{j}.Phase, col, 'LineWidth', 1);
            legend_labels{end+1} = lbl;
        end
        legend(legend_labels); xlabel('Time [s]'); ylabel('Phase [°]'); title(sprintf('Folder: %s', nome_png));
    end
    saveas(f, fullfile(folder_save, [nome_png, '_subplot_phase.png'])); close(f);
    
    % Single Magnitude Plots
    for i = 1:n
        segnali = cat_cell{i}; f = figure('Visible', 'off'); hold on;
        for j = 1:numel(segnali)
            if j == 1, col='b'; else, col='r'; end
            plot(segnali{j}.Time_s, segnali{j}.Magnitude, col);
        end
        saveas(f, fullfile(folder_save, sprintf('%s_%d_mag.png', nome_png, i))); close(f);
    end
    % Single Phase Plots
    for i = 1:n
        segnali = cat_cell{i}; f = figure('Visible', 'off'); hold on;
        for j = 1:numel(segnali)
            if j == 1, col='b'; else, col='r'; end
            plot(segnali{j}.Time_s, segnali{j}.Phase, col);
        end
        saveas(f, fullfile(folder_save, sprintf('%s_%d_phase.png', nome_png, i))); close(f);
    end
end

function plot_20min(cat_cell, folder_save, nome_png)
    n = numel(cat_cell);  
    for i = 1:n
        segnale = cat_cell{i};
        idx = segnale.Time_s >= 5;  
        f_mag = figure('Visible', 'off');
        plot(segnale.Time_s(idx), segnale.Magnitude(idx), 'b');
        xlabel('Time [s]'); ylabel('Mag [\Omega]'); title(sprintf('Stimulation %d - %s', i, nome_png));
        saveas(f_mag, fullfile(folder_save, sprintf('%s_Sample%d_magnitude.png', nome_png, i))); close(f_mag);
        
        f_phase = figure('Visible', 'off');
        plot(segnale.Time_s(idx), segnale.Phase(idx), 'r');
        xlabel('Time [s]'); ylabel('Phase [°]'); title(sprintf('Stimulation %d - %s', i, nome_png));
        saveas(f_phase, fullfile(folder_save, sprintf('%s_Sample%d_phase.png', nome_png, i))); close(f_phase);
    end
end

function [sorted_strings, sorted_idx] = sort_nat(cellstr_array)
    expr = '(\d+)';
    tokens = regexp(cellstr_array, expr, 'match');
    nums = cellfun(@(x) str2double(x{1}), tokens);
    [~, sorted_idx] = sort(nums);
    sorted_strings = cellstr_array(sorted_idx);
end
