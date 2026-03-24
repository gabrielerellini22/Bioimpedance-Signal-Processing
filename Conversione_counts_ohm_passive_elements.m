clear
close all
clc
% Folder containing all CSV files (both "a" and "b")
cartella = 'C:\Users\Gabriele Rellini\OneDrive\Desktop\Università\2. Magistrale\Tesi\Sperimentazione\Misure su resistori noti con configurazione per elettrodi';
% Find all CSV files in the folder
files_csv = dir(fullfile(cartella, '*.csv'));

% Filter calibration files (ending with "b") and measurement files (ending with "a")
files_a = files_csv(endsWith({files_csv.name}, 'a.csv'));
files_b = files_csv(endsWith({files_csv.name}, 'b.csv'));

% Check: files must be in pairs
if numel(files_a) ~= numel(files_b)
    error('Number of a.csv and b.csv files does not match!');
end
% Sort files alphabetically for safety
files_a = sort_nat({files_a.name});
files_b = sort_nat({files_b.name});
% Global variable for calibration
fattore_calibrazione_global = NaN;

% Loop through file pairs
for j = 1:numel(files_b)
    % Full paths for files a and b
    file_b_path = fullfile(cartella, files_b{j});
    file_a_path = fullfile(cartella, files_a{j});
    disp(['Processing pair: ', files_a{j}, ' + ', files_b{j}]);

    % === READING FILE B (calibration) ===
    opts_b = detectImportOptions(file_b_path);
    opts_b.VariableNamingRule = 'preserve';
    T = readtable(file_b_path, opts_b);
    % Extract the first row (text) for calibration coefficients
    fid_b = fopen(file_b_path, 'r');
    prima_riga = fgetl(fid_b);
    fclose(fid_b);
    pattern = ',(?=\D)|(?<=\D),';
    Calibrazione = regexp(prima_riga, pattern, 'split');
    Calibrazione = strrep(Calibrazione, ',', '.');

    % === READING FILE A (measurement) ===
    opts_a = detectImportOptions(file_a_path);
    opts_a.VariableNamingRule = 'preserve';
    C = readtable(file_a_path, opts_a);
    % Extract row 5: start time
    fid_a = fopen(file_a_path, 'r');
    for r = 1:5
        riga = fgetl(fid_a);
    end
    fclose(fid_a);
    Start_time = regexp(riga, pattern, 'split');
    Start_time = strrep(Start_time, ',', '.');

    % === BIOIMPEDANCE CALCULATION ===
    Iraw = table2array(C(:,3));
    Qraw = table2array(C(:,4));
    time_vector = ((table2array(C(:,1))) - str2double(Start_time{2})) / 1000;
    time_vector = time_vector(~isnan(time_vector));
    I_coef = str2double(Calibrazione{2});
    Q_coef = str2double(Calibrazione{4});
    I_phase_coef = str2double(Calibrazione{6});
    Q_phase_coef = str2double(Calibrazione{8});
    I_offset = str2double(Calibrazione{10});
    Q_offset = str2double(Calibrazione{12});
    I_cal_real = @(i)(((i - I_offset) / I_coef) * cosd(I_phase_coef));
    I_cal_imag = @(i)(((i - I_offset) / I_coef) * sind(I_phase_coef));
    Q_cal_real = @(q)(((q - Q_offset) / Q_coef) * sind(Q_phase_coef));
    Q_cal_imag = @(q)(((q - Q_offset) / Q_coef) * cosd(Q_phase_coef));
    I_cal_r = arrayfun(I_cal_real, Iraw);
    I_cal_i = arrayfun(I_cal_imag, Iraw);
    Q_cal_r = arrayfun(Q_cal_real, Qraw);
    Q_cal_i = arrayfun(Q_cal_imag, Qraw);
    Load_real = (I_cal_r - Q_cal_r) .* str2double(Calibrazione{14});
    Load_imag = (I_cal_i + Q_cal_i) .* str2double(Calibrazione{14});
    Load_mag = sqrt(Load_real.^2 + Load_imag.^2);
    Load_angle = atan2(Load_imag, Load_real);
    % Remove NaN
    Load_real = Load_real(~isnan(Load_real));
    Load_imag = Load_imag(~isnan(Load_imag));
    Load_mag = Load_mag(~isnan(Load_mag));
    Load_angle = Load_angle(~isnan(Load_angle));

    % === STIMULUS CURRENT CALCULATION SECTION ===
    fprintf('Calculating stimulus current...\n');
    VREF = 1.0; % ADC reference voltage (typical 1 V)
    % --- Automatic parameter extraction from filename ---
    nomefile = files_b{j};
    
    % === CAPACITOR VALUE EXTRACTION ===
    token_cond = regexp(nomefile, 'C\s*da\s*([\d,\.]+)\s*(uF|nF|pF|mF)', 'tokens', 'once');
    if ~isempty(token_cond)
        valore_C_str = strrep(token_cond{1}, ',', '.'); % number as string
        valore_C = str2double(valore_C_str);
        unita_C = token_cond{2};
        % Build readable string
        Capacita_descrizione = sprintf('%g %s', valore_C, unita_C);
        % Convert to Farad
        switch lower(unita_C)
            case 'uf'
                C_farad = valore_C * 1e-6;
            case 'nf'
                C_farad = valore_C * 1e-9;
            case 'pf'
                C_farad = valore_C * 1e-12;
            case 'mf'
                C_farad = valore_C * 1e-3;
            otherwise
                C_farad = NaN;
        end
    else
        Capacita_descrizione = 'N/A';
        C_farad = NaN;
        valore_C = NaN;
        unita_C = 'N/A';
    end
    
    % Frequency (e.g., "12,5 kHz" or "16 Hz")
    token_freq = regexp(nomefile, '([\d,\.]+)\s*(kHz|Hz)', 'tokens', 'once');
    if ~isempty(token_freq)
        valore = str2double(strrep(token_freq{1}, ',', '.'));
        unita = token_freq{2};
        if strcmpi(unita, 'kHz')
            frequenza_Hz = valore * 1e3;
        else
            frequenza_Hz = valore;
        end
    else
        frequenza_Hz = NaN;
    end
    
    % Gain (e.g., "2 Gain")
    token_gain = regexp(nomefile, '(\d+)\s*Gain', 'tokens', 'once');
    BIOZ_GAIN = str2double(token_gain{1});
    
    % === Nominal voltage, current and resistance table ===
    % Nominal RMS currents (in µArms)
    corr_nom = [ ...
        0.016, 0.032, 0.080, 0.160, ...
        0.320, 0.640, 1.6, 3.2, ... 
        6.4, 12.8, 32, 64,...
        128, 256, 640, 1280];
    % Corresponding nominal voltages (in mVrms)
    volt_nom = [ ...
        35.4, 70.7, 177, 354, ...
        35.4, 70.7, 177, 354, ...
        35.4, 70.7, 177, 354, ...
        35.4, 70.7, 177, 354];
    % Corresponding nominal resistances (in Ohm)
    res_nom = [ ...
        2210e3, 2210e3, 2210e3, 2210e3, ...
        110.5e3, 110.5e3, 110.5e3, 110.5e3, ...
        5.525e3, 5.525e3, 5.525e3, 5.525e3, ...
        276.25, 276.25, 276.25, 276.25];
    
    % Nominal current (e.g., "128 uArms")
    token_corrente = regexp(nomefile, '([\d,\.]+)\s*uArms', 'tokens', 'once');
    corrente_nominale_uA = str2double(strrep(token_corrente{1}, ',', '.'));
   
    % Find corresponding row in the table
    [~, idx_corr] = min(abs(corr_nom - corrente_nominale_uA));
    V_nominale_mV = volt_nom(idx_corr);
    R_nominale_ohm = res_nom(idx_corr);
   
    % === Automatic calculation of nominal voltage and resistance ===
    V_nominale_V = V_nominale_mV / 1000; % in Volt RMS
    I_nominale_A = corrente_nominale_uA * 1e-6; % in Ampere RMS
    R_calcolata = V_nominale_V / I_nominale_A;
    
    % Cut-off frequency fc and effective current
    if ~isnan(R_nominale_ohm) && ~isnan(C_farad) && ~isnan(frequenza_Hz)
        f_taglio = 1/(2*pi*R_nominale_ohm*C_farad); % Hz
        Xc = 1/(2*pi*frequenza_Hz*C_farad);
        Z_mod = sqrt(R_nominale_ohm.^2 + Xc.^2);
    else
        f_taglio = NaN;
        Z_mod = NaN;
    end
    Corrente_effettiva = ((V_nominale_mV/1000)/Z_mod)*1e6;
    
    % Known measurement resistor (e.g., "100 ohm")
    % === RESISTOR AND TOLERANCE EXTRACTION ===
    token_resistore = regexp(nomefile, '([\d.,]+)\s*ohm', 'tokens', 'once');
    Rmisura = str2double(strrep(token_resistore{1}, ',', '.'));
    % Known tolerance definitions
    if abs(Rmisura - 30.1) < 1e-3
        Tolleranza_tipo = '±1%';
        Tolleranza_valore = Rmisura * 0.01; % ohm
    elseif abs(Rmisura - 480) < 1e-3
        Tolleranza_tipo = '±5%';
        Tolleranza_valore = Rmisura * 0.05; % ohm
    elseif abs(Rmisura - 100) < 1e-3
        Tolleranza_tipo = '±5%';
        Tolleranza_valore = Rmisura * 0.05; % ohm
    else
        Tolleranza_tipo = 'N/A';
        Tolleranza_valore = NaN;
    end
    % Tolerance interval calculation
    limite_inferiore = Rmisura - Tolleranza_valore;
    limite_superiore = Rmisura + Tolleranza_valore;
    % Check if mean Magnitude is within tolerance
    if mean(Load_mag, 'omitnan') >= limite_inferiore && mean(Load_mag, 'omitnan') <= limite_superiore
        Entro_tolleranza = 'Yes';
    else
        Entro_tolleranza = 'No';
    end
    % Create combined string like "480 Ω ±5%"
    Resistore_descrizione = sprintf('%.3g Ω %s', Rmisura, Tolleranza_tipo);
    % Resistor used in calibration (e.g., "calibrated on 600")
    token_Rcal = regexp(nomefile, 'calibrato su\s*(\d+)', 'tokens', 'once');
    if ~isempty(token_Rcal)
        Rcalibrazione = str2double(token_Rcal{1});
    else
        Rcalibrazione = 600; % default value
    end
    
    % --- Calculate actual current using conversion coefficient ---
    Iraw_no_nan = Iraw(~isnan(Load_mag));
    Qraw_no_nan = Qraw(~isnan(Load_mag));
    ADC_counts_mag = sqrt(Iraw_no_nan.^2 + Qraw_no_nan.^2);
    I_calc_pk = VREF/(2^19*BIOZ_GAIN*(2/pi)*str2double(Calibrazione{14}));
    I_calc_uA_rms = I_calc_pk*1e6/sqrt(2); % Consistency with nominal current means software sets current correctly
    Vpk = (ADC_counts_mag * VREF)/(2^19*BIOZ_GAIN*(2/pi));
    Vpp = mean((Vpk*2)*1e3); %[mV]
    Corrente_misurata = ((Vpp*1e-3/(2*sqrt(2)))/Rmisura)*1e6;
    
    % --- Summary print ---
    fprintf('---------------------------------------------------\n');
    fprintf('File: %s\n', files_b{j});
    
    % === DRIVER COMPLIANCE CHECK (DRV_OOR) ===
    AVDD = 1.8; % [V] MAX30009 nominal supply
    VMID_TX = 0.81; % [V] Typical common mode voltage
    % --- DRV_OOR Internal Limits (simulating register 0x01 state) ---
    V_OOR_MIN = 0.2; % [V] Minimum threshold for DRV_OOR
    V_OOR_MAX = AVDD - 0.2; % [V] Maximum threshold for DRV_OOR
    % Calculate instantaneous voltage across measurement resistor
    Vload = (corrente_nominale_uA * 1e-6)*2*sqrt(2) * Rmisura; % [V p-p]
    Vload_PEAK = Vload / 2; % [V Peak]
    % Calculate signal instantaneous peaks around VMID_TX
    V_Istantaneo_Max = VMID_TX + Vload_PEAK;
    V_Istantaneo_Min = VMID_TX - Vload_PEAK;
    % Check if peaks exceed OOR limits
    DRV_OOR = (V_Istantaneo_Min < V_OOR_MIN) | (V_Istantaneo_Max > V_OOR_MAX);
    % --- Reporting ---
    if any(DRV_OOR)
        fprintf('⚠️ WARNING: driver out of compliance! ');
        fprintf('Instantaneous peaks [min-max] = [%.3f V – %.3f V] (limits %.2f–%.2f V)\n', ...
            V_Istantaneo_Min, V_Istantaneo_Max, V_OOR_MIN, V_OOR_MAX);
    else
        fprintf('✅ Driver within compliance window (%.2f–%.2f V)\n', V_OOR_MIN, V_OOR_MAX);
    end
    
    % === 1000 mV RECEPTION LIMIT CHECK ===
    VRX_LIMIT_mV = 1000; % [mV RMS]
    VRX_LIMIT_V = (VRX_LIMIT_mV / 1000) / BIOZ_GAIN;
    V_rx_theoretical = (corrente_nominale_uA * 1e-6) * Rcalibrazione; % [Theoretical V RMS]
    if V_rx_theoretical > VRX_LIMIT_V
        fprintf('⚠️ WARNING: product I_nominal * Rcalibration = %.3f V RMS > %.3f V RMS (limit)\n', ...
            V_rx_theoretical, VRX_LIMIT_V);
        RX_over = true;
    else
        fprintf('✅ Theoretical reception voltage within limits (%.3f V RMS ≤ %.3f V RMS)\n', ...
            V_rx_theoretical, VRX_LIMIT_V);
        RX_over = false;
    end
    fprintf('---------------------------------------------------\n\n');
    
    % --- Mean and standard deviation calculation for recorded magnitude and phase ---
    mag_clean = Load_mag;
    mag_clean(isoutlier(mag_clean, 'quartiles')) = NaN;
    phase_clean = Load_angle;
    phase_clean(isoutlier(phase_clean, 'quartiles')) = NaN;
    Magnitude_media = mean(mag_clean, 'omitnan');
    Magnitude_std = std(mag_clean, 'omitnan');
    Errore_assoluto = abs(Magnitude_media - Rmisura); % [ohm]
    Errore_percentuale = (Errore_assoluto / Rmisura) * 100; % [%]
    fprintf('Absolute Error on Magnitude: %.3f Ω\n', Errore_assoluto);
    fprintf('Percentage Error on Magnitude: %.2f %%\n', Errore_percentuale);
    Fase_media = mean(phase_clean*(180/pi), 'omitnan');
    Fase_std = std(phase_clean*(180/pi), 'omitnan');
    % --- Addition to summary table with standard deviation columns ---
    riga_tabella = table(frequenza_Hz, corrente_nominale_uA, Corrente_effettiva, I_calc_uA_rms, ...
    {Capacita_descrizione}, f_taglio, Fase_media, Fase_std, Magnitude_media, Magnitude_std, ...
    Errore_assoluto, Errore_percentuale, ...
    {Resistore_descrizione}, {Entro_tolleranza}, Rcalibrazione, ...
    'VariableNames', {'Frequency_Hz', 'Nominal_current_uArms','Effective_current_uArms','Current_calculated_from_C_uArms', ...
    'Capacitor_used','Cutoff_frequency', 'Mean_Phase_deg', 'Std_Phase_deg','Mean_Magnitude_ohm', 'Std_Magnitude_ohm', ...
    'Absolute_error_ohm','Percentage_error','Known_Resistor','Within_tolerance','Calibration_Resistor'});
    if j == 1
        Tabella_Risultati = riga_tabella;
    else
        Tabella_Risultati = [Tabella_Risultati; riga_tabella];
    end
    % Sort by calibration resistor value
    % Tabella_Risultati = sortrows(Tabella_Risultati, 'Condensatore usato', 'descend');
    % === TABLE CREATION ===
    Bioimpedance = table(time_vector, Load_real, Load_imag, Load_mag, Load_angle, ...
        'VariableNames', {'Time (s)', 'BiozI (ohm)', 'BioZQ (ohm)', 'Magnitude (ohm)', 'Phase (°)'});
    % === Save to workspace with dynamic name ===
    nome_var = sprintf('Bioimpedance_%d', j);
    assignin('base', nome_var, Bioimpedance);
end

% %% === BARPLOT COMPARISON OF CAPACITORS FOR EACH RESISTANCE (separate means and std) ===
% cartella_salvataggio = 'C:\Users\Gabriele Rellini\OneDrive\Desktop\Università\2. Magistrale\Tesi\Sperimentazione\Misure su resistori noti con configurazione per elettrodi\Grafici MATLAB';
% 
% res_uniche = unique(Tabella_Risultati.Known_Resistor);
% cond_uniche = unique(Tabella_Risultati.("Capacitor_used"));
% colori_cond = lines(length(cond_uniche));
% 
% for r = 1:length(res_uniche)
%     freq = unique(Tabella_Risultati.Frequency_Hz);
% 
%     mag_mean_mat = nan(length(freq), length(cond_uniche));
%     mag_std_mat  = nan(length(freq), length(cond_uniche));
%     phase_mean_mat = nan(length(freq), length(cond_uniche));
%     phase_std_mat  = nan(length(freq), length(cond_uniche));
% 
%     for c = 1:length(cond_uniche)
%         idx = strcmp(Tabella_Risultati.Known_Resistor, res_uniche{r}) & ...
%               strcmp(Tabella_Risultati.("Capacitor_used"), cond_uniche{c});
%         for f = 1:length(freq)
%             idx_f = idx & Tabella_Risultati.Frequency_Hz == freq(f);
%             if any(idx_f)
%                 mag_mean_mat(f,c)   = Tabella_Risultati.Mean_Magnitude_ohm(idx_f);
%                 mag_std_mat(f,c)    = Tabella_Risultati.Std_Magnitude_ohm(idx_f);
%                 phase_mean_mat(f,c) = Tabella_Risultati.Mean_Phase_deg(idx_f);
%                 phase_std_mat(f,c)  = Tabella_Risultati.Std_Phase_deg(idx_f);
%             end
%         end
%     end
% 
%     % --- Frequency labels with decimals always visible
%     freq_labels = arrayfun(@(x) sprintf('%.2f Hz', x), freq, 'UniformOutput', false);
% 
%     %% --- Mean Magnitude Barplot ---
%     fig = figure('Name', sprintf('Mean Magnitudes - R=%s', res_uniche{r}), 'NumberTitle','off');
%     b = bar(mag_mean_mat,'grouped');
%     for c = 1:length(cond_uniche)
%         b(c).FaceColor = colori_cond(c,:);
%     end
%     hold on;
%     [ngroups, nbars] = size(mag_mean_mat);
%     groupwidth = min(0.8, nbars/(nbars+1.5));
%     for i = 1:nbars
%         x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth/(2*nbars);
%         errorbar(x, mag_mean_mat(:,i), mag_std_mat(:,i), 'k.', 'LineWidth',1.5);
%     end
%     xticks(1:length(freq));
%     xticklabels(freq_labels);
%     xlabel('Frequency (Hz)'); ylabel('Mean Magnitude (Ω)');
%     title(sprintf('Mean Magnitude - R=%s', res_uniche{r}));
%     ylim([0 max(mag_mean_mat(:)+mag_std_mat(:))*1.1]);
%     legend(cond_uniche, 'Location','best'); grid on; box on; hold off;
% 
%     % --- Saving ---
%     saveas(fig, fullfile(cartella_salvataggio, sprintf('Mean_Magnitude_R%s.png', res_uniche{r})));
%     close(fig);
% 
%     %% --- Mean Phase Barplot ---
%     fig = figure('Name', sprintf('Mean Phases - R=%s', res_uniche{r}), 'NumberTitle','off');
%     b = bar(phase_mean_mat,'grouped');
%     for c = 1:length(cond_uniche)
%         b(c).FaceColor = colori_cond(c,:);
%     end
%     hold on;
%     [ngroups, nbars] = size(phase_mean_mat);
%     groupwidth = min(0.8, nbars/(nbars+1.5));
%     for i = 1:nbars
%         x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth/(2*nbars);
%         errorbar(x, phase_mean_mat(:,i), phase_std_mat(:,i), 'k.', 'LineWidth',1.5);
%     end
%     xticks(1:length(freq));
%     xticklabels(freq_labels);
%     xlabel('Frequency (Hz)'); 
%     ylabel('Mean Phase (°)'); 
%     title(sprintf('Mean Phase - R=%s', res_uniche{r}));
%     ylim([min(phase_mean_mat(:)-phase_std_mat(:))*1.2, max(phase_mean_mat(:)+phase_std_mat(:))*1.2]); 
%     legend(cond_uniche, 'Location','best'); grid on; box on; hold off;
% 
%     % --- Saving ---
%     saveas(fig, fullfile(cartella_salvataggio, sprintf('Mean_Phase_R%s.png', res_uniche{r})));
%     close(fig);
% 
%     %% --- Magnitude STD Barplot ---
%     fig = figure('Name', sprintf('Magnitude Variability - R=%s', res_uniche{r}), 'NumberTitle','off');
%     b = bar(mag_std_mat,'grouped');
%     for c = 1:length(cond_uniche)
%         b(c).FaceColor = colori_cond(c,:);
%     end
%     xticks(1:length(freq));
%     xticklabels(freq_labels);
%     xlabel('Frequency (Hz)'); ylabel('Magnitude Standard Deviation (Ω)');
%     title(sprintf('Magnitude Variability - R=%s', res_uniche{r}));
%     ylim([0 max(mag_std_mat(:))*1.2]);
%     legend(cond_uniche, 'Location','best'); grid on; box on;
% 
%     % --- Saving ---
%     saveas(fig, fullfile(cartella_salvataggio, sprintf('Magnitude_STD_R%s.png', res_uniche{r})));
%     close(fig);
% 
%     %% --- Phase STD Barplot ---
%     fig = figure('Name', sprintf('Phase Variability - R=%s', res_uniche{r}), 'NumberTitle','off');
%     b = bar(phase_std_mat,'grouped');
%     for c = 1:length(cond_uniche)
%         b(c).FaceColor = colori_cond(c,:);
%     end
%     xticks(1:length(freq));
%     xticklabels(freq_labels);
%     xlabel('Frequency (Hz)'); ylabel('Phase Standard Deviation (°)');
%     title(sprintf('Phase Variability - R=%s', res_uniche{r}));
%     ylim([0 max(phase_std_mat(:))*1.2]); 
%     legend(cond_uniche, 'Location','best'); grid on; box on;
% 
%     % --- Saving ---
%     saveas(fig, fullfile(cartella_salvataggio, sprintf('Phase_STD_R%s.png', res_uniche{r})));
%     close(fig);
% end
% 
% %% === PLOT MAGNITUDE AND PHASE TREND OVER TIME ===
% res_uniche = unique(Tabella_Risultati.Known_Resistor);
% cond_uniche = unique(Tabella_Risultati.("Capacitor_used"));
% colori_freq = lines(length(unique(Tabella_Risultati.Frequency_Hz)));
% 
% for r = 1:length(res_uniche)
%     for c = 1:length(cond_uniche)
%         % Find all file indices matching this combination
%         idx_files = find(strcmp(Tabella_Risultati.Known_Resistor, res_uniche{r}) & ...
%                          strcmp(Tabella_Risultati.("Capacitor_used"), cond_uniche{c}));
%         if isempty(idx_files)
%             continue
%         end
% 
%         % Magnitude Figure
%         fig_mag = figure('Name', sprintf('Magnitude vs Time - R=%s, C=%s', res_uniche{r}, cond_uniche{c}), 'NumberTitle', 'off'); hold on;
%         % Phase Figure
%         fig_phase = figure('Name', sprintf('Phase vs Time - R=%s, C=%s', res_uniche{r}, cond_uniche{c}), 'NumberTitle', 'off'); hold on;
% 
%         freq_uniche = Tabella_Risultati.Frequency_Hz(idx_files);
%         [freq_uniche, freq_sort_idx] = sort(freq_uniche);  % sort frequencies
%         idx_files = idx_files(freq_sort_idx);
% 
%         for k = 1:length(idx_files)
%             nome_var = sprintf('Bioimpedance_%d', idx_files(k));
%             Bioz = evalin('base', nome_var);
% 
%             % Remove Magnitude and Phase outliers
%             mag_clean = Bioz.("Magnitude (ohm)");
%             mag_clean(isoutlier(mag_clean, 'quartiles')) = NaN;
% 
%             phase_clean = rad2deg(Bioz.("Phase (°)"));
%             phase_clean(isoutlier(phase_clean, 'quartiles')) = NaN;
% 
%             % Magnitude over time
%             figure(fig_mag);
%             plot(Bioz.("Time (s)"), mag_clean, 'LineWidth', 1.5, 'Color', colori_freq(k,:));
% 
%             % Phase over time
%             figure(fig_phase);
%             plot(Bioz.("Time (s)"), phase_clean, 'LineWidth', 1.5, 'Color', colori_freq(k,:));
%         end
% 
%         % Setup Magnitude figure
%         figure(fig_mag);
%         xlabel('Time (s)'); ylabel('Magnitude (Ω)');
%         title(sprintf('Magnitude vs Time - R=%s, C=%s', res_uniche{r}, cond_uniche{c}));
%         legend(arrayfun(@(f) sprintf('%.2f Hz', f), freq_uniche, 'UniformOutput', false), 'Location', 'best');
%         grid on; box on;
%         saveas(fig_mag, fullfile(cartella_salvataggio, sprintf('Magnitude_Time_R%s_C%s.png', res_uniche{r}, cond_uniche{c})));
% 
%         % Setup Phase figure
%         figure(fig_phase);
%         xlabel('Time (s)'); ylabel('Phase (°)');
%         title(sprintf('Phase vs Time - R=%s, C=%s', res_uniche{r}, cond_uniche{c}));
%         legend(arrayfun(@(f) sprintf('%.2f Hz', f), freq_uniche, 'UniformOutput', false), 'Location', 'best');
%         grid on; box on;
%         saveas(fig_phase, fullfile(cartella_salvataggio, sprintf('Phase_Time_R%s_C%s.png', res_uniche{r}, cond_uniche{c})));
% 
%         close(fig_mag); close(fig_phase);  % close figures to free memory
%     end
% end
% Natural sorting function (e.g.: 1a, 2a, 10a, ...)

function sorted = sort_nat(files)
[~, idx] = sort_nat_helper(files);
sorted = files(idx);
end

function [sorted_strings, sorted_idx] = sort_nat_helper(cellstr_array)
expr = '(\d+)';
tokens = regexp(cellstr_array, expr, 'match');
nums = cellfun(@(x) str2double(x{1}), tokens);
[~, sorted_idx] = sort(nums);
sorted_strings = cellstr_array(sorted_idx);
end