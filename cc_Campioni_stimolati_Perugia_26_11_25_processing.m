clear
close all
clc

% Load raw data
dataPath = 'C:\Users\Gabriele Rellini\OneDrive\Desktop\Università\2. Magistrale\Tesi\Sperimentazione\4 - Misure Perugia 26-11-25\Elaborazione';
fileName = 'dati_Perugia_26_11.mat';
load(fullfile(dataPath, fileName), 'dati_per_cartella');

% Initialize struct for processed data
dati_proc = struct;

% Outlier detection parameters
mad_factor = 4; 

% FIR filter parameters
fc = 5;       % Cut-off frequency at 5 Hz
order = 200;  % FIR filter order

%% 1. SIGNAL PROCESSING
% Loop through all biological sample folders
salvaPath = 'C:\Users\Gabriele Rellini\OneDrive\Desktop\Università\2. Magistrale\Tesi\Sperimentazione\4 - Misure Perugia 26-11-25\Elaborazione\Comparative plots raw vs processed'; 
if ~exist(salvaPath, 'dir')
    mkdir(salvaPath);
end

cartelle = fieldnames(dati_per_cartella);
for k = 1:numel(cartelle)
    nome_cartella = cartelle{k};
    signals = dati_per_cartella.(nome_cartella);
    signals_proc = cell(size(signals)); 
    
    for j = 1:numel(signals)
        tab = signals{j}; 
       
        % 1.1 Sampling frequency
        if numel(signals) == 3
            if j == 2
                fs = 32;  % stimulation
            else
                fs = 25;  % pre/post
            end
        elseif numel(signals) == 2
            fs = 25;      % controls
        else
            error('Unexpected number of signals in %s', nome_cartella);
        end
        
        % 1.2 Outlier removal from Magnitude
        mag = tab.Magnitude;
        window_size = round(5 * fs); 
        
        mediane_locali = movmedian(mag, window_size);
        mad_locali     = movmad(mag, window_size, 1);
        soglia_sup = mediane_locali + mad_factor * mad_locali;
        soglia_inf = mediane_locali - mad_factor * mad_locali;
        outliers = (mag > soglia_sup) | (mag < soglia_inf);
        mag_no_out = mag;
        mag_no_out(outliers) = mediane_locali(outliers);
        
        % 1.3 Outlier removal from Phase
        phase = tab.Phase;
        mediane_locali = movmedian(phase, window_size);
        mad_locali     = movmad(phase, window_size, 1);
        soglia_sup = mediane_locali + mad_factor * mad_locali;
        soglia_inf = mediane_locali - mad_factor * mad_locali;
        outliers = (phase > soglia_sup) | (phase < soglia_inf);
        phase_no_out = phase;
        phase_no_out(outliers) = mediane_locali(outliers);
        
        % 1.4 FIR filter design
        Wn = fc/(fs/2); 
        b = fir1(order, Wn, 'low', hamming(order+1)); 
        
        % Zero-phase filtering
        mag_filt = filtfilt(b, 1, mag_no_out);
        phase_filt = filtfilt(b, 1, phase_no_out); 
        
        % Create processed table
        tab_proc = table(tab.Time_s, mag_filt, phase_filt, ...
            'VariableNames', {'Time_s','Magnitude','Phase'});
      
        % Save in dati_proc struct
        signals_proc{j}.data = tab_proc;
        signals_proc{j}.fs = fs;
        
        % Comparative plots
        figure('Visible', 'off'); 
        subplot(2,1,1);
        plot(tab.Time_s, mag_no_out, 'b-', 'LineWidth', 1); hold on;
        plot(tab.Time_s, mag_filt, 'r-', 'LineWidth', 1);
        xlabel('Time [s]'); ylabel('Magnitude [\Omega]');
        legend('Raw without outliers','Filtered');
        if numel(signals) == 3
            if j == 1
                titolo = 'Pre-treatment';
            elseif j == 2, titolo = 'Stimulation';
            else 
                titolo = 'Post-treatment'; end
        else
            titolo = sprintf('Control %d', j);
        end
        title(sprintf('%s - %s - Magnitude', nome_cartella, titolo));
        
        subplot(2,1,2);
        plot(tab.Time_s, phase_no_out, 'b-', 'LineWidth', 1); hold on;
        plot(tab.Time_s, phase_filt, 'r-', 'LineWidth', 1);
        xlabel('Time [s]'); ylabel('Phase [°]');
        legend('Raw without outliers','Filtered');
        title(sprintf('%s - %s - Phase', nome_cartella, titolo));
        
        filename = fullfile(salvaPath, sprintf('Comparative_%s_%s.png', nome_cartella, titolo));
        saveas(gcf, filename);
        close(gcf); 
    end
    dati_proc.(nome_cartella) = signals_proc;
end

salvaPath = 'C:\Users\Gabriele Rellini\OneDrive\Desktop\Università\2. Magistrale\Tesi\Sperimentazione\4 - Misure Perugia 26-11-25\Elaborazione';
if ~exist(salvaPath, 'dir')
    mkdir(salvaPath);
end
save(fullfile(salvaPath,'dati_Perugia_26_11_proc.mat'),'dati_proc');
disp('Processing completed: processed data and comparative plots generated');

%% 2. STATISTICS AND DECAY (SLOPE)
Statistica = struct;
win_len = 60;  % Window duration in seconds
% Stimulation Decay Plots folder
salvaPathDecay = fullfile(salvaPath, 'Stimulation Decay Plots');
if ~exist(salvaPathDecay, 'dir'), mkdir(salvaPathDecay); end

for k = 1:numel(cartelle)
    nome_cartella = cartelle{k};
    signals_proc = dati_proc.(nome_cartella);  
    n_signals = numel(signals_proc);
    Statistica.(nome_cartella) = struct;
    
   if n_signals == 3
        % -------------------- STIMULATION (20 min) --------------------
        sig_stim = signals_proc{2}.data;
        fs = signals_proc{2}.fs;
        N = length(sig_stim.Time_s);
        win_size = win_len * fs;  % points for 60s window
        n_windows = floor(N / win_size);
        stim_stats = cell(1, n_windows);
        for w = 1:n_windows
            idx_range = ( (w-1)*win_size + 1 ) : (w*win_size);
            stim_stats{w} = CalcolaStatistiche(sig_stim, idx_range);
        end
        Statistica.(nome_cartella).Stim = stim_stats;
        
        % --- Magnitude Slope Calculation ---
        time_rel = sig_stim.Time_s - sig_stim.Time_s(1); 
        p = polyfit(time_rel, sig_stim.Magnitude, 1); 
        slope = p(1); % slope in Ohm/s
        Statistica.(nome_cartella).Slope_Stim = slope;
        
        % --- Complete Decay Plot (Mag + Phase) for Thesis ---
        f = figure('Visible', 'off');
        titolo_pulito = getTitoloSintetico(nome_cartella);
        
        % Subplot 1: Magnitude
        subplot(2,1,1);
        plot(time_rel/60, sig_stim.Magnitude, 'b', 'LineWidth', 0.5); hold on;
        plot(time_rel/60, polyval(p, time_rel), 'r--', 'LineWidth', 2);
        ylabel('Magnitude [\Omega]', 'Interpreter', 'tex', 'FontSize', 10);
        legend('Measured Data', 'Linear Fit', 'Location', 'best');
        grid on;
        
        % Title on first panel
        riga1 = sprintf('\\bf Stimulation Decay: %s', titolo_pulito); 
        riga2 = sprintf('\\rm Slope Mag = %.4f \\Omega/s', slope);         
        title({riga1, riga2}, 'Interpreter', 'tex', 'FontSize', 11);
        
        % Subplot 2: Phase
        subplot(2,1,2);
        plot(time_rel/60, sig_stim.Phase, 'Color', [0 0.5 0], 'LineWidth', 0.5); % Dark Green
        ylabel('Phase [°]', 'Interpreter', 'tex', 'FontSize', 10);
        xlabel('Time [min]', 'Interpreter', 'tex', 'FontSize', 10);
        grid on;
        
        % Save
        saveas(f, fullfile(salvaPathDecay, sprintf('Decay_%s.png', nome_cartella)));
        close(f);
    else
        Statistica.(nome_cartella).Slope_Stim = NaN;
    end
    
    if n_signals == 3
        % -------------------- PRE --------------------
        sig_pre = signals_proc{1}.data;
        fs = signals_proc{1}.fs;
        N = length(sig_pre.Time_s);
        mid_idx = floor(N/2);
        half_win = round((win_len*fs)/2);
        idx_range = max(1,mid_idx-half_win):min(N,mid_idx+half_win-1);
        Statistica.(nome_cartella).Pre = CalcolaStatistiche(sig_pre, idx_range);
        % -------------------- POST --------------------
        sig_post = signals_proc{3}.data;
        fs = signals_proc{3}.fs;
        N = length(sig_post.Time_s);
        mid_idx = floor(N/2);
        half_win = round((win_len*fs)/2);
        idx_range = max(1,mid_idx-half_win):min(N,mid_idx+half_win-1);
        Statistica.(nome_cartella).Post = CalcolaStatistiche(sig_post, idx_range);
    elseif n_signals == 2
        % -------------------- PRE --------------------
        sig_pre = signals_proc{1}.data;
        fs = signals_proc{1}.fs;
        N = length(sig_pre.Time_s);
        mid_idx = floor(N/2);
        half_win = round((win_len*fs)/2);
        idx_range = max(1,mid_idx-half_win):min(N,mid_idx+half_win-1);
        Statistica.(nome_cartella).Pre = CalcolaStatistiche(sig_pre, idx_range);
        % -------------------- POST --------------------
        sig_post = signals_proc{2}.data;
        fs = signals_proc{2}.fs;
        N = length(sig_post.Time_s);
        mid_idx = floor(N/2);
        half_win = round((win_len*fs)/2);
        idx_range = max(1,mid_idx-half_win):min(N,mid_idx+half_win-1);
        Statistica.(nome_cartella).Post = CalcolaStatistiche(sig_post, idx_range);
    end
end
disp('Statistics created for all folders (using STD).');

%% 3. PRE/POST HISTOGRAMS (With STD)
salvaPath = 'C:\Users\Gabriele Rellini\OneDrive\Desktop\Università\2. Magistrale\Tesi\Sperimentazione\4 - Misure Perugia 26-11-25\Elaborazione\Pre and post signals bar charts';
if ~exist(salvaPath, 'dir')
    mkdir(salvaPath);
end

cartelle = fieldnames(Statistica);
% Grouping by conditions
gruppi.controlli = cartelle(contains(lower(cartelle),'controllo'));
gruppi.stim_45 = cartelle(contains(lower(cartelle),'4_5'));
gruppi.stim_169 = cartelle(contains(lower(cartelle),'16_9'));
campiGruppi = fieldnames(gruppi);

for g = 1:numel(campiGruppi)
    gruppo_nome = campiGruppi{g};
    lista_cartelle = gruppi.(gruppo_nome);
    
    mean_mag = [];
    std_mag  = [];
    mean_phase = [];
    std_phase  = [];
    labels = {};
    
    for k = 1:numel(lista_cartelle)
        nome_cartella = lista_cartelle{k};
        stats = Statistica.(nome_cartella);
        
        % Pre (Using Std_Mag)
        mean_mag(end+1)   = stats.Pre.Mean_Mag;
        std_mag(end+1)    = stats.Pre.Std_Mag; 
        mean_phase(end+1) = stats.Pre.Mean_Phase;
        std_phase(end+1)  = stats.Pre.Std_Phase;
        
        % Post (Using Std_Mag)
        mean_mag(end+1)   = stats.Post.Mean_Mag;
        std_mag(end+1)    = stats.Post.Std_Mag;
        mean_phase(end+1) = stats.Post.Mean_Phase;
        std_phase(end+1)  = stats.Post.Std_Phase;
        
        % Sample label
        labels{end+1} = getCampioneLabel(nome_cartella);
    end
    
    % Organize data: Pre = Blue, Post = Red
    nCampioni = numel(lista_cartelle);
    X = 1:nCampioni;
    
    % --- Magnitude ---
    figure('Visible','off'); 
    hold on;
    bar(X-0.15, mean_mag(1:2:end), 0.3, 'FaceColor','b'); % Pre
    bar(X+0.15, mean_mag(2:2:end), 0.3, 'FaceColor','r'); % Post
    set(gca,'XTick',X,'XTickLabel',labels);
    ylabel('Resistance [\Omega]'); 
    title(sprintf('%s - Resistance', gruppo_nome),'Interpreter','none');
    legend({'Pre','Post'});
    grid on;
    filename = fullfile(salvaPath, sprintf('Histogram_Mag_%s.png', gruppo_nome));
    saveas(gcf, filename); close(gcf);
    
    % --- Phase ---
    figure('Visible','off'); 
    hold on;
    bar(X-0.15, mean_phase(1:2:end), 0.3, 'FaceColor','b'); % Pre
    bar(X+0.15, mean_phase(2:2:end), 0.3, 'FaceColor','r'); % Post
    errorbar(X-0.15, mean_phase(1:2:end), std_phase(1:2:end), 'k.', 'LineWidth',1); % Using std_phase
    errorbar(X+0.15, mean_phase(2:2:end), std_phase(2:2:end), 'k.', 'LineWidth',1); % Using std_phase
    set(gca,'XTick',X,'XTickLabel',labels);
    ylabel('Phase [°] (Mean \pm STD)'); 
    title(sprintf('%s - Phase', gruppo_nome),'Interpreter','none');
    legend({'Pre','Post'});
    grid on;
    filename = fullfile(salvaPath, sprintf('Histogram_Phase_%s.png', gruppo_nome));
    saveas(gcf, filename); close(gcf);
end

% 20-minute Stimulation HISTOGRAMS (With STD)
salvaPath = 'C:\Users\Gabriele Rellini\OneDrive\Desktop\Università\2. Magistrale\Tesi\Sperimentazione\4 - Misure Perugia 26-11-25\Elaborazione\Stimulation signal bar charts';
if ~exist(salvaPath, 'dir')
    mkdir(salvaPath);
end

for g = 1:numel(campiGruppi)
    gruppo_nome = campiGruppi{g};
    lista_cartelle = gruppi.(gruppo_nome);
    for k = 1:numel(lista_cartelle)
        nome_cartella = lista_cartelle{k};
        if ~isfield(Statistica.(nome_cartella),'Stim')
            continue;
        end
        
        stim_stats = Statistica.(nome_cartella).Stim;
        nW = numel(stim_stats);
        
        mean_mag = zeros(1,nW);
        std_mag  = zeros(1,nW);
        mean_phase = zeros(1,nW);
        std_phase  = zeros(1,nW);
        
        for w = 1:nW
            mean_mag(w)   = stim_stats{w}.Mean_Mag;
            std_mag(w)    = stim_stats{w}.Std_Mag; 
            mean_phase(w) = stim_stats{w}.Mean_Phase;
            std_phase(w)  = stim_stats{w}.Std_Phase; 
        end
        
        % --- Magnitude ---
        figure('Visible','off');
        bar(1:nW, mean_mag, 'FaceColor',[0.2 0.6 0.8]);
        hold on;
        errorbar(1:nW, mean_mag, std_mag, 'k.', 'LineWidth',1);
        xlabel('Window (1 min each)');
        ylabel('Magnitude [\Omega] (Mean \pm STD)');
        title(sprintf('%s - 20 min Stimulation - Magnitude', nome_cartella),'Interpreter','none');
        grid on;
        filename = fullfile(salvaPath, sprintf('Stimulation_Mag_%s.png', nome_cartella));
        saveas(gcf, filename); close(gcf);
        
        % --- Phase ---
        figure('Visible','off');
        bar(1:nW, mean_phase, 'FaceColor',[0.8 0.4 0.4]);
        hold on;
        errorbar(1:nW, mean_phase, std_phase, 'k.', 'LineWidth',1);
        xlabel('Window (1 min each)');
        ylabel('Phase [°] (Mean \pm STD)');
        title(sprintf('%s - 20 min Stimulation - Phase', nome_cartella),'Interpreter','none');
        grid on;
        filename = fullfile(salvaPath, sprintf('Stimulation_Phase_%s.png', nome_cartella));
        saveas(gcf, filename); close(gcf);
    end
end

%% 4. NORMALIZATION (Post relative to Pre - using STD)
Normalizzato = struct;
for k = 1:numel(cartelle)
    nome_cartella = cartelle{k};
    stats = Statistica.(nome_cartella);
    if isfield(stats,'Pre') && isfield(stats,'Post')
        norm_stats = struct;
        
        % Pre normalized = 1
        norm_stats.Pre_Mag   = 1;
        norm_stats.Pre_Phase = 1;
        
        % Post normalized = Post/Pre
        norm_stats.Post_Mag   = stats.Post.Mean_Mag   / stats.Pre.Mean_Mag;
        norm_stats.Post_Phase = stats.Post.Mean_Phase / stats.Pre.Mean_Phase;
        
        % Also scale STANDARD DEVIATION relative to Pre (mean)
        norm_stats.STD_Pre_Mag   = stats.Pre.Std_Mag   / stats.Pre.Mean_Mag;
        norm_stats.STD_Post_Mag  = stats.Post.Std_Mag  / stats.Pre.Mean_Mag;
        norm_stats.STD_Pre_Phase = stats.Pre.Std_Phase / stats.Pre.Mean_Phase;
        norm_stats.STD_Post_Phase= stats.Post.Std_Phase/ stats.Pre.Mean_Phase;
        
        Normalizzato.(nome_cartella) = norm_stats;
    end
end
disp('Normalization completed');

% Normalized PLOTS (Pre=1, Post=Post/Pre)
salvaPath = 'C:\Users\Gabriele Rellini\OneDrive\Desktop\Università\2. Magistrale\Tesi\Sperimentazione\4 - Misure Perugia 26-11-25\Elaborazione\Normalization Plots'; 
if ~exist(salvaPath, 'dir')
    mkdir(salvaPath);
end

for g = 1:numel(campiGruppi)
    gruppo_nome = campiGruppi{g};
    lista_cartelle = gruppi.(gruppo_nome);
    
    mean_mag = [];
    std_mag  = [];
    mean_phase = [];
    std_phase  = [];
    labels = {};
    
    for k = 1:numel(lista_cartelle)
        nome_cartella = lista_cartelle{k};
        if isfield(Normalizzato, nome_cartella)
            norm_stats = Normalizzato.(nome_cartella);
            % Insert Pre then Post
            mean_mag(end+1)   = norm_stats.Pre_Mag;
            std_mag(end+1)    = norm_stats.STD_Pre_Mag;  
            mean_phase(end+1) = norm_stats.Pre_Phase;
            std_phase(end+1)  = norm_stats.STD_Pre_Phase;
            
            mean_mag(end+1)   = norm_stats.Post_Mag;
            std_mag(end+1)    = norm_stats.STD_Post_Mag;
            mean_phase(end+1) = norm_stats.Post_Phase;
            std_phase(end+1)  = norm_stats.STD_Post_Phase;
            
            % Labels
            labels{end+1} = getCampioneLabel(nome_cartella);
        end
    end
    
    nCampioni = numel(lista_cartelle);
    X = 1:nCampioni;
    
    % --- Magnitude ---
    figure('Visible','off');
    hold on;
    bar(X-0.15, mean_mag(1:2:end), 0.3, 'FaceColor','b'); % Pre = 1
    bar(X+0.15, mean_mag(2:2:end), 0.3, 'FaceColor','r'); % Post
    errorbar(X-0.15, mean_mag(1:2:end), std_mag(1:2:end), 'k.', 'LineWidth',1);
    errorbar(X+0.15, mean_mag(2:2:end), std_mag(2:2:end), 'k.', 'LineWidth',1);
    yline(1,'k--','LineWidth',1.5); 
    set(gca,'XTick',X,'XTickLabel',labels);
    ylabel('Normalized Magnitude (Pre=1) \pm STD');
    title(sprintf('%s - Normalized Magnitude', gruppo_nome),'Interpreter','none');
    legend({'Pre','Post'});
    grid on;
    filename = fullfile(salvaPath, sprintf('Normalized_Mag_%s.png', gruppo_nome));
    saveas(gcf, filename); close(gcf);
    
    % --- Phase ---
    figure('Visible','off');
    hold on;
    bar(X-0.15, mean_phase(1:2:end), 0.3, 'FaceColor','b'); % Pre = 1
    bar(X+0.15, mean_phase(2:2:end), 0.3, 'FaceColor','r'); % Post
    errorbar(X-0.15, mean_phase(1:2:end), std_phase(1:2:end), 'k.', 'LineWidth',1);
    errorbar(X+0.15, mean_phase(2:2:end), std_phase(2:2:end), 'k.', 'LineWidth',1);
    yline(1,'k--','LineWidth',1.5); 
    set(gca,'XTick',X,'XTickLabel',labels);
    ylabel('Normalized Phase (Pre=1) \pm STD');
    title(sprintf('%s - Normalized Phase', gruppo_nome),'Interpreter','none');
    legend({'Pre','Post'});
    grid on;
    filename = fullfile(salvaPath, sprintf('Normalized_Phase_%s.png', gruppo_nome));
    saveas(gcf, filename); close(gcf);
end
disp('Normalized plots completed.');

%% 5. DELTA (Median / MAD and Mean / STD) integrated in Statistica
for k = 1:numel(cartelle)
    nome_cartella = cartelle{k};
    stats = Statistica.(nome_cartella);
    if isfield(stats,'Pre') && isfield(stats,'Post')
        delta = struct;
        
        % --- Δ on Median ---
        delta.Delta_Mag   = stats.Post.Median_Mag   - stats.Pre.Median_Mag;
        delta.Delta_Phase = stats.Post.Median_Phase - stats.Pre.Median_Phase;
        delta.DeltaNorm_Mag   = delta.Delta_Mag   / stats.Pre.MAD_Mag;
        delta.DeltaNorm_Phase = delta.Delta_Phase / stats.Pre.MAD_Phase;
        
        % SE of the median delta (using MAD, unchanged)
        delta.SE_Delta_Mag   = sqrt(stats.Pre.MAD_Mag^2   + stats.Post.MAD_Mag^2);
        delta.SE_Delta_Phase = sqrt(stats.Pre.MAD_Phase^2 + stats.Post.MAD_Phase^2);
        
        % --- Δ on Mean ---
        delta.DeltaMean_Mag   = stats.Post.Mean_Mag   - stats.Pre.Mean_Mag;
        delta.DeltaMean_Phase = stats.Post.Mean_Phase - stats.Pre.Mean_Phase;
        
        % *** STD *** of the mean delta (error propagation using STD)
        delta.STD_DeltaMean_Mag   = sqrt(stats.Pre.Std_Mag^2   + stats.Post.Std_Mag^2);
        delta.STD_DeltaMean_Phase = sqrt(stats.Pre.Std_Phase^2 + stats.Post.Std_Phase^2);
        
        % Normalization by Pre STD
        if stats.Pre.Std_Mag > 0
            delta.DeltaMeanNorm_Mag = delta.DeltaMean_Mag / stats.Pre.Std_Mag;
        else
            delta.DeltaMeanNorm_Mag = NaN;
        end
        
        if stats.Pre.Std_Phase > 0
            delta.DeltaMeanNorm_Phase = delta.DeltaMean_Phase / stats.Pre.Std_Phase;
        else
            delta.DeltaMeanNorm_Phase = NaN;
        end
        
        % Save in main struct
        Statistica.(nome_cartella).Delta = delta;
    end
end

% --- Saving updated Statistica struct ---
salvaPath = 'C:\Users\Gabriele Rellini\OneDrive\Desktop\Università\2. Magistrale\Tesi\Sperimentazione\4 - Misure Perugia 26-11-25\Elaborazione';
if ~exist(salvaPath, 'dir'); mkdir(salvaPath); end
save(fullfile(salvaPath,'dati_Perugia_26_11_statistiche.mat'),'Statistica');
disp('Delta (median and mean) calculated and saved in Statistica struct.');

% --- NORMALIZED Δ PLOTS ---
% === 5a. Δ/MAD (Median) ===
salvaPath_MAD = fullfile(salvaPath,'Delta vs MAD plots');
if ~exist(salvaPath_MAD,'dir'); mkdir(salvaPath_MAD); end

for g = 1:numel(campiGruppi)
    gruppo_nome = campiGruppi{g};
    lista_cartelle = gruppi.(gruppo_nome);
    delta_mag = []; delta_phase = []; labels = {};
    for k = 1:numel(lista_cartelle)
        nome_cartella = lista_cartelle{k};
        if isfield(Statistica.(nome_cartella),'Delta')
            d = Statistica.(nome_cartella).Delta;
            delta_mag(end+1)   = d.DeltaNorm_Mag;
            delta_phase(end+1) = d.DeltaNorm_Phase;
            labels{end+1} = getCampioneLabel(nome_cartella);
        end
    end
    
    X = 1:numel(delta_mag);
    % --- Magnitude ---
    figure('Visible','off');
    hold on;
    for i = 1:numel(delta_mag)
        if abs(delta_mag(i)) < 1
            col = [0.2 0.8 0.2];
        elseif abs(delta_mag(i)) < 2
            col = [0.9 0.7 0.1];
        else
            col = [0.9 0.2 0.2];
        end
        bar(i, delta_mag(i), 0.5, 'FaceColor', col);
    end
    yline(0,'k--','LineWidth',1.5);
    set(gca,'XTick',X,'XTickLabel',labels);
    ylabel('Δ Median Magnitude / MAD_{Pre}');
    title(sprintf('%s - Δ Magnitude (median)', gruppo_nome),'Interpreter','none');
    grid on;
    h1 = bar(nan, nan, 'FaceColor', [0.2 0.8 0.2]);
    h2 = bar(nan, nan, 'FaceColor', [0.9 0.7 0.1]);
    h3 = bar(nan, nan, 'FaceColor', [0.9 0.2 0.2]);
    legend([h1 h2 h3], {'Weak (|Δ|<1 MAD)', 'Moderate (1–2 MAD)', 'Strong (≥2 MAD)'}, 'Location', 'best'); 
    saveas(gcf, fullfile(salvaPath_MAD, sprintf('Delta_Mag_%s.png', gruppo_nome)));
    close(gcf);
    
    % --- Phase ---
    figure('Visible','off');
    hold on;
    for i = 1:numel(delta_phase)
        if abs(delta_phase(i)) < 1
            col = [0.2 0.8 0.2];
        elseif abs(delta_phase(i)) < 2
            col = [0.9 0.7 0.1];
        else
            col = [0.9 0.2 0.2];
        end
        bar(i, delta_phase(i), 0.5, 'FaceColor', col);
    end
    yline(0,'k--','LineWidth',1.5);
    set(gca,'XTick',X,'XTickLabel',labels);
    ylabel('Δ Median Phase / MAD_{Pre}');
    title(sprintf('%s - Δ Phase (median)', gruppo_nome),'Interpreter','none');
    grid on;
    h1 = bar(nan, nan, 'FaceColor', [0.2 0.8 0.2]);
    h2 = bar(nan, nan, 'FaceColor', [0.9 0.7 0.1]);
    h3 = bar(nan, nan, 'FaceColor', [0.9 0.2 0.2]);
    legend([h1 h2 h3], {'Weak (|Δ|<1 MAD)', 'Moderate (1–2 MAD)', 'Strong (≥2 MAD)'}, 'Location', 'best');
    saveas(gcf, fullfile(salvaPath_MAD, sprintf('Delta_Phase_%s.png', gruppo_nome)));
    close(gcf);
end

% === 5b. Δ/STD (Mean) ===
salvaPath_STD = fullfile(salvaPath,'Delta vs STD plots');
if ~exist(salvaPath_STD,'dir'); mkdir(salvaPath_STD); end
for g = 1:numel(campiGruppi)
    gruppo_nome = campiGruppi{g};
    lista_cartelle = gruppi.(gruppo_nome);
    delta_mean_mag = []; delta_mean_phase = []; labels = {};
    for k = 1:numel(lista_cartelle)
        nome_cartella = lista_cartelle{k};
        if isfield(Statistica.(nome_cartella),'Delta')
            d = Statistica.(nome_cartella).Delta;
            delta_mean_mag(end+1)   = d.DeltaMeanNorm_Mag;
            delta_mean_phase(end+1) = d.DeltaMeanNorm_Phase;
            labels{end+1} = getCampioneLabel(nome_cartella);
        end
    end
    
    X = 1:numel(delta_mean_mag);
    % --- Magnitude ---
    figure('Visible','off');
    hold on;
    for i = 1:numel(delta_mean_mag)
        val = delta_mean_mag(i);
        if isnan(val)
            col = [0.6 0.6 0.6];
        elseif abs(val) < 1
            col = [0.2 0.8 0.2];
        elseif abs(val) < 2
            col = [0.9 0.7 0.1];
        else
            col = [0.9 0.2 0.2];
        end
        bar(i, val, 0.5, 'FaceColor', col);
    end
    yline(0,'k--','LineWidth',1.5);
    set(gca,'XTick',X,'XTickLabel',labels);
    ylabel('Δ Mean Magnitude / STD_{Pre}'); 
    title(sprintf('%s - Δ Magnitude (mean)', gruppo_nome),'Interpreter','none');
    grid on;
    h1 = bar(nan, nan, 'FaceColor', [0.2 0.8 0.2]);
    h2 = bar(nan, nan, 'FaceColor', [0.9 0.7 0.1]);
    h3 = bar(nan, nan, 'FaceColor', [0.9 0.2 0.2]);
    legend([h1 h2 h3], {'Weak (|Δ|<1 STD)', 'Moderate (1–2 STD)', 'Strong (≥2 STD)'}, 'Location', 'best');
    saveas(gcf, fullfile(salvaPath_STD, sprintf('DeltaMean_Mag_%s.png', gruppo_nome)));
    close(gcf);
    
    % --- Phase ---
    figure('Visible','off');
    hold on;
    for i = 1:numel(delta_mean_phase)
        val = delta_mean_phase(i);
        if isnan(val)
            col = [0.6 0.6 0.6];
        elseif abs(val) < 1
            col = [0.2 0.8 0.2];
        elseif abs(val) < 2
            col = [0.9 0.7 0.1];
        else
            col = [0.9 0.2 0.2];
        end
        bar(i, val, 0.5, 'FaceColor', col);
    end
    yline(0,'k--','LineWidth',1.5);
    set(gca,'XTick',X,'XTickLabel',labels);
    ylabel('Δ Mean Phase / STD_{Pre}');
    title(sprintf('%s - Δ Phase (mean)', gruppo_nome),'Interpreter','none');
    grid on;
    h1 = bar(nan, nan, 'FaceColor', [0.2 0.8 0.2]);
    h2 = bar(nan, nan, 'FaceColor', [0.9 0.7 0.1]);
    h3 = bar(nan, nan, 'FaceColor', [0.9 0.2 0.2]);
    legend([h1 h2 h3], {'Weak (|Δ|<1 STD)', 'Moderate (1–2 STD)', 'Strong (≥2 STD)'}, 'Location', 'best');
    saveas(gcf, fullfile(salvaPath_STD, sprintf('DeltaMean_Phase_%s.png', gruppo_nome)));
    close(gcf);
end

% NYQUIST PLOTS (Real vs Imaginary)
salvaPathNyquist = fullfile(salvaPath, 'Nyquist Plots');
if ~exist(salvaPathNyquist, 'dir'), mkdir(salvaPathNyquist); end

for g = 1:numel(campiGruppi)
    gruppo_nome = campiGruppi{g};
    lista_cartelle = gruppi.(gruppo_nome);
    figure('Visible','off'); hold on;
    for k = 1:numel(lista_cartelle)
        nome_cartella = lista_cartelle{k};
        stats = Statistica.(nome_cartella);
        % Reconstruct Real/Imaginary from mean values
        Z_pre_R = stats.Pre.Mean_Mag * cosd(stats.Pre.Mean_Phase);
        Z_pre_I = stats.Pre.Mean_Mag * sind(stats.Pre.Mean_Phase);
        Z_post_R = stats.Post.Mean_Mag * cosd(stats.Post.Mean_Phase);
        Z_post_I = stats.Post.Mean_Mag * sind(stats.Post.Mean_Phase);
        
        plot([Z_pre_R, Z_post_R], [Z_pre_I, Z_post_I], '-', 'Color', [0.6 0.6 0.6]);
        hPre = plot(Z_pre_R, Z_pre_I, 'bo', 'MarkerFaceColor','b');
        hPost = plot(Z_post_R, Z_post_I, 'rs', 'MarkerFaceColor','r');
        text(Z_pre_R, Z_pre_I, sprintf('  %d', k), 'FontSize',8);
    end
    xlabel('Real Part [\Omega]'); ylabel('Imaginary Part [\Omega]');
    title(sprintf('%s - Nyquist', gruppo_nome), 'Interpreter', 'none');
    legend([hPre, hPost], {'Pre', 'Post'}, 'Location', 'best');
    grid on; axis equal;
    saveas(gcf, fullfile(salvaPathNyquist, sprintf('Nyquist_%s.png', gruppo_nome))); close(gcf);
end
disp('Nyquist plots generated.');

%% 6. MEGA, MINI AND STIM SUMMARY TABLES
MegaTabella = table();
MiniTabella = table();
StimTabelle = table(); 

for g = 1:numel(campiGruppi)
    gruppo_nome = campiGruppi{g};
    lista_cartelle = gruppi.(gruppo_nome);
    for k = 1:numel(lista_cartelle)
        nome_cartella = lista_cartelle{k};
        
        % --- STIM TABLES ---
        if isfield(Statistica.(nome_cartella), 'Stim') && ~isempty(Statistica.(nome_cartella).Stim)
            stim_data = Statistica.(nome_cartella).Stim;
            n_windows = numel(stim_data);
            for w = 1:n_windows
                win_stats = stim_data{w};
                rigaStim = {
                    gruppo_nome, ...
                    getCampioneLabel(nome_cartella), ...
                    w, ... 
                    win_stats.Mean_Mag, ...
                    win_stats.Std_Mag, ...
                    win_stats.Mean_Phase, ...
                    win_stats.Std_Phase
                };
                StimTabelle = [StimTabelle; rigaStim];
            end
        end
        
        if ~isfield(Statistica.(nome_cartella),'Delta')
            continue;
        end
        stats = Statistica.(nome_cartella);
        delta = stats.Delta;
        
        if isfield(stats, 'Slope_Stim'), val_slope = stats.Slope_Stim; else, val_slope = NaN; end
        
        % Mean ± STD strings
        MeanStd_Pre_Mag  = sprintf('%.3f ± %.3f', stats.Pre.Mean_Mag,  stats.Pre.Std_Mag);
        MeanStd_Post_Mag = sprintf('%.3f ± %.3f', stats.Post.Mean_Mag, stats.Post.Std_Mag);
        MeanStd_Pre_Phase  = sprintf('%.3f ± %.3f', stats.Pre.Mean_Phase,  stats.Pre.Std_Phase);
        MeanStd_Post_Phase = sprintf('%.3f ± %.3f', stats.Post.Mean_Phase, stats.Post.Std_Phase);
        
        % Delta Mean ± STD
        DeltaMeanStd_Mag   = sprintf('%.3f ± %.3f', delta.DeltaMean_Mag,   delta.STD_DeltaMean_Mag);
        DeltaMeanStd_Phase = sprintf('%.3f ± %.3f', delta.DeltaMean_Phase, delta.STD_DeltaMean_Phase);
        
        riga1 = {
            gruppo_nome, ...
            getCampioneLabel(nome_cartella), ...
            val_slope, ... 
            stats.Pre.Mean_Mag, stats.Post.Mean_Mag, ...
            stats.Pre.Std_Mag, stats.Post.Std_Mag, ...
            stats.Pre.Median_Mag, stats.Post.Median_Mag, ...
            stats.Pre.MAD_Mag, stats.Post.MAD_Mag, ...
            delta.DeltaMean_Mag, delta.Delta_Mag, ...
            delta.DeltaMeanNorm_Mag, delta.DeltaNorm_Mag, ...
            stats.Pre.Mean_Phase, stats.Post.Mean_Phase, ...
            stats.Pre.Std_Phase, stats.Post.Std_Phase, ...
            stats.Pre.Median_Phase, stats.Post.Median_Phase, ...
            stats.Pre.MAD_Phase, stats.Post.MAD_Phase, ...
            delta.DeltaMean_Phase, delta.Delta_Phase, ...
            delta.DeltaMeanNorm_Phase, delta.DeltaNorm_Phase ...
            };
        
        riga2 = {
            gruppo_nome, ...
            getCampioneLabel(nome_cartella), ...
            val_slope, ... 
            MeanStd_Pre_Mag, MeanStd_Post_Mag, ...
            DeltaMeanStd_Mag,...
            MeanStd_Pre_Phase, MeanStd_Post_Phase, ...
            DeltaMeanStd_Phase,...
            };
            
        MegaTabella = [MegaTabella; riga1];
        MiniTabella = [MiniTabella; riga2];
    end
end

% MegaTabella Column names
MegaTabella.Properties.VariableNames = { ...
    'Group','Sample', ...
    'Stim_Slope_Ohm_s', ... 
    'Mean_Pre_Mag [Ω]','Mean_Post_Mag [Ω]','STD_Pre_Mag [Ω]','STD_Post_Mag [Ω]', ...
    'Median_Pre_Mag [Ω]','Median_Post_Mag [Ω]','MAD_Pre_Mag [Ω]','MAD_Post_Mag [Ω]', ...
    'DeltaMean_Mag [Ω]','DeltaMedian_Mag [Ω]','DeltaMeanNorm_Mag','DeltaMedianNorm_Mag', ...
    'Mean_Pre_Phase [°]','Mean_Post_Phase [°]','STD_Pre_Phase [°]','STD_Post_Phase [°]', ...
    'Median_Pre_Phase [°]','Median_Post_Phase [°]','MAD_Pre_Phase [°]','MAD_Post_Phase [°]', ...
    'DeltaMean_Phase [°]','DeltaMedian_Phase [°]','DeltaMeanNorm_Phase','DeltaMedianNorm_Phase'};

% MiniTabella Column names
MiniTabella.Properties.VariableNames = { ...
    'Group','Sample', ...
    'Stim_Slope_Ohm_s', ... 
    'Mean±STD_Pre_Mag [Ω]','Mean±STD_Post_Mag [Ω]', ...
    'DeltaMean±STD_Mag [Ω]',...
    'Mean±STD_Pre_Phase [°]','Mean±STD_Post_Phase [°]', ...
    'DeltaMean±STD_Phase [°]'};

% StimTabelle Column names
if ~isempty(StimTabelle)
    StimTabelle.Properties.VariableNames = { ...
        'Group', 'Sample', 'Window_Num', ...
        'Mean_Mag [Ω]', 'STD_Mag [Ω]', ...
        'Mean_Phase [°]', 'STD_Phase [°]'};
end

% Save to Excel and MAT files
salvaPath = 'C:\Users\Gabriele Rellini\OneDrive\Desktop\Università\2. Magistrale\Tesi\Sperimentazione\4 - Misure Perugia 26-11-25\Elaborazione';
save(fullfile(salvaPath,'MegaTabella.mat'),'MegaTabella');
writetable(MegaTabella, fullfile(salvaPath,'MegaTabella.xlsx'));
save(fullfile(salvaPath,'MiniTabella.mat'),'MiniTabella');
writetable(MiniTabella, fullfile(salvaPath,'MiniTabella.xlsx'));
save(fullfile(salvaPath,'StimTabelle.mat'),'StimTabelle'); 
writetable(StimTabelle, fullfile(salvaPath,'StimTabelle.xlsx'));
disp('Mega, Mini and Stim summary tables successfully created and saved.');

%% 7. FREQUENCY ANALYSIS (PSD 0–8 Hz)
salvaPath = 'C:\Users\Gabriele Rellini\OneDrive\Desktop\Università\2. Magistrale\Tesi\Sperimentazione\4 - Misure Perugia 26-11-25\Elaborazione\Frequency analysis plots';
if ~exist(salvaPath, 'dir')
    mkdir(salvaPath);
end

fs_target   = 16;          
maxFreqVis  = 8;           
bands = [0 0.1; 0.1 0.5; 0.5 2; 2 5; 5 8];
bandNames = {'ULF(0-0.1)','VLF(0.1-0.5)','LF(0.5-2)','MF(2-5)','HF(5-8)'};
PSD = struct;

for k = 1:numel(cartelle)
    nome_cartella = cartelle{k};
    titoloSint = getTitoloSintetico(nome_cartella); 
    signals_proc  = dati_proc.(nome_cartella);
    n_signals     = numel(signals_proc);
    
    if n_signals == 3
        legends = {'Pre','Stim','Post'};
    else
        legends = {'Pre','Post'};
    end
    
    % ---------- MAGNITUDE PSD ----------
    figure('Visible','off'); hold on;
    for j = 1:n_signals
        x   = signals_proc{j}.data.Magnitude;
        fs0 = signals_proc{j}.fs;
        if fs0 ~= fs_target
            x = resample(x, fs_target, fs0);
        end
        [f_mag, Pxx_mag] = stimaPSDWelch(x, fs_target);
        mask = (f_mag <= maxFreqVis);
        plot(f_mag(mask), 10*log10(Pxx_mag(mask)+eps), 'LineWidth', 1.2);
        PSD.(nome_cartella).Mag{j}.f   = f_mag(:);
        PSD.(nome_cartella).Mag{j}.Pxx = Pxx_mag(:);
        PSD.(nome_cartella).Mag{j}.BandPowers = bandPowers(f_mag, Pxx_mag, bands);
        PSD.(nome_cartella).Mag{j}.label = legends{j};
    end
    xlabel('Frequency [Hz]'); ylabel('PSD Magnitude [dB/Hz]');
    title('Magnitude PSD');
    legend(legends, 'Location','best'); grid on;
    saveas(gcf, fullfile(salvaPath, sprintf('PSD_Mag_%s.png', titoloSint))); close(gcf);
    
    % ---------- PHASE PSD ----------
    figure('Visible','off'); hold on;
    for j = 1:n_signals
        x   = signals_proc{j}.data.Phase;
        fs0 = signals_proc{j}.fs;
        if fs0 ~= fs_target
            x = resample(x, fs_target, fs0);
        end
        [f_ph, Pxx_ph] = stimaPSDWelch(x, fs_target);
        mask = (f_ph <= maxFreqVis);
        plot(f_ph(mask), 10*log10(Pxx_ph(mask)+eps), 'LineWidth', 1.2);
        PSD.(nome_cartella).Phase{j}.f   = f_ph(:);
        PSD.(nome_cartella).Phase{j}.Pxx = Pxx_ph(:);
        PSD.(nome_cartella).Phase{j}.BandPowers = bandPowers(f_ph, Pxx_ph, bands);
        PSD.(nome_cartella).Phase{j}.label = legends{j};
    end
    xlabel('Frequency [Hz]'); ylabel('PSD Phase [dB/Hz]');
    title('Phase PSD');
    legend(legends, 'Location','best'); grid on;
    saveas(gcf, fullfile(salvaPath, sprintf('PSD_Phase_%s.png', titoloSint))); close(gcf);
    
    % ---------- BARPLOT ----------
    if n_signals == 3
        % --- Pre vs Post (non-normalized) ---
        bpMag_pre  = PSD.(nome_cartella).Mag{1}.BandPowers;
        bpMag_post = PSD.(nome_cartella).Mag{3}.BandPowers;
        bpPh_pre   = PSD.(nome_cartella).Phase{1}.BandPowers;
        bpPh_post  = PSD.(nome_cartella).Phase{3}.BandPowers;
        
        figure('Visible','off');
        subplot(1,2,1);
        bar([bpMag_pre(:), bpMag_post(:)]);
        set(gca,'XTickLabel',bandNames);
        xlabel('Frequency Bands'); ylabel('Magnitude Power');
        legend({'Pre','Post'}); title('Magnitude'); grid on;
        
        subplot(1,2,2);
        bar([bpPh_pre(:), bpPh_post(:)]);
        set(gca,'XTickLabel',bandNames);
        xlabel('Frequency Bands'); ylabel('Phase Power');
        legend({'Pre','Post'}); title('Phase'); grid on;
        sgtitle(sprintf('%s - Pre vs Post', titoloSint),'Interpreter','none');
        saveas(gcf, fullfile(salvaPath, sprintf('Barplot_PrePost_%s.png', titoloSint))); close(gcf);
        
        % --- Pre vs Post (normalized to Pre = 1) ---
        bpMag_norm = [ones(size(bpMag_pre(:))), bpMag_post(:)./bpMag_pre(:)];
        bpPh_norm  = [ones(size(bpPh_pre(:))),  bpPh_post(:)./bpPh_pre(:)];
        
        figure('Visible','off');
        subplot(1,2,1);
        bar(bpMag_norm);
        set(gca,'XTickLabel',bandNames);
        xlabel('Frequency Bands'); ylabel('Normalized Power (Magnitude)');
        legend({'Pre=1','Post/Pre'}); title('Magnitude'); grid on;
        
        subplot(1,2,2);
        bar(bpPh_norm);
        set(gca,'XTickLabel',bandNames);
        xlabel('Frequency Bands'); ylabel('Normalized Power (Phase)');
        legend({'Pre=1','Post/Pre'}); title('Phase'); grid on;
        sgtitle(sprintf('%s - Pre vs Post NORMALIZED', titoloSint),'Interpreter','none');
        saveas(gcf, fullfile(salvaPath, sprintf('Barplot_PrePostNorm_%s.png', titoloSint))); close(gcf);
    end
end

% Save PSD variable
save(fullfile(salvaPath,'dati_Perugia_26_11_PSD.mat'),'PSD');

%% 8. NORMALIZED MEAN TREND PLOT (Pre -> Post -> 24h)
disp('--- Generating NORMALIZED mean trend plot (% of Pre) ---');
salvaPath24h = fullfile(salvaPath, '24h Trend Plots');
if ~exist(salvaPath24h, 'dir'), mkdir(salvaPath24h); end

% === 24H VALUES DEFINITION (Manual) ===
Valori_24h = struct();
Valori_24h.controlli.primo   = 222; 
Valori_24h.stim_45.primo     = 273; 
Valori_24h.stim_45.secondo   = 199;
Valori_24h.stim_45.terzo     = 190;
Valori_24h.stim_169.primo    = 196; 
Valori_24h.stim_169.secondo  = 217;
Valori_24h.stim_169.terzo    = 277;

% === AGGREGATE DATA BY GROUP ===
nomi_gruppi = fieldnames(gruppi);
colorMap = containers.Map({'controlli', 'stim_45', 'stim_169'}, {[0.5 0.5 0.5], [0.2 0.8 0.2], [0.2 0.6 0.8]});
labelMap = containers.Map({'controlli', 'stim_45', 'stim_169'}, {'Control', '4.5 \muA', '16.9 \muA'});

figure('Color','w', 'Position', [100 100 800 600]); 
hold on; grid on;

for g = 1:numel(nomi_gruppi)
    nome_gruppo = nomi_gruppi{g};
    lista_cartelle = gruppi.(nome_gruppo);
    if isempty(lista_cartelle), continue; end
    
    vec_pre_norm  = []; vec_post_norm = []; vec_24h_norm  = [];
    
    for k = 1:numel(lista_cartelle)
        nome_cartella = lista_cartelle{k};
        stats = Statistica.(nome_cartella);
        raw_pre  = stats.Pre.Mean_Mag;
        raw_post = stats.Post.Mean_Mag;
        raw_24h = NaN;
        key_g = lower(nome_gruppo);
        nome_lower = lower(nome_cartella);
        
        if contains(nome_lower, 'primo') || contains(nome_lower, 'a1') || contains(nome_lower, ' 1')
             if isfield(Valori_24h, key_g) && isfield(Valori_24h.(key_g),'primo')
                 raw_24h = Valori_24h.(key_g).primo; 
             end
        elseif contains(nome_lower, 'secondo') || contains(nome_lower, 'a3') || contains(nome_lower, ' 2')
             if isfield(Valori_24h, key_g) && isfield(Valori_24h.(key_g),'secondo')
                 raw_24h = Valori_24h.(key_g).secondo; 
             end
        elseif contains(nome_lower, 'terzo') || contains(nome_lower, ' 3')
             if isfield(Valori_24h, key_g) && isfield(Valori_24h.(key_g),'terzo')
                 raw_24h = Valori_24h.(key_g).terzo; 
             end
        end
        
        vec_pre_norm(end+1)  = 100;
        vec_post_norm(end+1) = (raw_post / raw_pre) * 100;
        vec_24h_norm(end+1)  = (raw_24h / raw_pre) * 100;
    end
    
    mean_pre  = mean(vec_pre_norm);
    mean_post = mean(vec_post_norm);
    mean_24h  = mean(vec_24h_norm, 'omitnan');

    if isnan(mean_24h) 
        continue; 
    end 
    
    sem_pre  = std(vec_pre_norm) / sqrt(length(vec_pre_norm));
    sem_post = std(vec_post_norm) / sqrt(length(vec_post_norm));
    sem_24h  = std(vec_24h_norm, 'omitnan') / sqrt(sum(~isnan(vec_24h_norm)));
    
    colore = colorMap(key_g); lab = labelMap(key_g);
    x_axis = [1, 2, 3]; y_axis = [mean_pre, mean_post, mean_24h]; err_bar = [sem_pre, sem_post, sem_24h];
    
    errorbar(x_axis, y_axis, err_bar, '-o', 'Color', colore, 'LineWidth', 2.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', colore, 'DisplayName', lab);
end

ylabel('Normalized TEER (% of Pre)');
title('Normalized Temporal Evolution (Pre = 100%)');
xticks([1 2 3]); xticklabels({'Pre-Stim', 'Post-Stim (20min)', 'Follow-up (24h)'});
xlim([0.8 3.2]); yline(100, '--k', 'DisplayName', 'Baseline (100%)');
legend('Location', 'best'); set(gca, 'FontSize', 12);
saveas(gcf, fullfile(salvaPath24h, 'Comparative_Normalized_Trend.png'));
disp('Normalized trend plot saved in "24h Trend Plots"');

%% Helper Functions
function [f, Pxx] = stimaPSDWelch(x, fs)
    x = detrend(x,'linear');
    N = numel(x);
    winLen = min(N, round(fs*60));
    winLen = max(winLen, min(N, 128));
    w = hamming(winLen);
    nover = floor(0.5*winLen);
    nfft  = max(1024, 2^nextpow2(winLen));
    [Pxx, f] = pwelch(x, w, nover, nfft, fs, 'onesided');
end

function bp = bandPowers(f, Pxx, bands)
    nb = size(bands,1); bp = zeros(1,nb);
    for b = 1:nb
        f1 = bands(b,1); f2 = bands(b,2);
        m  = (f >= f1) & (f <= f2);
        if nnz(m) >= 2
            bp(b) = trapz(f(m), Pxx(m));
        end
    end
end

function label = getCampioneLabel(nome)
    nome = lower(nome);
    if contains(nome,'primo')
        label = 'First';
    elseif contains(nome,'secondo')
        label = 'Second';
    elseif contains(nome,'terzo')
        label = 'Third';
    else, label = 'Sample'; 
    end
end

function stats = CalcolaStatistiche(tab, idx_range)
    mag = tab.Magnitude(idx_range);
    phase = tab.Phase(idx_range);
    N = numel(mag); 
    stats = struct;
    stats.Mean_Mag   = mean(mag);
    stats.Std_Mag    = std(mag);
    stats.SE_Mag     = stats.Std_Mag / sqrt(N);
    stats.Median_Mag = median(mag);
    stats.MAD_Mag    = mad(mag,1);
    stats.Mean_Phase   = mean(phase);
    stats.Std_Phase    = std(phase);
    stats.SE_Phase     = stats.Std_Phase / sqrt(N);
    stats.Median_Phase = median(phase);
    stats.MAD_Phase    = mad(phase,1);
end

function titolo = getTitoloSintetico(nome_cartella)
    nome = lower(nome_cartella);
    if contains(nome,'controllo')
        if contains(nome,'primo')
            titolo = 'Control - First';
        elseif contains(nome,'secondo')
            titolo = 'Control - Second'; 
        end
        return;
    end
    if contains(nome,'4_5')
        stim = '4.5 uApp';
    elseif contains(nome,'16_9')
        stim = '16.9 uApp'; 
    end

    if contains(nome,'primo')
        camp = 'First sample';
    elseif contains(nome,'secondo')
        camp = 'Second sample';
    elseif contains(nome,'terzo')
        camp = 'Third sample';
    else, camp = 'Sample'; 
    end
    
    titolo = sprintf('%s - %s', stim, camp);
end