clear all
clc
close all

% NOMINAL DATA
corr_nom = [ ...
    0.016, 0.032, 0.080, 0.160, ...
    0.320, 0.640, 1.6, 3.2, ...
    6.4, 12.8, 32, 64,...
    128, 256, 640, 1280]; % Nominal software µA RMS
volt_nom = [ ...
    35.4, 70.7, 177, 354, ...
    35.4, 70.7, 177, 354, ...
    35.4, 70.7, 177, 354, ...
    35.4, 70.7, 177, 354]; % mVrms
res_nom = [ ...
    2210e3, 2210e3, 2210e3, 2210e3, ...
    110.5e3, 110.5e3, 110.5e3, 110.5e3, ...
    5.525e3, 5.525e3, 5.525e3, 5.525e3, ...
    276.25, 276.25, 276.25, 276.25]; % Ohm

% TEST CONDITIONS
f = [12.5, 16, 50, 100, 500, 1000, 5000, 12500, 50000]; % Hz 
C_values = [47e-9, 4.7e-6]; % F

% STRUCT AND CALCULATION
dati = struct();
for iF = 1:length(f)
    freq = f(iF);
    freq_str = strrep(sprintf('%.3f', freq), '.', '_');
    field_name = ['f_' freq_str 'Hz'];
    dati.(field_name) = struct();
    for iC = 1:length(C_values)
        C = C_values(iC);
        Xc = 1 ./ (2*pi*freq*C);
        I_real = zeros(size(corr_nom));
        for i = 1:length(corr_nom)
            V = volt_nom(i) * 1e-3;
            R = res_nom(i);
            Z = sqrt(R.^2 + Xc.^2);
            I_real(i) = (V ./ Z) * 1e6; % µA RMS
        end
        C_str = strrep(sprintf('%.0e',C), '.', '_');
        C_str = strrep(C_str, '-', '_');
        C_field_name = ['C_' C_str];
        dati_tab = table( ...
            corr_nom', volt_nom', res_nom', ...
            repmat(C, length(corr_nom),1), ...
            repmat(freq, length(corr_nom),1), ...
            I_real', ...
            'VariableNames', {'I_nom_uA', 'V_nom_mV', 'R_Ohm', 'C_F', 'f_Hz', 'I_real_uA'});
        dati.(field_name).(C_field_name) = struct('table', dati_tab);
    end
end

% PLOT
for iC = 1:length(C_values)
    C = C_values(iC);
    figure('Name', sprintf('Actual current - C = %.2g F', C), 'NumberTitle', 'off');
    hold on
    grid on
    % Range and resistance definitions
    resistenze = [2.21e6, 110.5e3, 5.525e3, 276.25];
    intervalli = {[1 4], [5 8], [9 12], [13 16]}; % Indices in corr_nom vector
    colori = [
        1.0 0.7 0.7;   % light red
        0.7 1.0 0.7;   % light green
        0.7 0.8 1.0;   % light blue
        1.0 1.0 0.7    % light yellow
    ];
    nomi_legenda_R = { ...
        'R = 2.21 MΩ', ...
        'R = 110.5 kΩ', ...
        'R = 5.525 kΩ', ...
        'R = 276.25 Ω' ...
    };
    % Draw colored zones (before curves)
    for k = 1:length(resistenze)
        x_start = corr_nom(intervalli{k}(1));
        x_end = corr_nom(intervalli{k}(2));
        h_patch(k) = patch([x_start x_end x_end x_start], ...
            [1e-2 1e-2 1e4 1e4], colori(k,:), ...
            'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    % Plot frequency curves
    h_curve = gobjects(length(f), 1);
    for iF = 1:length(f)
        freq = f(iF);
        freq_str = strrep(sprintf('%.3f', freq), '.', '_');
        field_name = ['f_' freq_str 'Hz'];
        C_str = strrep(sprintf('%.0e',C), '.', '_');
        C_str = strrep(C_str, '-', '_');
        C_field_name = ['C_' C_str];
        dati_tab = dati.(field_name).(C_field_name).table;
        h_curve(iF) = plot(dati_tab.I_nom_uA, dati_tab.I_real_uA, '-o', ...
            'LineWidth', 1.5, 'DisplayName', sprintf('f = %.0f Hz', freq));
    end
    % Ideal line
    h_ideal = plot(corr_nom, corr_nom, '--k', 'LineWidth', 1.2, ...
        'DisplayName', 'Ideal (I_{actual} = I_{nom})');
    % Bring curves to front
    children = get(gca,'Children');
    set(gca,'Children',[children(end); children(1:end-1)]);
    % Combine legends (curves + resistances)
    lgd = legend([h_curve; h_ideal; h_patch'], ...
        [arrayfun(@(x) sprintf('f = %.0f Hz', f(x)), 1:length(f), 'UniformOutput', false), ...
        {'Ideal (I_{actual} = I_{nom})'}, nomi_legenda_R], ...
        'Location', 'northwest', 'FontSize', 8, ...
        'NumColumns', 2);
    % Graphical settings
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('Nominal current [\muA RMS]');
    ylabel('Actual current [\muA RMS]');
    title(sprintf('Nominal vs Actual Current Comparison\nC = %.2g F', C));
    grid on
    hold off
end