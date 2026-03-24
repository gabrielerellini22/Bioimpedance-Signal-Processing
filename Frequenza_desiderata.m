clear all
clc

% INPUT DESIRED FREQUENCY
F_BIOZ = input('Enter the desired stimulus frequency (Hz): ');

% SYSTEM CONSTANTS
REF_CLK     = 32768;        % Hz
PLL_MIN     = 14e6;         % Hz
PLL_MAX     = 28e6;         % Hz
ADC_CLK_MIN = 16e3;         % Hz
ADC_CLK_MAX = 36.375e3;     % Hz

% Acceptable relative error limit (es. 1% = 0.01)
MAX_FREQ_ERROR = 0.01;

% DISCRETE DIVISOR VALUES
KDIV_vals     = 2.^(0:13);            % 1..8192
DAC_OSR_vals  = [32, 64, 128, 256];
NDIV_vals     = [512, 1024];
ADC_OSR_vals  = 2.^(3:10);            % 8..1024
M_MIN = 427; M_MAX = 854;             % MDIV+1 range

found    = false;
bestPLL  = 0;
bestErr  = inf;

%% CASE A: Low frequencies (<54.688 Hz)
if F_BIOZ < 54688
    DAC_OSR_fixed = 256;
    Target_SYNTH_CLK = F_BIOZ * DAC_OSR_fixed;

    for KDIV = KDIV_vals
        PLL_target = Target_SYNTH_CLK * KDIV;
        if PLL_target < PLL_MIN || PLL_target > PLL_MAX
            continue;
        end

        M = round(PLL_target / REF_CLK);
        if M < M_MIN || M > M_MAX
            continue;
        end

        PLL_eff = REF_CLK * M;
        if PLL_eff < PLL_MIN || PLL_eff > PLL_MAX
            continue;
        end

        for NDIV = NDIV_vals
            ADC_CLK = PLL_eff / NDIV;
            if ADC_CLK < ADC_CLK_MIN || ADC_CLK > ADC_CLK_MAX
                continue;
            end

            for ADC_OSR = ADC_OSR_vals
                C_BIOZ = (NDIV * ADC_OSR) / (KDIV * DAC_OSR_fixed);

                % % Discard C_BIOZ = 1
                % if abs(C_BIOZ - 1) < 1e-6
                %     continue;
                % end

                % Must be an integer (or 0.5 if F_BIOZ=8 Hz)
                if F_BIOZ == 8
                    if ~(abs(C_BIOZ - round(C_BIOZ)) < 1e-6 || abs(C_BIOZ - 0.5) < 1e-6)
                        continue;
                    end
                else
                    if abs(C_BIOZ - round(C_BIOZ)) > 1e-6
                        continue;
                    end
                end

                F_BIOZ_eff = PLL_eff / (KDIV * DAC_OSR_fixed);
                relErr = abs(F_BIOZ_eff - F_BIOZ) / F_BIOZ;
                if relErr > MAX_FREQ_ERROR
                    continue;
                end

                % Selection: Highest PLL, then lowest error
                if (PLL_eff > bestPLL) || (abs(PLL_eff - bestPLL) < 1e-3 && relErr < bestErr)
                    found   = true;
                    bestPLL = PLL_eff; bestErr = relErr;
                    best.M = M; best.KDIV = KDIV; best.DAC_OSR = DAC_OSR_fixed;
                    best.NDIV = NDIV; best.ADC_OSR = ADC_OSR;
                    best.F_SAMPLE = PLL_eff / (NDIV * ADC_OSR);
                    best.F_BIOZ_eff = F_BIOZ_eff;
                    best.C_BIOZ = C_BIOZ;
                end
            end
        end
    end

else
    %% CASE B: High frequencies (≥54.688 Hz)
    KDIV_fixed = 1;

    for DAC_OSR = DAC_OSR_vals
        PLL_target = F_BIOZ * KDIV_fixed * DAC_OSR;
        if PLL_target < PLL_MIN || PLL_target > PLL_MAX
            continue;
        end

        M = round(PLL_target / REF_CLK);
        if M < M_MIN || M > M_MAX
            continue;
        end

        PLL_eff = REF_CLK * M;
        if PLL_eff < PLL_MIN || PLL_eff > PLL_MAX
            continue;
        end

        for NDIV = NDIV_vals
            ADC_CLK = PLL_eff / NDIV;
            if ADC_CLK < ADC_CLK_MIN || ADC_CLK > ADC_CLK_MAX
                continue;
            end

            for ADC_OSR = ADC_OSR_vals
                C_BIOZ = (NDIV * ADC_OSR) / (KDIV_fixed * DAC_OSR);

                % % Discard C_BIOZ = 1
                % if abs(C_BIOZ - 1) < 1e-6
                %     continue;
                % end

                % Must be an integer (or 0.5 if 8 Hz)
                if F_BIOZ == 8
                    if ~(abs(C_BIOZ - round(C_BIOZ)) < 1e-6 || abs(C_BIOZ - 0.5) < 1e-6)
                        continue;
                    end
                else
                    if abs(C_BIOZ - round(C_BIOZ)) > 1e-6
                        continue;
                    end
                end

                F_BIOZ_eff = PLL_eff / (KDIV_fixed * DAC_OSR);
                relErr = abs(F_BIOZ_eff - F_BIOZ) / F_BIOZ;
                if relErr > MAX_FREQ_ERROR
                    continue;
                end

                if (PLL_eff > bestPLL) || (abs(PLL_eff - bestPLL) < 1e-3 && relErr < bestErr)
                    found   = true;
                    bestPLL = PLL_eff; bestErr = relErr;
                    best.M = M; best.KDIV = KDIV_fixed; best.DAC_OSR = DAC_OSR;
                    best.NDIV = NDIV; best.ADC_OSR = ADC_OSR;
                    best.F_SAMPLE = PLL_eff / (NDIV * ADC_OSR);
                    best.F_BIOZ_eff = F_BIOZ_eff;
                    best.C_BIOZ = C_BIOZ;
                end
            end
        end
    end
end

% FINAL RESULTS
if ~found
    error('❌ No valid combination found. Increase MAX_FREQ_ERROR or vary the frequency.');
end

fprintf('\n=== Main Parameters ===\n');
fprintf('Actual PLL: %.3f MHz\n', bestPLL / 1e6);
fprintf('Actual F_BIOZ: %.6f Hz (error %.4f%%)\n', best.F_BIOZ_eff, bestErr * 100);
fprintf('C_BIOZ: %.3f\n', best.C_BIOZ);
fprintf('M Divider:   %d\n', best.M);
fprintf('K Divider:   %d\n', best.KDIV);
fprintf('DAC OSR:     %d\n', best.DAC_OSR);
fprintf('N Divider:   %d\n', best.NDIV);
fprintf('ADC OSR:     %d\n', best.ADC_OSR);
fprintf('Sample rate: %.3f Hz\n', best.F_SAMPLE);
