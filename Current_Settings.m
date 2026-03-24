%% Function
function [Veff, Reff, Ieff, flag_R] = Current_Settings(I, f)
% Function to calculate actual voltage, resistance, and current
% considering a capacitor in series with the resistor

% RVAR datasheet data
RvarMin = 10;
RvarMax = 5e3;
RvarNtaps = 128;

% Tolerance on the desired current value
max_current_error_perc = 0.005; 
max_current_error = max_current_error_perc * I;

% Voltage and resistance data
V = [354e-3, 177e-3, 70.7e-3, 35.4e-3]; % Nominal voltages [V]
Rembedded = [276.25, 5525, 110.5e3, 2210e3]; % Internal resistances [Ohm]
Rvar = linspace(RvarMin, RvarMax, RvarNtaps); % Variable resistances
R = [Rembedded, Rvar]; % All available resistances

% Series capacitance
C = 4.7e-6;
% C = 47e-9; % Farad

% Initialization
Veff = 0; Ieff = NaN; Reff = NaN; bestIndex = [nan,nan];
flag_R = 0; minDelta = Inf;

% Search for optimal combination
for i = 1:length(V)
    Vc = V(i);
    for j = 1:length(R)
        Rc = R(j);
        Xc = 1 / (2 * pi * f * C);             % Capacitive reactance
        Z = sqrt(Rc^2 + Xc^2);                 % Impedance magnitude
        Icalculated = Vc / Z;                  % Actual current
        
        newDelta = abs(Icalculated - I);
        
        % If current is within tolerance
        if(newDelta <= max_current_error)
            
            % Priority conditions
            is_internal = j <= length(Rembedded);
            best_is_internal = bestIndex(2) <= length(Rembedded);
            
            % Criterion 1: higher voltage
            better_voltage = (Vc > Veff);
            
            % Criterion 2: for same voltage, prefer internal
            better_internal = (Vc == Veff) && is_internal && ~best_is_internal;
            
            % Criterion 3: if same type and voltage, lower current error
            better_current = (Vc == Veff) && (is_internal == best_is_internal) && (newDelta < minDelta);
            
            % If any of the three improvement conditions are met
            if(better_voltage || better_internal || better_current)
                minDelta = newDelta;
                bestIndex = [i, j];
                Veff = Vc;
                Reff = Rc;
                Ieff = Icalculated;
            end
        end
    end
end

% Handle special cases
if(isnan(bestIndex(2)))
    Vc = max(V);
    Xc = 1 / (2 * pi * f * C);
    Z_target = Vc / I;
    if(Z_target <= Xc)
        Reff = RvarMin;
        flag_R = -2; % Current too high -> R below minimum
    else
        Reff = sqrt(Z_target^2 - Xc^2);
        flag_R = -1; % Resistance not available among discrete values
    end
    Veff = Vc;
    Ieff = Veff / sqrt(Reff^2 + Xc^2);
end
% External R identification
if(bestIndex(2) > length(Rembedded))
    flag_R = 1; % External R
end
end