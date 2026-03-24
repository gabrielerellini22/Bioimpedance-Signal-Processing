%% MAIN
clear all
clc

I = input('Enter the desired current in Ampere RMS: ');  % Desired current
f = input('Enter the frequency in Hz: ');               % Signal frequency
[Veff, Reff, Ieff, flag_Rext] = Current_Settings(I, f);
fprintf('Parameters found:\n');
fprintf('V: %.4f mV\n', Veff*1e3);
fprintf('Actual I: %.6f uArms\n', Ieff*1e6);

if(flag_Rext)
    fprintf('External R: %.2f Ohm\n', Reff);
else
    fprintf('Internal R: %.2f Ohm\n', Reff);
end