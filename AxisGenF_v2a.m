function [E_t,freq_t,lambda_t] = AxisGenF_v2a(t,isPadding,numpad,isFFTshift,t_stepsize,ref_freq,Undersample_win,flag)
%Generates frequency axis for FTSA etc. using padding or not.   
speedC = 2.99709e+5; %(nm/ps), speed of light in air.
planck = 4.135667662e-3; %(eV*ps), from NIST.
S3_flag = flag(1);
Zero_flag = flag(2);
if isPadding
    tpad = 2e+3/speedC*(0:t_stepsize:(numpad-1)*t_stepsize)'; %ps
else
    tpad = t; %ps
end
Fs = speedC / ( 2e+3*t_stepsize ); %THz
NumPnts = length(tpad);
if isFFTshift | Zero_flag == 1
    freq_t = (-floor(NumPnts/2):ceil(NumPnts/2)-1)*Fs/NumPnts; %THz. See MatLab help examples for FFTshift. Also, note that this should work for both odd-length and even-length arrays, unlike the MatLab example. 
else
    freq_t = (0:NumPnts-1)*Fs/NumPnts;
end 
if(S3_flag == 0 && Zero_flag == 0) %S1 and S2
    freq_t = freq_t + ref_freq - Undersample_win*Fs; %THz
    lambda_t = speedC./freq_t;  %Units in nm
    E_t = freq_t*planck*(10^3); %Units in meV
elseif(S3_flag == 1) %S3
    freq_t = freq_t + 2 * (ref_freq - Undersample_win * Fs); %THz
    lambda_t = speedC./freq_t;  %Units in nm
    E_t = freq_t*planck*(10^3); %Units in meV
elseif(Zero_flag == 1) %Zero Quantum
    freq_t = (freq_t - Undersample_win*Fs); %THz
    lambda_t = speedC./freq_t;  %Units in nm
    E_t = freq_t*planck*(10^3); %Units in meV
end     
end