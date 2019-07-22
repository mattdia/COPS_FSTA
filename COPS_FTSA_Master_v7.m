%% Updated for Mac from v2a by Chris Smallwood starting 2015-12-02.
%Revised from version3b by Chris Smallwood on 2017-03-16.
%Version 5: revised 2017-04-27 by Chris Smallwood so that ArbAxisPlot is not a separate function.
    %Plots are no longer contour plots.
    %Setting the frequency axis of the FFT used to be wrong, by a pixel or two, but is now corrected. (Might have even been wrong by more than that, in fact...)
    %Separate frequency axes are defined for S1 and for S2. This is necessary because of the asymmetric conventions for fftshift.
    %Now includes an option to start at t != 0, for local oscillator measurements.
%Version 5a:
    %Modified to allow limited segments of t and tau.
    %Modified to allow easier view of 3d scans.
    %Modified for better view of zero-quantum scans -> FFTshift is always on for zero-quantum.
    %Modified to allow viewing simultaneous functions of S1 and S2.
    %Improved axis labels.
    %Added a folder to set dir_path locally, for when the script is saved locally as well.
%Version 5b/5c:
    %Implemented an option for windowing using a Tukey function.
    %Improved the methodology for correcting phase using time0.
%Version 5d:
    %Implemented an option to implement a Gaussian window function for photon echos.
%Version 6: Revised 2017-08-16 by Chris Smallwood.
    %Previous version was backwards in phase conventions!!
    %All fft's are changed to iffts now. Quadratures for PrepDataF_v9 are swapped by making Y quadrature negative.
    %To accomodate this, sign of t0 correction has been changed.
    %Version 6 also includes the ability to plot as a function of frequency instead of just energy.
%Version 6a: Revised 2017-10-25 by Chris Smallwood.
    %Generalize so that program is capable of extracting MD data sets of (for example) 1 time dimension and 1 frequency dimension.

%Version6c: Revised 2017-04-03 by Matt Day
    % Fixed photon echo windowing so the windowing occurs along the proper
    % rotated coordinates in the time-time domain
    % Changed code to save processing parameters
    
 %Version 6e: Revised 2018-08-28 by Matt Day
    %Partial fix of data loading, as prepdata spits out a 7D matrix, but
    %this program only looks at the first 6 dimensions of the data spit out
    %of PrepData, so added the 7th dimension. The last four dimensions
    %still need to be rearranged to make more sense and mesh better with
    %the labview outputs.
    
% Changes in v7:
%   - Changed to FindParameters2D_v6 (for 2D-COPS_v6.vi) (3/12/19)
%   - Changed to PrepDataF_v11 (3/12/19)

%% Variables

clear all; clc; %clf;% Clear variables, close MuPad engine, clear command window.
speedC = 2.99709e+5; % nm/ps, speed of light in air.
speedCvac = 2.99792458e+5; % nm/ps, speed of light in vacuum. For wavemeter measurements.
planck = 4.135667662e-3;  % eV*ps, or eV/THz, from NIST. Uncertainty is in the last 2 digits.

year     = '2019';
month    = '07';
day      = '11';
scan_num = '19';
fig_num  = str2num(scan_num);
% fig_num  = 1;

% ref_freq = speedCvac/736.57082;
% ref_freq = speedCvac/739.93416;
ref_freq = speedCvac/743.32835;

% pre_file_path = 'E:/Data';                    % Lab computer
% pre_file_path = '/Users/Chris2/Desktop/Data'; % Chris (outdated?)
pre_file_path = '/Volumes/cundiff/COPS/Data'; % Mac
% pre_file_path = 'R:/COPS/Data';               % Windows

file_path = [pre_file_path '/' year '/' year '_' month '/' year '_' month '_' day '/scan' scan_num '/'];
% file_path = ['/Volumes/cundiff/COPS/Data/2019/2019_07/2019_07_11/scan19/'];
% file_path = pwd;
% file_path = '.';

Delay_t0_um = 0; %um. Use this for Local oscillator measurement.
isFFTshift = 1;
isPadding = 2; %Pad with zeros up to numpad if set to 1. Pad by factor of 2 if set to 2.
numpad = 1024;  %fft prefers 2^n points
Undersample_win = 0;
isContourPlot = 0;
NbContours=10;  %Sets the number of contours if using contour plots.
CrtlFlags = [0,3,2,0,0,0];
    %Flags correspond to [tau,T,t,V,aux2,aux1],  Flags used to correspond to [tau,T,t,V,aux,pwr] - CLS, 2017-10-25.
    %Value of 0 means do nothing
    %Value of 1 means plot time domain
    %Value of 2 means plot frequency domain for S1/S2.
    %Value of 3 means plot S3 (only for T)
    %Value of 4 means plot ZeroQuantum (only for T)
PlotIndx = [1,1,1,1,1,1]; %Flags correspond to the slice number extracted for elements of CrtlFlags that are not plotted.
StepLimit = [0,0,0]; %Step limit for [tau, T, t]. Entering 0 leaves them at full length.
isCorrectOverallPhase = 1; %Enter 1 to correct everything by the Tstep specified by PhaseCorrectionIndx, 2 to correct each Tstep independently, 0 for no correction.
PhaseCorrectionIndx = 1;
isS1andS2 = 0; %Enter 1 if both S1 and S2 data sets were collected, 0 if only S1.
isFrequencyUnits = 1; %Enter 1 for frequency units (THz). Enter 0 for energy units (meV).
isWindowFunction_tau = 0; %Enter 1 to window along the tau axis.
isWindowFunction_T = 0; %Enter 1 to window along the T axis.
isWindowFunction_t = 0; %Enter 1 to window along the t axis.
isWindowPhotonEcho = 0; %Enter 1 for photon echo windowing across AND along the t/tau diagonal, enter 2 for across.
TukeyAlpha_tau = .5;     % Select a decimal between 0 (no window) and 1 (Hanning window).
TukeyAlpha_T =.8;     % Select a decimal between 0 (no window) and 1 (Hanning window).
TukeyAlpha_t = .5;     % Select a decimal between 0 (no window) and 1 (Hanning window).

sigma_diag = 1000000; %px
sigma_cross_diag = 40; %px
%norm = 1/(sqrt(2*pi)*sigma);
x_offset = 0; %(right: +, left:-) (along t_axis in pixels)
%stdev_window_time = 2.5; %in ps, t axis;
%time_slope = 1; %in ps/ps
%time_offset = -.5; %in ps/ps
isSaveProcessedData =0; %Set to 1 to save processed data.
isRemoveCW = 0;
% Eliminate the dialog box below in favor of hard-coding the values.
% isub = [d(:).isdir];
% nameFolds = {d(isub).name}';
% b = nameFolds{end};
% [m , ~] = size(b);
% default_num = b(5:end);
% prompt = {'Enter Scan Number','Enter Undersampling Ratio','S1 Demodulator','S2 Demodulator'};
% INPUT = inputdlg(prompt,'Input',1,{default_num,'0','1','4'});
% scan_num = INPUT{1};
% Undersample_win = str2num(INPUT{2});
% S1_Demod = INPUT{3};
% S2_Demod = INPUT{4};

%% Prepping

%Load Laser Spectrum
%    laserspec_path = [file_path '00\absFFT_Z.txt'];
%    laserwave_path = [file_path '00\Wavelength.txt'];
%     laser_spec = load(laserspec_path);
%     laser_wav = load(laserwave_path);

parameters_path = [file_path 'MD_parameters.txt'];
% power_path = [file_path 'MD_Power_measured.txt'];
% pwr = load(power_path); % Power comes in from the bs
% [a,b] = size(pwr);

%Call data from PrepDataF_v9 reading in all demodulators at once.
[MatrixX1,MatrixY1,MatrixX2,MatrixY2,MatrixX3,MatrixY3,MatrixX4,MatrixY4,MatrixX5,MatrixY5,MatrixX6,MatrixY6] = PrepDataF_v11(file_path);
parameters = FindParameters2D_v6(parameters_path); %NB! This also gets executed internally in PrepData above. Not sure if could be faster. I think PrepData does not need to be a separate function.
NumSteps_pwr = 1;
%Define Number of Steps, allow for both scan limiting and interrupted scans
%in time axes and scan limiting in all six axes.
NumSteps3d = ones([1,7]);
for (i= 1:7)
   if i<=3
        if StepLimit(i)==0
            NumSteps3d(i) = size(MatrixX1,i);
        elseif StepLimit(i)~=0
            NumSteps3d(i) = StepLimit(i);
        elseif ScanInterrupted(i)~=0
            NumSteps3d(i) = size(MatrixX1,i);
        end
   elseif i>3 && size(MatrixX1,i)~=1
       NumSteps3d(i) = size(MatrixX1,i); %NB! The last two dimensions are swapped relative to what you would think they should look like in LabView. Gross! - CLS, 2017-10-25
%    elseif i>3
%       NumSteps3d(4) =  parameters(28,:);
%       NumSteps3d(5) =  parameters(34,:);
%       NumSteps3d(6) =  parameters(31,:);
   end
end
NumSteps_tau = NumSteps3d(1);
NumSteps_T = NumSteps3d(2);
NumSteps_t = NumSteps3d(3);
NumSteps_V = NumSteps3d(4);
NumSteps_LC = NumSteps3d(5);
NumSteps_aux = NumSteps3d(6); %NB! The last two dimensions are swapped relative to what you would think they should look like in LabView. Gross! - CLS, 2017-10-25
NumSteps_aux2 = NumSteps3d(7);
%NumSteps_tau = parameters(7,:);
%NumSteps_T = parameters(8,:);
%NumSteps_t = parameters(9,:);
%NumSteps_V = parameters(28,:);
%NumSteps_aux = parameters(31,:);
%NumSteps_aux2 = parameters(34,:);
%Define StepSize
tau_stepsize = abs(parameters(4,1));
T_stepsize = abs(parameters(5,1));
t_stepsize = abs(parameters(6,1));
V_stepsize = parameters(27,:);
aux_stepsize = parameters(30,:); %NB! The last two dimensions are swapped relative to what you would think they should look like in LabView. Gross! - CLS, 2017-10-25
aux2_stepsize = parameters(33,:);
%Other Parameters needed
V_init = parameters(26,1);
aux_init = parameters(29,1); %NB! The last two dimensions are swapped relative to what you would think they should look like in LabView. Gross! - CLS, 2017-10-25
aux2_init = parameters(32,1);
%Cut data to fit the step limit

%Define StepMatrix
%StepMatrix = [NumSteps_tau,NumSteps_T,NumSteps_t,NumSteps_V,NumSteps_aux,NumSteps_pwr];
StepMatrix = [NumSteps_tau,NumSteps_T,NumSteps_t,NumSteps_V,NumSteps_LC,NumSteps_aux,NumSteps_aux2]; %NB! The last two dimensions are swapped relative to what you would think they should look like in LabView. Gross! - CLS, 2017-10-25

%Photon echo windowing parameters. These need to be integers.
%pix_slope = time_slope*(t_stepsize/tau_stepsize);
%tau_stepsize_ps = 2e3*tau_stepsize/speedC;
%t_stepsize_ps = 2e3*t_stepsize/speedC;
%pix_offs = round(time_offset/tau_stepsize_ps);
%stdev_window= stdev_window_time/t_stepsize_ps;

%old phase correction. Discarded in favor of phase correction which uses
%correct time zero
% ZS1_m = complex(MatrixX1(1:NumSteps_tau,1:NumSteps_T,1:NumSteps_t,:,:,:),MatrixY1(1:NumSteps_tau,1:NumSteps_T,1:NumSteps_t,:,:,:));
% ZS4_m = complex(MatrixX4(1:NumSteps_tau,1:NumSteps_T,1:NumSteps_t,:,:,:),MatrixY4(1:NumSteps_tau,1:NumSteps_T,1:NumSteps_t,:,:,:));

% ZS1_phase = angle(ZS1_m);
% ZS4_phase = angle(ZS4_m);
% for(l=1:1:NumSteps_aux2)
% for(m=1:1:NumSteps_pwr)
% for(j=1:1:NumSteps_aux)
% for(i=1:1:NumSteps_V)
%     ZS1_phase(:,:,:,i,m,j,l) = mod(ZS1_phase(:,:,:,i,m,j,l) - ZS1_phase(1,1,3,i,m,j,l)+3*pi,2*pi) - pi;
%     ZS4_phase(:,:,:,i,m,j,l) = mod(ZS4_phase(:,:,:,i,m,j,l) - ZS4_phase(1,1,3,i,m,j,l)+3*pi,2*pi) - pi;
% end
% end
% end
% end

%ZS1_m = abs(ZS1_m).*exp(complex(0,1)*ZS1_phase);
%ZS4_m = abs(ZS4_m).*exp(complex(0,1)*ZS4_phase);

%Re-insert correct phase
% i=1;
% while (i<=numel(ZS1_m))
%     if ZS1_m(i) ~= 0
%        ZS1_m(i) = abs(ZS1_m(i)).*exp(complex(0,1).* ZS1_phase(i));
%     end
%     if ZS4_m(i) ~= 0
%        ZS4_m(i) = abs(ZS4_m(i)).*exp(complex(0,1).* ZS4_phase(i));
%     end
%     i=i+1;
% end

%%Define Time/freq etc Axis assuming steps stepped correctly;
StepSizeMatrix = [tau_stepsize,T_stepsize,t_stepsize,V_stepsize,aux_stepsize];
t = (10^3)*2/speedC*(-Delay_t0_um:t_stepsize:(NumSteps_t-1)*t_stepsize-Delay_t0_um); %ps. %Funny conventions for sign of Delay_t0_um because negative delay = positive time.
tau = (10^3)*2/speedC*(0:tau_stepsize:(NumSteps_tau-1)*tau_stepsize);
T = (10^3)*2/speedC*(0:T_stepsize:(NumSteps_T-1)*T_stepsize);
bias = transpose(((V_init:V_stepsize:((NumSteps_V-1)*V_stepsize+V_init))));
aux = ((aux_init:aux_stepsize:((NumSteps_aux-1)* aux_stepsize+aux_init)));
aux2 = ((aux2_init:aux2_stepsize:((NumSteps_aux2-1)* aux2_stepsize+aux2_init)));
[m ,n] = size(bias);
if(m==0)
   bias = V_init;
end

%% Set global phase
%Set global phase using phase at time zero by taking the input data into
%polar coordinates and finding the phase of tau=t=0.
[d1_theta,d1_r] = cart2pol(MatrixX1,MatrixY1);
[d4_theta,d4_r] = cart2pol(MatrixX4,MatrixY4);
[t_zero,idx_t_zero] = min(abs(t));
[tau_zero,idx_tau_zero] = min(abs(tau));
if isCorrectOverallPhase == 2
    for q=1:StepMatrix(2)
        d1_phase_offset = d1_theta(idx_tau_zero,q,idx_t_zero); %Correct each value of T independently.
        d1_theta(:,q,:) = d1_theta(:,q,:)-d1_phase_offset;
        d4_phase_offset = d4_theta(idx_tau_zero,q,idx_t_zero);
    end
elseif isCorrectOverallPhase == 1
    d1_phase_offset = d1_theta(idx_tau_zero,PhaseCorrectionIndx,idx_t_zero); %Correct by the specified index point. Apply globally.
    disp(['Phase offset correction for D1: ',num2str(d1_phase_offset),' radians.'])
    d1_theta = d1_theta-d1_phase_offset;
    d4_phase_offset = d4_theta(idx_tau_zero,PhaseCorrectionIndx,idx_t_zero);
    disp(['Phase offset correction for D4: ',num2str(d4_phase_offset),' radians.'])
    d4_theta = d4_theta-d4_phase_offset;
end
[MatrixX1,MatrixY1]=pol2cart(d1_theta,d1_r);
[MatrixX4,MatrixY4]=pol2cart(d4_theta,d4_r);

ZS1_m = complex(MatrixX1(1:NumSteps_tau,1:NumSteps_T,1:NumSteps_t,:,:,:,:),MatrixY1(1:NumSteps_tau,1:NumSteps_T,1:NumSteps_t,:,:,:,:));
ZS4_m = complex(MatrixX4(1:NumSteps_tau,1:NumSteps_T,1:NumSteps_t,:,:,:,:),MatrixY4(1:NumSteps_tau,1:NumSteps_T,1:NumSteps_t,:,:,:,:));

%% Remove CW

if isRemoveCW == 1
  ZS1_mr = real(ZS1_m);
  ZS1_mi = imag(ZS1_m);
  ZS4_mr = real(ZS4_m);
  ZS4_mi = imag(ZS4_m);
  
  ZS1_mr = ZS1_mr - mean(mean(ZS1_mr));
  ZS1_mi = ZS1_mi - mean(mean(ZS1_mi));
  
  ZS4_mr = ZS4_mr - mean(mean(ZS4_mr));
  ZS4_mr = ZS4_mr - mean(mean(ZS4_mr));
 
  ZS1_m  = ZS1_mr + complex(0,1)*ZS1_mi;
  ZS4_m  = ZS4_mr + complex(0,1)*ZS4_mi;
end

%% Photon Echo
if isWindowPhotonEcho ~= 0
    slope = abs(t_stepsize/tau_stepsize);
    mask = zeros(size(ZS1_m));
    norm = 1;
    for i = 1:size(mask,1)
        i_offs = i+x_offset;
        for j = 1:size(mask,3)
            j_offs= slope*j;
            if isWindowPhotonEcho ==1
            mask(i,:,j) = norm *exp(-.5*((i_offs-j_offs)/(sigma_cross_diag))^2)*...
                    exp(-.5*((i_offs+j_offs)/(sigma_diag))^2);
            elseif isWindowPhotonEcho == 2
                mask(i,:,j) = norm *exp(-.5*((i_offs-j_offs)/(sigma_cross_diag))^2);
            end
            ZS1_m(i,:,j,:,:,:,:) = mask(i,:,j)*ZS1_m(i,:,j,:,:,:,:) ;
    end
end
%     for k = 1:NumSteps_T
%         for j = 1:NumSteps_t
%             offs(j) = round(pix_slope*j)+pix_offs;
%             for i =1:NumSteps_tau
%                 ZS1_m(i,k,j) = ZS1_m(i,k,j)*exp(-(i-offs(j))^2/(2*stdev_window^2));
%             end
%         end
%     end
clear offs i j k j_offs i_offs norm
end

%% Remaking of the function ArbAxisPlot (v4).

pwr = 1;
ZS1_m1 = ZS1_m(:,:,:,:,:,:,:);
ZS4_m1 = ZS4_m(:,:,:,:,:,:,:);

if isWindowFunction_tau
    alpha = TukeyAlpha_tau;
    alphaStepMatrix = round(alpha*StepMatrix(1));
    WindowFunc_tau(1:alphaStepMatrix) = 0.5*(1+cos(pi*(2*(0:alphaStepMatrix-1)/2/alphaStepMatrix-1)));
    WindowFunc_tau(alphaStepMatrix+1:StepMatrix(1)) = 1;
    WindowFunc_tau = flipud(transpose(WindowFunc_tau));
else
    WindowFunc_tau = ones(StepMatrix(1),1);
end
if isWindowFunction_T
    alpha = TukeyAlpha_T;
    alphaStepMatrix = round(alpha*StepMatrix(2));
    WindowFunc_T(1:alphaStepMatrix) = 0.5*(1+cos(pi*(2*(0:alphaStepMatrix-1)/2/alphaStepMatrix-1)));
    WindowFunc_T(alphaStepMatrix+1:StepMatrix(2)) = 1;
    WindowFunc_T = flipud(transpose(WindowFunc_T));
else
    WindowFunc_T = ones(StepMatrix(2),1);
end
if isWindowFunction_t
    alpha = TukeyAlpha_t;
    alphaStepMatrix = round(alpha*StepMatrix(3));
    WindowFunc_t(1:alphaStepMatrix) = 0.5*(1+cos(pi*(2*(0:alphaStepMatrix-1)/2/alphaStepMatrix-1)));
    WindowFunc_t(alphaStepMatrix+1:StepMatrix(3)) = 1;
    WindowFunc_t = flipud(transpose(WindowFunc_t));
else
    WindowFunc_t = ones(StepMatrix(3),1);
end
pp = 1:StepMatrix(1);
qq = 1:StepMatrix(2);
rr = 1:StepMatrix(3);

WindowFunc = ones(size(ZS1_m1,1),size(ZS1_m1,2),size(ZS1_m1,3),size(ZS1_m1,4),size(ZS1_m1,5),size(ZS1_m1,6),size(ZS1_m1,7));
for(pp=1:1:StepMatrix(1));
for(qq=1:1:StepMatrix(2));
for(rr=1:1:StepMatrix(3));
    WindowFunc(pp,qq,rr,:,:,:,:) = WindowFunc_tau(pp)*WindowFunc_T(qq)*WindowFunc_t(rr);
end
end
end
ZS1_m1 = ZS1_m1 .* WindowFunc;
ZS4_m1 = ZS4_m1 .* WindowFunc;

[p,o,u,~,~,~,~] = size(ZS1_m1);
Pad = zeros(StepMatrix(1),StepMatrix(2),StepMatrix(3),StepMatrix(4),StepMatrix(5),StepMatrix(6),StepMatrix(7));
if isPadding
    %Determine if tau should be padded
    if(CrtlFlags(1) == 2)
        if isPadding == 1
            a = numpad;
        elseif isPadding == 2
            a = p*2;
        end
    else
        a = p;
    end
    %Determine if T should be padded
    if(CrtlFlags(2) ==2 | CrtlFlags(2) == 3 | CrtlFlags(2) == 4)
        if isPadding == 1
            b = numpad;
        elseif isPadding == 2
            b = o*2;
        end
    else
        b = o;
    end
    %Determine if t should be padded
    if(CrtlFlags(3) == 2)
        if isPadding == 1
            m = numpad;
        elseif isPadding == 2
            m = u*2;
        end
    else
        m = u;
    end
    Pad = zeros(a,b,m,StepMatrix(4),StepMatrix(5),StepMatrix(6),StepMatrix(7));
end
ZS1 = Pad;
ZS4 = Pad;
ZS1(1:p,1:o,1:u,1:StepMatrix(4),1:StepMatrix(5),1:StepMatrix(6),1:StepMatrix(7))  = ZS1_m1;
ZS4(1:StepMatrix(1),1:o,1:StepMatrix(3),1:StepMatrix(4),1:StepMatrix(5),1:StepMatrix(6),1:StepMatrix(7))  = ZS4_m1;
% ZS1(1:StepMatrix(1),1:56,1:60,1:StepMatrix(4),1:StepMatrix(5),1:StepMatrix(6))  = ZS1_m1;
% ZS4(1:StepMatrix(1),1:StepMatrix(2),1:StepMatrix(3),1:StepMatrix(4),1:StepMatrix(5),1:StepMatrix(6))  = ZS4_m1;
S3_flag = (CrtlFlags(2)==3);
Zero_flag = (CrtlFlags(2)==4);
axisflag= [S3_flag,Zero_flag];

%Generate plot axes.
%The separate S1 and S2 conventions seem to be necessary because of MatLab's annoying FFTshift conventions.
if isPadding == 2
    numpad = 2*length(t);
end
[ E_t,freq_t,lambda_t] = AxisGenF_v2a(t,isPadding,numpad,isFFTshift,StepSizeMatrix(3),ref_freq,Undersample_win,[0,0]);
if isPadding == 2
    numpad = 2*length(T);
end
[ E_T,freq_T,lambda_T] = AxisGenF_v2a(T,isPadding,numpad,isFFTshift,StepSizeMatrix(2),ref_freq,Undersample_win,axisflag);
if isPadding == 2
    numpad = 2*length(tau);
end
if isFFTshift
    [ E_tauS1,freq_tauS1,lambda_tauS1] = AxisGenF_v2a(tau,isPadding,numpad,isFFTshift,StepSizeMatrix(1),-ref_freq,Undersample_win,[0,0]);
else
    [ E_tauS1,freq_tauS1,lambda_tauS1] = AxisGenF_v2a(tau,isPadding,numpad,isFFTshift,StepSizeMatrix(1),-ref_freq-speedC/(2e+3*tau_stepsize),Undersample_win,[0,0]);
end
[ E_tauS2,freq_tauS2,lambda_tauS2] = AxisGenF_v2a(tau,isPadding,numpad,isFFTshift,StepSizeMatrix(1),ref_freq,Undersample_win,[0,0]);

%% Determine Axis to use

i=1;
if(CrtlFlags(1) == 1)
    axis{1} = tau;
    axislabel{1} = '\tau (ps)';
    i = i+1;
elseif(CrtlFlags(1) == 2)
    if isFrequencyUnits == 1
        axis{1} = freq_tauS1;
        axis{3} = freq_tauS2;
        %axislabel{1} = 'E_\tau (meV)';
        axislabel{1} = 'Excitation frequency (THz)';
    else
        axis{1} = E_tauS1;
        axis{3} = E_tauS2;
        %axislabel{1} = 'E_\tau (meV)';
        axislabel{1} = 'Excitation frequency (meV)';
    end
    if isFFTshift
        ZS1 = fftshift(ifft(ZS1,[],1),1);
        ZS4 = fftshift(ifft(ZS4,[],1),1);
    else
        ZS1 = ifft(ZS1,[],1);
        ZS4 = ifft(ZS4,[],1);
    end
    i = i+1;
end
if(CrtlFlags(2) == 1)
    axis{i} = T;
    axislabel{i} = 'T (ps)';
    i=i+1;
elseif(CrtlFlags(2) == 2 | CrtlFlags(2) == 3)
    if isFrequencyUnits == 1
        axis{i} = freq_T;
        axislabel{i} = 'f_T (THz)';
    else
        axis{i} = E_T;
        axislabel{i} = 'E_T (meV)';
    end
    if isFFTshift
        ZS1 = fftshift(ifft(ZS1,[],2),2);
        ZS4 = fftshift(ifft(ZS4,[],2),2);
    else
        ZS1 = ifft(ZS1,[],2);
        ZS4 = ifft(ZS4,[],2);
    end
    i=i+1;
elseif(CrtlFlags(2) == 4) %Always FFTshift for zero-quantum.
    if isFrequencyUnits == 1
        axis{i} = freq_T;
        axislabel{i} = 'Zero-quantum frequency (THz)';
    else
        axis{i} = E_T;
        axislabel{i} = 'Zero-quantum frequency (meV)';
    end
    ZS1 = fftshift(ifft(ZS1,[],2),2);
    ZS4 = fftshift(ifft(ZS4,[],2),2);
    i=i+1;
end
if((CrtlFlags(3) == 1) & (i < 3))
    axis{i} = t;
    axislabel{i} = 't (ps)';
    i=i+1;
elseif((CrtlFlags(3) == 2) & (i < 3))
    if isFrequencyUnits == 1
        axis{i} = freq_t;
        axislabel{i} = 'Detection frequency (THz)';
    else
        axis{i} = E_t;
        axislabel{i} = 'Detection frequency (meV)';
    end
    ZS1 = ifft(ZS1,[],3);
    ZS4 = ifft(ZS4,[],3);
    Delay_t0 = Delay_t0_um * 2e+3/speedC;
    ReducedFreq = freq_t - ref_freq;
    PhaseAdjustment_t0 = exp( -complex(0,1)*2*pi*ReducedFreq*Delay_t0 );
    for r=1:length(PhaseAdjustment_t0)
        ZS1(:,:,r) = ZS1(:,:,r) * PhaseAdjustment_t0(r);
        ZS4(:,:,r) = ZS4(:,:,r) * PhaseAdjustment_t0(r);
    end
    if isFFTshift
        ZS1 = fftshift(ZS1,3);
        ZS4 = fftshift(ZS4,3);
    end
    i=i+1;
end
if((CrtlFlags(4) == 1) & (i < 3))
    axis{i} = bias;
    axislabel{i} = 'Bias';
    i=i+1;
end
if((CrtlFlags(5) == 1) & (i < 3))
%     %Older version that may not be relevant. Calibration for angle axis after spectra. First degree should be 14.1842 degrees.
%     axis{i} = atan( (aux) / (.05*10^6))*(180/pi) +14.1842-6.788;
%     axislabel{i} = 'angle';
    axis{i} = aux2;
    axislabel{i} = 'aux2';
    i=i+1;
end
if((CrtlFlags(6) == 1) & (i < 3))
    axis{i} = aux;
    axislabel{i} = 'aux';
%     %Older version that may not be relevant.
%     axis{i} = pwr;
%     axislabel{i} = 'Power';
    i=i+1;
end

%% Plot the figure.

if (CrtlFlags(1) ~= 0) & (CrtlFlags(3) ~= 0) %Prefer tau vs. t
    Z1plot = ZS1(:,PlotIndx(2),:,PlotIndx(4),PlotIndx(5),PlotIndx(6));
    Z4plot = ZS4(:,PlotIndx(2),:,PlotIndx(4),PlotIndx(5),PlotIndx(6));
elseif (CrtlFlags(2) ~= 0) & (CrtlFlags(3) ~= 0) %Then T vs. t
    Z1plot = ZS1(PlotIndx(1),:,:,PlotIndx(4),PlotIndx(5),PlotIndx(6));
    Z4plot = ZS4(PlotIndx(1),:,:,PlotIndx(4),PlotIndx(5),PlotIndx(6));
elseif (CrtlFlags(1) ~= 0) & (CrtlFlags(2) ~= 0) %Then tau vs. T
    Z1plot = ZS1(:,:,PlotIndx(3),PlotIndx(4),PlotIndx(5),PlotIndx(6));
    Z4plot = ZS4(:,:,PlotIndx(3),PlotIndx(4),PlotIndx(5),PlotIndx(6));
elseif (CrtlFlags(3) ~= 0) & (CrtlFlags(6) ~= 0) %Then t vs. aux
    Z1plot = ZS1(PlotIndx(1),PlotIndx(2),:,PlotIndx(4),PlotIndx(5),:);
    Z4plot = ZS4(PlotIndx(1),PlotIndx(2),PlotIndx(3),PlotIndx(4),PlotIndx(5),PlotIndx(6));
elseif (CrtlFlags(3) ~= 0) & (CrtlFlags(5) ~= 0) %Then t vs. aux2
    Z1plot = ZS1(PlotIndx(1),PlotIndx(2),:,PlotIndx(4),:,PlotIndx(6));
    Z4plot = ZS4(PlotIndx(1),PlotIndx(2),:,PlotIndx(4),:,PlotIndx(6));
elseif (CrtlFlags(1) ~= 0) & (CrtlFlags(6) ~= 0) %Then tau vs. aux
    Z1plot = ZS1(:,PlotIndx(2),PlotIndx(3),PlotIndx(4),PlotIndx(5),:);
    Z4plot = ZS4(:,PlotIndx(2),PlotIndx(3),PlotIndx(4),PlotIndx(5),:);
elseif (CrtlFlags(1) ~= 0) & (CrtlFlags(5) ~= 0) %Then tau vs. aux2
    Z1plot = ZS1(:,PlotIndx(2),PlotIndx(3),PlotIndx(4),:,PlotIndx(6));
    Z4plot = ZS4(:,PlotIndx(2),PlotIndx(3),PlotIndx(4),:,PlotIndx(6));
elseif (CrtlFlags(5) ~= 0) & (CrtlFlags(6) ~= 0) %Then aux2 vs. aux
    Z1plot = ZS1(PlotIndx(1),PlotIndx(2),PlotIndx(3),PlotIndx(4),:,:);
    Z4plot = ZS4(PlotIndx(1),PlotIndx(2),PlotIndx(3),PlotIndx(4),:,:);
elseif (CrtlFlags(3) ~= 0) & (CrtlFlags(4) ~= 0) %Then t vs. V
    Z1plot = ZS1(PlotIndx(1),PlotIndx(2),:,:,PlotIndx(5),PlotIndx(6));
    Z4plot = ZS4(PlotIndx(1),PlotIndx(2),:,:,PlotIndx(5),PlotIndx(6));
end
Z1plot = squeeze(Z1plot);
Z4plot = squeeze(Z4plot);
Z1procc{1} = Z1plot;
Z4procc{1} = Z4plot;
VmaxZ1 = max(max(abs(Z1plot)));
VminZ1 = min(min(abs(Z1plot)));
VmaxZ4 = max(max(abs(Z4plot)));
VminZ4 = min(min(abs(Z4plot)));
axis1 = axis{1};
axis1label = axislabel{1};
axis2 = axis{2};
axis2label = axislabel{2};
if CrtlFlags(1) == 2
    axis1S2 = axis{3};
elseif CrtlFlags(1) == 1
    axis1S2 = axis{1};
end
[m,n] = size(Z1plot);

figure(fig_num);
if isS1andS2
    set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
        'Position', [0.1,0.1,.8,.8] ) ;
else
    set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
        'Position', [0.1,0.1,.8,.4] ) ;
end

if isS1andS2
    ax1 = subplot(2,2,1);
else
    ax1 = subplot(1,2,1);
end
xlim_min = 1;
xlim_max = n;
ylim_min = 1;
ylim_max = m;
% xlim_min = 5;
% xlim_max = 44;
% ylim_min = 40;
% ylim_max = 80;
if(isContourPlot)
    hFig = contourf(axis2(1:n),axis1(1:m),abs(Z1plot),linspace(0,VmaxZ1,NbContours),'linestyle','none');
else
    %hFig = imagesc(axis2(ylim_min:ylim_max),axis1(xlim_min:xlim_max),abs(Z1plot(ylim_min:ylim_max,xlim_min:xlim_max))); set(gca,'Ydir','Normal');
    hFig = imagesc(axis2(xlim_min:xlim_max),axis1(ylim_min:ylim_max),abs(Z1plot(ylim_min:ylim_max,xlim_min:xlim_max)));
    set(gca,'Ydir','Normal');
    %hFig = surf(axis2(xlim_min:xlim_max),axis1(ylim_min:ylim_max),abs(Z1plot(ylim_min:ylim_max,xlim_min:xlim_max)),'EdgeColor','none'); set(gca,'Ydir','Normal');
end
title('S1 Absolute Value')
colormap(sunset)
%colormap(parula)
%colormap(viridis)
%colormap(beach)
%colormap(flipud(bone))
%colormap(gray)
%colormap(flipud(gray))
if (CrtlFlags(1) == 2) & (CrtlFlags(3) == 2)
   % x = linspace(axis2(1),axis2(end),20); y = -x; line(x,y,'Color','White')%,'LineStyle', ':','MarkerSize',16)
end
colorbar();
colormap(sunset)
% ylabel('${\hbar\omega_{\tau}}$', 'interpreter','latex','FontSize',18)
% xlabel('${\hbar\omega_{t}}$', 'interpreter','latex','FontSize',18)
ylabel(axis1label, 'FontSize',12)
xlabel(axis2label,'FontSize',12)

if isS1andS2
    ax2 = subplot(2,2,2);
else
    ax2 = subplot(1,2,2);
end
if(isContourPlot)
    contourf(axis2(1:n),axis1(1:m),real(Z1plot)/VmaxZ1,linspace(-1,1,NbContours),'linestyle','none');
else
    %hFigReal = imagesc(axis2(xlim_min:xlim_max),axis1(ylim_min:ylim_max),angle(Z1plot(ylim_min:ylim_max,xlim_min:xlim_max))); %set(gca,'Ydir','Normal');
    hFigReal = imagesc(axis2(xlim_min:xlim_max),axis1(ylim_min:ylim_max),real(Z1plot(ylim_min:ylim_max,xlim_min:xlim_max)),[-VmaxZ1,VmaxZ1]); set(gca,'Ydir','Normal');
end
title('S1 Real part')
if (CrtlFlags(1) == 2) & (CrtlFlags(3) == 2)
    x = linspace(axis2(1),axis2(end),20); y = -x; line(x,y,'Color','Black','LineStyle', ':')%,'MarkerSize',16)
end
colorbar();
% xlim([1450,1480])
% ylim([-1480,-1450])
ylabel(axis1label, 'FontSize',12)
xlabel(axis2label,'FontSize',12)

if isS1andS2
    subplot(2,2,3)
    if(isContourPlot)
        hFig = contourf(axis2(1:n),axis1S2(1:m),abs(Z4plot),linspace(0,VmaxZ4,NbContours),'linestyle','none');
    else
        %hFig = imagesc(axis2(ylim_min:ylim_max),axis1(xlim_min:xlim_max),abs(Z1plot(ylim_min:ylim_max,xlim_min:xlim_max))); set(gca,'Ydir','Normal');
        hFig = imagesc(axis2(xlim_min:xlim_max),axis1S2(ylim_min:ylim_max),abs(Z4plot(ylim_min:ylim_max,xlim_min:xlim_max))); set(gca,'Ydir','Normal');
    end
    title('S2 Absolute Value')
    x = linspace(axis2(1),axis2(end),20); y = x; line(x,y,'Color','White')%,'LineStyle', ':','MarkerSize',16)
    colorbar();
    ylabel(axis1label, 'FontSize',12)
    xlabel(axis2label,'FontSize',12)

    subplot(2,2,4)
    if(isContourPlot)
        contourf(axis2(1:n),axis1S2(1:m),real(Z4plot)/VmaxZ4,linspace(-1,1,NbContours),'linestyle','none');
    else
        hFigReal = imagesc(axis2(xlim_min:xlim_max),axis1S2(ylim_min:ylim_max),real(Z4plot(ylim_min:ylim_max,xlim_min:xlim_max)),[-VmaxZ4,VmaxZ4]); set(gca,'Ydir','Normal');
    end
    title('S2 Real Part')
    x = linspace(axis2(1),axis2(end),20); y = x; line(x,y,'Color','Black','LineStyle', ':')%,'MarkerSize',16)
    colorbar();
    ylabel(axis1label, 'FontSize',12)
    xlabel(axis2label,'FontSize',12)
end

%% Extra linear plots

% subplot(3,2,5)
% if(isContourPlot)
%     %contourf(axis2(1:n),-axis1(1:m),flipud(abs(Z4plot))/VmaxZ4,linspace(-1,1,NbContours),'linestyle','none')
%     contourf(axis2(1:n),axis1(1:m),(abs(Z4plot))/VmaxZ4,linspace(-1,1,NbContours),'linestyle','none')
% else
% end
% colorbar();
% % ylim([1450,1480])
% % %xlim([1450,1480])
% ylabel('${\hbar\omega_{\tau}}$', 'interpreter','latex','FontSize',18)
% xlabel('${\hbar\omega_{t}}$', 'interpreter','latex','FontSize',18)
% title('linear abs')
% subplot(3,2,6)
% if(isContourPlot)
%     %contourf(axis2(1:n),-axis1(1:m),flipud(real(Z4plot))/VmaxZ4,linspace(-1,1,NbContours),'linestyle','none')
%     contourf(axis2(1:n),axis1(1:m),(real(Z4plot))/VmaxZ4,linspace(-1,1,NbContours),'linestyle','none')
% else
% end
% colorbar();
% % ylim([1450,1480])
% % %xlim([1450,1480])
% ylabel('${\hbar\omega_{\tau}}$', 'interpreter','latex','FontSize',18)
% xlabel('${\hbar\omega_{t}}$', 'interpreter','latex','FontSize',18)
% title('linear re')

%% Save Processed data

OutDataPath = strcat(file_path, 'Processed_Output/');
if isSaveProcessedData


    if ~isdir(OutDataPath) mkdir(file_path, 'Processed_Output'); end
%     dlmwrite([OutDataPath 'ZS1Real.txt'],real(ZS1));
%     dlmwrite([OutDataPath 'ZS1Imag.txt'],imag(ZS1));
    dlmwrite([OutDataPath 'ZS1Real.txt'],squeeze(real(ZS1)));
    dlmwrite([OutDataPath 'ZS1Imag.txt'],squeeze(imag(ZS1)));
%     dlmwrite([OutDataPath 'ZS4Real.txt'],real(ZS4));
%     dlmwrite([OutDataPath 'ZS4Imag.txt'],imag(ZS4));
    dlmwrite([OutDataPath 'ZS4Real.txt'],squeeze(real(ZS4)));
    dlmwrite([OutDataPath 'ZS4Imag.txt'],squeeze(imag(ZS4)));
    dlmwrite([OutDataPath 'axis1.txt'], axis1','precision',8);
    dlmwrite([OutDataPath 'axis2.txt'], axis2','precision',8);
    if CrtlFlags(1) == 1 || CrtlFlags(1) == 2
        dlmwrite([OutDataPath 'axis1S2.txt'], axis1S2');
    end
fid = fopen(strcat(file_path,'Processed_Output/','processing_parameters.txt'),'wt');
fprintf(fid, '%s\n', ['ref_freq = ' num2str(ref_freq)], ['Delay_t0_um = ' num2str(Delay_t0_um)], ...
    ['isFFTshift = ' num2str(isFFTshift)], ['isPadding = ' num2str(isPadding)], ['numpad = ' num2str( numpad)], ...
    ['Undersample_win = ' num2str( Undersample_win)], ['isContourPlot = ' num2str(isContourPlot)], ['NbContours = ' num2str( NbContours)], ...
    ['CrtlFlags = ' mat2str(CrtlFlags)], ['PlotIndx = ' mat2str(PlotIndx)], ['StepLimit = ' mat2str(StepLimit)], ...
    ['isCorrectOverallPhase = ' num2str(isCorrectOverallPhase)], ['PhaseCorrectionIndx = ' mat2str(PhaseCorrectionIndx)], ...
    ['isS1andS2 = ' num2str(isS1andS2)], ['isFrequencyUnits = ' num2str(isFrequencyUnits)], ['isWindowFunction_tau = ' num2str( isWindowFunction_tau)], ...
    ['isWindowFunction_T = ' num2str(isWindowFunction_T )], ['isWindowFunction_t = ' num2str(isWindowFunction_t)], ...
    ['isWindowPhotonEcho = ' num2str(isWindowPhotonEcho)], ['TukeyAlpha_tau = ' num2str(TukeyAlpha_tau)], ...
    ['TukeyAlpha_T = ' num2str(TukeyAlpha_T)], ['TukeyAlpha_t = ' num2str(TukeyAlpha_t)], ['sigma_diag = ' num2str(sigma_diag)], ...
    ['sigma_cross_diag = ' num2str(sigma_cross_diag)], ['x_offset = ' num2str(x_offset)]);
fclose(fid);

if CrtlFlags(1) == 1 || CrtlFlags(1) == 2
    save([OutDataPath 'pData.mat'],'ZS1','ZS4','axis1','axis2','axis1S2', ...
        'ref_freq','Delay_t0_um', ...
        'isFFTshift','isPadding','numpad',...
        'Undersample_win','isContourPlot','NbContours',...
        'CrtlFlags','PlotIndx','StepLimit',...
        'isCorrectOverallPhase','PhaseCorrectionIndx',...
        'isS1andS2','isFrequencyUnits','isWindowFunction_tau',...
        'isWindowFunction_T','isWindowFunction_t',...
        'isWindowPhotonEcho','TukeyAlpha_tau',...
        'TukeyAlpha_T','TukeyAlpha_t','sigma_diag',...
        'sigma_cross_diag','x_offset');
else
    save([OutDataPath 'pData.mat'],'ZS1','ZS4','axis1','axis2', ...
        'ref_freq','Delay_t0_um', ...
        'isFFTshift','isPadding','numpad',...
        'Undersample_win','isContourPlot','NbContours',...
        'CrtlFlags','PlotIndx','StepLimit',...
        'isCorrectOverallPhase','PhaseCorrectionIndx',...
        'isS1andS2','isFrequencyUnits','isWindowFunction_tau',...
        'isWindowFunction_T','isWindowFunction_t',...
        'isWindowPhotonEcho','TukeyAlpha_tau',...
        'TukeyAlpha_T','TukeyAlpha_t','sigma_diag',...
        'sigma_cross_diag','x_offset');
end
end
