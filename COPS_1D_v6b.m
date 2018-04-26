% Updated by Chris Smallwood starting 2015-12-02.
% Changes from v3c:
%   - Using speed of light in air, not vacuum.
%   - Deleted many things that had been commented out.
%   - Updated paths and file saving to be compatible with Mac.
%   - Removed undersampling options. These never got used anyway.
%   - Lots of changes to make the fft phase extraction better.
% Changes in v6:
%   - Use ifft instead of fft.
%   - Corrected frequency sign starting in version 6b to go with ifft.

clear all;
%clf;
%file_path = ['E:/Data/2018/2018_04/2018_04_09/scan22/'];
%file_path = ['/Users/Chris2/Desktop/Data/2015/2015_12/2015_12_01/scan13/'];
file_path = ['/Volumes/cundiff/COPS/Data/2018/2018_04/2018_04_09/scan22/'];
%file_path = ['/Volumes/cundiff/COPS/Data/2018/2017_01/2017_01_02/scan00'];
%file_path = ['/Volumes/cundiff/COPS/Data/2017/2017_11/2017_11_10 inc/scan05/'];
data_path = [file_path '1D_output.txt'];
parameters_path = [file_path '1D_parameters.txt'];
Data = load(data_path);

prompt = {'Demodulator','Measuring tau? (change to 0 if measuring T or t)','Phase gradient option (choose 1, 2, or 3)'};
INPUT = inputdlg(prompt,'Input',1,{'1','0','1'});
Demod = str2num(INPUT{1});
Measuring_tau = str2num(INPUT{2});
phase_gradient_option = str2num(INPUT{3});

speedC = 2.99709e+5; %(nm/ps), speed of light in air. This value is from Wolfram Alpha. 
speedCvac = 2.99792458e+5; % nm/ps, speed of light in vacuum. For wavemeter measurements.
%lambda_ref = 738.49; %lambda reference beam in nm
lambda_ref = 737.81961; %lambda reference beam in nm
ref_freq = speedCvac/lambda_ref; % in THz
%ref_freq = speedC/lambda_ref; % in THz

%% Notes on corrections of phase  
% There are two phase corrections in the code.
  
%One is from the phase gradient (spectral domain) induced by not starting the scan at the
%zero delay. we remove the phase gradient in the spectral domain by either 
%the "reference wave" method or by substracting a linear phase component
%in the spectrum
  
%The other is a constant phase offset. We remove it by setting the phase
%(or Y) to be zero at zero delay in the time domain
  
%phase_gradient_option=3
% phase_gradient_option=1 : Correction of phase gradient using a gaussian
% fit. It is limited by the resolution of the zero delay fit.
      
% phase_gradient_option=2 : Correction of the phase gradient under the
% assumption that zero delay is actually given by the position where the
% array crosses zero.

% phase_gradient_option=3 : Like option 2, but prompts an option to choose
% delay = 0 manually.
        
% no correction of phase gradient for other numbers

%% extract partial data from data matrix

[a b]=size(Data);
data_init=0; %if 0 is choosen, starts from first data point
data_fin=0;  %if 0 is choosen, ends at last data point

if data_fin
    Data_extract=Data(1:data_fin,:);
    Data=Data_extract;
end

if data_init
    Data_extract=Data(data_init:a,:);
    Data=Data_extract;
end

%% retrieving data (6 demodulators)

X1 = Data(:,1);
Y1 = -Data(:,2);
X2 = Data(:,3);
Y2 = -Data(:,4);
X3 = Data(:,5);
Y3 = -Data(:,6);
X4 = Data(:,7);
Y4 = -Data(:,8);
X5 = Data(:,9);
Y5 = -Data(:,10);
X6 = Data(:,11);
Y6 = -Data(:,12);
RefFreq1= Data(:,13); %applying to Demodulators 1-3
RefFreq2= Data(:,14); %applying to Demodulators 4-6
AuxIn0= Data(:,15);
AuxIn1= Data(:,16);

if (Demod == 1)
    X_ = X1;
    Y_ = Y1;
elseif (Demod == 2)
    X_ = X2;
    Y_ = Y2;
elseif (Demod == 3)
    X_ = X3;
    Y_ = Y3;
elseif (Demod == 4)
    X_ = X4;
    Y_ = Y4;
elseif (Demod == 5)
    X_ = X5;
    Y_ = Y5;
elseif (Demod == 6)
    X_ = X6;
    Y_ = Y6;
else
    msgbox('Valid demodulator numbers: 1 to 6');
end

Z_ = complex(X_,Y_);
R_ = abs(Z_);
theta = unwrap(angle(Z_)) ;%(in radians) [-pi,pi];

position_measured = Data(:,12);
position_measured = 1000*position_measured; % Converts mm to um

Parameters = FindParameters1D_v3(parameters_path);
%default values of parameters are 0
%global InitialPosition StepSize NumberOfSteps %These probably don't need to be global variables...
InitialPosition = Parameters(2,1);
StepSize = Parameters(3,1);
[NumberOfSteps n] = size(R_);
p = transpose(0:1:NumberOfSteps-1); %Step index
  
position_calculated =InitialPosition+ p*StepSize;  %Calculated Expected Steps
  
%% choose if using measured or calculated positions

%position=position_measured;
position=position_calculated;

%% fitting time domain signal with a gaussian to find 0 delay

fit_gaussian = fittype('gauss1');
gaussian_fit = fit(position,R_,fit_gaussian);
[W_coef] = coeffvalues(gaussian_fit); %First element is amplitude, second is x0, third is 1/e half-width.
FittedPositionOfZeroDelay = W_coef(2)

%% scaling data into delays [ps], frequency [THz], wavelength [nm]

if(Measuring_tau == 1)
    delay_ps = (position)*( 10^3 / speedC)*2; 
%     if( position(NumberOfSteps) < position(1) ) %Correct wave orientation
%     if measured in the wrong direction. Actually, I think this is not appropriate. Gets taken care of by Fs.
%         X_ = flipud(X_);
%         Y_ = flipud(Y_);
%         Z_ = flipud(Z_);
%         R_ = flipud(R_);
%         theta = flipud(theta);
%         position = flipud(position);
%         delay_ps = flipup(delay_ps);
%     end
else
    delay_ps = -(position)*( 10^3 / speedC)*2; % negative for T/t since stage is mounted backward -> negative position gives positive delay
end

Fs = 1/(delay_ps(2)-delay_ps(1));
NumPnts = length(delay_ps);
ReducedFrequency_THz = (-floor(NumPnts/2):ceil(NumPnts/2)-1)*Fs/NumPnts; %Does the FFTshift by default.
Frequency_THz = ReducedFrequency_THz + ref_freq;
lambda_t = (speedC./Frequency_THz);
Wavelength_nm = (lambda_t);

%% Fourier transform and phase correction

FFT_Z = ifft(Z_);
FFT_Z = fftshift(FFT_Z);

if (phase_gradient_option==1 || phase_gradient_option==2 || phase_gradient_option==3) % Removal of phase gradient in the FFT assuming that zero is calibrated   

        PositionOfZeroDelay = FittedPositionOfZeroDelay;
    if (phase_gradient_option==2)
        PositionOfZeroDelay = 0;
    elseif phase_gradient_option==3
        prompt = {'Position of zero delay in um?'};
        INPUT = inputdlg(prompt,'Input',1,{'0'});
        PositionOfZeroDelay = str2num(INPUT{1});
    end
    [placeholder,PositionOfZeroDelayIndex] = min(abs(position - PositionOfZeroDelay));

    if(Measuring_tau == 1)
        EffectiveTime0_ps = (PositionOfZeroDelay-position(1))*( 10^3 / speedC)*2; 
    else
        EffectiveTime0_ps = -(PositionOfZeroDelay-position(1))*( 10^3 / speedC)*2; % negative for t and T since stage is mounted backward -> negative position gives positive delay
    end
   
    %Correct for zero delay not being at the beginning of the data set...
    FFT_Z = FFT_Z .* exp(-complex(0,1)*2*pi*ReducedFrequency_THz'*EffectiveTime0_ps);
 
    %and for the overall offset...
    phaseoffset2 = theta(ceil(PositionOfZeroDelayIndex));
    if ceil(PositionOfZeroDelayIndex) ~= 1;
        phaseoffset1 = theta(ceil(PositionOfZeroDelayIndex)-1);
    else
        phaseoffset1 = phaseoffset2;
    end
    phaseoffset = mod((phaseoffset1*abs(PositionOfZeroDelayIndex-ceil(PositionOfZeroDelayIndex))+phaseoffset2*abs(PositionOfZeroDelayIndex-ceil(PositionOfZeroDelayIndex)+1)),2*pi);
    FFT_Z = FFT_Z .* exp(-complex(0,1)*phaseoffset);  
   
end

%% Look for the difference between true zero-delay and the closest available step for a given stage
  
n=0;  %Number of  Digits in after (1um) Actual Position
q = 10^n;
Apos = round(FittedPositionOfZeroDelay*q)/q
%pos_dif = (PositionOfZeroDelay - Apos)

%%

fig1 = figure(8);
plot(position,R_,'k')
hold on
 plot(gaussian_fit,'c')
plot(position,real(Z_),'r')
  xlabel('Position (um)')
  ylabel('Volts (V)')
 %     ylim([0 .002])
tdomainr(:,1) = position;
tdomainr(:,2) = R_;
savezpath = [file_path 'TimeDomainR.txt'];
dlmwrite(savezpath,tdomainr);
tdomainx(:,1) = position;
tdomainx(:,2) = X_;
savezpath = [file_path 'TimeDomainX.txt'];
dlmwrite(savezpath,tdomainx);
tdomainy(:,1) = position;
tdomainy(:,2) = Y_;
savezpath = [file_path 'TimeDomainY.txt'];
dlmwrite(savezpath,tdomainy);

%%

fig2 = figure(9);
[m,n] = size(Wavelength_nm);
%figg= plotyy(Wavelength_nm,abs(FFT_Z(1:n)),Wavelength_nm,Z_(angle(FFT_Z(1:n)))/pi);
figg= plotyy(Wavelength_nm,abs(FFT_Z(1:n)),Wavelength_nm,angle(FFT_Z(1:n)));
xlabel('Wavelength (nm)')
ylabel('FFT Ab. Units')
saveas(fig2,[file_path 'spectrum.png'])
%set(figg,'xlim',[750,850]);
%plot(Wavelength_nm,abs(FFT_Z(1:end-1)))
savezpath = [file_path 'absFFT_Z.txt'];
%dlmwrite(savezpath,abs(FFT_Z));
dlmwrite(savezpath,flipud(abs(FFT_Z)));
Wavelength_nm = transpose(Wavelength_nm);
savezpath = [file_path 'Wavelength.txt'];
dlmwrite(savezpath,Wavelength_nm);

dataz(:,1) = Wavelength_nm;
%dataz(:,2) = abs(FFT_Z);
dataz(:,2) = flipud(abs(FFT_Z));
savezpath = [file_path 'data.txt'];
dlmwrite(savezpath,dataz);