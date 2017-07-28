planck = 4.135667662e-3; %eV/THz

M_dir_path = ['R:\MONSTR1\2015\2015_07_14 - phased CG153, offres CG156\2D4\Output'];
filename = strcat(M_dir_path,'\DTMatched.dat');

raw_dat = dlmread(filename);

pump_probe = [raw_dat(:,1), raw_dat(:,3)];

freq_dat = pump_probe(:,1);
sig_dat = pump_probe(:,2);

COPS_dir_path = ['R:/COPS/Data/2017/2017_05/2017_05_10 DQW 5nm/scan26 - high stats S1 3uW/Processed_Output/'];

Z1_phase_corr_man = -1*pi;


Z1_real = dlmread(strcat(COPS_dir_path,'/ZS1Real.txt'));
Z1_imag = dlmread(strcat(COPS_dir_path,'/ZS1Imag.txt'));

Z4_real = dlmread(strcat(COPS_dir_path,'/ZS4Real.txt'));
Z4_imag = dlmread(strcat(COPS_dir_path,'/ZS4Imag.txt'));

t_axis_en = dlmread(strcat(COPS_dir_path,'/axis2.txt'));
det_axis_THz = (t_axis_en/1000)/planck;

Z1_real = Z1_real(1,:);
Z1_imag = Z1_imag(1,:);

[cops_peaks ~]= find(det_axis_THz >= 353 & det_axis_THz <= 355);
[monstr_peaks ~] = find(freq_dat >= 353 & freq_dat <= 355);

real_interp = interp(Z1_real(cops_peaks),round(length(sig_dat(monstr_peaks))/length(Z1_real(cops_peaks))));
imag_interp = interp(Z1_imag(cops_peaks),round(length(sig_dat(monstr_peaks))/length(Z1_real(cops_peaks))));

Z1 = complex(real_interp,imag_interp);
%Z1 = Z1/max(abs(Z1));

sig_peaks = sig_dat(monstr_peaks);
sig_peaks = sig_peaks/(max(abs(sig_peaks)));
sig_peaks = sig_peaks';


   [z1_theta,z1_r]= cart2pol(real_interp/max(abs(real_interp)),imag_interp/max(abs(imag_interp)));
%     
    for i = 1:360
        [z1_theta,z1_r]= cart2pol(real_interp/max(abs(real_interp)),imag_interp/max(abs(imag_interp))); 
        phase_offs = (i*2*pi)/360;
        
        z1_theta = z1_theta+phase_offs;
% 
        [z1x,z1y] = pol2cart(z1_theta,z1_r);
        diff = abs(z1x) - abs(sig_peaks(1:170));
        diffs(i,:)=  diff;
        


    end

%     


    
figure()
surf(diffs,'EdgeColor','none')
grid off
% figure()
% plot(freq_dat,sig_dat)
% 
% figure()
% plot(det_axis_THz,Z1_real)

