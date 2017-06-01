%Gaussian window function for photon echo. Assumes you've run the master
%FTSA file

stdev_window = 5;
norm = 1/(sqrt(2*pi)*stdev_window);
pix_slope = (19/20);
pix_offs = -4;
windowed = [];
for k = 1:NumSteps_T
    for j = 1:NumSteps_t
        offs(j) = j-4;
        for i =1:NumSteps_tau
            dummy_filt(i,k,j) = ZS1(i,k,j)*exp(-(i-offs(j))^2/(2*stdev_window^2));
        end
    end
end