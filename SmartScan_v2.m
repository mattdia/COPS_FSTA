%SmartScan_v2
%Author Matthew Day


%This file generates a 7 column long matrix which is read into LabView
%and which can set the scan parameters. This file is to replace a slower
%version of SmartScan which lives in the COPS old software folder.

%% Flags
inhom_flag = 1; %Purely inhomogeneous mask
hom_flag = 0; %Purely homogeneous mask
other_flag = 0; %Change to 1 to load custom scan mask

%% Scan Parameters

NumPnts_tau=10;
NumPnts_T=1;
NumPnts_t=10;
stepsize_tau=45;
stepsize_T=-30;
stepsize_t=-45;
tau_init=0;
T_init=-30;
t_init=0;
t_offset = 0;

tau_cutoff_index=10; %Number of points on either side of the 
%diagonal to take for a purely inhomogeneous scan

%Misc scan parameters
V_init=.2;
stepsize_V=.2;
NumPnts_V=1;
NumPnts_LCVolt=1;
% LCVolt = [60,86,150,500,1000,2000,3500];
LCVolt_init = 5;
LCVolt = [];

%Microscope stages NEED TO AGREE WITH LABVIEW
aux_init = -24259; %x (µm)
aux2_init = -25586; %y (µm)

stepsize_aux = 0;
stepsize_aux2 = 0;

NumPnts_aux = 1;
NumPnts_aux2 = 1;

 
%Variable initialization
num_pnts_total = NumPnts_tau*NumPnts_T*NumPnts_t;
position_matrix = zeros(NumPnts_t,NumPnts_tau);
t_position_matrix = [];
tau_position_vector = [];
t_abs_position = [];
t_position = [];

mask = zeros(NumPnts_t,NumPnts_tau);

%% Full Scan, no masking
if inhom_flag== 0 && hom_flag == 0 && other_flag == 0
    for i = 1:NumPnts_t
        for j = 1:NumPnts_tau
            t_position_matrix(i,j) = j*stepsize_t;
        end
    end
t_position_vector = reshape(t_position_matrix,[],1);
for i = 1:NumPnts_tau
    for j = 1:NumPnts_t
        tau_position_matrix(i,j) = j*stepsize_tau;
    end
end
tau_position_vector = reshape(tau_position_matrix,[],1);


%% Inhomogeneous Masking
elseif inhom_flag==1 && hom_flag==0 && other_flag==0   
mask = zeros(NumPnts_tau,NumPnts_tau+tau_cutoff_index);


for j= 1:NumPnts_tau
    i_low= j-tau_cutoff_index;
    i_high = j+tau_cutoff_index;
    if i_low <= 0
        for i = 1:i_high
            mask(i,j) = (i)*stepsize_tau;
            t_position_matrix(i,j) = (i)*stepsize_t;
            
        end
    else
        for i = i_low:i_high
            mask(i,j) = (i)*stepsize_tau;
            t_position_matrix(i,j) = (i)*stepsize_t;
        end
    end

end


    t_position_matrix = t_position_matrix'; %needs to be transposed because t is the slow axis.
    
 %These matricies have nonzero offsets because the masking program cuts out zero elements, 
 %so we need to get rid of offsets to start at zero.
    tau_position_vector = mask(mask ~=0)-stepsize_tau;  
    t_position_vector = t_position_matrix(t_position_matrix ~= 0)-stepsize_t;
    disp('t/tau mask done')
end
%% Pure homogeneous masking
if hom_flag ==1
    mask = ones(NumPnts_tau,NumPnts_t);
    md_mask = ones(NumPnts_tau,NumPnts_t,NumPnts_T);
    mask_idx = NumPnts_tau+NumPnts_t;
    for i = 1:NumPnts_tau
        for j = 1:NumPnts_t
            if (i-1 + j-1) >= mask_idx/2
                mask(i,j) = 0;
            end
        end
    end
    
    if NumPnts_T ~=1  
        for k = 1:NumPnts_T
            md_mask(:,:,k) = mask;
        end
    
    
        [row, col, page] = ind2sub(size(md_mask),find(md_mask>0)); %find nonzero mask elements
        tau_coordinate_vector = reshape(row,[],1);%make coordniate vectors
        t_coordinate_vector = reshape(col,[],1);
        T_coordinate_vector = reshape(page,[],1);
        T_position_vector = T_coordinate_vector; %Force T to be a column vector, as it was a special case for some reason

        for i = 1:numel(tau_coordinate_vector) %make position vectors, all positions should start from zero
            tau_position_vector(i) = (tau_coordinate_vector(i)-1)*stepsize_tau;
            t_position_vector(i) = (t_coordinate_vector(i)-1)*stepsize_t - t_offset;
            T_position_vector(i) = (T_coordinate_vector(i)-1)*stepsize_t;
        end 
    elseif NumPnts_T ==1 
    [row, col] = find(mask>0);
    tau_coordinate_vector = reshape(row,[],1);%make coordniate vectors
    t_coordinate_vector = reshape(col,[],1);
    
    
        for i = 1:numel(tau_coordinate_vector) %make position vectors
            tau_position_vector(i,1) = (tau_coordinate_vector(i)-1)*stepsize_tau;
            t_position_vector(i,1) = (t_coordinate_vector(i)-1)*stepsize_t - t_offset;

        end   
    end
    disp('done creating hom. mask')
end

%% Other Mask
if other_flag ==1 %takes mask (for now just a 3D matrix with tau,t,T coordinates)
%and converts into proper labview inputfiles
    mask_file = uigetdir %find the mask file
    mask = dlmread(maskfile,'\t'); %read in mask file
    if NumPnts_T ~=1
        [row, col, page] = ind2sub(size(md_mask),find(md_mask>0)); %find nonzero mask elements
        tau_coordinate_vector = reshape(row,[],1);%make coordniate vectors
        t_coordinate_vector = reshape(col,[],1);
        T_coordinate_vector = reshape(page,[],1);
    
        for i = 1:numel(tau_coordinate_vector) %make position vectors
            tau_position_vector(i) = (tau_coordinate_vector(i)-1)*stepsize_tau;
            t_position_vector(i) = (t_coordinate_vector(i)-1)*stepsize_t - t_offset;
            T_position_vector(i) = (T_coordinate_vector(i)-1)*stepsize_t;
        end
    end 
    if NumPnts_T == 1
        [row, col] = find(mask>0);
        tau_coordinate_vector = reshape(row,[],1);%make coordniate vectors
        t_coordinate_vector = reshape(col,[],1);

    
        for i = 1:numel(tau_coordinate_vector) %make position vectors
                tau_position_vector(i) = stepsize_tau*(tau_coordinate_vector(i)-1);
                t_position_vector(i) = stepsize_t*(t_coordinate_vector(i)-1)-t_ofset;

        end
    end
        
end


%% Creating Correct T for Inhomogeneous scans (if a 3D scan is desired)
size_tau = size(tau_position_vector,1);

if NumPnts_T ~=1 && other_flag == 0
    for k = 1:NumPnts_T
        T_position_matrix(:,k) = [ones(size_tau,1)]*(k-1)*stepsize_T;
        tau_composite_matrix(:,k) = tau_position_vector;
        t_composite_matrix(:,k) = t_position_vector;
    end
T_position_vector = reshape(T_position_matrix,[],1);
t_position_vector = reshape(t_composite_matrix,[],1);
tau_position_vector = reshape(tau_composite_matrix,[],1);
disp('all three masks done')
end

%% Case Structure for making vectors of unused dimensions
%if we don't want any T (or other) points, we can just create a 1D vector
%which just repeats in value for the mask

if NumPnts_T ==1 
    for i = 1:size(t_position_vector)
        T_position_vector(i) = T_init;
        T_coordinate_vector(i) = 1;
    end
    T_position_vector=reshape(T_position_vector,[],1);
    T_coordinate_vector=reshape(T_coordinate_vector,[],1);
end

if NumPnts_V == 1
    for i = 1:size(t_position_vector)
        V_position_vector(i) = V_init;
        V_coordinate_vector(i) = 1;
    end
    V_position_vector=reshape(V_position_vector,[],1);
    V_coordinate_vector=reshape(V_coordinate_vector,[],1);
end
           
if NumPnts_aux == 1
    for i = 1:size(t_position_vector)
        aux_position_vector(i) = aux_init;
        aux_coordinate_vector(i) = 1;
    end
    aux_position_vector=reshape(aux_position_vector,[],1);
    aux_coordinate_vector=reshape(aux_coordinate_vector,[],1);
end
if NumPnts_aux2 == 1
    for i = 1:size(t_position_vector)
        aux2_position_vector(i) = aux2_init;
        aux2_coordinate_vector(i) = 1;
    end
    aux2_position_vector=reshape(aux2_position_vector,[],1);
    aux2_coordinate_vector=reshape(aux2_coordinate_vector,[],1);
end

if NumPnts_LCVolt == 1
    for i = 1:size(t_position_vector)
        LCVolt_position_vector(i) = LCVolt_init;
        LCVolt_coordinate_vector(i) = 1;
    end
    LCVolt_position_vector = reshape(LCVolt_position_vector,[],1);
    LCVolt_coordinate_vector = reshape(LCVolt_coordinate_vector,[],1);
end
disp('all position vectors done')
%% Building Global Position Array

% Column 1 is tau, 2 is T, 3 is t, 4 is V, 5 is LCVolt, 6 is aux, 7 is aux2

global_position = zeros(size(t_position_vector,1),7);
global_position(:,1) = tau_position_vector;
global_position(:,2) = T_position_vector;
global_position(:,3) = t_position_vector-t_offset;
global_position(:,4) = V_position_vector;
global_position(:,5) = LCVolt_position_vector;
global_position(:,6) = aux_position_vector;
global_position(:,7) = aux2_position_vector;


global_coordinate = zeros(size(t_position_vector,1),7);

global_coordinate(:,1) = tau_position_vector/stepsize_tau;
global_coordinate(:,2) = abs(T_position_vector/stepsize_T);
global_coordinate(:,3) = abs(t_position_vector/stepsize_t);
global_coordinate(:,4) = V_position_vector/V_init;
global_coordinate(:,5) = LCVolt_position_vector/LCVolt_init;
global_coordinate(:,6) = aux_position_vector/aux_init;
global_coordinate(:,7) = aux2_position_vector/aux2_init;

% disp('creating files')
% mask_file = strcat('MD_SmartScan_Mask.txt');
% dlmwrite(mask_file,global_position,'\t');
% % 
% position_file = strcat('MD_Calculated_Positions.txt');
% dlmwrite(position_file,global_position,'\t');
% % 
% coordinate_file = strcat('MD_Calculated_coordinates.txt');
% dlmwrite(coordinate_file,global_coordinate,'\t');
% disp('done')




