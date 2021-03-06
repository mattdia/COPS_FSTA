%SmartScan_v4
%Author Matthew Day

clear
%This file generates a 7 column long matrix which is read into LabView
%and which can set the scan parameters. This file is to replace a slower
%version of SmartScan which lives in the COPS old software folder.
%Changes after v3: 
%{
v4a: updated the Inhomogeneous windowing to work for an arbitrary stepsize
in t/tau. 

things to do: 
?update such that we can take multiple inhomogeneously
windowed T steps for 3D scans
?update such that we can take multiple steps in the aux directions

%}
   


%% Flags
inhom_flag = 1; %Purely inhomogeneous mask
hom_flag = 0; %Purely homogeneous mask
other_flag = 0; %Change to 1 to load custom scan mask
cont_scan = 1; %Change to 1 for continuous scanning

%% Write parameters to File
is_writetoFile = 1;
%% Scan Parameters

NumPnts_tau = 400;
NumPnts_T = 1;
NumPnts_t = 800;
stepsize_tau= 240;
stepsize_T = -30;
stepsize_t= -120;
tau_init= 0;
T_init= -90;
t_init=0;
t_offset = 0;

tau_cutoff_index = 50; %Number of points on either side of the 
%diagonal to take for a purely inhomogeneous scan

%Misc scan parameters
V_init = -5;
stepsize_V= .2;
NumPnts_V= 1;
NumPnts_LCVolt= 1;
% LCVolt = [60,86,150,500,1000,2000,3500];
LCVolt_init = 5;
LCVolt = [];

%Microscope stages NEED TO AGREE WITH LABVIEW
aux_init = 20577; %x (�m)
aux2_init = 21948; %y (�m)

stepsize_aux = 0;
stepsize_aux2 = 100;

NumPnts_aux = 1;
NumPnts_aux2 = 1;

 
%Variable initialization
num_pnts_total = NumPnts_tau*NumPnts_T*NumPnts_t;
position_matrix = zeros(NumPnts_t,NumPnts_tau);
t_position_matrix = [];
tau_position_matrix = [];
T_position_matrix = [];
t_composite_matrix = [];
tau_composite_matrix = [];
t_position_matrix = [];
tau_position_vector = [];
t_position_vector = [];
t_abs_position = [];
t_position = [];

mask = zeros(NumPnts_t,NumPnts_tau);

%% no mask 
if inhom_flag==0 && hom_flag==0 &&other_flag==0

    for i = 1:NumPnts_tau
        for j = 1:NumPnts_t
            t_position_matrix(i,j) = j*stepsize_t;
        end
    end
t_position_vector = reshape(t_position_matrix,[],1);


    for i = 1:NumPnts_tau
        for j = 1:NumPnts_t
            tau_position_matrix(i,j) = i*stepsize_tau;
        end
    end

tau_position_vector = reshape(tau_position_matrix,[],1');
t_position_vector = reshape(t_position_matrix,[],1)';
end
clear i j

%% Inhom masking (along diagonal pixel, plus/minus index cutoff in t direction)
if inhom_flag == 1    
    %this part of the code works for photon echo masks of unequal t,tau
    %stepsize. 
        
        for i = 1:NumPnts_t
            
                tau_center = floor(abs((i-1)*stepsize_t)/abs(stepsize_tau))*stepsize_tau;
                tau_center_non_integer = (abs((i-1)*stepsize_t)/abs(stepsize_tau))*stepsize_tau;
                tau_upper = tau_center+abs(tau_cutoff_index*stepsize_tau);
                tau_lower = tau_center-abs(tau_cutoff_index*stepsize_tau);
               
                if tau_cutoff_index ==0
                    vec = tau_lower;
                else 
                   for j = 1:(2*tau_cutoff_index)+1
                        vec(j) = tau_lower + (j-1)*stepsize_tau;
                   end
                end
                   
                   
                 vec(vec/stepsize_tau < 0) = NaN;
                 vec(abs(vec)> abs(NumPnts_t * stepsize_t)) = NaN;
                 
                 t_vec(1:size(vec,2)) = (i-1)*stepsize_t;
                 
                t_position_vector = [t_position_vector, t_vec];
                tau_position_vector = [tau_position_vector,vec];
                
                
        
        end
        corrected_tau_idx = find(isnan(tau_position_vector)==0);
        tau_position_vector = tau_position_vector(corrected_tau_idx);
        t_position_vector = t_position_vector(corrected_tau_idx);
end       
        


        tau_position_vector = tau_position_vector';
        t_position_vector = t_position_vector';
    %   t_position_matrix = t_position_matrix'; %needs to be transposed because t is the slow axis.

     %These matricies have nonzero offsets because the masking program cuts out zero elements, 
     %so we need to get rid of offsets to start at zero.

        disp('t/tau mask done')

%% Pure Homogeneous windowing
%want to scan over window defined by a right triangle (on the time-time
%plane) subtended by the tau and t axes.


if hom_flag==1
    mask = ones(NumPnts_tau,NumPnts_t);
for i = 1:NumPnts_tau
    for j = 1:NumPnts_t
        if (j-1) + round(NumPnts_tau/NumPnts_t)*(i-1) >= NumPnts_tau 
            mask(i,j)=0;
        end
    end
end

    md_mask = [];
if NumPnts_T ~=1
    mask_threed = zeros(NumPnts_tau,NumPnts_t,NumPnts_T);
    
    for k = 1:NumPnts_T
        md_mask(:,:,k) = mask;
    end
    
    [row, col, page] = ind2sub(size(md_mask),find(md_mask>0));
    tau_coordinate_vector = reshape(row,[],1);%make coordniate vectors
    t_coordinate_vector = reshape(col,[],1);
    T_coordinate_vector = reshape(page,[],1);
    
    for i = 1:numel(tau_coordinate_vector) %make position vectors
        tau_position_vector = (i-1)*tau_coordinate_vector(i);
        t_position_vector = (i-1)*t_coordinate_vector(i);
        T_position_vector = (i-1)*T_coordinate_vector(i);
    end
else
    
    [row, col,~] = find(md_mask>0);
    tau_coordinate_vector = reshape(row,[],1);%make coordniate vectors
    t_coordinate_vector = reshape(col,[],1);   
    for i = 1:numel(tau_coordinate_vector) %make position vectors
        tau_position_vector = (i-1)*tau_coordinate_vector(i);
        t_position_vector = (i-1)*t_coordinate_vector(i);
        T_position_vector = (i-1)*T_coordinate_vector(i);
    end
end
end

%% Loading a custom mask (NEEDS TO BE 3D MATRIX STRUCTURED IN tau,t,T FORMAT)
if other_flag ==1 %takes mask (for now just a 3D matrix with tau,t,T coordinates)
%and converts into proper labview inputfiles
    mask_file = uigetfile %find the mask file
    mask = dlmread(maskfile,'\t'); %read in mask file
    [row, col, page] = ind2sub(size(md_mask),find(md_mask>0)); %find nonzero mask elements
    tau_coordinate_vector = reshape(row,[],1);%make coordniate vectors
    t_coordinate_vector = reshape(col,[],1);
    T_coordinate_vector = reshape(page,[],1);
    
    for i = 1:numel(tau_coordinate_vector) %make position vectors
        tau_position_vector = (i-1)*tau_coordinate_vector(i);
        t_position_vector = (i-1)*t_coordinate_vector(i);
        T_position_vector = (i-1)*T_coordinate_vector(i);
    end
disp('mask loaded')
end


%% Creating Correct T for Inhomogeneous scans (if a 3D scan is desired)
size_tau = size(tau_position_vector,1);
clear k i j mask
if NumPnts_T ~=1 
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

%% Creating Aux2 vector

size_tau = size(tau_position_vector,1);
clear k i j mask
if NumPnts_aux2 ~=1 
    for k = 1:NumPnts_aux2
        aux2_position_matrix(:,k) = [ones(size_tau,1)]*(k-1)*stepsize_aux2 + aux2_init;
        aux2_coordinate_matrix(:,k) = [ones(size_tau,1)]*(k);
        tau_composite_matrix(:,k) = tau_position_vector;
        t_composite_matrix(:,k) = t_position_vector;
    end
aux2_position_vector = reshape(aux2_position_matrix,[],1);
aux2_coordinate_vector = reshape(aux2_coordinate_matrix,[],1);
t_position_vector = reshape(t_composite_matrix,[],1);
tau_position_vector = reshape(tau_composite_matrix,[],1);
end

%% Case Structure for making vectors of unused dimensions
%if we don't want any T (or other) points, we can just create a 1D vector
%which just repeats in value for the mask. WILL BREAK UNLESS UPDATED FOR
%UNUSED DIMENSIONS

if NumPnts_T ==1 
    for i = 1:size(t_position_vector)
        T_position_vector(i) = T_init;
        T_coordinate_vector(i) = 1;
    end
    T_position_vector = reshape(T_position_vector,[],1);
    T_coordinate_vector=reshape(T_coordinate_vector,[],1);
end

if NumPnts_V == 1
    V_position_vector = [];
    V_coordinate_vector=[];
    for i = 1:size(t_position_vector)
        V_position_vector(i) = V_init;
        V_coordinate_vector(i) = 1;
    end
    V_position_vector=reshape(V_position_vector,[],1);
    V_coordinate_vector=reshape(V_coordinate_vector,[],1);
end
           
if NumPnts_aux == 1
    aux_position_vector =[];
    aux_coordinate_vector = [];
    
    for i = 1:size(t_position_vector)
        aux_position_vector(i) = aux_init;
        aux_coordinate_vector(i) = 1;
    end
    aux_position_vector=reshape(aux_position_vector,[],1);
    aux_coordinate_vector=reshape(aux_coordinate_vector,[],1);
end
if NumPnts_aux2 == 1
    aux2_position_vector=[];
    aux2_coordinate_vector = [];
    for i = 1:size(t_position_vector)
        aux2_position_vector(i) = aux2_init;
        aux2_coordinate_vector(i) = 1;
    end
    aux2_position_vector=reshape(aux2_position_vector,[],1);
    aux2_coordinate_vector=reshape(aux2_coordinate_vector,[],1);
end

if NumPnts_LCVolt == 1
    LCVolt_position_vector =[];
    LCVolt_coordinate_vector= [];
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
global_position(:,3) = t_position_vector+t_offset;
global_position(:,4) = V_position_vector;
global_position(:,5) = LCVolt_position_vector;
global_position(:,6) = aux_position_vector;
global_position(:,7) = aux2_position_vector;

if cont_scan
    locs = find([1; diff(tau_position_vector)<0; 1]);
    cont_position = zeros(length(locs) - 1,8);
    for i = 1:length(locs) - 1
        cont_position(i,:) = [global_position(locs(i),1),global_position(locs(i + 1) - 1,1),global_position(locs(i),2:7)];
    end
end

global_coordinate = zeros(size(t_position_vector,1),7);

global_coordinate(:,1) = tau_position_vector/stepsize_tau + 1;
global_coordinate(:,2) = abs(T_coordinate_vector);
global_coordinate(:,3) = abs(t_position_vector/stepsize_t) + 1;
global_coordinate(:,4) = V_position_vector/V_init;
global_coordinate(:,5) = LCVolt_position_vector/LCVolt_init;
global_coordinate(:,6) = aux_position_vector/aux_init;
global_coordinate(:,7) = aux2_coordinate_vector;

if is_writetoFile ==1
    disp('creating files')
    
    mask_file = strcat('MD_SmartScan_Mask.txt');
    dlmwrite(mask_file,global_position,'\t');

    position_file = strcat('MD_Calculated_Positions.txt');
    dlmwrite(position_file,global_position,'\t');
    
    cont_file = strcat('MD_ContinuousScan_Positions.txt');
    dlmwrite(cont_file,cont_position,'\t');

    coordinate_file = strcat('MD_Calculated_coordinates.txt');
    dlmwrite(coordinate_file,global_coordinate,'\t');
    disp('done')
end




