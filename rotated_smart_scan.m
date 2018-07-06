tic;
clear
steps = struct;
stepsize = struct;
initial_pos=struct;
is_writetoFile =1;
%set diagonal scan parameters

steps.t_prime= 200;
stepsize.t_prime = 100; %um, must be positive; t and T position vectors will be created with proper sign at the end
stepsize.tau_prime = 20;
cutoff_idx = 0;
cut_neg_t =1;
cut_neg_tau = 1;
shortside = 2*cutoff_idx+1;

%set other scan parameters

steps.T = 1;
stepsize.T = -45;
initial_pos.T= -90;

%Misc scan parameters
initial_pos.V = -5;
stepsize.V= .2;
steps.V= 1;
steps.LCVolt= 1;
% LCVolt = [60,86,150,500,1000,2000,3500];
initialpos.LCVolt = 5;
LCVolt = [];

%Microscope stages NEED TO AGREE WITH LABVIEW
initial_pos.aux1 = -23057; %x (�m)
initial_pos.aux2 = -24722; %y (�m)

stepsize_aux = 0;
stepsize_aux2 = 0;

steps.aux1 = 1;
steps.aux2 = 1;



%creating mask IN ROTATED FRAME
mask = ones(shortside,steps.t_prime);
%constructing coordinates
coords = struct;
[row,col] = find(mask>0);
coords.tau_prime = row-(cutoff_idx+1);
coords.t_prime = col-1;

%constructing primed positions
positions = struct;
positions.tau_prime = (coords.tau_prime)*stepsize.tau_prime;
positions.t_prime = (coords.t_prime)*abs(stepsize.t_prime);
%constructing unprimed using transfomr t=t'+tau', tau = t-tau, sqrts of two
%will divide out since, for example, we have stepsize_tau' =stepsize_tau*sqrt(2), and the transfomation will have a
%1/sqrt(2) out front

positions.t = positions.t_prime+positions.tau_prime;
positions.tau=positions.t_prime-positions.tau_prime;



%cut points before time zero
if cut_neg_tau==1
    cut_tau = find(positions.tau>=0);
else
    cut_tau = 1:length(positions.tau);
end
%cut_tau = find(positions.tau>=0);
positions.tau = positions.tau(cut_tau);
positions.t=  positions.t(cut_tau);

if cut_neg_t==1
    cut_t = find(positions.t>=0);
else
    cut_t = 1:length(positions.t);
end

if cutoff_idx ~=0
    pos = [positions.tau(cut_t),positions.t(cut_t)];
else
    pos = [positions.tau(cut_t)',positions.t(cut_t)'];
end
%calculate coordinate matrix; coordinate = (position/stepsize)+1 of the
%proper axis
stepsize.tau = abs(abs(positions.tau(4)) -abs(positions.tau(3))) ;
stepsize.t =abs(abs(positions.t(4)) - abs(positions.t(3)))
coordinates = [(pos(:,1)/stepsize.tau)+1,(pos(:,2)/abs(stepsize.t))+1];


%% Case Structure for making vectors of unused dimensions
%if we don't want any T (or other) points, we can just create a 1D vector
%which just repeats in value for the mask. WILL BREAK UNLESS UPDATED FOR
%UNUSED DIMENSIONS WHEN REQUIRED (for example a position scan)

if steps.T ==1
    for i = 1:size(pos,1)
        T_position_vector(i) = initial_pos.T;
        T_coordinate_vector(i) = 1;
    end
    T_position_vector = reshape(T_position_vector,[],1);
    T_coordinate_vector=reshape(T_coordinate_vector,[],1);
end

if steps.V == 1
    V_position_vector = [];
    V_coordinate_vector=[];
    for i = 1:size(pos,1)
        V_position_vector(i) = initial_pos.V;
        V_coordinate_vector(i) = 1;
    end
    V_position_vector=reshape(V_position_vector,[],1);
    V_coordinate_vector=reshape(V_coordinate_vector,[],1);
end

if steps.aux1 == 1
    aux_position_vector =[];
    aux_coordinate_vector = [];

    for i = 1:size(pos,1)
        aux_position_vector(i) = initial_pos.aux1;
        aux_coordinate_vector(i) = 1;
    end
    aux_position_vector=reshape(aux_position_vector,[],1);
    aux_coordinate_vector=reshape(aux_coordinate_vector,[],1);
end

if steps.aux2 == 1
    aux2_position_vector=[];
    aux2_coordinate_vector = [];
    for i = 1:size(pos,1)
        aux2_position_vector(i) = initial_pos.aux2;
        aux2_coordinate_vector(i) = 1;
    end
    aux2_position_vector=reshape(aux2_position_vector,[],1);
    aux2_coordinate_vector=reshape(aux2_coordinate_vector,[],1);
end

if steps.LCVolt == 1
    LCVolt_position_vector =[];
    LCVolt_coordinate_vector= [];
    for i = 1:size(pos,1)
        LCVolt_position_vector(i) = initialpos.LCVolt;
        LCVolt_coordinate_vector(i) = 1;
    end
    LCVolt_position_vector = reshape(LCVolt_position_vector,[],1);
    LCVolt_coordinate_vector = reshape(LCVolt_coordinate_vector,[],1);
end
disp('all position vectors done')


%% Building Global Position Array

% Column 1 is tau, 2 is T, 3 is t, 4 is V, 5 is LCVolt, 6 is aux, 7 is aux2

global_position = zeros(size(pos,1),7);
global_position(:,1) = pos(:,1);
global_position(:,2) = T_position_vector;
global_position(:,3) = -pos(:,2);
global_position(:,4) = V_position_vector;
global_position(:,5) = LCVolt_position_vector;
global_position(:,6) = aux_position_vector;
global_position(:,7) = aux2_position_vector;


global_coordinate = zeros(size(pos,1),7);

global_coordinate(:,1) = coordinates(:,1);
global_coordinate(:,2) = abs(T_coordinate_vector);
global_coordinate(:,3) = coordinates(:,2);
global_coordinate(:,4) = V_coordinate_vector;
global_coordinate(:,5) = LCVolt_coordinate_vector;
global_coordinate(:,6) = aux_coordinate_vector;
global_coordinate(:,7) = aux2_coordinate_vector;

if is_writetoFile ==1
    disp('creating files')
    mask_file = strcat('MD_SmartScan_Mask.txt');
    dlmwrite(mask_file,global_position,'\t');

    position_file = strcat('MD_Calculated_Positions.txt');
    dlmwrite(position_file,global_position,'\t');

    coordinate_file = strcat('MD_Calculated_coordinates.txt');
    dlmwrite(coordinate_file,global_coordinate,'\t');
    disp('done')
end
toc;
