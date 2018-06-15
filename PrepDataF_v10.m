%PrepDataF_v10: copied from _v9 by Matt Day, 06,06,2018


function  [MatrixX1,MatrixY1,MatrixX2,MatrixY2,MatrixX3,MatrixY3,MatrixX4,MatrixY4,MatrixX5,MatrixY5,MatrixX6,MatrixY6] = PrepDataF_v10(file_path)

Data_path = [file_path 'MD_output.txt'];
parameters_path = [file_path 'MD_parameters.txt'];
coordinate_path = [file_path 'MD_Calculated_coordinates.txt'];

parameters = FindParameters2D_v5(parameters_path);

Data = dlmread(Data_path);
power_path = [file_path 'MD_Power_measured.txt'];
power_data = load(power_path);
[a,b] = size(power_data);

coordinates = load(coordinate_path);

%These parameters define the size of the DataL1
NumSteps_t = parameters(9,:); 
NumSteps_tau = parameters(7,:); 
NumSteps_T = parameters(8,:); 
NumSteps_V =  parameters(28,1);
NumSteps_aux = parameters(31,1);
NumSteps_aux2 = parameters(34,1);
NumSteps_pwr = a;
NumSteps_scan = parameters(12,:);

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
RefFreq1= Data(:,13); %applying to Demodulators 1-3, depricated
RefFreq2= Data(:,14); %applying to Demodulators 4-6, depricated
AuxIn0= Data(:,15);
AuxIn1= Data(:,16);


%inserting zeros for incomplete scans
[m n] = size(X1);
if(m <= NumSteps_t*NumSteps_tau*NumSteps_T*NumSteps_V*NumSteps_aux*NumSteps_pwr*NumSteps_aux2)
    X1(m+1:(NumSteps_scan)) = 0;
    X2(m+1:(NumSteps_scan)) = 0;
    X3(m+1:(NumSteps_scan)) = 0;
    X4(m+1:(NumSteps_scan)) = 0;
    X5(m+1:(NumSteps_scan)) = 0;   
    X6(m+1:(NumSteps_scan)) = 0;
    Y1(m+1:(NumSteps_scan)) = 0;
    Y2(m+1:(NumSteps_scan)) = 0;
    Y3(m+1:(NumSteps_scan)) = 0;
    Y4(m+1:(NumSteps_scan)) = 0;
    Y5(m+1:(NumSteps_scan)) = 0;   
    Y6(m+1:(NumSteps_scan)) = 0;
    RefFreq1(m+1:(NumSteps_scan)) = 0;
    RefFreq2(m+1:(NumSteps_scan)) = 0;
    AuxIn0(m+1:(NumSteps_scan)) = 0;
    AuxIn1(m+1:(NumSteps_scan)) = 0;
end

%Each quadrant of each demodulator gets its own 7-dimensional matrix, set
%by the vectors (X1...Y6), but is no longer than needed incase scan is
%incomplete
for i = 1:NumSteps_scan
    if coordinates(i,2) ~= round(coordinates(i,2))
        coordinates(i,2) = 1;
    end
    MatrixX1(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = X1(i);
    MatrixY1(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = Y1(i);
    MatrixX2(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = X2(i);
    MatrixY2(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = Y2(i);
    MatrixX3(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = X3(i);
    MatrixY3(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = Y3(i);
    MatrixX4(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = X4(i);
    MatrixY4(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = Y4(i);
    MatrixX5(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = X5(i);
    MatrixY5(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = Y5(i);
    MatrixX6(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = X6(i);
    MatrixY6(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = Y6(i);
    if i>= m
            if i<NumSteps_scan
                disp('Incomplete scan detected.')
            end
       break
    end
    i=i+1;

end

end

