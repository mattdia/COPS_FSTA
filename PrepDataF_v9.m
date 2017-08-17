function [MatrixX1,MatrixY1,MatrixX2,MatrixY2,MatrixX3,MatrixY3,MatrixX4,MatrixY4,MatrixX5,MatrixY5,MatrixX6,MatrixY6] = PrepDataF_v9(file_path)
% upgraded for single loop (SmartScan)
% upgraded for Zurich Instrument Lock-In data  - Gael 5/6/2014
% upgraded (version 9) for more logical phase conventions - Chris Smallwood 8/16/2017.
% Modified to return data and prevent FTSA from loading it again.

%Script to load and format Data correctly

%The correct format of ZI lock in outputs is
% X1 = Data(:,1);
% Y1 = Data(:,2);
% X2 = Data(:,3);
% Y2 = Data(:,4);
% X3 = Data(:,5);
% Y3 = Data(:,6);
% X4 = Data(:,7);
% Y4 = Data(:,8);
% X5 = Data(:,9);
% Y5 = Data(:,10);
% X6 = Data(:,11);
% Y6 = Data(:,12);
% RefFreq1= Data(:,13); %applying to Demodulators 1-3
% RefFreq2= Data(:,14); %applying to Demodulators 4-6
% AuxIn0= Data(:,15);
% AuxIn1= Data(:,16);

a=1;

% scan_num = '7';
% file_path = ['E:\Data\2013\blobtest\Scan' scan_num '\'];
% file_path = ['E:\Data\2014\New folder\testaux\']
Data_path = [file_path 'MD_output.txt'];
parameters_path = [file_path 'MD_parameters.txt'];
coordinate_path = [file_path 'MD_Calculated_coordinates.txt'];

parameters = FindParameters2D_v5(parameters_path);

Data = load(Data_path);
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

%% retrieving data (6 demodulators)
X1 = Data(:,1);
Y1 = Data(:,2);
X2 = Data(:,3);
Y2 = Data(:,4);
X3 = Data(:,5);
Y3 = Data(:,6);
X4 = Data(:,7);
Y4 = Data(:,8);
X5 = Data(:,9);
Y5 = Data(:,10);
X6 = Data(:,11);
Y6 = Data(:,12);
RefFreq1= Data(:,13); %applying to Demodulators 1-3
RefFreq2= Data(:,14); %applying to Demodulators 4-6
AuxIn0= Data(:,15);
AuxIn1= Data(:,16);

%All the DataL1 should be the same size even if aborted.   To enforce that
%everything is the right size we pad the end with zeros for all DataL1 we use if it did
%not go to the final parameters (i.e. if it stoppped at tau =4 instead of
%tau =10
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

%Here we define our Data matrix

ZMatrixX1 = zeros(NumSteps_tau,NumSteps_T,NumSteps_t,NumSteps_V,NumSteps_pwr,NumSteps_aux,NumSteps_aux2);
ZMatrixY1 = zeros(NumSteps_tau,NumSteps_T,NumSteps_t,NumSteps_V,NumSteps_pwr,NumSteps_aux,NumSteps_aux2);

ZMatrixX2 = zeros(NumSteps_tau,NumSteps_T,NumSteps_t,NumSteps_V,NumSteps_pwr,NumSteps_aux,NumSteps_aux2);
ZMatrixY2 = zeros(NumSteps_tau,NumSteps_T,NumSteps_t,NumSteps_V,NumSteps_pwr,NumSteps_aux,NumSteps_aux2);

ZMatrixX3 = zeros(NumSteps_tau,NumSteps_T,NumSteps_t,NumSteps_V,NumSteps_pwr,NumSteps_aux,NumSteps_aux2);
ZMatrixY3 = zeros(NumSteps_tau,NumSteps_T,NumSteps_t,NumSteps_V,NumSteps_pwr,NumSteps_aux,NumSteps_aux2);

ZMatrixX4 = zeros(NumSteps_tau,NumSteps_T,NumSteps_t,NumSteps_V,NumSteps_pwr,NumSteps_aux,NumSteps_aux2);
ZMatrixY4 = zeros(NumSteps_tau,NumSteps_T,NumSteps_t,NumSteps_V,NumSteps_pwr,NumSteps_aux,NumSteps_aux2);

ZMatrixX5 = zeros(NumSteps_tau,NumSteps_T,NumSteps_t,NumSteps_V,NumSteps_pwr,NumSteps_aux,NumSteps_aux2);
ZMatrixY5 = zeros(NumSteps_tau,NumSteps_T,NumSteps_t,NumSteps_V,NumSteps_pwr,NumSteps_aux,NumSteps_aux2);

ZMatrixX6 = zeros(NumSteps_tau,NumSteps_T,NumSteps_t,NumSteps_V,NumSteps_pwr,NumSteps_aux,NumSteps_aux2);
ZMatrixY6 = zeros(NumSteps_tau,NumSteps_T,NumSteps_t,NumSteps_V,NumSteps_pwr,NumSteps_aux,NumSteps_aux2);

DATA_SIZE_tau_T_t=size(ZMatrixX1);

%Here I define the Data matrices that will be used.
i=1;
while (i<= NumSteps_scan)
    MatrixX1(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = X1(i);
    MatrixY1(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = -Y1(i);
    MatrixX2(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = X2(i);
    MatrixY2(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = -Y2(i);
    MatrixX3(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = X3(i);
    MatrixY3(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = -Y3(i);
    MatrixX4(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = X4(i);
    MatrixY4(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = -Y4(i);
    MatrixX5(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = X5(i);
    MatrixY5(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = -Y5(i);
    MatrixX6(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = X6(i);
    MatrixY6(coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5),coordinates(i,6),coordinates(i,7)) = -Y6(i);
    if i>= m;
        
            if i<NumSteps_scan
                disp('Incomplete scan detected.')
            end
       break
    end
    i=i+1;

end

%End User has to load in the parameters to reshape the data matrices using
%the reshape() command.
 
%save_varstr = {MatrixX1,MatrixY1,MatrixX2,MatrixY2,MatrixX3,MatrixY3,MatrixX4,MatrixY4,MatrixX5,MatrixY5,MatrixX6,MatrixY6};
%save_string = {[file_path 'MatrixX1.txt'],[file_path 'MatrixY1.txt'],[file_path 'MatrixX2.txt'],[file_path 'MatrixY2.txt'],[file_path 'MatrixX3.txt'],[file_path 'MatrixY3.txt'],[file_path 'MatrixX4.txt'],[file_path 'MatrixY4.txt'],[file_path 'MatrixX5.txt'],[file_path 'MatrixY5.txt'],[file_path 'MatrixX6.txt'],[file_path 'MatrixY6.txt']};

%[m n] = size(save_varstr);
%for(i=1:1:n)
%    dlmwrite(save_string{i},save_varstr{i})
%end
end
