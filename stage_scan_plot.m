folder = 'E:\Data\2018\2018_05\2018_05_02\scan08';
outputfile = strcat(folder,'\MD_output.txt');
output = dlmread(outputfile);
X1 = output(:,1);
Y1 = output(:,2);

Z1 = X1 + i*Y1;

