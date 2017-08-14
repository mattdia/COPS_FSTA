folder = 'E:\Data\2017\2017_07\2017_07_20\scan07';
outputfile = strcat(folder,'\MD_output.txt');
output = dlmread(outputfile);
X1 = output(:,1);
Y1 = output(:,2);

Z1 = X1 + i*Y1;

imagesc(abs(Z1))