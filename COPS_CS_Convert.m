% This takes continuous scannning data and converts it to a format that can
% be read by COPS_FTSA_Master_v7.

% Inputs:
%   filePathIn:       Directory of scan (e.g. '.../scan03/')
%   stepSize:         Size of steps along scanned axis

% Note: does not properly deal with demodulators 2 through 6.

function [res] = COPS_CS_Convert(filePathIn,stepSize)
%% Import Data

scanPos = importdata([filePathIn 'CS_scanning_positions.txt'],'\t',2);
scanAx = scanPos.textdata{1}(12:end);
sweepsPerLoc = str2num(scanPos.textdata{2}(16:end));
scanPos = scanPos.data;

if strcmp(scanAx,'tau')
    axNum = 1;
elseif strcmp(scanAx,'T')
    axNum = 2;
elseif strcmp(scanAx,'t')
    axNum = 3;
end

scanCrd = importdata([filePathIn 'CS_scanning_coordinates.txt']);

recPos = importdata([filePathIn 'CS_Recorded_Positions.txt']);

data = importdata([filePathIn 'CS_output.txt']);

zero = importdata([filePathIn 'CS_zero.txt']);

zeroTheta = cart2pol(zero(1),zero(2));

%% Output Files

filePathOut = [filePathIn 'Converted/'];

if ~exist(filePathOut, 'dir')
    mkdir(filePathIn,'Converted');
end

outputFile = [filePathOut 'MD_output.txt'];
FID = fopen(outputFile,'w');
fclose(FID);

positionsFile = [filePathOut 'MD_Calculated_Positions.txt'];
FID = fopen(positionsFile,'w');
fclose(FID);

maskFile = [filePathOut 'MD_SmartScan_Mask.txt'];
FID = fopen(maskFile,'w');
fclose(FID);

coordinatesFile = [filePathOut 'MD_Calculated_coordinates.txt'];
FID = fopen(coordinatesFile,'w');
fclose(FID);

FID = fopen([filePathOut 'MD_Power_measured.txt'],'w');
fclose(FID);

%% Organize and Interpolate

if scanPos(1,axNum) < scanPos(1,7)
    first = min(scanPos(:,axNum));
    last = max(scanPos(:,7));
    
    newPos = first:stepSize:last;
    newCrd = 1:length(newPos);
    
    sweepInd = find(diff(recPos(:,1)) < -1);
else
    first = max(scanPos(:,axNum));
    last = min(scanPos(:,7));
    
    newPos = first:-stepSize:last;
    newCrd = 1:length(newPos);
    
    sweepInd = find(diff(recPos(:,1)) > 1);
end

sweepInd = [[1; sweepInd + 1] [sweepInd; length(recPos)]];
tooSmall = find((sweepInd(:,2) - sweepInd(:,1)) < 2);
if ~isempty(tooSmall)
    disp('Removing short sweeps');
    sweepInd(tooSmall:end,:) = [];
end

[theta,r] = cart2pol(data(:,1),data(:,2));
data(:,[1 2]) = [unwrap(theta) r];
data(:,1) = data(:,1) - zeroTheta;

numSweeps = length(sweepInd(:,1));

numScanPts = 0;

toWriteDat = [];
toWritePos = [];
toWriteCrd = [];

for i = 1:floor(numSweeps / sweepsPerLoc)
    if mod(i,25) == 0
        disp([num2str(i * sweepsPerLoc/numSweeps * 100) '%']);
    end
    
    start = scanPos(i,axNum);
    finish = scanPos(i,7);
    
    [~,startCrd] = min(abs(newPos - start));
    [~,finishCrd] = min(abs(newPos - finish));
    
    newPosCurr = newPos(startCrd:finishCrd);
    newCrdCurr = newCrd(startCrd:finishCrd);
    numPts = length(newPosCurr);
    numScanPts = numScanPts + numPts;
    newDat = zeros(numPts,16,sweepsPerLoc);
    
    for j = 1:sweepsPerLoc
        currPos = recPos(sweepInd(i + j - 1,1):sweepInd(i + j - 1,2),1);
        currDat = data(sweepInd(i + j - 1,1):sweepInd(i + j - 1,2),:);
        
        [currPos,index] = unique(currPos);
        currDat = currDat(index,:);
        
        newDat(:,:,j) = interp1(currPos,currDat,newPosCurr);
        
        [x,y] = pol2cart(newDat(:,1,j),newDat(:,2,j));
        newDat(:,[1 2],j) = [x y];
    end
    
    newDatMean = nanmean(newDat,3);
    newDatMean(isnan(newDatMean)) = 0;
    
    newPosAll = ones(numPts,7) .* scanPos(i,:);
    newPosAll(:,7) = newPosAll(:,4);
    newPosAll(:,axNum) = newPosCurr;
    
    newCrdAll = ones(numPts,6) .* scanCrd(i,:);
    newCrdAll(:,7) = newCrdAll(:,4);
    newCrdAll(:,axNum) = newCrdCurr;
    
    toWriteDat = [toWriteDat; newDatMean];
    toWritePos = [toWritePos; newPosAll];
    toWriteCrd = [toWriteCrd; newCrdAll];
    
%     dlmwrite(outputFile,newDatMean,'-append','delimiter','\t','precision',10);
%     dlmwrite(positionsFile,newPosAll,'-append','delimiter','\t');
%     dlmwrite(maskFile,newPosAll,'-append','delimiter','\t');
%     dlmwrite(coordinatesFile,newCrdAll,'-append','delimiter','\t');
end

dlmwrite(outputFile,toWriteDat,'-append','delimiter','\t','precision',10);
dlmwrite(positionsFile,toWritePos,'-append','delimiter','\t');
dlmwrite(maskFile,toWritePos,'-append','delimiter','\t');
dlmwrite(coordinatesFile,toWriteCrd,'-append','delimiter','\t');

parametersIn=fopen([filePathIn 'CS_parameters.txt']);
parametersOut = fopen([filePathOut 'MD_parameters.txt'],'w');

tline = fgetl(parametersIn);

while tline ~= -1
    if contains(tline, ['Axis ' num2str(axNum) ' initial position [microns]'])
        fprintf(parametersOut, ['Axis ' num2str(axNum) ' initial position [microns]\t' num2str(first) '\n']);
    elseif contains(tline, ['Axis ' num2str(axNum) ' step size [microns]'])
        fprintf(parametersOut, ['Axis ' num2str(axNum) ' step size [microns]\t' num2str(stepSize) '\n']);
    elseif contains(tline, ['Axis ' num2str(axNum) ' number of steps'])
        fprintf(parametersOut, ['Axis ' num2str(axNum) ' number of steps\t' num2str(length(newPos)) '\n']);
    elseif contains(tline,'Number of scan points')
        fprintf(parametersOut, ['Number of scan points ' num2str(numScanPts) '\n']);
    else
        fprintf(parametersOut,[strrep(tline,'\','\\') '\n']);
    end
    
    tline = fgetl(parametersIn);
end
    
fclose(parametersIn);

res = 1;

end