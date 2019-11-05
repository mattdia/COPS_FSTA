% This takes continuous scannning data and converts it to a format that can
% be read by COPS_FTSA_Master_v7.

% Inputs:
%   filePathIn:       Directory of scan (e.g. '.../scan03/')
%   stepSize:         Size of steps along scanned axis
%   polarInterpolate: 0 for interpolation along x and y, 1 for
%                     interpolation along theta and r.

% Note: does not properly deal with demodulators 2 through 6.

function [res] = COPS_CS_Convert(filePathIn,stepSize,polarInterpolate)
%% Import Data

scanPos = importdata([filePathIn 'CS_scanning_positions.txt'],'\t',1);
scanAx = scanPos.textdata{1}(12:end);
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

% copyfile([filePathIn 'CS_parameters.txt'], [filePathOut 'MD_parameters.txt']);

%% Organize and Interpolate

if scanPos(1,axNum) < scanPos(1,7)
    first = min(scanPos(:,axNum));
    last = max(scanPos(:,7));
    
    newPos = first:stepSize:last;
    newCrd = 1:length(newPos);
    
    dir = 1;
else
    first = max(scanPos(:,axNum));
    last = min(scanPos(:,7));
    
    newPos = first:-stepSize:last;
    newCrd = 1:length(newPos);
    
    dir = 0;
end

numSweeps = length(scanPos);

sweepIndCheck = 1:6;
sweepIndCheck(axNum) = [];

sweepInd = find(ismember(recPos(:,sweepIndCheck), scanPos(:,sweepIndCheck),'rows') == 1);
recPos(sweepInd,:) = [];

sweepInd = sweepInd - (0:(numSweeps - 1))';
sweepInd = [sweepInd [sweepInd(2:end) - 1; length(recPos)]];

numScanPts = 0;

for i = 1:numSweeps
    dat = data(sweepInd(i,1):sweepInd(i,2),:);
    pos = recPos(sweepInd(i,1):sweepInd(i,2),1);
    
    start = scanPos(i,axNum);
    finish = scanPos(i,7);
    
    repInd = find(-diff(pos) > finish - start);
    repInd = [[1; repInd + 1] [repInd; length(pos)]];
    numReps = length(repInd);
    
    [~,startCrd] = min(abs(newPos - start));
    [~,finishCrd] = min(abs(newPos - finish));
    
    newPosCurr = newPos(startCrd:finishCrd);
    newCrdCurr = newCrd(startCrd:finishCrd);
    numPts = length(newPosCurr);
    numScanPts = numScanPts + numPts;
    newDat = zeros(numPts,16,numReps);
    
    for j = 1:numReps
        currPos = pos(repInd(j,1):repInd(j,2));
        currDat = dat(repInd(j,1):repInd(j,2),:);    
        
        if polarInterpolate
            [theta,r] = cart2pol(currDat(:,1),currDat(:,2));
            currDat(:,[1 2]) = [theta r];
        end
        
        newDat(:,:,j) = interp1(currPos,currDat,newPosCurr);
        
        if ~polarInterpolate
            [theta,r] = cart2pol(newDat(:,1,j),newDat(:,2,j));
            newDat(:,[1 2],j) = [theta r];
        end
        
        newDat(:,1,j) = newDat(:,1,j) - zeroTheta;
        
        [x,y] = pol2cart(newDat(:,1,j),newDat(:,2,j));
        newDat(:,[1 2],j) = [x y];
    end
    
    newDatMean = mean(newDat,3);
    
    newPosAll = ones(numPts,7) .* scanPos(i,:);
    newPosAll(:,7) = newPosAll(:,4);
    newPosAll(:,axNum) = newPosCurr;
    
    newCrdAll = ones(numPts,6) .* scanCrd(i,:);
    newCrdAll(:,7) = newCrdAll(:,4);
    newCrdAll(:,axNum) = newCrdCurr;
    
    dlmwrite(outputFile,newDatMean,'-append','delimiter','\t','precision',10);
    dlmwrite(positionsFile,newPosAll,'-append','delimiter','\t');
    dlmwrite(maskFile,newPosAll,'-append','delimiter','\t');
    dlmwrite(coordinatesFile,newCrdAll,'-append','delimiter','\t');
end

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