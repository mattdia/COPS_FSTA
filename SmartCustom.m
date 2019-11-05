% Now with fewer for loops!

% Creates SmartScans that can be read and imported by LabVIEW. Save three
% output files in desired scan folder.

% A plot will also be created to help visualize generated SmartScan.

% See inputs below.

% Written by Kelsey Bates.

numSteps = [371 1 741 1 1 1];
initial = [0 -90 0 0 0 0];
stepSize = [240 1 -120 1 1 1];

inhom = 1;            % Purely inhomogeneous mask
iso = 1;              % Isosceles mask (requires inhom)
hom = 0;              % Purely homogeneous mask
homOpt = 2;           % Options for homogeneous scan
                      %   1: Basic, slope of 1, but has odd behavior for unequal axes
                      %   2: Custom slope and rectangle of data along ax2
                      %   3: Triangle defined by end points
plus = 0;             % Plus (+) mask

cutoffIndex = 60;     % Inhomogeneous width (defined by ax1) (for inhom)
endIndex = 5;         % Inhomogeneous width at end of scan (for isosceles)

homWidth = 3;         % Width of rectangle along ax2 (for homOpt = 2)
homSlope = 2;         % Slope of homogeneous window (2 for 2Q) (for homOpt = 2)

plusWidth = 5;        % Width of arm of plus mask (for plus)

ax1 = 1;              % First SmartScan axis (typically 1)
ax2 = 3;              % Second SmartScan axis (typically 3)

final = initial + (numSteps-1) .* stepSize;

%% Create tau, t

% Creates vectors of the coordinates and positions for tau and t.
vecC1 = 1:numSteps(ax1);
vecC2 = 1:numSteps(ax2);
vec1 = initial(ax1):stepSize(ax1):final(ax1);
vec2 = initial(ax2):stepSize(ax2):final(ax2);
 
% Creates matrices of the coordinates and positions for tau and t.
stageCoordinates = zeros(prod(numSteps([ax1 ax2])),2);
stagePositions = zeros(prod(numSteps([ax1 ax2])),2);
stageCoordinates(:,1) = repmat(vecC1',numSteps(ax2),1);
stageCoordinates(:,2) = repelem(vecC2,numSteps(ax1))';
stagePositions(:,1) = repmat(vec1',numSteps(ax2),1);
stagePositions(:,2) = repelem(vec2,numSteps(ax1))';

% Removes points not in an inhomogeneous scan.
if inhom
    if iso == 1
        m = abs(stepSize(ax1)) * cutoffIndex;
        n = abs(stepSize(ax1)) * endIndex;
        d = abs(final(ax2));
        inRange = (...
            ((sign(stepSize(ax2)) * stagePositions(:,2) - m) / (d - m) ...
            > (sign(stepSize(ax1)) * stagePositions(:,1)) / (d - n)) ...
            | ...
            ((sign(stepSize(ax2)) * stagePositions(:,2)) / (d - n) ...
            < (sign(stepSize(ax1)) * stagePositions(:,1) - m) / (d - m)));
    else
        inRange = abs(stagePositions(:,1) - sign(sign(stepSize(1) * stepSize(ax2))) * stagePositions(:,2)) ...
            > abs(stepSize(ax1) * cutoffIndex);
    end
    stageCoordinates(inRange,:) = [];
    stagePositions(inRange,:) = [];
end

% Removes points not in a homogeneous scan.
if hom
    if homOpt == 1
        inRange = abs(stagePositions(:,1)) + abs(stagePositions(:,2)) ...
            > max(abs(final(ax1)),abs(final(ax2)));
    elseif homOpt == 2
        inRange = ((abs(stagePositions(:,1)) + abs(stagePositions(:,2))/homSlope) ...
            > max(abs(stagePositions(:,1)))) ...
            & (abs(stagePositions(:,1)) > abs(stepSize(ax1)) * homWidth);
    elseif homOpt == 3
        inRange = abs(stagePositions(:,1)) / max(abs(final(ax1))) ...
            + abs(stagePositions(:,2)) / max(abs(final(ax2))) ...
            > 1;
    end
    stageCoordinates(inRange,:) = [];
    stagePositions(inRange,:) = [];
end

% Removes points not in a plus scan.
if plus
    inRange = (...
            (abs(stagePositions(:,1)) >= stepSize(ax1) * plusWidth) ...
            & ...
            (abs(stagePositions(:,2)) >= stepSize(ax1) * plusWidth));
    stageCoordinates(inRange,:) = [];
    stagePositions(inRange,:) = [];
end

%% Other parameters

otherax = 1:3;
otherax([ax1 ax2]) = [];

stagePosLen = length(stagePositions);
coordinates = ones(stagePosLen * prod(numSteps([otherax 4:6])),7);
positions = zeros(stagePosLen * prod(numSteps([otherax 4:6])),7);

% Turns the coordinate and position matrices for the stages into the final desired matrices.
coordinates(:,[ax1 ax2]) = repmat(stageCoordinates,prod(numSteps([otherax 4:6])),1);
positions(:,[ax1 ax2]) = repmat(stagePositions,prod(numSteps([otherax 4:6])),1);

% Adds the T coordinates and positions.
coordinates(:,otherax) = repmat(repelem(1:numSteps(otherax), ...
    stagePosLen)', ...
    prod(numSteps(4:end)),1);
positions(:,otherax) = repmat(repelem(initial(otherax):stepSize(otherax):final(otherax), ...
    stagePosLen)', ...
    prod(numSteps(4:end)),1);

% Adds the aux coordinates and positions.
for i = 4:6
    coordinates(:,i) = repmat(repelem(1:numSteps(i), ...
        stagePosLen * prod(numSteps([otherax 4:i-1])))', ...
        prod(numSteps(i+1:end)),1);
    positions(:,i) = repmat(repelem(initial(i):stepSize(i):final(i), ...
        stagePosLen * prod(numSteps([otherax 4:i-1])))', ...
        prod(numSteps(i+1:end)),1);
end

coordinates(:,5:7) = coordinates(:,4:6);
positions(:,5:7) = positions(:,4:6);

% Finds the total number of points to scan.
numPoints = length(positions);

%% Save files

mask_file = 'MD_SmartScan_Mask.txt';
position_file = 'MD_Calculated_Positions.txt';
dlmwrite(mask_file,positions,'\t');
dlmwrite(position_file,positions,'\t');

coordinate_file = 'MD_Calculated_coordinates.txt';
dlmwrite(coordinate_file,coordinates,'\t');

figure(43);clf;plot(-stagePositions(:,2),-stagePositions(:,1),'.');