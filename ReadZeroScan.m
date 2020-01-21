function [theta1,theta4] = ReadZeroScan(filePathIn)

dat = importdata([filePathIn 'MD_zero.txt']);

x1 = dat(:,1);
y1 = dat(:,2);
x4 = dat(:,7);
y4 = dat(:,8);

if all(x1 == 0) && all(y1 == 0) && all(x4 == 0) && all(y4 == 0)
    disp('Possibly missing phase data')
end

thetaAll1 = cart2pol(x1,y1);
thetaAll4 = cart2pol(x4,y4);

theta1 = mean(unwrap(thetaAll1));
theta4 = mean(unwrap(thetaAll4));

end