vars = load('vars');
newVars = load('vars_afterFF');

CCSDAT = vars.CCSDAT;
WALL = vars.WALL;
SENSOR_NPRB = vars.SENSOR_NPRB;
SENSOR_FLXLP = vars.SENSOR_FLXLP;
SENSOR_TPRB = vars.SENSOR_TPRB;

figure()
plot(WALL.REV, WALL.ZEV)
hold on
plot(WALL.RES, WALL.ZES)

sum(CCSDAT.NCCN)
sum(CCSDAT.NCCN)
sum(WALL.KNN)
sum(WALL.KSN)


NMAX = SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM + sum(CCSDAT.NCCN)
JMAX = sum(CCSDAT.NCCN) + sum(CCSDAT.NCCN) + sum(WALL.KNN) + sum(WALL.KSN)

% 90 x 87 センサー数90　未知数87
AA = zeros(NMAX, JMAX);

newAA = newVars.AA;
newFF = newVars.FF;

figure()
plot(1:length(FF), FF)
hold on 
plot(1:length(newFF), newFF);
plot(newVars.FC)
legend('FF', 'newFF', 'FC')
[x y] = size(newAA);
figure() 
contour(newAA,'ShowText', 'on')

