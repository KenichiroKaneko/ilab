vars = load("vars");
close all;
PARAM = vars.PARAM;
NE = PARAM.NE;
SENSOR_FLXLP = vars.SENSOR_FLXLP;
SENSOR_NPRB = vars.SENSOR_NPRB;
SENSOR_TPRB = vars.SENSOR_TPRB;
NAPB = SENSOR_TPRB.NUM + SENSOR_NPRB.NUM;
NCCN = vars.CCSDAT.NCCN;
NFLX = SENSOR_FLXLP.NUM;
KNN = vars.WALL.KNN;

% c = 1
% for III = 1:2
%     for I = 1:SENSOR_FLXLP.NUM
%         for k = 1:NE(III)
%             for JJ = 1:3

%                 KK = 3*(k-1) + JJ;
%                 KKK = KK + 3 * sum(NE(1:III - 1));
                
%                 % disp([c III I k JJ KK KKK])
%                 disp([I + NAPB KK + 3 * sum(NE(1:III - 1)) KK + 3 * sum(NE(1:III - 1)) + sum(NCCN)])
%                 c =c+1;
%             end
%         end
%     end
% end

% c = 1
% for I = 1:NAPB
%     for K = 1:NE(III)
%         for JJ = 1:3
%             KK = 3*(K-1)+JJ;
%             KKK = KK + 3 * sum(NE(1:III - 1));
%             % disp([c I K JJ KK KKK])
%             disp([I KK + 3 * sum(NE(1:III - 1)) KK + 3 * sum(NE(1:III - 1)) + sum(NCCN)])
%         end
%     end
% end

% for I = 1:NCCN(III)
%     for K = 1:NE(III)
%         for JJ = 1:3
%             KK = 3 * (K - 1) + JJ;
%             disp([I + sum(NCCN(1:III - 1)) + NAPB + NFLX KK + 3 * sum(NE(1:III - 1)) KK + 3 * sum(NE(1:III - 1)) + sum(NCCN)])
%         end
%     end
%     for J = 1:NCCN(III)
%     end
%     disp([I + sum(NCCN(1:III - 1)) + NAPB + NFLX I + sum(NCCN(1:III - 1)) + sum(NCCN)])
% end


% III = 2;
% I = NCCN(III);
% K = KNN - 1;
% JJ = 3;
% I + sum(NCCN(1:III - 1)) + NAPB + NFLX
% sum(NCCN) * 2 + 1
% sum(NCCN) * 2 + K - 1 + JJ

I = NFLX;
K = vars.WALL.KNE;
JJ= 3;
KK = 3 * (K - 1) + JJ;
disp([I + NAPB sum(NCCN) * 2 + KK])


III = 1;
I = 1;
K = 1;
JJ = 1;
KK = 3 * (K - 1) + JJ;
I
KK + 3 * sum(NE(1:III - 1)) + sum(NCCN)

I + sum(NCCN(1:III - 1)) + NAPB + NFLX
KK + 3 * sum(NE(1:III - 1)) + sum(NCCN)

