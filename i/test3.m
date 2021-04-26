% suffix = "test"
% save_dir = "GSPLOT_OUTPUT";

% fp = fopen(save_dir+'\jeddy_'+suffix+'.txt', 'w')
% fprintf(fp, 'a');


SEN0 = dlmread('i/sensorCoordinate0.txt');
r = SEN0(:, 1);
z = SEN0(:, 2);
figure()
plot(r, z, 'o')
