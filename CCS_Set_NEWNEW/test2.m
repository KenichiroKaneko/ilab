
sensordata_B0 = fileread('CCS_temporary/@UTST_SenPos.txt');
sensordata_B =  strsplit(sensordata_B0,{'\n',' ','\r'});
sensornum_B =   length(sensordata_B)/2 -1 ;

for i=1:sensornum_B
    R(i) = str2double(sensordata_B{i});
    Z(i) = str2double(sensordata_B{i+1});
end
figure()
plot(R, Z);
R = 0;
Z = 0;
