
fid0 = fopen("helicalTrajectoryDeltaT=1.000e-13.tsv");
O = textscan(fid0, '%f%f%f%f','EndOfLine','\r\n');
fclose(fid0);

figure
TitleText = "Electron Trajectory - pitch angle 87 degrees";

x = O{1}; y = O{2}; z = O{3};

plot3(x,y,z)
view(10,20)
title(TitleText);
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
ax =  gca;
grid on;

fileSave = strcat('idealCyclotronMotion.png'); 
fig = gcf;            
saveas(fig, fileSave);