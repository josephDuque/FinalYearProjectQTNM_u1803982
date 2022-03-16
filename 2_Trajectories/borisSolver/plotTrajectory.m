fid0 = fopen("trajectoriesBoris.tsv");
O = textscan(fid0, '%f%f%f%f%f%f%f%f%f%f','EndOfLine','\r\n');
fclose(fid0);

figure
TitleText = "Electron Trajectory";

x = O{2}; y = O{3}; z = O{4};

plot3(x,y,z)
view(10,20)
title(TitleText);
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
ax =  gca;
grid on;

fileSave = strcat('borisTrajectory.png'); 
fig = gcf;            
saveas(fig, fileSave);