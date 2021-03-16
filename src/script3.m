% load("/Volumes/KODAK/data/Taxi_070220/compress0.40121.dat",'-mat','x0','y0','x1','y1');
% xl = x1;
% yl = y1;
% load("/Volumes/KODAK/data/Taxi_070220/compress0.53677.dat",'-mat','x0','y0','x1','y1');
% xh = x1;
% yh = y1;
% 
% 
% subplot(1,2,1);
% plot(x0,y0,xl,yl,xh,yh);
% 
% subplot(1,2,2);
load("/Volumes/KODAK/data/tianchi/VIS/hy_round2_train_20200225/compress0.43064.dat",'-mat','x0','y0','x1','y1');
xl = x1;
yl = y1;
load("/Volumes/KODAK/data/tianchi/VIS/hy_round2_train_20200225/compress0.51879.dat",'-mat','x0','y0','x1','y1');
xh = x1;
yh = y1;
plot(y0,x0,yl,xl,yh,xh);
legend('raw','0.43064','0.51879')
