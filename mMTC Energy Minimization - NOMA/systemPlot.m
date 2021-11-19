function [] = systemPlot(Xm,Ym,Xu,Yu,radius)

figure;
scatter(Xm,Ym,'filled','k')

hold on
scatter(Xu,Yu,'filled','b')

xc = 0;
yc = 0;
theta = linspace(0,2*pi);

hold on
x = radius*cos(theta) + xc;
y = radius*sin(theta) + yc;
plot(x,y,'--b');
axis equal

hold on
plot(xc,yc,'rx','LineWidth',2,'MarkerSize',10);


legend('mMTC users','URLLC users','','Base Station')
% grid on
set(gca,'XColor', 'none','YColor','none')
hold off