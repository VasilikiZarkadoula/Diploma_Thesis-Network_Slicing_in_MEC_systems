function [] = systemPlot(Xm,Ym,Xu,Yu,radius_m,radius_u)

figure;
scatter(Xm,Ym,'filled','k')

hold on
scatter(Xu,Yu,'filled','b')

hold on
xc = 0;
yc = 0;
theta = linspace(0,2*pi);
x = radius_m*cos(theta) + xc;
y = radius_m*sin(theta) + yc;
plot(x,y,'--k');
axis equal

hold on
xu = radius_u*cos(theta) + xc;
yu = radius_u*sin(theta) + yc;
plot(xu,yu,'--b');
axis equal

hold on
plot(xc,yc,'rx','LineWidth',2,'MarkerSize',10);

legend('mMTC users','URLLC users','','','Base Station')
% grid on
set(gca,'XColor', 'none','YColor','none')
hold off