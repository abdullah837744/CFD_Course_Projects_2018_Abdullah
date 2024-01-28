%%%%%

load u20.mat;
load u40.mat;
load u60.mat;
load u80.mat;
load u100.mat;
load E.mat;
load div.mat;
load ut.mat;
load Pt.mat;

x = linspace(1,68,68);
y1 = u20(:,68);
y2 = u40(:,68);
y3 = u60(:,68);
y4 = u80(:,68);
y5 = u100(:,68);


figure(1);
plot(y1,x,y2,x,y3,x,y4,x,y5,x);
ylim([1 68]);
legend('t=20%','t=40%','t=60%','t=80%','t=100%');
xlabel('Vx');
ylabel('Y');
h1 = figure(1);
saveas(h1,'steady state velocity.png');

figure(2);
plot(E);
xlabel('Iteration');
ylabel('Energy');
h2 = figure(2);
saveas(h2,'Energy Plot.png');

figure(3);
plot(div);
xlabel('Iteration');
ylabel('Divergence');
h3 = figure(3);
saveas(h3,'Divergence.png');

figure(4);
contourf(ut,25);
colorbar
xlabel('X -->');
ylabel('Y -->');
h4 = figure(4);
saveas(h4,'velocity contour.png');

figure(5);
contourf(Pt,25);
colorbar;
xlabel('X -->');
ylabel('Y -->');
h5 = figure(5);
saveas(h5,'Pressure contour.png');

figure(6);
plot(Pt');
ylim([100 100.16]);
xlim([1 138]);
xlabel('X');
ylabel('P');
h6 = figure(6);
saveas(h6,'steady state Pressure.png');



