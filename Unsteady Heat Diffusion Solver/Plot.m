clear 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Amplification factor for Explicit scheme
gamma1 = 0.15;
gamma2= 0.3;
gamma3 = 1.5;
phi=0:pi/10:pi;

G1_Exp = (1-8*gamma1*sin(phi/2).^2);
G2_Exp = (1-8*gamma2*sin(phi/2).^2);
G3_Exp = (1-8*gamma3*sin(phi/2).^2);
figure(1)
plot(phi,G1_Exp,phi,G2_Exp,phi,G3_Exp);


legend('gamma=0.15','gamma=0.3','gamma=1.5');
xlabel('Phi');
ylabel('Amplification factor, G');
ylim([-1 1]);
h1=figure(1)
set(h1, 'Visible', 'off');
saveas(h1,'G_Explicit.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Amplification factor for Implicit scheme
% gamma1 = 0.15;
% gamma2= 0.3;
% gamma3 = 1.5;
% phi=0:pi/10:pi;
% 
G1_Imp = 1./(1+8*gamma1*sin(phi/2).^2);
G2_Imp = 1./(1+8*gamma2*sin(phi/2).^2);
G3_Imp = 1./(1+8*gamma3*sin(phi/2).^2);
figure(2)
plot(phi,G1_Imp,phi,G2_Imp,phi,G3_Imp);


legend('gamma=0.15','gamma=0.3','gamma=1.5');
xlabel('Phi');
ylabel('Amplification factor, G');
ylim([-1 1]);
h2=figure(2)
set(h2, 'Visible', 'off');
saveas(h2,'G_Implicit.png'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Amplification factor for Exact solution
% gamma1 = 0.15;
% gamma2= 0.3;
% gamma3 = 1.5;
% phi=0:pi/10:pi;

G1_Exact = exp(-gamma1*phi.^2);
G2_Exact = exp(-gamma2*phi.^2);
G3_Exact = exp(-gamma3*phi.^2);
figure(3)
plot(phi,G1_Exact,phi,G2_Exact,phi,G3_Exact);


legend('gamma=0.15','gamma=0.3','gamma=1.5');
xlabel('Phi');
ylabel('Amplification factor, G');
ylim([-1 1]);
h3=figure(3)
set(h3, 'Visible', 'off');
saveas(h3,'G_Exact.png'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Error Plot

% load T.mat;
% load T_final.mat;
% 
% for i = 1:100+1
%     for j= 1:100+1
%         Err(i,j)=abs(T_final(i,j)-T(i,j));
%     end
% end
% save('Err.mat','Err');
% x = 1:100+1;
% y = 1:100+1;
% [X,Y] = meshgrid(x,y);
% Z = Err;
% 
% 
% figure(4)
% contourf(X,Y,Z,50);
% colorbar
% h4= figure(4)
% saveas(h4,'Err.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Impact of time-step size (various gamma)

% m=100;
% 
% load T_Imp_gamma1.mat;
% load T_Imp_gamma2.mat;
% load T_Imp_gamma3.mat;
% 
% load T_Exp_gamma1.mat;
% load T_Exp_gamma2.mat;
% 
% 
% 
% T_Imp_half_gamma1 = T_Imp_gamma1(m/2+1,:); % T at y=0.5
% T_Imp_half_gamma2 = T_Imp_gamma2(m/2+1,:); % T at y=0.5
% T_Imp_half_gamma3 = T_Imp_gamma3(m/2+1,:); % T at y=0.5
% 
% T_Exp_half_gamma1 = T_Exp_gamma1(m/2+1,:); % T at y=0.5
% T_Exp_half_gamma2 = T_Exp_gamma2(m/2+1,:); % T at y=0.5
% 
% x = 1:m+1;
% 
% figure(5)
% plot(x,T_Imp_half_gamma1,x,T_Imp_half_gamma2,x,T_Imp_half_gamma3);
% 
% legend('gamma=0.15','gamma=0.3','gamma=1.5','Location','northwest');
% xlabel('X at Y=0.5');
% ylabel('T');
% 
% h5= figure(5)
% saveas(h5,'100x100 Implicit for different gamma at t=0.1.png');
% 
% figure(6)
% plot(x,T_Exp_half_gamma1,x,T_Exp_half_gamma2,x,T_Imp_half_gamma1);
% 
% legend('gamma=0.15','gamma=0.3','Exact','Location','northwest');
% xlabel('X at Y=0.5');
% ylabel('T');
% 
% h6= figure(6)
% saveas(h6,'100x100 Explicit for different gamma at t=0.1.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impact of spatial resolution (various mesh size)


% load T_Imp_Hundred.mat;
% load T_Imp_Fifty.mat;
% load T_Imp_Twenty.mat;
% 
% load T_Exp_Hundred.mat;
% load T_Exp_Fifty.mat;
% load T_Exp_Twenty.mat;
% 
% x1=1:101
% x2=1:51
% x3=1:21
% 
% T_Imp_half_Hundred = T_Imp_Hundred(51,:); % T at y=0.5
% T_Imp_half_Fifty = T_Imp_Fifty(26,:); % T at y=0.5
% T_Imp_half_Twenty = T_Imp_Twenty(11,:); % T at y=0.5
% 
% T_Exp_half_Hundred = T_Exp_Hundred(51,:); % T at y=0.5
% T_Exp_half_Fifty = T_Exp_Fifty(26,:); % T at y=0.5
% T_Exp_half_Twenty = T_Exp_Twenty(11,:); % T at y=0.5
% 
% figure(7)
% plot(x1,T_Imp_half_Hundred,x2,T_Imp_half_Fifty,x3,T_Imp_half_Twenty);
% 
% legend('100x100','50x50','20x20','Location','southeast');
% xlabel('X at Y=0.5');
% ylabel('T');
% ylim([0 1.6]);
% 
% h7= figure(7)
% saveas(h7,'Implicit for different mesh size at gamma=0.15 at t=0.1.png');
% 
% figure(8)
% plot(x1,T_Exp_half_Hundred,x2,T_Exp_half_Fifty,x3,T_Exp_half_Twenty);
% 
% legend('100x100','50x50','20x20','Location','southeast');
% xlabel('X at Y=0.5');
% ylabel('T');
% ylim([0 1.6]);
% h8= figure(8)
% saveas(h8,'Explicit for different mesh size at gamma=0.15 at t=0.1.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impact of numerical method (various time)

% m=100;
% 
% load T_Imp_t1.mat;
% load T_Imp_t2.mat;
% load T_Imp_t3.mat;
% 
% load T_Exp_t1.mat;
% load T_Exp_t2.mat;
% load T_Exp_t3.mat;
% 
% 
% 
% T_Imp_half_t1 = T_Imp_t1(m/2+1,:); % T at y=0.5
% T_Imp_half_t2 = T_Imp_t2(m/2+1,:); % T at y=0.5
% T_Imp_half_t3 = T_Imp_t3(m/2+1,:); % T at y=0.5
% 
% T_Exp_half_t1 = T_Exp_t1(m/2+1,:); % T at y=0.5
% T_Exp_half_t2 = T_Exp_t2(m/2+1,:); % T at y=0.5
% T_Exp_half_t3 = T_Exp_t3(m/2+1,:); % T at y=0.5
% 
% x = 1:m+1;
% 
% figure(9)
% plot(x,T_Imp_half_t1,x,T_Exp_half_t1);
% 
% legend('Implicit','Explicit','Location','northwest');
% xlabel('X at Y=0.5');
% ylabel('T');
% title('at t=0.1');
% ylim([0 1.6]);
% h9= figure(9)
% saveas(h9,'100x100 Implicit_Explicit for at gamma=0.15 at t=0.1.png');
% 
% figure(10)
% plot(x,T_Imp_half_t2,x,T_Exp_half_t2);
% 
% legend('Implicit','Explicit','Location','northwest');
% xlabel('X at Y=0.5');
% ylabel('T');
% title('at t=0.15');
% ylim([0 1.6]);
% h10= figure(10)
% saveas(h10,'100x100 Implicit_Explicit for at gamma=0.15 at t=0.15.png');
% 
% figure(11)
% plot(x,T_Imp_half_t3,x,T_Exp_half_t3);
% 
% legend('Implicit','Explicit','Location','northwest');
% xlabel('X at Y=0.5');
% ylabel('T');
% title('at t=0.2');
% ylim([0 1.6]);
% h11= figure(11)
% saveas(h11,'100x100 Implicit_Explicit for at gamma=0.15 at t=0.2.png');

