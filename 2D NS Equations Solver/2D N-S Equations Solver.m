%%%%%%% 2D N-S equations Solver by Abdullah %%%%%%

Re=100;
Ly=1;
Lx=2;
neu=1e-2;
ro=1;
Vpeak = (Re*neu)/Ly ;
Pout = 100;
Q = Vpeak*Ly;
meu = neu*ro;
dp = -0.16 ;
Pin = Pout - dp ;

MI = 0;           % MI = 1 for with momentum interpolation and 0 for without

gamma = 0.2;
Pe = 1.5;
delx = (Pe*neu)/Vpeak;
dely = delx ;
delt = (gamma*delx^2)/neu ;

nx = round(Lx/delx) + 1 ;
ny = round(Ly/dely) + 1 ;

u = zeros(nx,ny);
v = zeros(nx,ny);
P = zeros(nx,ny);
Hx = zeros(nx,ny);
Hy = zeros(nx,ny);
E = zeros (2447,1);
div = zeros (2447,1);

for i = 1:nx
    for j = 1:ny
        P(1,j)= Pin ;
        P(nx,j)= Pout;
        end
end
for i = 2:nx-1
    for j = 1:ny
P(i,j) = Pin + dp*(i-1)/(nx-1);
    end
end
P(nx+1,:)=Pout + dp/(nx-1);
P(nx+2,:)=Pout + (2*dp/(nx-1));

P1 = zeros(ny,1);
P1(:,1)= Pin - dp/(nx-1);


norm_u=100;
norm_v=100;
tol_vel = 1e-4;
time_iter = 0;

while (norm_u > tol_vel || norm_v > tol_vel)
time_iter = time_iter+1;

u(1,:) = u(nx,:);
u(nx+1,:)=u(2,:);
v(1,:) = v(nx,:);
v(nx+1,:)=v(2,:);
Hx(1,:) = Hx(nx,:);
Hx(nx+1,:)=Hx(2,:);
Hy(1,:) = Hy(nx,:);
Hy(nx+1,:)=Hy(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 2:nx
    for j = 2:ny-1

Hx(i,j) = -(u(i,j)*((u(i,j)-u(i-1,j))/delx)+v(i,j)*((u(i,j)-u(i,j-1))/dely))...
    +(1/Re)*((u(i+1,j)-2*u(i,j)+u(i-1,j))/delx^2 ...
    +(u(i,j+1)-2*u(i,j)+u(i,j-1))/dely^2) ;

Hy(i,j) = -(u(i,j)*((v(i,j)-v(i-1,j))/delx)+v(i,j)*((v(i,j)-v(i,j-1))/dely))...
    +(1/Re)*((v(i+1,j)-2*v(i,j)+v(i-1,j))/delx^2 ...
    +(v(i,j+1)-2*v(i,j)+v(i,j-1))/dely^2) ;

    end
end
err = 0;
norm = 100;
tol = 1e-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Without Momentum %%%%%%%%
if (MI == 0)

while (norm > tol)

Pold = P;
for i = 2:nx-1
    for j = 2:ny-1

P(i,j)=((Hx(i+1,j)-Hx(i-1,j))/(2*delx)+(Hy(i,j+1)-Hy(i,j-1))/(2*dely)...
    -(P(i+1,j)+P(i-1,j))/delx^2 - (P(i,j+1)+P(i,j-1))/dely^2)*...
      (-1/((2/delx^2)+(2/dely^2))) ;

  err = err + abs(Pold(i,j)-P(i,j));

  Pold(i,j)=P(i,j);
    end
end
norm = err/((nx-2)*(ny-2));


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_u = 0;
err_v = 0;
u_new = zeros(nx,ny);
v_new = zeros(nx,ny);

for i = 2:nx
    for j = 2:ny-1

        u_new(i,j)=u(i,j)+delt*(Hx(i,j)-(P(i+1,j)-P(i-1,j))/2*delx);

        v_new(i,j)=v(i,j)+delt*(Hy(i,j)-(P(i,j+1)-P(i,j-1))/2*dely);

        err_u = err_u + abs(u(i,j)-u_new(i,j));
        err_v = err_v + abs(v(i,j)-v_new(i,j));

    end
end

norm_u = err_u/(nx-2)*(ny-2);
norm_v = err_v/(nx-2)*(ny-2);


u(1:nx,1:ny) = u_new(1:nx,1:ny);
u(1,:)= u(nx,:);
v(1:nx,1:ny) = v_new(1:nx,1:ny);
v(1,:)= v(nx,:);
disp(['Time iter ', int2str(time_iter)]);

%%% Energy & Mass%%%%
e = 0;
ut=u';
vt=v';
dsum=0;
for i = 2:nx-1
    for j=2:ny-1
        e= e + (ro * u(i,j)^2);

    end
end

for i=2:ny-1
    for j=2:nx-1
        da = ut(i,j)*(ut(i+1,j)-ut(i-1,j))/(2*delx);
        db = vt(i,j)*(ut(i,j+1)-ut(i,j-1))/(2*dely);
        dc = ut(i,j)*(vt(i+1,j)-vt(i-1,j))/(2*delx);
        dd = vt(i,j)*(vt(i,j+1)-vt(i,j-1))/(2*dely);
        dsum = dsum + da+db+dc+dd ;
    end
end

%%%%%%%%%%

    E(time_iter,1) = e;
    div(time_iter,1) = dsum;


    if time_iter==490
       u20 = ut;
       save u20.mat;
    elseif time_iter==980
        u40 = ut;
         save u40.mat;
    elseif time_iter==1470
        u60 = ut;
        save u60.mat;
    elseif time_iter==1960
        u80 = ut;
        save u80.mat;
    elseif time_iter==2447
        u100 = ut;
         save u100.mat;
    end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% With Momentum %%%%%%%%%%%
if (MI == 1)

     P = [P1' ; P];

    for i=1:nx
        for j = 2:ny-1
            u(i,j) = (u(i,j)+(delt*Hx(i,j)))/delx - (P(i+2,j)-P(i,j))/2 ;
            v(i,j) = (v(i,j)+(delt*Hy(i,j)))/dely ;
        end
    end

u(nx+1,:)=u(2,:);

%     while (norm > tol)
for k = 1:1000
Piter = Piter+1;
Pold = P;
for i = 3:nx
    for j = 2:ny-1

P(i,j)= ((u(i-2,j)-u(i,j))*2/delx - P(i+2,j)+2*P(i+1,j)+4*P(i-1,j)+P(i-2,j))/6 ;

  err = err + abs(Pold(i,j)-P(i,j));

  Pold(i,j)=P(i,j);
    end
end
norm = err/((nx-2)*(ny-2));

% P(1:nx,1)=P(1:nx,2);
% P(1:nx,ny)=P(1:nx,ny-1);
disp(['PPE iter ', int2str(Piter)]);
    end


    err_u = 0;
err_v = 0;
u_new = zeros(nx,ny);
v_new = zeros(nx,ny);

for i = 2:nx
    for j = 2:ny-1

        u_new(i,j)=u(i,j)+delt*(Hx(i,j)-(P(i+1,j)-P(i-1,j))/2*delx);

        v_new(i,j)=v(i,j)+delt*(Hy(i,j)-(P(i,j+1)-P(i,j-1))/2*dely);

        err_u = err_u + abs(u(i,j)-u_new(i,j));
        err_v = err_v + abs(v(i,j)-v_new(i,j));

    end
end

norm_u = err_u/(nx-2)*(ny-2);
norm_v = err_v/(nx-2)*(ny-2);


u(1:nx,1:ny) = u_new(1:nx,1:ny);
u(1,:)= u(nx,:);
v(1:nx,1:ny) = v_new(1:nx,1:ny);
v(1,:)= v(nx,:);
disp(['Time iter ', int2str(time_iter)]);

end

end


Pt = P';
save E.mat;
save div.mat;
save ut.mat;
save Pt.mat;
