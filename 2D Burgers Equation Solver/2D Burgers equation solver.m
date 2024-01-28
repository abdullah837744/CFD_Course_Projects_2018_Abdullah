%%% 2D Burgers equation solver by Abdullah %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

Method = 1; % 1 for Explicit and 2 for Implicit

X_length = 2;
Y_length = 1;

Re = 40 ; % Global Reynolds number
Pe = 1 ; % Local Peclet number
gamma = 0.1 ;
V = 1; % Upper wall velocity

neu = V*Y_length/Re ;
delx = Pe*neu/V ;
dely = delx ;
delt = (gamma*delx^2)/neu ;

nx = X_length/delx +1 ; % Total nodes in x direction
ny = Y_length/dely +1 ; % Total nodes in y direction

u = zeros(nx,ny);
v = zeros(nx,ny);

% Initial condition
for i = 2:nx-1
    for j = 2:ny-1
            u(i,j)=0;
            v(i,j)=0;
    end
 end

% Boundary condition
u(:,1) = 0;
u(:,ny) = 1;
v(:,1) = 0;
v(:,ny) = 0;

u(nx,2:ny-1) = 0;
v(nx,2:ny-1) = 0;

if (Method == 1)

tol = 1.0e-6;
iter = 0;
norm1 = 100;
norm2 = 100;

while (norm1 > tol || norm2 > tol)
    iter = iter+1;
    diff1 = 0;
    diff2 = 0;

    for j = 2:ny-1
            u(1,j)= u(nx,j);
            v(1,j)= v(nx,j);
            u(nx+1,j)= u(2,j);
            v(nx+1,j)= v(2,j);

    end

    u_old=u;
    v_old=v;


    Ex_advective = 1;  % 1 for central difference and 2 for First Order Upwinding

    if (Ex_advective == 1)

    for i=2:nx
        for j=2:ny-1

    u(i,j)=u_old(i,j) - u_old(i,j)*(delt/2*delx)*(u_old(i+1,j)-u_old(i-1,j))...
        - v_old(i,j)*(delt/2*dely)*(u_old(i,j+1)-u_old(i,j-1))...
        +gamma*(u_old(i+1,j)-2*u_old(i,j)+u_old(i-1,j))...
        +gamma*(u_old(i,j+1)-2*u_old(i,j)+u_old(i,j-1)) ;

    v(i,j)=v_old(i,j) - u_old(i,j)*(delt/2*delx)*(v_old(i+1,j)-v_old(i-1,j))...
        - v_old(i,j)*(delt/2*dely)*(v_old(i,j+1)-v_old(i,j-1))...
        +gamma*(v_old(i+1,j)-2*v_old(i,j)+v_old(i-1,j))...
        +gamma*(v_old(i,j+1)-2*v_old(i,j)+v_old(i,j-1)) ;

    diff1 = diff1 + abs(u(i,j) - u_old(i,j));
    diff2 = diff2 + abs(v(i,j) - v_old(i,j));

        end
    end
    end

    if (Ex_advective == 2)

    for i=2:nx
        for j=2:ny-1

    if (u(i,j)>= 0)  % For flow going left to right

    u(i,j)=u_old(i,j) - u_old(i,j)*(delt/delx)*(u_old(i,j)-u_old(i-1,j))...
        - v_old(i,j)*(delt/dely)*(u_old(i,j)-u_old(i,j-1))...
        +gamma*(u_old(i+1,j)-2*u_old(i,j)+u_old(i-1,j))...
        +gamma*(u_old(i,j+1)-2*u_old(i,j)+u_old(i,j-1)) ;

    v(i,j)=v_old(i,j) - u_old(i,j)*(delt/delx)*(v_old(i,j)-v_old(i-1,j))...
        - v_old(i,j)*(delt/dely)*(v_old(i,j)-v_old(i,j-1))...
        +gamma*(v_old(i+1,j)-2*v_old(i,j)+v_old(i-1,j))...
        +gamma*(v_old(i,j+1)-2*v_old(i,j)+v_old(i,j-1)) ;

    diff1 = diff1 + abs(u(i,j) - u_old(i,j));
    diff2 = diff2 + abs(v(i,j) - v_old(i,j));

    elseif (u(i,j)< 0)    % For flow going right to left

    u(i,j)=u_old(i,j) - u_old(i,j)*(delt/delx)*(u_old(i+1,j)-u_old(i,j))...
        - v_old(i,j)*(delt/dely)*(u_old(i,j+1)-u_old(i,j))...
        +gamma*(u_old(i+1,j)-2*u_old(i,j)+u_old(i-1,j))...
        +gamma*(u_old(i,j+1)-2*u_old(i,j)+u_old(i,j-1)) ;

    v(i,j)=v_old(i,j) - u_old(i,j)*(delt/delx)*(v_old(i+1,j)-v_old(i,j))...
        - v_old(i,j)*(delt/dely)*(v_old(i,j+1)-v_old(i,j))...
        +gamma*(v_old(i+1,j)-2*v_old(i,j)+v_old(i-1,j))...
        +gamma*(v_old(i,j+1)-2*v_old(i,j)+v_old(i,j-1)) ;

    diff1 = diff1 + abs(u(i,j) - u_old(i,j));
    diff2 = diff2 + abs(v(i,j) - v_old(i,j));
    end

        end
    end
    end


    norm1 = diff1/((nx-1)*(ny-2)) ;
    norm2 = diff2/((nx-1)*(ny-2)) ;

    Q = [' Now running ',num2str(iter),' number iteration'];
    disp(Q)

end

u_final = zeros(ny,nx+1);
for i=1:nx+1
    for j = 1:ny
        u_final(j,i)=u(i,j);
    end
end


u_final = u_final(:,1:nx);

contourf(u_final,25);
colormap jet

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (Method==2)

    c = gamma ;
    a = delt/(2*delx);
    b = delt/delx;

% Time Iteration
tol_t = 1.0e-6;
iter = 0;
norm_t_1 = 100;
norm_t_2 = 100;

while (norm_t_1 > tol_t || norm_t_2 > tol_t)
    iter = iter+1;
    diff_t_1 = 0;
    diff_t_2 = 0;

    for j = 2:ny-1
            u(1,j)= u(nx,j);
            v(1,j)= v(nx,j);
            u(nx+1,j)= u(2,j);
            v(nx+1,j)= v(2,j);

    end

    u_old_t = u;
    v_old_t = v;

    u_old=u;
    v_old=v;

    Im_advective = 2; % 1 for central difference and 2 for First Order Upwinding

    % Gauss-Seidel Iteration

tol_g = 1.0e-4;
norm_g_1 = 100;
norm_g_2 = 100;

    while (norm_g_1 > tol_t || norm_g_2 > tol_g)
    diff_g_1 = 0;
    diff_g_2 = 0;

    if (Im_advective == 1)

    for i = 2:nx
        for j = 2:ny-1
   u(i,j) = (u_old(i,j)+(a*v(i,j)+c)*u(i,j-1)+(a*u(i,j)+c)*u(i-1,j)...
       -(a*u(i,j)-c)*u(i+1,j)-(a*v(i,j)-c)*u(i,j+1))/(1+4*c);

    v(i,j) = (v_old(i,j)+(a*v(i,j)+c)*v(i,j-1)+(a*u(i,j)+c)*v(i-1,j)...
        -(a*u(i,j)-c)*v(i+1,j)-(a*v(i,j)-c)*v(i,j+1))/(1+4*c);


    diff_g_1 = diff_g_1 + abs(u(i,j)-u_old(i,j));
    diff_g_2 = diff_g_2 + abs(v(i,j)-v_old(i,j));

    u_old(i,j) = u(i,j);
    v_old(i,j) = v(i,j);

        end
    end
    end

    if (Im_advective == 2)

    for i = 2:nx
        for j = 2:ny-1

            if (u(i,j)>= 0)  % For flow going left to right

   u(i,j) = (u_old(i,j)+ c*(u(i+1,j)+u(i,j+1))+ (b*u(i,j)+c)*u(i-1,j)...
       +(b*v(i,j)+c)*u(i,j-1))/(1+4*c+b*u(i,j)+b*v(i,j));

    v(i,j) = (v_old(i,j)+ c*(v(i+1,j)+v(i,j+1))+ (b*u(i,j)+c)*v(i-1,j)...
       +(b*v(i,j)+c)*v(i,j-1))/(1+4*c+b*u(i,j)+b*v(i,j));


    diff_g_1 = diff_g_1 + abs(u(i,j)-u_old(i,j));
    diff_g_2 = diff_g_2 + abs(v(i,j)-v_old(i,j));

    u_old(i,j) = u(i,j);
    v_old(i,j) = v(i,j);

            elseif (u(i,j)< 0)  % For flow going right to left

    u(i,j) = (u_old(i,j)+ (c-b*u(i,j))*u(i+1,j)+ (c-b*v(i,j))*u(i,j+1)...
        + c*(u(i-1,j)+u(i,j-1)))/(1+4*c-b*u(i,j)-b*v(i,j));

    v(i,j) = (v_old(i,j)+ (c-b*u(i,j))*v(i+1,j)+ (c-b*v(i,j))*v(i,j+1)...
        + c*(v(i-1,j)+v(i,j-1)))/(1+4*c-b*u(i,j)-b*v(i,j));


            end
        end
    end
    end

    norm_g_1 = diff_g_1 / ((nx-1)*(ny-2));
    norm_g_2 = diff_g_2 / ((nx-1)*(ny-2));



    end

    for i = 2:nx
        for j = 2:ny-1
  diff_t_1 = diff_t_1 + abs(u(i,j)-u_old_t(i,j));
  diff_t_2 = diff_t_2 + abs(v(i,j)-v_old_t(i,j));
        end
    end

    norm_t_1 = diff_t_1 / ((nx-1)*(ny-2));
    norm_t_2 = diff_t_2 / ((nx-1)*(ny-2));

    Q = [' Now running ',num2str(iter),' number iteration'];
    disp(Q)
end

u_final = zeros(ny,nx+1);
for i=1:nx+1
    for j = 1:ny
        u_final(j,i)=u(i,j);
    end
end


u_final = u_final(:,1:nx);

contourf(u_final,25);
colormap jet

end
