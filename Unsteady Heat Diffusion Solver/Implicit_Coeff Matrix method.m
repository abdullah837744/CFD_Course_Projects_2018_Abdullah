%%%% Implicit Scheme by Abdullah %%%%%%%%
% Let LxL be the size of the mesh
% n is the number of nodes in x or y direction to be computed
% N is the length of coefficient matrix or number of equation in system
clear
L = 1; %length of the system in each direction
m = 20 ; % So we have a mxm sized mesh
n = m-1; % Nodes in each direction to compute except Boundary
N = n*n; % Size of the coefficient matrix
gamma = 0.15;
delx = L/m ;
dely = delx ;
delt = gamma*(delx^2) ;
Endtime = 0.1;   % time level upto which go iteration
a = gamma ; %a = gamma = (delt/delx^2);
b = gamma ; %b = gamma = (delt/dely^2);
c = 1+2*a+2*b ; %c = 1+2*gamma+2*gamma;

%Create the Co-efficient Matrix

A = zeros (N);
for i = 1:N
for j = i
A(i,j)= c ;
end
end

for i = n+1:N
    for j = i-n
    A(i,j)= -b;
    end
end

for i = 1:N-n
for j = n+i
A(i,j)= -b;
end
end

for i = 2:N
    for j = i-1
    A(i,j) = -a;
    end
end
for i = n+1:n:N
    for j = n:n:N-1
    A(i,j) = 0;
    end
end


for i = 1:N-1
    for j = i+1
    A(i,j) = -a;
    end
end
for i = n:n:N-1
    for j = n+1:n:N
    A(i,j) = 0;
    end
end



% Vector coefficient which are known from Boundary and Initial condition
% At t=0, T(x,y)=0 for 0<x<1 & 0<y<1 ; Otherwise T(x,y)=x+y

aT = zeros(N,1);
bT = zeros(N,1);

Tfa1 = zeros(n,1);
Tfa2 = zeros(n,1);

for i = 1:n

       Tfa1(i)= Tfa1(i)+ delx*i;

end
for i=1:n
        Tfa2(i)= Tfa2(i)+delx*i+L;
end

for i=1:n
      for j=1
          aT(n*(i-1)+1,j)= a*Tfa1(i,j);
      end
end
for i=1:n
    for j=1
        aT(n*i,j)= a*Tfa2(i,j);
    end
end

 Tfb1 = Tfa1;
 Tfb2 = Tfa2;

 for i = 1:n
     for j=1
         bT(i,j)= b*Tfb1(i,j);
     end
 end
 for i = 1:n
     for j=1
         bT(n*(n-1)+i,j)= b*Tfb2(i,j);
     end
 end


% Create the matrix B

B = zeros(N,1);



X = zeros (N,1);


% %Iteration for time progress

Itnumber = round(Endtime/delt);
Time_iter = 0;


for i = 1:Itnumber


for k = 1:N

       B(k) = X(k)  + aT(k) + bT(k);

end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Gauss-Seidel
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
X = zeros (N,1); %initial guess
% alpha = 1.5;
tol = 1.0e-4 ;
normVal = 100 ;
while normVal>tol
    X_old=X;

    for u=1:N

        sigma=0;

        for j=1:u-1
                sigma=sigma+A(u,j)*X(j);
        end

        for j=u+1:N
                sigma=sigma+A(u,j)*X_old(j);
        end

        X(u)=(1/A(u,u))*(B(u)-sigma);
    end

%      X_guess = X ;
%
%      for i=1:N
%          X(i)= alpha*X_guess(i) + (1-alpha)*X_old(i);
%      end


err = 0;

    for p = 1:N
       err = err + abs(X_old(p)-X(p));

    end

    normVal = err/N;

end

Time_iter = Time_iter +1;
Q = [' Now running ',num2str(Time_iter),' number iteration'];
disp(Q)


end

dlmwrite('Result_X.out',X);

Xmat = reshape(X,[n,n]); % Creating 2D matrix from vector



Xmat_whole = zeros(m+1,m+1);
%Creating a full mesh size matrix with boundary conditions
for i = 1:m+1
    for j = 1:m+1
        Xmat_whole(i,j) = Xmat_whole(i,j)+ (i-1)*delx + (j-1)*dely;
    end
end
% Inserting the boundary values to get full mesh data
for i = 1:n
    for j = 1:n
        Xmat_whole(i+1,j+1) = Xmat(i,j);
    end
end

dlmwrite('Xmat_whole',Xmat_whole,'delimiter','\t','precision',8);
save('Xmat_whole.mat','Xmat_whole');

x = 1:m+1;
y = 1:m+1;
[X,Y] = meshgrid(x,y);
Z = Xmat_whole;


figure
contourf(Y,X,Z,50);
