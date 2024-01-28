%%% Explicit Scheme by Abdullah %%%

m = 100 ; % mesh size m x m
L = 1; % Length in each direction
n = m-1; % Grids without boubdary condition
N = n*n;
gamma = 0.15;
delx = L/m ;
dely = delx ;
delt = gamma*(delx^2) ;
Endtime = 0.1;   % time level upto which go iteration
a = gamma ; %a = gamma = (delt/delx^2);
b = gamma ; %b = gamma = (delt/dely^2);
c = 1-2*a-2*b ; %c = 1-4*gamma;


% Matrix with boundary values
T = zeros (n+2,n+2);

for i = 1:n+2
    for j = 1:n+2
        T (i,j) = (i+j-2)*delx;
    end
end
for j = 1:n+2
    for i = 1:n+2
        T(i,j) = (i+j-2)*delx;
    end
end
for i=2:n+1
    for j = 2:n+1
         T(i,j)=0;
     end
end


Itnumber = round((Endtime-delt)/delt);

% Compute

for k = 1:Itnumber

    for i = 1:m+1
        for j = 1:m+1
            T(i,j)=T(i,j);
        end
    end

for i = 1+1:n+1
    for j = 1+1:n+1
        T(i,j) = a*(T(i-1,j)+T(i+1,j)+T(i,j-1)+T(i,j+1)) + c*(T(i,j));
    end
end

end

dlmwrite('T',T,'delimiter','\t','precision',8);
save('T.mat','T');

T_Exp_gamma1 = T;
save('T_Exp_gamma1.mat','T_Exp_gamma1');

T_Exp_Hundred = T;
save('T_Exp_Hundred.mat','T_Exp_Hundred');

T_Exp_t1 = T;
save('T_Exp_t1.mat','T_Exp_t1');

x = 1:m+1;
y = 1:m+1;
[X,Y] = meshgrid(x,y);
Z = T;

figure(1)
contourf(Y,X,Z,50);
colorbar
xlabel('X');
ylabel('Y');
h1 = figure(1)
saveas(h1,'Explicit_100_0.15_0.1.png');
