% Author: Nicholas Piercy
% Date : 9/30/2021
% Numerical Methods, Homework #3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%% Problem 2
clear all
clc
x2_guess = 0;
A2 = [3 2 1;
      2 3 1;
      1 2 3];
b2 = [1; 6; 10];
hard_stop = 50; % max number of iterations to perform
epsilon_stop = 10^(-6); % norm of the residual error stopping criteria
norm = 1;
[x2a,k_J] = Jacobi(A2,b2,x2_guess,hard_stop,epsilon_stop,norm);
[x2b,k_GS] = GS(A2,b2,x2_guess,hard_stop,epsilon_stop,norm);

%% Problem 3
clear all
clc
 x3_guess = 0;
 A3 = [9 5 1;
       3 0 -5
      -4 11 2];
 b3 = [8; -13; 41];
 % Performing row swapping
%  A3_swap = A3(2,:); A3(2,:) = A3(3,:); A3(3,:) = A3_swap;
%  b3_swap = b3(2); b3(2) = b3(3); b3(3) = b3_swap;
 
 hard_stop = 10000000;
 epsilon_stop = [10^(-4), 10^(-8)];
 norm = 1;
 for i = 1:2
%      [x3_J{i},k_J(i),epsilon_J(i)] = Jacobi(A3,b3,x3_guess,hard_stop,epsilon_stop(i),norm);
     [x3_GS{i},k_GS(i),epsilon_GS(i)] = GS(A3,b3,x3_guess,hard_stop,epsilon_stop(i),norm);
 end
 
 %% Problem 4
 clear all
 clc
 format long
 n = 9;
 A4 = zeros(n,n);
 b4 = (1/16)*(ones(n,1));
 x4_guess = 0;
 hard_stops = [2, 20, 50];
 epsilon_stop = 0;
 % omega > 1 ==> over relaxation, omega <1 ==> under relaxation
 omega = [1.1716, 1.67]; 
 norm = 2;
 for i = 1:n
    A4(i,i) = 4;
    if i > 1
        A4(i,i-1) = -1;
        A4(i-1,i) = -1;
    end
    if i > 3
        A4(i,i-3) = -1;
        A4(i-3,i) = -1;
    end
 end
   
 for j = 1:length(hard_stops)
     hard_stop = hard_stops(j);
     [x4_GS{j},k_GS(j),epsilon_GS(j)] = GS(A4,b4,x4_guess,hard_stop,epsilon_stop,norm);
     [x4_J{j},k_J(j),epsilon_J(j)] = Jacobi(A4,b4,x4_guess,hard_stop,epsilon_stop,norm);
     for i = 1:length(omega)
     [x4_SOR{j}(:,i),k_SOR(j,i),epsilon_SOR_GS(j,i)] = SOR_GS(A4,...
         b4,x4_guess,hard_stop,epsilon_stop,omega(i),norm);
     end
 end
 
%% problem 5
clear all
clc

n = 9;
A5 = zeros(n,n);
b5 = (1/16)*(ones(n,1));
x5_guess = 0;
for i = 1:n
   A5(i,i) = 4;
   if i > 1
       A5(i,i-1) = -1;
       A5(i-1,i) = -1;
   end
   if i > 3
       A5(i,i-3) = -1;
       A5(i-3,i) = -1;
   end
end  
maxits = [2, 20, 50];
format long
for j = 1:length(maxits)
    
    [x5_PCG{j},flag(j),relres(j),iter(j),resvec{j}] = pcg(A5,b5,[],maxits(j),[],[],[]);
end

%% Problem 5 Plotting
clear all
clc

n = 9;
A5 = zeros(n,n);
b5 = (1/16)*(ones(n,1));
x5_guess = 0;
epsilon_stop = 0;
norm = 2;
omega = [1.1716, 1.67]; 
for i = 1:n
   A5(i,i) = 4;
   if i > 1
       A5(i,i-1) = -1;
       A5(i-1,i) = -1;
   end
   if i > 3
       A5(i,i-3) = -1;
       A5(i-3,i) = -1;
   end
end  
maxits = 50;
format long

for j = 1:maxits
     
     [x5_GS,k_GS,epsilon_GS(j)] = GS(A5,b5,x5_guess,j,epsilon_stop,norm);
     [x5_J,k_J,epsilon_J(j)] = Jacobi(A5,b5,x5_guess,j,epsilon_stop,norm);
     
     [x4_SOR1,k_SOR1,epsilon_SOR_GS1(j)] = SOR_GS(A5,...
         b5,x5_guess,j,epsilon_stop,omega(1),norm);
     [x4_SOR2,k_SOR2,epsilon_SOR_GS2(j)] = SOR_GS(A5,...
         b5,x5_guess,j,epsilon_stop,omega(2),norm);
     
     [x5_PCG,flag,relres(j),iter,resvec] = pcg(A5,b5,[],j,[],[],[]);
end
    
X = (1:1:maxits);
     semilogy(X,epsilon_J(:),'linewidth',2)
     grid on
     hold on
     semilogy(X,epsilon_GS(:),'linewidth',2)
     semilogy(X,epsilon_SOR_GS1(:),'linewidth',2)
     semilogy(X,epsilon_SOR_GS2(:),'linewidth',2)
     semilogy(X,relres(:),'linewidth',2)
 xlabel('Number of Iterations','interpreter','latex','fontsize',20)
 ylabel('$||r||_2$','interpreter','latex','fontsize',20)
 legend('Jacobi','Gauss-Seidel','SOR GS $\omega = 1.1716$','SOR GS $\omega = 1.67$'...
     ,'PCG','interpreter','latex','fontsize',15,'linewidth',2)
 hold off


%% Problem 6
clear all
clc
A6 = [3.02 3.06 1.99;
     1.27 4.16 -1.23;
     0.987 -4.81 9.34];
b6 = [1; 1; 1];

A6_perturbed = A6;
A6_perturbed(1,1) = 2.99;
A6_perturbed(3,1) = 0.990;

[A6_new,b6_new] = GEPP(A6,b6);
[A6_perturbed_new,b6_perturbed_new] = GEPP(A6_perturbed,b6);

[x6,check] = BackwardSub(A6_new,b6_new);
[x6_perturbed,check] = BackwardSub(A6_perturbed_new,b6_perturbed_new);

x6_checker = A6\b6;
x6_perturbed_checker = A6_perturbed\b6;

for i = 1:length(x6) % *100 ==> percent difference
    relativeDifference(i) = abs((x6(i)- x6_perturbed(i)))*100/abs(x6(i)); 
end

relativeDiff_system = (norm((A6-A6_perturbed),2)/norm(A6,2))*100;
relativeDiff_solution = (norm((x6-x6_perturbed),2)/norm(x6,2))*100;



