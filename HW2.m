% Author: Nicholas Piercy
% Date : 9/14/2021
% Numerical Methods, Homework #2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%% Problem 1
U1a = [3 3 1;
     0 -4 -3;
     0 0 2];
b1a = [12; -10; 4];
[x1a,Check1a] = BackwardSub(U1a,b1a); 

U1b = [3 0 0;
      2 -1 0;
      5 1 -2];
b1b = [15; 10; 5];
[x1b,Check1b] = ForwardSub(U1b,b1b);

T = table(x1a(:),Check1a(:),x1b(:),Check1b(:),'VariableNames',{'x1a','Check1a','x1b','Check1b'});
disp(T)


%% Problem 2
format long

A2 = [10^(-13) 1;
         1     1];
b2 = [10^(-13)+1; 2];
[U2a,b2a_new] = GE(A2,b2);
[x2a] = BackwardSub(U2a,b2a_new);

[U2b,b2b_new] = GEPP(A2,b2);
[x2b] = BackwardSub(U2b,b2b_new);

T = table(x2a(:,:),x2b(:,:),'VariableNames',{'Without Pivoting', 'With Partial Pivoting'});
disp(T)

%% Problem 3
A3 = [-2 1 -1;
    -3 4 -6;
    2 7 15];

[A3_Inverse,check] = GaussJordanInverseV3(A3);

T1 = table(A3_Inverse(:,:),'VariableNames',{'A3_Inverse'});
T2 = table(check(:,:),'VariableNames',{'MATLAB "inv"'});
disp(T1)
disp(T2)

%% Problem 4
A4 = [8 1 0;
     1 25 3;
     0 3 9];
b4 = [7; -18; 15];

[H4,H4T] = Cholesky(A4); %H4 == Lower triang, H4T == Upper triang

[y4] = ForwardSub(H4,b4);
[x4] = BackwardSub(H4T,y4);

 Check4 = A4*x4;
 
 T1 = table(H4(:,:),'VariableNames',{'H'});
 T2 = table(H4T(:,:),'VariableNames',{'Htranspose'});
 T3 = table(x4,b4,Check4,'VariableNames',{'x (solution)','b (RHS)','check'});
 disp(T1)
 disp(T2)
 disp(T3)
 
 %% Problem 5
n = 10;
A5 = zeros(n,n);
b5 = zeros(n,1);
for i = 1:n
    A5(i,i) = 3*i;
    if i<n
        A5(i,i+1) = -(i+1);
        A5(i+1,i) = -i;
    end
    b5(i) = i;
end
[Thom_U,Thom_L] = ThomasAlgorithm(A5);

[y5] = ForwardSub(Thom_L,b5);
[x5] = BackwardSub(Thom_U,y5);

check5 = A5*x5;

T = table(x5(:),b5(:),check5(:),'VariableNames',{'x (solution)','b (RHS)','check'});
disp(T)

%% Problem 6
 A6 = [1 -1 0;
     -1 2 -1;
     0 -1 1];
 [L,U,P] = lu(A6);
 
 T = table(L(:,:),U(:,:),P(:,:),'VariableNames',{'L','U','Permutations'});
 disp(T)
 
 %% Problem 7
alpha = sqrt(2)/2;
A7 = zeros(14,14);
% Linear system given in computer problem 2.3
A7(2,2) = 1; A7(2,6) = -1;
A7(3,3) = 1;
A7(4,1) = alpha; A7(4,4) = -1; A7(4,5) = -alpha;
A7(5,1) = alpha; A7(5,3) = 1; A7(5,5) = alpha;
A7(6,4) = 1; A7(6,8) = -1;
A7(7,7) = 1;
A7(8,5) = alpha; A7(8,6) = 1; A7(8,9) = -alpha; A7(8,10) = -1; 
A7(9,5) = alpha; A7(9,7) = 1; A7(9,9) = alpha;
A7(10,10) = 1; A7(10,13) = -1;
A7(11,11) = 1;
A7(12,8) = 1; A7(12,9) = alpha; A7(12,12) = -alpha ;
A7(13,9) = alpha; A7(13,11) = 1; A7(13,12) = alpha;
% A7(14,12) = 1; A7(14,13) = alpha; % Incorrect ordering
A7(14,12) = alpha; A7(14,13) = 1; %Corrected ordering

RHS = [0 10 0 0 0 0 0 15 0 20 0 0 0]';
A7 = A7(2:end,:);

f = A7\RHS;

T = table(f(:),'VariableNames',{'forces, f, (solution)'});
disp(T)
 



