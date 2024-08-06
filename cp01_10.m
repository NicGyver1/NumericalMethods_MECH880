function cp01_10 % robust quadratic equation solver
a = [ 6 6*10^154 0 1	1	10^(-155)];
b = [ 5 5*10^154 1 -10^5 -4 	-10^155];
c  = [-4 -4*10^154 1 1	3.999999 10^155];
root1  = zeros(size(a));
root2  = zeros(size(a));
for i=1:length(a)    
    [root1(i),root2(i)] = QuadSolve(a(i),b(i),c(i));
    enddisp('	a	b	c 	Root1	Root2')
    fprintf('%10.7g	%10.7g	%10.7g	%10.6g	%10.7g\n',...    
    [a; b; c; root1; root2]);
end



function [root1, root2] = QuadSolve(unscaled_a,unscaled_b,unscaled_c)% Returns Nan for roots that don’t exist or can’t be computed.
root1  = NaN;
root2  = NaN; % initialize to NaN
factor = max(abs([unscaled_a unscaled_b unscaled_c]));
a = unscaled_a/factor; 
b = unscaled_b/factor; 
c = unscaled_c/factor;
disc  = b*b-4*a*c;	% compute discriminant

% handle degenerate cases first
if abs(a) < eps	% linear case (a is zero)    
    if abs(b) >= eps, 
        root1  = -c/b;
    end
elseif disc <0	% imaginary roots    
    fprintf('imaginary roots for %13.6g,%13.6g,%13.6g\n',...        
    unscaled_a,unscaled_b,unscaled_c);
elseif sqrt(disc) < eps	% repeated root    
    root1  = -b/(2*a); 
    root2  = root1;
else
    temp1 = sqrt(disc);	    
    % choose formula that avoids cancellation    
    if(b>0), 
        root1  = (2*c)/(-b-temp1); 
        root2  = (-b-temp1)/(2*a);    
    else
        root1  = (-b+temp1)/(2*a); 
        root2  = (2*c)/(-b+temp1); 
    end
end
end
end