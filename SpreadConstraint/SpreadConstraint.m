function SpreadConstraint()
% Note: If you wrote: 
% function y = SpreadConstraint()
% It will always print the value of y at the end of this 
% function execution  with or without the ";" 
% Therefore, do not write it so won't mess up commandline with value of 
% y 

% Plot the dots by itself 
% Define elements
x = [-5 -4 -3 0 3]; % 5 points by itself 
% Submit element-by-element to the function
for ix = 1 : length(x)
    y(ix) = piecewise1(x(ix));
end
% Plot discrete values
% plot(x, y, 'ro');
% Define x and y ranges to display
axis([-8 4 -2 10]);


% Plot lines on the points
x = -8 : .01 : 4; % A lot of points continuous from -8 to 4 
for ix = 1 : length(x)
    y(ix) = piecewise1(x(ix));
end
%figure; 
% plot(x, y);
axis([-8 4 -2 9]);

% Working with symbolic variables 
% symbolicVariables(); 
% Uncomment above to check for solving for q0 from DMax = 0 calculation

% Generation of convex graph of OPT(q) 
axis([4 20 0 20]);
x = 4 : .01 : 20; % A lot of points continuous from 4 to 20 
for ix = 1 : length(x)
    y(ix) = piecewiseOpt(x(ix));
end 
%figure; 
 plot(x, y);
 
 
% Generation of concave graph of dmax(q) 
axis([4 20 0 20]);
x = 4 : .01 : 20; % A lot of points continuous from 4 to 20 
for ix = 1 : length(x)
    y(ix) = piecewiseDMax(x(ix), 1, 50)
end 
figure; 
plot(x, y); 
 
end


function y = piecewiseDMax(q,xmin,optMax)
v = 0; 
m = 0; 
opt = 0; 
a = 1; 
b = -2; 
c = 0; 
if 6 <= q & q <= 18
    if  6 <= q & q <= 7
        v = (q-5)/1; 
        m = 1;
    elseif 7 <= q & q <= 9
        v = (q-3)/2; 
        m = 2; 
    elseif 9 <= q  & q <= 15
        v = (q-3)/2; 
        m = 2; 
    elseif 15 <= q & q <= 18
        v = (q-9)/1; 
        m = 1; 
    else
        v = (q-5)/1; 
        m = 1;
    end
    n = 3; 
    opt = piecewiseOpt(q); 
    a = 1+ (1/m); 
    b = xmin - v; 
    c = opt - n*((optMax)^2); 
    y = (-b + (((b^2) - a*c)^(1/2)))/a;
else 
    y = 0; 
end 
end 


function y = piecewiseOpt(q)
if  6 <= q & q <= 7
    y = 13 + ((q-5)^2)/1 - (q^2)/3;
elseif 7 <= q & q <= 9
    y = 9 + ((q-3)^2)/2 - (q^2)/3;
elseif 9 <= q  & q <= 15
    y = 9 + ((q-3)^2)/2 - (q^2)/3;
elseif 15 <= q & q <= 18
    y = 45 + ((q-9)^2)/1 - (q^2)/3;
else
    y = 0; 

end
end 

function y = piecewise1(x)
if x <= -4
    y = 3;
elseif -4 < x & x <= -3
    y = -4*x - 13;
elseif -3 < x & x <= 0
    y = x^2 + 6*x + 8;
else
    y = 8; 
    
end 
end %Note: Double click this end to see which function it closes 


function y = symbolicVariables()

syms a b c e x
  solve(a*x^2 + b*x + c == 0, a) % To see output on commandline, don't type ";"
  solve(a*x^2 + b*x + c == 0, b) % To see output on commandline, don't type ";"
  solve((2*e*((a*x^2 + b*x + c)^(1/2)) + 2*a*x + b) == 0,x)


syms d m Xjmin OptMax q n ES C

solve ((1.0+1.0/m)*d^2 + (2.0*(Xjmin - ((q-ES)/m)))*d + (C + (((q-ES)^2)/m) - ((q^2)/n) - OptMax^2)  ,d)

% Above line outputs this as the positive solution 
% -(ES - q + Xjmin*m - m*((OptMax^2*n - ES^2*n - C*n + m*q^2 - n*q^2 + q^2 + 2*ES*Xjmin*n - C*m*n + 2*ES*n*q - 2*Xjmin*n*q + OptMax^2*m*n + Xjmin^2*m*n)/(m*n))^(1/2))/(m + 1)

% Have to retype it so that its a*b^2 NOT b^2*a
 -(ES - q + Xjmin*m - m*(((n*OptMax^2 - n*ES^2 - C*n + m*q^2 - n*q^2 + q^2 + 2*ES*Xjmin*n - C*m*n + 2*n*q*ES- 2*n*q*Xjmin + m*n*OptMax^2 + m*n*Xjmin^2)/(m*n))^(1/2)))/(m + 1)


diff(-(ES - q + Xjmin*m - m*(((n*OptMax^2 - n*ES^2 - C*n + m*q^2 - n*q^2 + q^2 + 2*ES*Xjmin*n - C*m*n + 2*n*q*ES- 2*n*q*Xjmin + m*n*OptMax^2 + m*n*Xjmin^2)/(m*n))^(1/2)))/(m + 1), q)

%Differentiate the line above with respect to q 
% Which gives you 
((2*q + 2*ES*n - 2*Xjmin*n + 2*m*q - 2*n*q)/(2*n*((OptMax^2*n - ES^2*n - C*n + m*q^2 - n*q^2 + q^2 + 2*ES*Xjmin*n - C*m*n + 2*ES*n*q - 2*Xjmin*n*q + OptMax^2*m*n + Xjmin^2*m*n)/(m*n))^(1/2)) + 1)/(m + 1)

% Solve for q with above equation == 0 
solve ((((2*q + 2*ES*n - 2*Xjmin*n + 2*m*q - 2*n*q)/(2*n*((OptMax^2*n - ES^2*n - C*n + m*q^2 - n*q^2 + q^2 + 2*ES*Xjmin*n - C*m*n + 2*ES*n*q - 2*Xjmin*n*q + OptMax^2*m*n + Xjmin^2*m*n)/(m*n))^(1/2)) + 1)/(m + 1)) == 0, q)

% The 2 output solutions are: 
% Solution 1: 
% (ES*n^2 - Xjmin*n^2 + n*(-(m - n)*(C + C*m - C*n - OptMax^2*m + OptMax^2*n - Xjmin^2*m + Xjmin^2*n + ES^2 - OptMax^2 - 2*ES*Xjmin))^(1/2) - ES*m*n + Xjmin*m*n)/(m^2 - 2*m*n + m + n^2 - n)

% Solution 2: 
% -(Xjmin*n^2 - ES*n^2 + n*(-(m - n)*(C + C*m - C*n - OptMax^2*m + OptMax^2*n - Xjmin^2*m + Xjmin^2*n + ES^2 - OptMax^2 - 2*ES*Xjmin))^(1/2) + ES*m*n - Xjmin*m*n)/(m^2 - 2*m*n + m + n^2 - n)

end 


