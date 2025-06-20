% Lets see the anomaly for part1 of this question, that tells us to evaluate 1/(y1+y2) 
% wrt y1 from limits 0 to 15.

%% FIRSTLY WE WILL PLOT THE FUNCTION
% The given expression is worked out in the written solution.
y1 = linspace(0,15,100);
f1 = @(t) (t+9)./(t.^2+19*t-10);
figure;
plot(y1,f1(y1));
grid on;

% As observed from the plot, we can see, that the graph abruptly changes near 0.6,
% (as it is a POLE) and as it reaches 15, it tends to zero, without touching the axis. 

%Lets try to evaluate the sum
sum = integral(f1,0,15); 
%It proves that the integral do not exist due to a very high bound-on
%error(~2.1e+00) as it is NOT BOUNDED!!
