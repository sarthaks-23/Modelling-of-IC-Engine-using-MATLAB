%% Explanation of the problem given
% Given: y1 = 2*x+1, y2 = 2*x/(1+0.2x)
% 1. Integral of 1/(y1+y2) wrt y1 cannot be calculated between 0 to 15 (SUM NOT BOUNDED IN THE RANGE)
% So, for part 1 code is written as 'Anomaly_in_part_1.m' attached herewith.
% Now we define the function f2 = 1/(y1-y2)
% It has been given that the integration of f2 has to be done wrt y1.
% Therefore, expressing f2 as a function of y1, by manipulation.
% The substitutions done to get final expression:
% 1. x = (y1-1)/2
% 2. Substituting the value of x in y2, to get y2 in terms of y1.
% 3. Thereby, getting f in terms of y1.
% Then we used the integral (@func, lower_limit, upper_limit) to calculate the result.
y1 = linspace(0,15,100);
f2 = @(y1) (y1+9)./(y1.^2-y1+10);
sum = integral(f2,0,15);
 %% NOTE: ACTUAL IMPLEMENTATION OF F2 is in WRITTEN SOLUTION provided with the code.