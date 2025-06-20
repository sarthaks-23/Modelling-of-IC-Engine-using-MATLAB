function eq = dydx(s)
eq(1) = 2*s(1)/(1+s(1)) - s(2);
eq(2) = 1+ 2*s(1)/(1+s(1))^2 - s(2);
end
