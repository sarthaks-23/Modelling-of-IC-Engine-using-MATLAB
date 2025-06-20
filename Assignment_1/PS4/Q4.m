s0 = [0,0]; %first guess point
t0 = [-0.5,0]; %second guess point

%Calling function dydx as f in this code
f = @dydx;

%point of tangency for first guess point
pt1=fsolve(f,s0); 
xt1 = pt1(1);
yt1 = pt1(2);

%point of tangency for second guess point
pt2 = fsolve(f,t0);
xt2 = pt2(1);
yt2 = pt2(2);
%SLOPE OF TANGENTS for both points
m1 = (yt1-1)/xt1;
m2 = (yt2-1)/xt2;
%% Equation of curve and tangent for both points
curve = @(x) (2*x)./(1+x);
% Equation of tangents: y = mx + c, as it passes through A(0,1), c = 1.
c = 1;
y1 = @(x) m1.*x + c; %tangent 1 passing through (xt1,yt1)
y2 = @(x) m2.*x + c; %tangent 2 passing through (xt2,yt2)

x = linspace(-1,5,1000);

figure;
hold on;
plot(x,curve(x),'r','LineWidth',1.5);% Equation of curve: y=2*x/(1+x)
plot(x,y1(x),'b--','LineWidth',1.5); % Equation of tangent to first point of tangency
plot(x,y2(x),'--','LineWidth',1.5); % Equation of tangent to second point of tangency
plot(xt1,yt1,'ko','MarkerSize',8,'MarkerFaceColor', 'k'); %First point of tangency
plot(xt2,yt2,'o','MarkerSize',8,'MarkerFaceColor','auto'); %Second point of tangency
plot(0,1,'go','MarkerSize',8,'MarkerFaceColor', 'g'); %The point A(0,1)
xlabel('x');
ylabel('y');
legend('equation of curve: y = (2*x)/(1+x)','equation of tangent: y=0.1716x+1','equation of tangent: y=5.8284x+1','1^st point of tangency (x_t_1,y_t_1)','2^nd point of tangency (x_t_2,y_t_2)','point A(0,1)','Location','best');
grid on;
hold off;





