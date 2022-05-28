%% Exercise 1
% 1
y = 30:-1:1
% 2
y(1:2:30) = sin(y(1:2:30))
% 3
y(2:2:30) = sort(y(2:2:30))

%% Exercise 2
% 1
x = 0:0.25:4
% 2
y1 = 2*exp(-x.^2)+2*sin(0.67*x+0.1)
% 3
y2 = 2.2159+1.2430*x-2.6002*x.^2+1.7223*x.^3-0.4683*x.^4+0.0437*x.^5
% 4
hold
plot(x,y1)
plot(x,y2,'--')
legend('True values y','Approximate values yhat')
% 5
e = y1-y2
figure
plot(x,e)
% 6
error = 1/length(e)*sum(e.^2)
error = round(error,4)
string = sprintf('MSE = %f',error)
string = string(1:length(string)-2)
title(string)