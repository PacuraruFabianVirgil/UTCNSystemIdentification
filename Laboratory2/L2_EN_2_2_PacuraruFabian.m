% Laboratory 2
% Pacuraru Fabian-Virgil

%% first-order function

load lab2_order1_4.mat
plot(data)
figure
plot(data(1:100))
clc
y0 = data.y(1);
yss = sum(data.y((length(t)/5-10):(length(t)/5)))/11;
u0 = 0;
uss = data.u(length(t)/5);
K = (yss-y0)/(uss-u0);
fprintf('For the first order transfer function of index 4\n\n');
fprintf('The gain (K) is: %f\n', K);
valueOfyOfT = y0 + 0.632*(yss-y0);
for i=1:length(t)/5
    if data.y(i)>=valueOfyOfT
        indexForTime = i;
        break
    end
end
T = t(indexForTime);
fprintf('The time constant (T) is: %f\n', T);
H = tf(K,[T 1]);
fprintf('The transfer function (H) is:');
H
yhat = lsim(H,data.u,t);
validationY = data.y(2*length(data.y)/5:end);
validationYhat = yhat(2*length(yhat)/5:end);
validationTime = t(2*length(t)/5:end);
MSE = sum((validationYhat-validationY).^2)/length(validationTime);
figure
hold on
plot(validationTime,validationYhat)
plot(validationTime,validationY)
legend('Yhat (Model)','Y (System)')
error = sprintf('MSE = %f',MSE);
title(error)
hold off

%% second-order function

load lab2_order2_4.mat
plot(data)
figure
plot(data(1:100))
clc
y0 = data.y(1);
yss = sum(data.y((length(t)/5-10):(length(t)/5)))/11;
u0 = 0;
uss = data.u(length(t)/5);
K = (yss-y0)/(uss-u0);
fprintf('For the second order transfer function of index 4\n\n');
fprintf('The gain (K) is: %f\n', K);
[~,indexOfFirstMax] = max(data.y(1:length(data.y)/5));
[~,indexOfFirstMin] = min(data.y(indexOfFirstMax:length(data.y)/5));
indexOfFirstMin = indexOfFirstMax + indexOfFirstMin - 1;
[~,indexOfSecondMax] = max(data.y(indexOfFirstMin:length(data.y)/5));
indexOfSecondMax = indexOfFirstMin + indexOfSecondMax - 1;
M = (data.y(indexOfFirstMax)-yss)/(yss-y0);
fprintf('The overshoot (M) is: %f\n', M);
E = log(1/M)/sqrt(pi^2+log(M)^2);
fprintf('The damping factor (E) is: %f\n', E);
To = t(indexOfSecondMax)-t(indexOfFirstMax);
fprintf('The oscillation period (To) is: %f\n', To);
Wn = (2*pi)/(To*sqrt(1-E^2));
fprintf('The natural frequency (Wn) is: %f\n', Wn);
H = tf(K*Wn^2,[1 2*E*Wn Wn^2]);
fprintf('The transfer function (H) is:');
H
yhat = lsim(H,data.u,t);
validationY = data.y(2*length(data.y)/5:end);
validationYhat = yhat(2*length(yhat)/5:end);
validationTime = t(2*length(t)/5:end);
MSE = sum((validationYhat-validationY).^2)/length(validationTime);
figure
hold on
plot(validationTime,validationYhat)
plot(validationTime,validationY)
legend('Yhat (Model)','Y (System)')
error = sprintf('MSE = %f',MSE);
title(error)
hold off