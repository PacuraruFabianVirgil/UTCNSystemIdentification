% Laboratory 3
% Pacuraru Fabian-Virgil

%% first-order function

load lab3_order1_7.mat
plot(data)
figure
plot(data(31:130))
clc
y0 = sum(data.y(120:130))/11; yss = y0;
u0 = data.u(1); uss = u0;
K = yss/uss;
[ymax,i1]=max(data.y(31:130));
i1 = i1+30;
t1 = t(i1);
fprintf('For the first order transfer function of index 7\n\n');
fprintf('The gain (K) is: %f\n', K);
valueOfyOfT = y0 + 0.368*(ymax-y0);
for i=i1:130
    if data.y(i)<=valueOfyOfT
        indexForTime = i;
        break
    end
end
t2 = t(indexForTime);
T = t2-t1;
fprintf('The time constant (T) is: %f\n', T);
H = tf(K,[T 1]);
fprintf('The transfer function (H) is:');
H
A = -1/T;
B = K/T;
C = 1;
D = 0;
Hss = ss(A, B, C, D);
Hss
yhat = lsim(Hss,data.u,t);
validationY = data.y(131:end);
validationYhat = yhat(131:end);
validationTime = t(131:end);
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

load lab3_order2_7.mat
plot(data)
figure
plot(data(31:130))
clc
y0 = sum(data.y(120:130))/11; yss = y0;
u0 = data.u(1); uss = u0;
K = yss/uss;
fprintf('For the second order transfer function of index 7\n\n');
fprintf('The gain (K) is: %f\n', K);

[~,indexOfFirstMax] = max(data.y(31:130));
indexOfFirstMax = indexOfFirstMax + 30;
[~,indexOfFirstMin] = min(data.y(indexOfFirstMax:130));
indexOfFirstMin = indexOfFirstMax + indexOfFirstMin - 1;
[~,indexOfSecondMax] = max(data.y(indexOfFirstMin:130));
indexOfSecondMax = indexOfFirstMin + indexOfSecondMax - 1;
indexesFor0Values = zeros(1,3);
for i=indexOfFirstMax:-1:1
    if data.y(i)<=yss
        indexesFor0Values(1) = i;
        break
    end
end
for i=indexOfFirstMax:indexOfFirstMin
    if data.y(i)>=yss && data.y(i+1)<yss
        indexesFor0Values(2) = i;
        break
    end
end
for i=indexOfFirstMin:indexOfSecondMax
    if data.y(i)<=yss && data.y(i+1)>yss
        indexesFor0Values(3) = i;
        break
    end
end
Ts = t(2)-t(1);
fprintf('The sampling time (Ts) is: %f\n', Ts);
Aplus = Ts*sum(data.y(indexesFor0Values(1):indexesFor0Values(2))-y0);
fprintf('The positive area (Aplus) is: %f\n', Aplus);
Aminus = Ts*sum(y0-data.y(indexesFor0Values(2):indexesFor0Values(3)));
fprintf('The negative area (Aminus) is: %f\n', Aminus);
M = Aminus/Aplus;
fprintf('The overshoot (M) is: %f\n', M);
E = log(1/M)/sqrt(pi^2+log(M)^2);
fprintf('The damping factor (E) is: %f\n', E);
To = 2*(t(indexOfFirstMin)-t(indexOfFirstMax));
fprintf('The oscillation period (To) is: %f\n', To);
Wn = (2*pi)/(To*sqrt(1-E^2));
fprintf('The natural frequency (Wn) is: %f\n', Wn);
H = tf(K*Wn^2,[1 2*E*Wn Wn^2]);
fprintf('The transfer function (H) is:');
H
A = [0 1; -Wn^2 -2*E*Wn];
B = [0; K*Wn^2];
C = [1 0];
D = 0;
Hss = ss(A, B, C, D);
Hss
yhat = lsim(Hss,data.u,t);
validationY = data.y(131:end);
validationYhat = yhat(131:end);
validationTime = t(131:end);
MSE = sum((validationYhat-validationY).^2)/length(validationTime);
figure
hold on
plot(validationTime,validationYhat)
plot(validationTime,validationY)
legend('Yhat (Model)','Y (System)')
error = sprintf('MSE = %f',MSE);
title(error)
hold off