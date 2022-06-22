% Laboratory 5
% Pacuraru Fabian-Virgil

load lab5_3.mat
clc
detrendedInput = detrend(id.u);
detrendedOutput = detrend(id.y);
fprintf('The mean of the original input is: %f\n',mean(id.u));
fprintf('The mean of the detrended input is: %e\n',mean(detrendedInput));
fprintf('The mean of the original output is: %f\n',mean(id.y));
fprintf('The mean of the detrended output is: %e\n',mean(detrendedOutput));
figure('Name','Identification Data');
subplot(2,2,1)
plot(tid,id.u)
title('Original Input')
subplot(2,2,2)
plot(tid,detrendedInput)
title('Detrended Input')
subplot(2,2,3)
plot(tid,id.y)
title('Original Output')
subplot(2,2,4)
plot(tid,detrendedOutput)
title('Detrended Output')
ryu = zeros(1,length(detrendedInput));
ru = zeros(1,length(detrendedInput));
for tau = 0:(length(detrendedOutput)-1)
    u_calc1 = detrendedInput(1:(end-tau));
    u_calc2 = detrendedInput((tau+1):end);
    u_calc2 = u_calc2.';
    y_calc = detrendedOutput((tau+1):end);
    y_calc = y_calc.';
    ryu(tau+1) = sum(y_calc*u_calc1)/length(u_calc1);
    ru(tau+1) = sum(u_calc2*u_calc1)/length(u_calc1);
end
T = 1000;
Ryu = ryu(1:T).';
MSE = zeros(1,length(100:100:T));
for M = 100:100:T
    H = calculateH(T,M,Ryu,ru);
    val_y_hat = approxOutput(val.u,H,M);
    MSE(M/100) = sum((val_y_hat-val.y).^2)/length(val.y);
end
[~,bestM] = min(MSE);
bestM = bestM*100;
H = calculateH(T,bestM,Ryu,ru);
id_y_hat = approxOutput(id.u,H,bestM);
val_y_hat = approxOutput(val.u,H,bestM);
figure
subplot(2,1,1)
hold on
plot(tid,id.y)
plot(tid,id_y_hat)
hold off
title('Identification data')
legend('Real signal','Approximated signal')
subplot(2,1,2)
hold on
plot(tval,val.y)
plot(tval,val_y_hat)
hold off
title('Validation data')
legend('Real signal','Approximated signal')
figure
plot(100:100:T,MSE)
xlabel('M')
ylabel('MSE')


function H = calculateH(T,M,Ryu,ru)
    Ru = zeros(T,M);
    for i = 1:T
        for j = 1:M
            Ru(i,j) = ru(abs(j-i)+1);
        end
    end
    H = Ru\Ryu;
    H = H.';
end


function y_hat = approxOutput (u,h,M)
    y_hat = zeros(length(u),1);
    for i = 1:length(u)
        y_hat(i) = sum(h(1:min(i,M))*u(i:-1:max(1,i-M+1)));
    end
end