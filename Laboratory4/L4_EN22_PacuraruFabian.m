% Laboratory 4
% Pacuraru Fabian-Virgil

load lab4_04.mat
n_max = 20;
v = 1:n_max;
MSEid = zeros(1,n_max);
MSEval = zeros(1,n_max);
for n=1:n_max
    coeffVector = computePolyCoeff(n,id);
    figure;
    subplot(2,1,1);
    str = 'Identification';
    MSEid(n) = computeMSE(n,coeffVector,id,str);
    subplot(2,1,2);
    str = 'Validation';
    MSEval(n) = computeMSE(n,coeffVector,val,str);
end
figure;
subplot(2,1,1);
plot(v,MSEid);
title('Identification set');
xlabel('polynomial degree');
ylabel('mean squared error');
subplot(2,1,2);
plot(v,MSEval);
title('Validation set');
xlabel('polynomial degree');
ylabel('mean squared error');
clc
[~,perfectPolDegree] = min(MSEval);
fprintf('The polynomial that most acurately represents the ');
fprintf('validation set output is of degree %d',perfectPolDegree);


function [d] = computePolyCoeff(n,id)
    A = zeros(length(id.X),n);
    for i=1:length(id.X)
        for j=1:n
            A(i,j) = id.X(i)^(j-1);
        end
    end
    id.Y=id.Y.';
    d = A\id.Y;
end


function [MSE] = computeMSE(n,d,set,str)
    my_y = zeros(1,length(set.Y));
    for i=1:n
        my_y = my_y+d(i).*set.X.^(i-1);
    end
    hold on
    plot(set.X,my_y);
    plot(set.X,set.Y);
    hold off
    MSE = sum((my_y-set.Y).^2)/length(set.Y);
    error = sprintf('%s set (n = %d)\nMSE = %f',str,n,MSE);
    title(error);
    legend('Aproximated output','Real output');
end