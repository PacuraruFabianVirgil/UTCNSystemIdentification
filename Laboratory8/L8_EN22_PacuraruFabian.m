load lab8_4.mat
warning('off','all'); clc;
theta = calculateTheta([1; 1], 0.1, 0, id);
myModelId = idpoly(1, [0 theta(2)], 1, 1, [1 theta(1)], 0, id.Ts);
myModelVal = idpoly(1, [0 theta(2)], 1, 1, [1 theta(1)], 0, val.Ts);
figure('Name','Identification data'); compare(id, myModelId);
figure('Name','Validation data'); compare(val, myModelVal);



function theta = calculateTheta(thetaNew, alfa, l, set)
    while(l<=200)
        thetaOld = thetaNew;
        error = calculateError(set.u, set.y, thetaOld);
        derivativeOfB = calculateDerivativeOfB(set.u, thetaOld(1));
        derivativeOfF = calculateDerivativeOfF(set.y, error, thetaOld(1));
        derivative = [derivativeOfF, derivativeOfB].';
        gradient = zeros(2,1);
        hassian = zeros(2,2);
        for k = 1:length(error)
            gradient = gradient+(error(k)*derivative(:,k));
            hassian = hassian+(derivative(:,k)*(derivative(:,k).'));
        end
        gradient = gradient.*(2/length(error));
        hassian = hassian.*(2/length(error));
        thetaNew = thetaOld-alfa.*(hassian\gradient);
        l = l+1;
    end
    theta = thetaNew;
end

function e = calculateError(u, y, theta)
    e = zeros(length(u),1);
    thetaError = [1 theta(1) -theta(1) -theta(2)].';
    e(1) = y(1);
    for k = 2:length(u)
        parameters = [y(k) y(k-1) e(k-1) u(k-1)];
        e(k) = parameters*thetaError;
    end
end


function eB = calculateDerivativeOfB(u, f)
    eB = zeros(length(u),1);
    thetaB = [-1 -f].';
    for k = 2:length(u)
        parameters = [u(k-1) eB(k-1)];
        eB(k) = parameters*thetaB;
    end
end


function eF = calculateDerivativeOfF(y, e, f)
    eF = zeros(length(y),1);
    thetaF = [1 -1 -f].';
    for k = 2:length(y)
        parameters = [y(k-1) e(k-1) eF(k-1)];
        eF(k) = parameters*thetaF;
    end
end