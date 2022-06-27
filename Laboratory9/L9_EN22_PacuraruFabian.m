load lab9_4.mat
warning('off','all'); clc;
na = n; nb = n;
thetaARX = calculateThetaARX(id, na, nb);
y_id_hatARXpred = zeros(length(id.u),1);
y_val_hatARXpred = zeros(length(val.u),1);
for k = 1:length(id.u)
    philine = calculatePhiLine(id.u,id.y,k,na,nb);
    y_id_hatARXpred(k) = philine*thetaARX;
end
for k = 1:length(val.u)
    philine = calculatePhiLine(val.u,val.y,k,na,nb);
    y_val_hatARXpred(k) = philine*thetaARX;
end
y_id_hatARXsym = calculateY_hat(id.u, thetaARX, na, nb);
y_val_hatARXsym = calculateY_hat(val.u, thetaARX, na, nb);
idModelARXpred = iddata(y_id_hatARXpred, id.u, id.Ts);
valModelARXpred = iddata(y_val_hatARXpred, val.u, val.Ts);
idModelARXsym = iddata(y_id_hatARXsym, id.u, id.Ts);
valModelARXsym = iddata(y_val_hatARXsym, val.u, val.Ts);
thetaIV = calculateThetaIV(id, idModelARXsym, na, nb);
y_id_hatIVpred = zeros(length(id.u),1);
y_val_hatIVpred = zeros(length(val.u),1);
for k = 1:length(id.u)
    philine = calculatePhiLine(id.u,id.y,k,na,nb);
    y_id_hatIVpred(k) = philine*thetaIV;
end
for k = 1:length(val.u)
    philine = calculatePhiLine(val.u,val.y,k,na,nb);
    y_val_hatIVpred(k) = philine*thetaIV;
end
y_id_hatIVsym = calculateY_hat(id.u, thetaIV, na, nb);
y_val_hatIVsym = calculateY_hat(val.u, thetaIV, na, nb);
idModelIVpred = iddata(y_id_hatIVpred, id.u, id.Ts);
valModelIVpred = iddata(y_val_hatIVpred, val.u, val.Ts);
idModelIVsym = iddata(y_id_hatIVsym, id.u, id.Ts);
valModelIVsym = iddata(y_val_hatIVsym, val.u, val.Ts);
figure; compare(id, idModelARXpred, idModelIVpred);
figure; compare(val, valModelARXpred, valModelIVpred);
figure; compare(id, idModelARXsym, idModelIVsym);
figure; compare(val, valModelARXsym, valModelIVsym);


function theta = calculateThetaARX(id, na, nb)
    matrix = zeros(length(id.u), na+nb);
    for k = 1:length(id.u)
        phi_k = calculatePhiLine(id.u, id.y, k, na, nb);
        matrix(k,:) = phi_k;
    end
    Phi = matrix;
    theta = Phi\id.y;
end


function theta = calculateThetaIV(id, idModelARX, na, nb)
    Phi = zeros(na+nb, na+nb);
    Y = zeros(na+nb, 1);
    for k = 1:length(id.u)
        Z = (calculatePhiLine(idModelARX.u, idModelARX.y, k, na, nb)).';
        f = (calculatePhiLine(id.u, id.y, k, na, nb)).';
        Phi = Phi + Z*(f.');
        Y = Y + Z.*id.y(k);
    end
    Phi = Phi./length(id.u);
    Y = Y./length(id.u);
    theta = Phi\Y;
end


function y_hat = calculateY_hat(u, theta, na, nb)
    y_hat = zeros(length(u),1);
    for k = 1:length(u)
        phi_k = calculatePhiLine(u, y_hat, k, na, nb);
        y_hat(k) = phi_k*theta;
    end
end


function phi_k = calculatePhiLine(u, y, k, na, nb)
    phi_k = zeros(1, na+nb);
    for i = 1:na
       if(i<k)
           phi_k(i) = -y(k-i);
       end
    end
    for i = 1:nb
        if(i<k)
            phi_k(i+na) = u(k-i);
        end
    end
end