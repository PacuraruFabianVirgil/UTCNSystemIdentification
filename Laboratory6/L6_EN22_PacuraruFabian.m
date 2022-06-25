load lab6_4.mat
warning('off','all'); clc;
figure('Name','Identification data'); plot(id)
figure('Name','Validation data'); plot(val)
na = 7;
nb = 7;
Phi_id = calculatePhi(id, na, nb);
theta = Phi_id\id.y;
y_id_hat = calculateY_hat(id.u, theta, na, nb);
y_val_hat = calculateY_hat(val.u, theta, na, nb);
model_id_hat = iddata(y_id_hat, id.u, id.Ts);
model_val_hat = iddata(y_val_hat, val.u, val.Ts);
figure;
subplot(2,1,1); compare(id, model_id_hat);
subplot(2,1,2); compare(val, model_val_hat);
model_id_arx = arx(id, [na nb 1]);
model_val_arx = arx(val, [na nb 1]);
figure;
subplot(2,1,1); compare(model_id_arx, model_id_hat);
subplot(2,1,2); compare(model_val_arx, model_val_hat);


function Phi = calculatePhi(set, na, nb)
    matrix = zeros(length(set.u), na+nb);
    for k = 1:length(set.u)
        phi_k = calculatePhiLine(set.u, set.y, k, na, nb);
        matrix(k,:) = phi_k;
    end
    Phi = matrix;
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
           phi_k(i) = y(k-i);
       end
    end
    for i = 1:nb
        if(i<k)
            phi_k(i+na) = u(k-i);
        end
    end
end