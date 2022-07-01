%System Identification Project - Part 2
%Nonlinear ARX

%Students: Nemeth Raymond, Pacuraru Fabian-Virgil, Zah Elena
%Group indexes: 07/15

warning('off','all');
clc; clearvars; close all;
load iddata-15.mat

% Initializing MSE matrices
% Rows - m; Columns - na=nb;
MSE_prediction = zeros(5,3);
MSE_simulation = zeros(5,3);

%% Running the algorithms for different orders of the approximator (m) and
% different number of parameters (na, nb)
for m = 1:5
    for na = 1:3
        nb = na;
        matrixOfPowers = powers(m, na, nb);
        phi_id = calculatePhi(id, na, nb, matrixOfPowers);
        phi_val  = calculatePhi(val, na, nb, matrixOfPowers);
        theta = phi_id\id.y; clc
        y_val_prediction = phi_val*theta;
        y_val_simulation = calculateY_hat(val.u, theta, na, nb, matrixOfPowers);
        
        % saving the MSE
        MSE_prediction(m, na) = calculateMSE(val.y, y_val_prediction)
        MSE_simulation(m, na) = calculateMSE(val.y, y_val_simulation)
    end
end
%% Plotting Mean Squared Errors w.r.t. m, na, nb
plotMSE(MSE_prediction, 'MSE on prediction model');
plotMSE_log(MSE_prediction, 'MSE on prediction model, logaritmic scale');
plotMSE(MSE_simulation, 'MSE on simulation model');
plotMSE_log(MSE_simulation, 'MSE on simulation model, logaritmic scale');

%% Plotting the best approximations
[bestPred, bestSim] = plotBestModel(MSE_prediction, MSE_simulation, id, val);
fprintf(['\nSmallest MSE on predicted validation data is %.2d for m=%.2d and na=nb=%.2d.\nMSE on the identification data is %.2d.\n'],bestPred);
fprintf(['\nSmallest MSE on simulated validation data is %.2d for m=%.2d and na=nb=%.2d.\nMSE on the identification data is %.2d.\n'],bestSim);


%% Functions
function powersMat = powers(m,na,nb)
%   Returns a matrix of powers.
%   Each parameter of phi is computed raising the delayed signals to the
%   powers contained in one line of powersMat.
    dim = 1;
    for i=1:(na+nb)
        dim = dim*(m+i)/i;
    end
    p = zeros(1,na+nb);
    powersMat = zeros(dim,na+nb);
    k = 1;
    while(k<=dim)
        powersMat(k,:)=p;
        if(sum(p)<m)
            c = na+nb;
        else
            if(c<(na+nb))
                p(c:end) = 0;
            else
                p(end) = 0;
            end
            c = max(c-1,1);
        end
        p(c) = p(c)+1;
        k = k+1;
    end
end

function Phi = calculatePhi(set, na, nb, pMat)
%   Returns the matrix of parameters phi. Uses calculatePhiLine to build
%   phi line-by-line.
    matrix = zeros(length(set.y), length(pMat));
    for k = 1:length(set.y)
        phi_k = calculatePhiLine(set.u, set.y, k, na, nb, pMat);
        matrix(k,:) = phi_k;
    end
    Phi = matrix;
end

function y_hat = calculateY_hat(u, theta, na, nb, pMat)
%   Returns the vector of simulated outputs.
    y_hat = zeros(length(u),1);
    for k = 1:length(u)
        phi_k = calculatePhiLine(u, y_hat, k, na, nb, pMat);
        y_hat(k) = phi_k*theta;
    end
end

function phi_k = calculatePhiLine(u, y, k, na, nb, pMat)
%   Returns all the combinations of the elements of the delayed signals
%   vector raised at the corresponding powers.
    phi_k = ones(1, length(pMat));
    param = parameters(u,y,k,na,nb);
    for i = 1:length(pMat)
        for j = 1:(na+nb)
            prod = param(j)^pMat(i,j);
            phi_k(i) = phi_k(i)*prod;
        end
    end
end

function param = parameters(u, y, k, na, nb)
%   Returns the vector of delayed inputs and outputs.
    param = zeros(1, na+nb);
    for i = 1:nb
       if(i<k)
           param(i) = u(k-i);
       end
    end
    for i = 1:na
        if(i<k)
            param(i+nb) = y(k-i);
        end
    end   
end

function [bestP, bestS] = plotBestModel(MSE_pred, MSE_sim, id, val)
%   Recalculates and plots the best models, based the given MSE vectors.
%   Two figures will be made, containing:
%       1. Identification outputs: predicted, simulated, real
%       2. Validation outputs: predicted, simulated, real
%   The function returns two vectors, each containing minMSE, m, na, MSEid 
%   for the predicted, respectively the simulated outputs

    %[minMse_val,m,na, MSE_id] 
    bestP=[min(MSE_pred(:)), 0, 0, 0];
    bestS=[min(MSE_sim(:)), 0, 0, 0];

    [bestP(2), bestP(3)] = find(bestP(1) == MSE_pred);
    [bestS(2), bestS(3)] = find(bestS(1) == MSE_sim);
    
    matrixOfPowers_P = powers(bestP(2), bestP(3), bestP(3));
    phi_id_p = calculatePhi(id, bestP(3), bestP(3), matrixOfPowers_P);
    phi_val_p  = calculatePhi(val, bestP(3), bestP(3), matrixOfPowers_P);
    theta_p = phi_id_p\id.y;
    
    matrixOfPowers_S = powers(bestS(2), bestS(3), bestS(3));
    phi_id_s = calculatePhi(id, bestS(3), bestS(3), matrixOfPowers_S);
    theta_s = phi_id_s\id.y;
    
    %Ploting identification outputs for prediction, simulation and real
    y_id_prediction = phi_id_p*theta_p;
    y_id_simulation = calculateY_hat(id.u, theta_s, bestS(3), bestS(3), matrixOfPowers_S);
    figure
    plot(y_id_prediction, 'LineWidth', 2); hold on;
    plot(y_id_simulation, 'r--', 'LineWidth', 2);
    plot(id.y, ':', 'color', '#EDB120', 'LineWidth', 2); hold off;
    legend('Prediction', 'Simulation', 'Real'); title('Identification');
    
    %Calculating MSE for identification data
    bestP(4)=calculateMSE(id.y, y_id_prediction);
    bestS(4)=calculateMSE(id.y, y_id_simulation);
    
    %Ploting validation outputs for prediction, simulation and real
    y_val_prediction = phi_val_p*theta_p;
    y_val_simulation = calculateY_hat(val.u, theta_s, bestS(3), bestS(3), matrixOfPowers_S);
    figure;
    plot(y_val_prediction, 'LineWidth', 2); hold on;
    plot(y_val_simulation, 'r--', 'LineWidth', 2);
    plot(val.y, ':', 'color', '#EDB120', 'LineWidth', 2); hold off;
    legend('Prediction', 'Simulation', 'Real'); title('Validation');
end

function MSE = calculateMSE(y_real, y_approx)
    % Returns MSE between two given vectors
    MSE = sum((y_approx-y_real).^2)/length(y_approx);
end

function plotMSE(MSE, text)
%   Plots MSE w.r.t m, na, nb
    figure;
    bar(MSE);
    title(text);
    legend('na=nb=1', 'na=nb=2', 'na=nb=3');
    xlabel('m'); ylabel('MSE');
end

function plotMSE_log(MSE, text)
%   Plots MSE w.r.t m, na, nb, on logarithmic scale
    figure;
    bar(MSE);
    set(gca,'YScale','log')
    title(text);
    legend('na=nb=1', 'na=nb=2', 'na=nb=3');
    xlabel('m'); ylabel('MSE');
end