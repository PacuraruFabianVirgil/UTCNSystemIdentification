%System Identification Project - Part 1
%Fitting an unknown function

%Students: Nemeth Raymond, Pacuraru Fabian-Virgil, Zah Elena
%Group indexes: 07/15

%% Reading and formatting input data
load ("proj_fit_07.mat")
x1_id = id.X{1, 1};
x2_id = id.X{2, 1};
x1_val = val.X{1, 1};
x2_val = val.X{2, 1};
y_id = id.Y(:);
y_val = val.Y(:);
clc; warning ('off', 'all');

%% Plotting input data
figure;
subplot(2, 1, 1); mesh(x2_id, x2_id, id.Y')
colorbar; title('Input data for identification')
subplot(2, 1, 2); mesh(x1_val, x2_val, val.Y');
colorbar; title('Input data for validation')

%% Computing lowest MSE by varying the degree of the approximating polynomial
maxM = 20; % maximum degree of the polynomial to try
MSEid = zeros(1, maxM);
MSEval = zeros(1, maxM);
for m=1:maxM
    phi_id = calculatePhi(m, id, x1_id, x2_id);
    phi_val = calculatePhi(m, val, x1_val, x2_val);
    theta = phi_id\y_id;
    y_id_approx = phi_id*theta;
    y_val_approx = phi_val*theta;
    MSEid(m) = sum((y_id_approx-y_id).^2)/id.dims(1)/id.dims(2);
    MSEval(m) = sum((y_val_approx-y_val).^2)/val.dims(1)/val.dims(2);
end

%% Plotting the obtained MSEs with respect to the polynomial degree
figure; plot(2:maxM, MSEid(2:maxM), 2:maxM, MSEval(2:maxM))
xlabel("Degree of the polynomial"); ylabel("MSE");
legend("Identification dataset", "Validation dataset")
[~, bestM] = min(MSEval);
fprintf(['\nThe smallest mean squared error on the validation data' ...
    ' was obtained for degree %.2d.\n'],bestM)

%% Recalculating the results for the best MSE 
phi_id = calculatePhi(bestM, id, x1_id, x2_id);
phi_val = calculatePhi(bestM, val, x1_val, x2_val);
theta = phi_id\y_id;
y_id_approx = phi_id*theta;
y_val_approx = phi_val*theta;
y_id_approx = reshape(y_id_approx, id.dims);
y_val_approx = reshape(y_val_approx, val.dims);

%% Plotting results compared to the desired outputs
doPlot(x1_id, x2_id, id.Y, y_id_approx, 'Identification');
doPlot(x1_val, x2_val, val.Y, y_val_approx, 'Validation');

%% Function for computing phi
function phi = calculatePhi(m, set, x1, x2)
    % Preallocation of memory for the matrix phi and coeffs vector is 
    % necessary for reducing the running time
    phi = zeros(set.dims(1)*set.dims(2), (m+1)*(m+2)/2);
    coeffs = zeros(1, (m+1)*(m+2)/2);
    
    % i and j are used for iterating over every x1 and x2
    for i = 1:set.dims(1)
        for j = 1:set.dims(2)
            count = 1;
            % p1 and p2 represent the powers for x1 and x2
            for p1 = 0:m
                for p2 = 0:m
                    if p1+p2<=m
                        coeffs(count) = x1(i)^(p1)*x2(j)^(p2);
                        count = count+1;
                    end
                end
            end
            % for every pair of x1(i)-x2(j) a vector of regressors is added
            % to the matrix phi
            phi((i-1)*set.dims(2)+j,:) = coeffs;
        end
    end
end

%% Function for plotting the real and approximated output on the same figure
function f = doPlot(x1, x2, y, y_approx, whichSet)
    f = figure('Name',strcat(whichSet, ' Data'));
    subplot(2, 1, 1); mesh(x1, x2, y);
    colorbar; title('Real Output')
    subplot(2, 1, 2); mesh(x1, x2, y_approx);
    colorbar; title('Approximated Output')
end