load uval.mat
index = 5;
val = system_simulator(index,u);
arxdata = arx(val,[15 15 1]);
figure; compare(val,arxdata); title('Validation data');
u_1 = prbsgenerator(300,3,0.5,1);
u_2 = prbsgenerator(300,10,0.5,1);
id_1 = system_simulator(index,u_1);
arxdata_1 = arx(id_1,[15 15 1]);
figure; compare(id_1,arxdata_1); title('m = 3');
id_2 = system_simulator(index,u_2);
arxdata_2 = arx(id_2,[15 15 1]);
figure; compare(id_2,arxdata_2); title('m = 10');


function [u] = prbsgenerator(N,m,b,c)
    a = zeros(1,m);
    switch m
        case 3
            a(1) = 1;
            a(3) = 1;
        case 4
            a(1) = 1;
            a(4) = 1;
        case 5
            a(2) = 1;
            a(5) = 1;
        case 6
            a(1) = 1;
            a(6) = 1;
        case 7
            a(1) = 1;
            a(7) = 1;
        case 8
            a(1) = 1;
            a(2) = 1;
            a(7) = 1;
            a(8) = 1;
        case 9
            a(4) = 1;
            a(9) = 1;
        case 10
            a(3) = 1;
            a(10) = 1;
        otherwise
            return;
    end
    X = ones(N, m);
    for i=2:N
        for k = 1:length(a)
            X(i, 1) = mod(X(i,1), a(k)*X(i-1,k));
        end
        for j = 2:m
            X(i, j) = X(i-1, j-1);
        end
    end
    u = X(:, end);
    u = b + (c-b).*u;
end