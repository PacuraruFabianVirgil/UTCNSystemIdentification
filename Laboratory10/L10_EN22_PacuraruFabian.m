load lab10_6.mat
warning('off','all'); clc;
na = 3*n; nb = 3*n;
s = 1e-4;
theta = calculateTheta(id.u,id.y,s,na,nb).';
A = ones(1,na+1);
B = zeros(1,nb+1);
for i = 1:na
    A(i+1) = theta(end,i);
end
for i = 1:nb
    B(i+1) = theta(end,na+i);
end
model = idpoly(A,B,[],[],[],0,id.Ts);
figure; compare(id,model);
figure; compare(val,model);
X = rarx([id.y id.u],[na nb 1],'ff',1);
figure;
for i = 1:na
    subplot(3,n,i); hold on; plot(theta(:,i));
    plot(X(:,i),'--'); hold off; title(strcat('a',num2str(i)));
end
figure;
for i = 1:nb
    subplot(3,n,i); hold on; plot(theta(:,i+na));
    plot(X(:,i+na),'--'); hold off; title(strcat('b',num2str(i)));
end


function theta = calculateTheta(u,y,s,na,nb)
    Pinv = 1/s.*eye(na+nb);
    theta = zeros(na+nb,length(y));
    for k = 2:length(y)
        phi = calculatePhiLine(u,y,k,na,nb).';
        e = y(k)-(phi.')*theta(:,k-1);
        Pinv = Pinv-Pinv*phi*(phi.')*Pinv./(1+(phi.')*Pinv*phi);
        W = Pinv*phi;
        theta(:,k) = theta(:,k-1)+W.*e;
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