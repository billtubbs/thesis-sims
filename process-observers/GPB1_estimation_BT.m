% MATLAB Code of the GPB1
% Adapted from:
% -  Bayesian State Estimations for Markovian Jump Systems, 
%    by Shunyi Zhao, Choon Ki Ahn, Peng Shi, Yuriy S. Shmaliy, 
%    and Fei Liu, 2019.
%

function [x,P,u] = GPB1_estimation(x, P, y, Model, u)
    TP = Model.TP;
    M = size(TP,1);
    N = size(Model.A,1);
    A = Model.A;
    B = Model.B;
    C = Model.C;
    D = Model.D;
    Q = Model.Q;
    R = Model.R;

    for j = 1:M
        Pred_u = TP(:,j) .* u;
        Pred_x = A{j} * x;
        Pred_P = A{j} * P * A{j}' + Q{j};
        S = C{j} * Pred_P * C{j}' + R{j};
        K = Pred_P * C{j}' / S;
        liklihood_u(j) = sum(Pred_u) * mvnpdf(y, C{j} * Pred_x, S);
        updx(:,j) = Pred_x + K *(y - C{j} * Pred_x);
        updP(:,:,j) = Pred_P - K * C{j} * Pred_P;
    end

    u = liklihood_u / sum(liklihood_u);
    x = sum(updx .* repmat(u,N,1), 2);
    P = zeros(N,N);
    for i = 1:M
        summP = u(i) * (updP(:,:,i) + (updx(:,i) - x) * (updx(:,i) - x)');
        P = P + summP;
    end
    u = u';
end
