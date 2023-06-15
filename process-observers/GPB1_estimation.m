% MATLAB Code of the GPB1 estimation algorthm
% Adapted from code provided by Shunyi Zhao in private email
%
% See this publication:
% -  Bayesian State Estimations for Markovian Jump Systems, 
%    by Shunyi Zhao, Choon Ki Ahn, Peng Shi, Yuriy S. Shmaliy, 
%    and Fei Liu, IEEE Systems, Man, & Cybernetics magazine, 2019.
%

function [x,P,mu] = GPB1_estimation(x,P,y,u,Model,mu)

    TP = Model.TP;  % transition probability matrix
    M = size(TP,1);  % number of modes (models)
    N = size(Model.A{1},1);  % dimension of state
    A = Model.A;
    Bu = Model.Bu;  % added a known input matrix
    Bw = Model.Bw;  % disturbance input matrix
    C = Model.C;
    D = Model.D;
    Q = Model.Q;
    R = Model.R;

    liklihood_u = nan(1, M);
    updx = nan(N,M);
    updP = nan(N,N,M);
    for j = 1:M

        % M Kalman filter predictions
        Pred_u = mu' * TP(:,j);
        Pred_x = A{j} * x + Bu{j} * u;
        Pred_P = A{j} * P * A{j}' + Bw{j} * Q{j} * Bw{j}';
        S = C{j} * Pred_P * C{j}'+ D{j} * R{j} * D{j}';

        % Update the conditional likelihood given the data
        liklihood_u(j) = Pred_u * mvnpdf(y, C{j}*Pred_x, S);

        % Update the M KF estimates and error covariances
        K = Pred_P * C{j}' / S;  % or * S^-1;
        updx(:,j)   = Pred_x + K * (y - C{j} * Pred_x);
        %updP(:,:,j) = Pred_P - K * C{j} * Pred_P;

        % Alternative method to update P(k).
        % This is the 'Joseph Stabilized' version of
        % the update equation which guarantees that Pk is positive
        % semi-definite in the presence of roundoff error.
        % (See p73 of Lewis etal. Book Optimal and Robust Estimation).
        z = eye(size(Pred_P)) - K * C{j};
        updP(:,:,j) = z * Pred_P * z' + K * D{j} * R{j} * D{j}' * K';

    end

    % Normalize conditional likelihood
    mu = liklihood_u' / sum(liklihood_u);

    % Updated state estimates and error covariance
    x = sum(updx .* mu', 2);      % mixed
    % Above is same as Out_x = sum(updx .* repmat(mu, N, 1), 2)
    P = zeros(N,N);
    for i = 1:M
        summP = mu(i) * (updP(:,:,i) + (updx(:,i) - x) * (updx(:,i) - x)');
        P = P + summP;
    end

end