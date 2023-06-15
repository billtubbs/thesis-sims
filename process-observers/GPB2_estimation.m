% MATLAB Code of the GPB2 estimation algorthm
% Adapted from code provided by Shunyi Zhao in private email
%
% See this publication:
% -  Bayesian State Estimations for Markovian Jump Systems, 
%    by Shunyi Zhao, Choon Ki Ahn, Peng Shi, Yuriy S. Shmaliy, 
%    and Fei Liu, IEEE Systems, Man, & Cybernetics magazine, 2019.
%

function [x,P,mu,Out_x,Out_P] = GPB2_estimation(x,P,y,u,Model,mu)

    TP = Model.TP;  % transition probability matrix
    M = size(TP,1);  % number of modes (models)
    N = size(Model.A{1}, 1);  % dimension of state
    A = Model.A;
    Bu = Model.Bu;  % added a known input matrix
    Bw = Model.Bw;  % disturbance input matrix
    C = Model.C;
    D = Model.D;
    Q = Model.Q;
    R = Model.R;

    Mix_mu = nan(M,M);
    Mix_x = nan(N,M);
    Mix_P = nan(N,N,M);
    for j = 1:M

        Updx = nan(N,M);
        UpdP = nan(N,N,M);
        for i = 1:M

            % M^2 Kalman filter predictions
            Pred_x = A{j} * x(:,i) + Bu{j} * u;
            Pred_P = A{j} * P(:,:,i) * A{j}' + Bw{j} * Q{j} * Bw{j}';

            % Calculate the mixing probabilities
            S = C{j} * Pred_P * C{j}'+ D{j} * R{j} * D{j}';
            Mix_mu(i,j) = TP(i,j) * mu(i) * mvnpdf(y, C{j} * Pred_x, S);
            %assert(sum(Mix_mu(i,j) ~= 0),"Numerical underflow")

            % Update the M^2 KF estimates and error covariances
            K = Pred_P * C{j}' / S;  % or * S^-1;
            Updx(:,i) = Pred_x + K *(y - C{j} * Pred_x);
            %UpdP(:,:,i) = Pred_P - K * C{j} * Pred_P;

            % Alternative method to update P(k).
            % This is the 'Joseph Stabilized' version of
            % the update equation which guarantees that Pk is positive
            % semi-definite in the presence of roundoff error.
            % (See p73 of Lewis etal. Book Optimal and Robust Estimation).
            z = eye(size(Pred_P)) - K * C{j};
            UpdP(:,:,i) = z * Pred_P * z' + K * D{j} * R{j} * D{j}' * K';

        end

        % Normalize the mode probabilities
        Upd_u = Mix_mu(:,j) / sum(Mix_mu(:,j));

        % Mix the estimates
        Mix_x(:,j) = sum(Updx .* repmat(Upd_u', N, 1), 2);
        Mix_P(:,:,j) = zeros(N, N);
        for i = 1:M
            summP = Upd_u(i) * (UpdP(:,:,i) + ...
                (Updx(:,i) - Mix_x(:,j)) * (Updx(:,i) - Mix_x(:,j))');
            Mix_P(:,:,j) =  Mix_P(:,:,j) + summP;
        end

    end

    % Updated mixed state estimates and error covariances
    x = Mix_x;
    P = Mix_P;

    % Normalize the mixing probabilities
    mu = sum(Mix_mu, 1)' ./ sum(Mix_mu, [1 2]);
    % This is the same as having Stormu(j) = sum(Mix_mu(:,j)); in
    % the loop above and then mu = Stormu / sum(Stormu);

    % Outputs
    Out_x = sum(x .* mu', 2);
    % Above is same as Out_x = sum(x .* repmat(mu, N, 1), 2)
    Out_P = zeros(N, N);
    for i = 1:M
        summP = mu(i) * (P(:,:,i) + (x(:,i) - Out_x) * (x(:,i) - Out_x)');
        Out_P = Out_P + summP;
    end

end