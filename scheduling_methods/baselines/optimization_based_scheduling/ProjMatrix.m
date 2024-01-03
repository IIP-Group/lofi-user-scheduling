% =========================================================================
% -- Projection to satisfy row and column constraint
% -------------------------------------------------------------------------
% Reference:
%           V. Palhares and C. Studer, "An Optimization-Based User
%           Scheduling Framework for mmWave Massive MU-MIMO Systems," 2022
%           IEEE 23rd International Workshop on Signal Processing Advances
%           in Wireless Communication (SPAWC), Oulu, Finland, 2022,
%           pp. 1-5, doi: 10.1109/SPAWC51304.2022.9833954.
%
% Last Updated: 20/12/2023
%
% -- (c) 2023 Victoria Palhares, Christoph Studer
% -- e-mails: <palhares@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================

function final_X = ProjMatrix(par,Z)

beta = 1; % tuning factor of Douglas-Rachford Splitting
tol = 10^(-10);

% beta = 1/(1+theta); % set to 1 for good measure
%
% f(X) = beta*0.5*|A-X|^2 + Chi1(X)
% g(X) = beta*0.5*|A-X|^2 + Chi2(X)
%
% Chi1(X) is indicator function for X'*1 = 1*N/M
% Chi2(X) is indicator function for X*1 = 1

err_2 = inf;

maxiter = 2000;

size_Z = size(Z);

Yk = zeros(size_Z);

for ii=1:maxiter
    % COLUMN CONSTRAINT: CHi1(X)
    Xk = Proximal_Algorithm(((beta*Z+Yk)/(beta+1)),par.min_sum_columns,par.max_sum_columns); % Xk = prox_f(Yk) = argmin beta*0.5*|A-X|^2 + 0.5*|X-Yk|^2
    sum_rows = sum(Xk,2);
    % CHECK ROW CONSTRAINT: lower bound <= num <= upper bound
    temp_1 = find(sum_rows < par.min_sum_rows);
    temp_2 = find(sum_rows > par.max_sum_rows);
    if(~isempty(temp_1)) % if at least one of the sums is below the lower bound
        err_1 = norm(sum_rows(temp_1)-(ones(size(temp_1,1),1)*par.min_sum_rows),'inf'); % infinity norm -> max(abs(xi)) for i=1...n
    elseif(~isempty(temp_2)) % if at least one of the sums is above the upper bound
        err_1 = norm(sum_rows(temp_2)-(ones(size(temp_2,1),1)*par.max_sum_rows),'inf'); % infinity norm -> max(abs(xi)) for i=1...n
    else % if the column constraint is satisfied
        err_1 = 0;
    end
    if (err_1<tol && err_2<tol)
        break;
    end

    % ROW CONSTRAINT: Chi2(X)
    Yk = (Yk-Xk+Proximal_Algorithm(((beta*Z+(2*Xk-Yk))/(beta+1)).',par.min_sum_rows,par.max_sum_rows).');
    sum_columns = sum(Xk,1);
    % CHECK COLUMN CONSTRAINT: lower bound <= num <= upper bound
    temp_1 = find(sum_columns < par.min_sum_columns);
    temp_2 = find(sum_columns > par.max_sum_columns);
    if(~isempty(temp_1)) % if at least one of the sums is below the lower bound
        err_2 = norm(sum_columns(temp_1)-(ones(1,size(temp_1,2))*par.min_sum_columns),'inf'); % infinity norm -> max(abs(xi)) for i=1...n
    elseif(~isempty(temp_2)) % if at least one of the sums is above the upper bound
        err_2 = norm(sum_columns(temp_2)-(ones(1,size(temp_2,2))*par.max_sum_columns),'inf'); % infinity norm -> max(abs(xi)) for i=1...n
    else % if the column constraint is satisfied
        err_2 = 0;
    end
    if (err_1<tol && err_2<tol)
        break;
    end

end
if ii==maxiter
    warning('maximum iterations of %i reached',maxiter);
end

final_X = Xk;

end






