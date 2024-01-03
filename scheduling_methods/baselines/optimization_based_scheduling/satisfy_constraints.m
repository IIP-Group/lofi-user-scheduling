% =========================================================================
% -- Constraint checker
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
function [status_C,constraints] = satisfy_constraints(par,C)

[U] = size(C);

tol = 10^(-10);

sum_columns = sum(C,1);
sum_rows = sum(C,2);

for i=1:length(sum_columns)
    if(sum_columns(i) < par.min_sum_columns)
        if(abs(sum_columns(i) - par.min_sum_columns) < tol)
            sum_columns(i) = par.min_sum_columns;
        end
    elseif(sum_columns(i) > par.max_sum_columns)
        if(abs(sum_columns(i) - par.max_sum_columns) < tol)
            sum_columns(i) = par.max_sum_columns;
        end
    end
end

for i=1:length(sum_rows)
    if(sum_rows(i) < par.min_sum_rows)
        if(abs(sum_rows(i) - par.min_sum_rows) < tol)
            sum_rows(i) = par.min_sum_rows;
        end
    elseif(sum_rows(i) > par.max_sum_rows)
        if(abs(sum_rows(i) - par.max_sum_rows) < tol)
            sum_rows(i) = par.max_sum_rows;
        end
    end
end

constraints = zeros(3,1);

constraints(1) = (sum(sum((0 <= C) & (C <= 1))) == (U(1)*U(2))); % check all elements
constraints(2) = (sum(sum_rows >= par.min_sum_rows & sum_rows <= par.max_sum_rows)== U(1)); % check rows
constraints(3) = (sum(sum_columns >= par.min_sum_columns & sum_columns <= par.max_sum_columns)== U(2)); % check columns

if (sum(constraints) == length(constraints))
    status_C = 1;
else
    status_C = 0;
end

end