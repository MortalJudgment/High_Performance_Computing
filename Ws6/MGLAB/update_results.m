%========================
% Update array of results
%========================

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function results = update_results (results,method,iter,rn)

%jr results = [results; [iter, flops, toc, rn]];
results = [results; [iter, 0, toc, rn]];
disp (sprintf('%s iteration %g  residual %g',method, iter,rn));



