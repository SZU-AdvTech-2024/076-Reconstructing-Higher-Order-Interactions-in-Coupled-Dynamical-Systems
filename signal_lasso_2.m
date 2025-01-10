% function of signal lasso method using cvx package
% Output is adjacency matrix of network

function t=signal_lasso_2(y,Phi,SIZE,Lambda1,Lambda2)
    %XSIZE=SIZE*SIZE;
    cvx_clear
    cvx_begin quiet
        variable x(SIZE) nonnegative;
        minimize(Lambda1*norm(x,1)+Lambda2*norm(x-1,1)+square_pos(norm(y-Phi*x,2))/2 );
    cvx_end
    t = x;
    clear x;
end

 