% function of NNLS method using cvx package
% Output is adjacency matrix of network

function t = NNLS_cvx(y,Phi,SIZE)
    %XSIZE=SIZE*SIZE;
    cvx_clear
    cvx_begin quiet
        variable x(SIZE) nonnegative;
        minimize(square_pos(norm(y-Phi*x,2))/2 );
    cvx_end
    t = x;
    clear x;
end
