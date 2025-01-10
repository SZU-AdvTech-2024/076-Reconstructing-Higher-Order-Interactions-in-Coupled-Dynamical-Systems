% Function of LASSO code using cvx package

function t=lasso_2(y,Fai,SIZE,lambda)
  %XSIZE=SIZE*SIZE;
  cvx_begin
      variable x(SIZE) nonnegative;
      minimize(lambda*norm(x,1)+square_pos(norm(y-Fai*x,2)));        
  cvx_end
  t = x;
  clear x;
end

