function [Q,R] = chol_qr_it(A)
  i=0;  
  cn = 200;  
  Q = A;
  G = Q'*Q;
  n = size(A,2);   
  R = eye(n);

  while cn > 100,
    i = i + 1
    [u,s,v]=svd(G);
    [q,r]=qr(sqrt(s)*v');
    R = r * R;
    cn = sqrt(cond(s));
    Q = Q * inv(r);
    if cn>100
      G = Q'*Q;
    end;
  end;  
return
  
