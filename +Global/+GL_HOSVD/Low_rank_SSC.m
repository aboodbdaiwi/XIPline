function  [X W]   =   Low_rank_SSC( Y1,Y2, nsig)

[Sigma2 U1] = hosvd2(full(Y1),full(Y2));                    % perform HOSVD

Sigma2(abs(Sigma2) < nsig*sqrt( 2*log(length(Y1(:))) ) )=0; % apply klocal thresholding to singular values
r   =   sum( abs(Sigma2(:))>0 );                            % sum of singular values
X   =   tprod(Sigma2, U1);                                  % get U only within kept singular values
wei =   1/(1+r);                                            % weight singular values
W   =   wei*ones( size(X) );                                % weight matrix
X   =   X*wei;
return;