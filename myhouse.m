function [A,err,Piv]=myhouse(Psi,np);
% [A,err,Piv]=myhouse(Psi,np);
% Matrix Computations 2nd Ed. pg 212
% Given matrix Psi (m,n), this function finds Q such that
% Q'*Psi=V is top-triangular. The top triangular part of A
% is replaced by the top-triangular part of V.
% 
% It is assumed that the last column of Psi is the vector
% of observations, y(k).
% np is the number of regressors chosen to composed the model
% err is a vector of np values that indicate the error reduction ratio
%	of each of the chosen regressors
% Piv is a vector that contains the chosen regressor indices, that is,
% 	such indices indicate the columns that were used in pivoting matrix
%	Psi np times.

[m,n]=size(Psi);
A=Psi;
yy=Psi(:,n)'*Psi(:,n);
piv=1:n-1;

for j=1:np	
   
   for k=j:n-1 % ate completar o numero de termos candidatos
   	c(k)=((A(j:m,k)'*A(j:m,n))^2)/((A(j:m,k)'*A(j:m,k))*yy);  % err of regressor k
	end;
   
   [ans aux]=max(c(j:n-1));
   jm=j+aux-1;
   err(j)=ans;
   aux=A(:,jm); % column of regressor with greatest err
	A(:,jm)=A(:,j);
   A(:,j)=aux;
   aux=piv(jm); % index of regressor with greatest err
	piv(jm)=piv(j);
	piv(j)=aux;
   
   
   x=A(j:m,j);
   % v=house(x)
	% Matrix Computations 2nd Ed. pg 196
	% given vector x, returns vector v such that
	% (I-2vv'/v'v)x is zero but for the first element

	nx=length(x);
	u=norm(x,2);
	v=x;
	if u ~= 0
	  b=x(1) + sign(x(1))*u;
	  v(2:nx) = v(2:nx)/b;
	end;
	v(1)=1;
   % end house(x)

   a=A(j:m,j:n);
   
	% a=rowhouse(a,v)
	% Matrix Computations 2n Ed. pg 197
	% Given a matrix A (m,n), and a vector v of length m 
	% whose first element is 1, this algorithm replaces
	% A with P*A, where P=I-2vv'/v'v

	b=-2/(v'*v);
	w=b*a'*v;
   a=a+v*w';
   % endrowhouse(a,v)
   
   A(j:m,j:n)=a;   

end;
% end myhouse(A)
Piv=piv(1:np);

