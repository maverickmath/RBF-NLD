function [A, B]=mtrsres(nraic, pivlin)
%Esta funcao monta a matriz com valor de restricoes A e os coeficientes de restricoes B 
%para o caso de imposicao de simetria nos pesos e nos centros da rede RBF, de forma que 
%B=A*Teta. Os parametros de entrada sao o numero de centros simetricos nraic, o vetor 
%pivlin obtido na funcao CTRSIM. Para o caso de imposicao de simetria somente nos centros 
%nao e necessario chamar esta funcao.
   
B=zeros((nraic/2),1);
A=zeros((nraic/2),nraic+length(pivlin));
for r=1:(nraic/2)
   A(r,r)=1;
   A(r,r+(nraic/2))=1;
end
   
   
  