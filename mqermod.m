function x=mqermod(Psim,X0,A,B)
%Esta funcao e uma modificacao do algoritmo de Minimos Quadrados com Restricao implementado
%por Marcio Barroso do grupo MACSIN. A funcao retorna o vetor simetrico de pesos da rede. Os
%parametros de entrada sao a matriz de regressores obtida com os centros simetricos Psim, o
%vetor de pesos estimado pelo algoritmo de Minimos Quadrados X0, a matriz com valores de 
%restricoes A e o vetor com os coeficientes de restricoes B.    


%Coloca o problema na forma de Lagrange
%------------------------------------------------
x_corr = inv(Psim'*Psim)*A'*inv(A*inv(Psim'*Psim)*A')*(A*X0-B);
%------------------------------------------------

%Resultado da correção
%------------------------------------------------
x = X0 - x_corr;
%------------------------------------------------
