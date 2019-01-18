function [Caic, pivlin]=ctrsim(Piv, tipo, nr, X)
%Esta funcao monta uma matriz que contem centros simetricos, em que cada linha corresponde 
%a um centro e o numero de colunas da matriz e igual a nu+ny. A segunda metade dos centros
%e igual aa primeira metade multiplicada por -1. A matriz de centros simetricos e montada a 
%partir do vetor Piv, retornado pela funcao MYHOUSE e do numero de regressores (que incluem 
%centros e/ou termos lineares) que se deseja ter, passado no parametro tipo. A funcao retorna
%a matriz de centros simetricos e um vetor pivlin, que contem indices dos termos lineares, 
%apresentando valores que podem estar contidos na faixa de 1 ate nul+nyl.

pivlin=[];
Caic=[];
%Determina tipo
piv=sort(Piv(1:tipo));
for i=1:tipo
   %Monta matriz de centros sem simetria ainda
   if (piv(i)<=nr+1) & (piv(i)~=1)
      Caic=[Caic; X(piv(i)-1,:)];
   end
   %Verifica inclusao de termos lineares
   if piv(i)>nr+1
      pivlin=[pivlin piv(i)-(nr+1)];
   end   
end
   
Caic=[Caic; -Caic];


      


