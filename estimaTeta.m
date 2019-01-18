%ESTIMATETA    Estima��o de par�metros de uma rede RBF.

%Esta rotina estima os par�metros de uma rede RBF, que s�o recebidos no vetor Teta.

%As vari�veis de entrada no argumento da fun��o ESTIMATETA s�o:

%c - Matriz contendo os centros. Esta matriz � composta de forma que cada linha contenha
%um centro. Assim, ela possui dimens�o nr (n�mero de centros) X nu + ny (soma dos m�ximos 
%atrasos na entrada u(k) e na sa�da y(k)).
%nr - N�mero de centros.
%dp - Desvio padr�o da fun��o de base (nesta rotina utiliza-se como fun��o de base a 
%gaussiana).
%ut - Vetor de entradas para treinamento.
%Nt - Comprimento do vetor ut ou yt.
%yt - Vetor de sa�das para treinamento.
%nu, ny - M�ximos atrasos na sa�da e na entrada, respectivamente.
%linear - Se o valor de linear for igual a 1, a rotina estimar� par�metros para os termos
%lineares tamb�m. Sen�o, isto n�o ocorrer� e o modelo obtido n�o ter� termos lineares.

function Teta=estimaTeta(c, nr, dp, ut, Nt, yt, nu, ny, linear)
%Primeiramente, verifica-se qual dos m�ximos atrasos � maior, para, � medida que vai se
%pegando dados nos vetores ut e yt, n�o extrapolar o limite destes.
if nu>=ny
   %Nas linhas abaixo � montada a matriz Psi dos regressores.
   Psi=ones(Nt-nu,1);
   X=ones(Nt-nu,1);
   for i=1:nu
      X=[X ut(nu+1-i:Nt-i)];
   end
   for l=1:ny
      X=[X yt(nu+1-l:Nt-l)];
   end
   X(:,1)=[];
   if nr>0 %Esta linha testa a possibilidade de uma rede com somente os termos lineares.
      for k=1:Nt-nu
         for j=1:nr
            V(k,j)=((X(k,:)-c(j,:))*((X(k,:)-c(j,:))'))^0.5;
            %if V(k,j)==0
	         %   V(k,j)=.0000001;   
	         %end   
         end
      end
      %Fi=V.*V.*(log10(V));
      Fi=exp(-1*V.*V*1/(dp^2));
      Psi=[Psi Fi];
   end
   %Nas tr�s linhas abaixo � testado o valor de linear e, se acaso esta vari�vel contiver 
   %o valor 1, s�o inclusos regressores dos termos lineares na matriz dos regressores.
   if linear==1
      Psi=[Psi X];
   end
   %Finalizada a montagem da matriz Psi dos regressores, o vetor Teta � estimado a partir 
   %da solu��o do algoritmo dos M�nimos Quadrados.
   Teta=((inv((Psi')*Psi))*(Psi'))*yt(nu+1:Nt);
%A partir daqui � seguido o mesmo procedimento, somente quando ny > nu.   
else 
   %Montagem da matriz Psi de regressores.
   Psi=ones(Nt-ny,1);
   X=ones(Nt-ny,1);
   for i=1:nu
      X=[X ut(ny+1-i:Nt-i)];
   end
   for l=1:ny
      X=[X yt(ny+1-l:Nt-l)];
   end
   X(:,1)=[];
   if nr>0 %Teste para inclus�o dos termos que comp�em a forma generalizada de uma RBF.
      for k=1:Nt-ny
         for j=1:nr
            V(k,j)=((X(k,:)-c(j,:))*((X(k,:)-c(j,:))'))^0.5;
         end
      end
      Fi=exp(-1*V.*V*1/(dp^2));     
      Psi=[Psi Fi];
   end
   %Teste para inclus�o dos regressores dos termos lineares na matriz dos regressores.
   if linear==1
      Psi=[Psi X];
   end
   %Ap�s finalizada a montagem da matriz Psi dos regressores, os par�metros s�o estimados 
   %a partir do algoritmo dos M�nimos Quadrados e recebidos no vetor Teta.
   Teta=((inv((Psi')*Psi))*(Psi'))*yt(ny+1:Nt);
end            
