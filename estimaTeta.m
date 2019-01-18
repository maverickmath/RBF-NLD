%ESTIMATETA    Estimação de parâmetros de uma rede RBF.

%Esta rotina estima os parâmetros de uma rede RBF, que são recebidos no vetor Teta.

%As variáveis de entrada no argumento da função ESTIMATETA são:

%c - Matriz contendo os centros. Esta matriz é composta de forma que cada linha contenha
%um centro. Assim, ela possui dimensão nr (número de centros) X nu + ny (soma dos máximos 
%atrasos na entrada u(k) e na saída y(k)).
%nr - Número de centros.
%dp - Desvio padrão da função de base (nesta rotina utiliza-se como função de base a 
%gaussiana).
%ut - Vetor de entradas para treinamento.
%Nt - Comprimento do vetor ut ou yt.
%yt - Vetor de saídas para treinamento.
%nu, ny - Máximos atrasos na saída e na entrada, respectivamente.
%linear - Se o valor de linear for igual a 1, a rotina estimará parâmetros para os termos
%lineares também. Senão, isto não ocorrerá e o modelo obtido não terá termos lineares.

function Teta=estimaTeta(c, nr, dp, ut, Nt, yt, nu, ny, linear)
%Primeiramente, verifica-se qual dos máximos atrasos é maior, para, à medida que vai se
%pegando dados nos vetores ut e yt, não extrapolar o limite destes.
if nu>=ny
   %Nas linhas abaixo é montada a matriz Psi dos regressores.
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
   %Nas três linhas abaixo é testado o valor de linear e, se acaso esta variável contiver 
   %o valor 1, são inclusos regressores dos termos lineares na matriz dos regressores.
   if linear==1
      Psi=[Psi X];
   end
   %Finalizada a montagem da matriz Psi dos regressores, o vetor Teta é estimado a partir 
   %da solução do algoritmo dos Mínimos Quadrados.
   Teta=((inv((Psi')*Psi))*(Psi'))*yt(nu+1:Nt);
%A partir daqui é seguido o mesmo procedimento, somente quando ny > nu.   
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
   if nr>0 %Teste para inclusão dos termos que compõem a forma generalizada de uma RBF.
      for k=1:Nt-ny
         for j=1:nr
            V(k,j)=((X(k,:)-c(j,:))*((X(k,:)-c(j,:))'))^0.5;
         end
      end
      Fi=exp(-1*V.*V*1/(dp^2));     
      Psi=[Psi Fi];
   end
   %Teste para inclusão dos regressores dos termos lineares na matriz dos regressores.
   if linear==1
      Psi=[Psi X];
   end
   %Após finalizada a montagem da matriz Psi dos regressores, os parâmetros são estimados 
   %a partir do algoritmo dos Mínimos Quadrados e recebidos no vetor Teta.
   Teta=((inv((Psi')*Psi))*(Psi'))*yt(ny+1:Nt);
end            
