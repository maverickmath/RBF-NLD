%SIMULA��O    Simula��o livre de uma rede RBF com restri��es de estrutura.

%Esta fun��o faz a simula��o livre de uma rede RBF. O vetor y recebe os valores calculados 
%pela rede na ordem temporal. Assim, para verificar graficamente a simula��o, basta somente
%executar um comando plot(y) no prompt do MATLAB.

%As vari�veis de entrada no argumento da fun��o SIMULA��O s�o:

%c - Matriz contendo os centros. Esta matriz � composta de forma que cada linha contenha
%um centro. Assim, ela possui dimens�o nr (n�mero de centros) X nu + ny (soma dos m�ximos 
%atrasos na entrada u(k) e na sa�da y(k)).
%nr - N�mero de centros.
%nu, ny - M�ximos atrasos na sa�da e na entrada da parte linear, respectivamente.
%nurbf, nyrbf - M�ximos atrasos na sa�da e na entrada das RBFs, respectivamente.
%N - N�mero de pontos que pretende-se simular. Deve ser igual ao comprimento do vetor us.
%us - Vetor de entradas contendo os valores das entradas na janela de simula��o.
%yi - Vetor contendo os valores iniciais da sa�da. Este vetor deve ter comprimento ny. 
%Teta - Vetor de par�metros da rede.
%linear - Se o valor de linear for igual a 1, a rede possui termos lineares, sen�o, n�o 
%possui.
%dp - Desvio padr�o da fun��o de base (nesta rotina utiliza-se como fun��o de base a 
%gaussiana).
%const - 1 para uma rede com termo de offset, 0 para uma rede sem offset


function y=simulacao_r(c,nr,nu,ny,nurbf,nyrbf,const,N,us,yi,Teta,linear,dp)
%Primeiramente � testado qual dos m�ximos atrasos � maior para n�o extrapolar os limites
%dos vetores de entrada e/ou sa�da.
if ny>=nu
   y(1:ny)=yi;
   %Em cada itera��o do loop principal nas linhas abaixo � simulado um valor pela rede.
   for t=ny+1:N
      %Abaixo, monta-se o vetor de regressores Psi que � atualizado a cada itera��o do 
      %loop principal.
      Psi(1)=1;
      for k=1:nu
         X(k)=us(t-k);
      end 
      for l=1:ny
         X(l+nu)=y(t-l);
      end
      if nurbf>0
         Xrbf=X(1);
         if nurbf>1
            for o=2:nurbf
               Xrbf=[Xrbf X(o)];
            end
         end
      end
      if nyrbf>0 & nurbf==0
         Xrbf=X(nu+1);
         if nyrbf>1
            for p=nu+2:nu+1+nyrbf
               Xrbf=[Xrbf X(p)];
            end
         end
      end
      if nyrbf>0 & nurbf>0
         Xrbf=[Xrbf X(nu+1)];
         if nyrbf>1
            for p=nu+2:nu+1+nyrbf
               Xrbf=[Xrbf X(p)];
            end
         end
      end
      if nr>0 %Nesta linha testa-se a exist�ncia dos termos n�o-lineares da forma 
              %generalizada de uma rede RBF.
         for m=1:nr
            V=Xrbf-c(m,:);
            %norma=(V*(V'))^.5;
            %if norma==0
            %   norma=.000001;
            %end;
            %w=(norma^2)*log10(norma);
            w=exp(-1*(((V*(V'))^0.5)^2)/(dp^2));
            %if m == 1
            %   w=exp(-1*(((V*(V'))^0.5)^2)/(.9^2));
				%end;
            Psi(m+1)=w;
         end
      end
      %Nas tr�s linhas abaixo � testado se a rede possui termos lineares. Se linear for 
      %igual a 1, a rede possui termos lineares, sen�o n�o possui.
      if linear==1
         Psi=[Psi X];
      end
      if const==0
         Psi(:,1)=[];
      end
      %Montado o vetor dos regressores, o valor simulado na itera��o � obtido pelo produto
      %interno do vetor Psi dos regressores pelo vetor Teta dos par�metros.
      y(t)=Psi*Teta;
      %Nas tr�s linhas abaixo s�o limpas as posi��es do vetor de regressores ocupadas pelos
      %regressores dos termos lineares. Estas s�o limpas unicamente se a rede possuir termos
      %lineares.
      if const==1 & linear==1
         Psi(nr+2:nr+1+nu+ny)=[];
      end
      if const==0 & linear==1
         Psi(nr+1:nr+nu+ny)=[];
      end
   end
%Abaixo � seguido o mesmo procedimento de simula��o, somente para o caso em que nu � maior
%que ny.
else  
   y(1:nu)=yi;
   %Loop principal come�a na linha abaixo. 
   for t=nu+1:N
      %Nas linhas abaixo monta-se o vetor de regressores Psi.
      Psi(1)=1;
      for k=1:nu
         X(k)=us(t-k);
      end 
      for l=1:ny
         X(l+nu)=y(t-l);
      end
      if nurbf>0
         Xrbf=X(1);
         if nurbf>1
            for o=2:nurbf
               Xrbf=[Xrbf X(o)];
            end
         end
      end
      if nyrbf>0 & nurbf==0
         Xrbf=X(nu+1);
         if nyrbf>1
            for p=nu+2:nu+1+nyrbf
               Xrbf=[Xrbf X(p)];
            end
         end
      end
      if nyrbf>0 & nurbf>0
         Xrbf=[Xrbf X(nu+1)];
         if nyrbf>1
            for p=nu+2:nu+1+nyrbf
               Xrbf=[Xrbf X(p)];
            end
         end
      end
      if nr>0 %Nesta linha testa-se a inclus�o dos termos n�o-lineares.
         for m=1:nr
            V=Xrbf-c(m,:);
            %if m = 1
            %   w=exp(-1*(((V*(V'))^0.5)^2)/(.2^2));
				%end;
            w=exp(-1*(((V*(V'))^0.5)^2)/(dp^2));
            %if m == 1
            %   w=exp(-1*(((V*(V'))^0.5)^2)/(.9^2));
				%end;
            Psi(m+1)=w;
         end
      end   
      %Nas tr�s linhas abaixo s�o acrescentados regressores dos termos lineares, somente
      %se a rede possuir estes termos.
      if linear==1
         Psi=[Psi X];
      end
      if const==0
         Psi(:,1)=[];
      end
      %Montado o vetor dos regressores, o valor simulado na itera��o � obtido pelo produto
      %interno do vetor Psi dos regressores pelo vetor Teta dos par�metros.
      y(t)=Psi*Teta;
      %Nas tr�s linhas abaixo s�o limpas as posi��es do vetor de regressores ocupadas pelos
      %regressores dos termos lineares. Estas s�o limpas unicamente se a rede possuir termos
      %lineares.
      if const==1 & linear==1
         Psi(nr+2:nr+1+nu+ny)=[];
      end
      if const==0 & linear==1
         Psi(nr+1:nr+nu+ny)=[];
      end
   end
end
   
