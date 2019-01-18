%Simula��o livre de um modelo RBF com estrutura selecionada segundo o crit�rio de taxa de redu��o de erro

function y=simulacao_errn(c,nr,nu,ny,N,us,yi,Teta,linear,dp,const,u1,u2,y1,y2,y3)

%Primeiramente � testado qual dos m�ximos atrasos � maior para n�o extrapolar os limites
%dos vetores de entrada e/ou sa�da.
if ny>=nu
   y(1:ny)=yi;
   %Em cada itera��o do loop principal nas linhas abaixo � simulado um valor pela rede.
   for t=ny+1:N
      %Abaixo, monta-se o vetor de regressores Psi que � atualizado a cada itera��o do 
      %loop principal.
      if const == 1,
          Psiu(1)=1; %Termo constante
      end;
      %Psiu
      %pause
      for k=1:nu
         X(k)=us(t-k);  %X -> vetor de entrada para os termos nao-lineares
      end 
      for l=1:ny
         X(l+nu)=y(t-l);
      end
      %X
      %pause
      if nr>0 %Nesta linha testa-se a exist�ncia dos termos n�o-lineares da forma 
              %generalizada de uma rede RBF.
         for m=1:nr
            V=X(nu+1:end)-c(m,:);
            %V=((X(nu+1:end)-c(m,:))*((X(nu+1:end)-c(m,:))'))^0.5;
            %w=((V.*V)+(dp^2)).^0.5;	%multiquadr�tica
            w=exp(-1*(((V*(V'))^0.5)^2)/(dp^2));
            %if V==0
            %    V=.0000001;   
            %end
            %w=V.*V.*(log10(V));  %Fun�ao: thin-plate spline
            %w=(1-exp(-V))./(1+exp(-V));
            %pause
            %w=V.^3;
            if const == 1,
                Psiu(m+1)=w;
            else Psiu(m)=w;
            end;
        end
      end
      %whos Psiu
      %pause
      %Nas tr�s linhas abaixo � testado se a rede possui termos lineares. Se linear for 
      %igual a 1, a rede possui termos lineares, sen�o n�o possui.
      %Psiu
      %pause
      if linear==1
          if u1 == 1
              Psiu=[Psiu X(1)];
              %Psiu
              %X(1)
              %pause
          end;
          if u2 == 1
              Psiu=[Psiu X(2)];
          end;
          %if u3 == 1
          %    Psiu=[Psiu X(3)];
          %end;
          %X
          %y1
          %pause;
          if y1 == 1
              Psiu=[Psiu X(1)];
              %Psiu
              %X(3)
              %pause
          end;
          if y2 == 1
              Psiu=[Psiu X(2)];
              %Psiu
              %X(4)
              %pause
          end;
          if y3 == 1
             Psiu=[Psiu X(3)];
          end;
      end
      %Montado o vetor dos regressores, o valor simulado na itera��o � obtido pelo produto
      %interno do vetor Psi dos regressores pelo vetor Teta dos par�metros.
      %length(Psiu)
      %Psiu
      %pause
      %length(Teta)
%       if (t==478)
%           Psi
%           pause
%       end;
%whos Psiu Teta;
y(t)=Psiu*Teta;
%pause
      
      %----------------------------------------------------------------------
      %Nas tr�s linhas abaixo s�o limpas as posi��es do vetor de regressores ocupadas pelos
      %regressores dos termos lineares. Estas s�o limpas unicamente se a rede possuir termos
      %lineares.
      if linear==1
         Psiu(nr+const+1:length(Psiu))=[];
         %tamanho=length(Psiu)
      end
   end
end
   
   %----------------------------------------------------------------------
%Abaixo � seguido o mesmo procedimento de simula��o, somente para o caso em que nu � maior
%que ny.
% else  
%    y(1:nu)=yi;
   %Loop principal come�a na linha abaixo. 
%    for t=nu+1:N
      %Nas linhas abaixo monta-se o vetor de regressores Psi.
%       Psi(1)=1;
%       for k=1:nu
%          X(k)=us(t-k);
%       end 
%       for l=1:ny
%          X(l+nu)=y(t-l);
%       end
%       if nr>0 %Nesta linha testa-se a inclus�o dos termos n�o-lineares.
%          for m=1:nr
%             V=X-c(m,:);
%             w=exp(-1*(((V*(V'))^0.5)^2)/(dp^2));
%             Psi(m+1)=w;
%          end
%       end   
      %Nas tr�s linhas abaixo s�o acrescentados regressores dos termos lineares, somente
      %se a rede possuir estes termos.
%       if linear==1
%          Psi=[Psi X];
%       end
      %Montado o vetor dos regressores, o valor simulado na itera��o � obtido pelo produto
      %interno do vetor Psi dos regressores pelo vetor Teta dos par�metros.
%      y(t)=Psi*Teta;
      %Nas tr�s linhas abaixo s�o limpas as posi��es do vetor de regressores ocupadas pelos
      %regressores dos termos lineares. Estas s�o limpas unicamente se a rede possuir termos
      %lineares.
%       if linear==1
%          Psi(nr+2:nr+1+nu+ny)=[];
%       end
%    end
% end
%    
