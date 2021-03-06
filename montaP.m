function P=montaP(ut,yt,nu,ny,nul,nyl,linear,C,Nt,nr,dp)
%Esta funcao monta a matriz de regressores para redes RBF a partir dos vetores de entrada, ut,
%e de saida, yt. Os atrasos passados em nu e ny definem a dimensao do vetor de entrada aas 
%funcoes de base radial. Se a rede tambem for composta por termos lineares, o valor passado
%em linear deve ser igual a 1, senao deve ser colocado qualquer outro valor. Os atrasos dos 
%termos lineares sao passados em nul e nyl, para entrada e saida, respectivamente. Em C sao 
%passados os centros, sendo que cada centro deve ser uma linha na matriz C, que deve possuir 
%o numero de colunas igual a nu+ny. Nt e o tamanho de ut ou yt e nr e o numero de linhas de C.
%O retorno da funcao e, portanto, a matriz de regressores, que posteriormente pode ser passada 
%a um algoritmo como o myhouse, que calcula o ERR das colunas da matriz de regressores. Entao
%pode-se estimar os pesos da rede a partir do algoritmo dos Minimos Quadrados. Inicialmente
%esta funcao utiliza a funcao thin-plate spline como funcao de base, mas a funcao de base a 
%ser utilizada pode ser facilmente configurada (na verdade ja esta implementada, mas em 
%comentarios)

atr=max([nu ny nul nyl]);

P=ones(Nt-atr,1);	%Modelo sem o termo constante
%P = [];
X=[];
Xl=[];
if nu>=1
   for l=1:nu
      X=[X ut(atr+1-l:Nt-l)];   
   end   
end
if ny>=1
   for m=1:ny
      X=[X yt(atr+1-m:Nt-m)];
   end   
end
if nul>=1
   for ll=1:nul
      Xl=[Xl ut(atr+1-ll:Nt-ll)];   
   end   
end
if nyl>=1
   for ml=1:nyl
      Xl=[Xl yt(atr+1-ml:Nt-ml)];   
   end   
end
if nr>0
   for n=1:Nt-atr
      for o=1:nr
         V(n,o)=((X(n,:)-C(o,:))*((X(n,:)-C(o,:))'))^0.5;
         %A linha de baixo testa se a distancia euclidiana entre a entrada e um centro e nula
         %se isto for verdade, ela substitui o valor calculado acima por um valor muito pequeno,
         %para que seja viavel o calculo de seu logaritmo. Este artificio nao influi na 
         %identificacaco de sistemas, pois se utiliza de metodos estatisticos.
         %if V(n,o)==0
         %   V(n,o)=.0000001;   
         %end   
      end   
   end
   PP=(-1*V.*V)*(1/(dp^2));
   Fi=exp(PP);
   %Fi=V.*V.*(log10(V));
   %Fi=V.^3;
   %Fi=(1-exp(-V))./(1+exp(-V));   
   P=[P Fi];   
end

if linear==1
   P=[P Xl];   
end


