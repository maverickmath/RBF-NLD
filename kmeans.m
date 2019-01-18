function [cn,vc] = kmeans(k,ci,x)
% Funcao que calcula os k centros das classes de dados x, utilizando
% o algoritmo Kmeans.
% Entrada:
%      k  = número de classes a ser classificado o conjunto de dados
%      ci = matriz com k centros iniciais, centros em colunas
%      x = matriz do conjunto de dados, amostras em colunas
% Retorno:
%      cn = matriz de centros das k classes, em colunas
%      vc = vetor de classificação de cada amostra em k classes

troca = 1;
% Repete a iteracao
while troca==1,

	troca = 0;
	% Calculo das distancias de cada x aos centros
	for i=1:(size(x,2)),
		for j=1:k,
			d(j,i) = norm( x(:,i) - ci(:,j) );
		end;
	end;

	% Encontra o centro de menor distancia a x e associa a ele a classe
	[dm, vc] = min(d);

	% Calcula novos valores para os centros pela media dos pontos x
	% da mesma classe

	for j=1:k,
		ps(j) = 0;               % armazena numero de pontos encontrados
		cn(:,j) = zeros(size(x,1),1);      % limpa matriz de centros
		for i=1:(size(x,2)),    % varre todas as amostras de x
			if vc(i) == j    % procura a amostra de classe j
				% Calcula a soma dos x de mesma classe
				cn(:,j) = cn(:,j) + x(:,i);
				% Calcula a quantidade de x de mesma classe
				ps(j) = ps(j) + 1;
			end;
		end;
		cn(:,j) = cn(:,j)./ps(j);  % Calcula a nova media
		if cn(:,j) ~= ci(:,j)       % Compara com a media anterior
			ci(:,j) = cn(:,j);
			troca = 1;
		end;
	end;


end;