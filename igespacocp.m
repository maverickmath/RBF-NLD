%Define centros igualmente espaçados [no tempo] de acordo com o espaço de entrada de modelos RBF

function c=igespacocp(nr, nu, ny, Nt, ut, yt)
if nu>=ny
   cont=1;
   for i=nu+1:fix((Nt-nu)/nr)+1:Nt
      if nu>=1
        for j=1:nu
           c(cont,j)=ut(i-j);
        end
      end
      if ny>=1
        for k=1:ny
            c(cont,nu+k)=yt(i-k);
        end
      end
        cont=cont+1;
   end
else
   cont=1;
   for i=ny+1:fix((Nt-ny)/nr)+1:Nt
      if nu>=1
        for j=1:nu
           c(cont,j)=ut(i-j);
        end
      end
      if ny>=1
        for k=1:ny
           c(cont,nu+k)=yt(i-k);
        end
      end
      cont=cont+1;
   end
end