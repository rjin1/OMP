function [Xsparse,reconerr,indexpic]=OrthMatchPur (ycompress, Atran,threshold,sparsity1)
%__________________Orthogonal Matching Pursuit____________________
%1.Upadate whole sparse representation in each iteration
%2.Pick principle:unit norm inner product.
%3.Try non-sparse case, check iter, balance threshold and sparsity
%____________________2/28/2019-coded by Rui_______________________
k=size(Atran,2);
resi=ycompress;
A=normc(Atran);
indexpic=zeros(k,1);
%indexrep=zeros(k,1);%For testing repeated column uses
%indrep=0;%For testing repeated column uses
Xsparse=zeros(k,1);
%Pick most aligned atom:
for iter=1:k
y=normc(resi);
Apikcol=y'*A;
[~,indexpic(iter,1)]=max(abs(Apikcol));
%Repeated atom use situation check, maximum precision in MATLAB is E-16
induniq=unique(indexpic(1:iter,1),'stable');
if length(induniq)~=nnz(indexpic(1:iter,1))
    %indrep=indrep+1;%For testing repeated column uses
    display(strcat(num2str(indexpic(iter,1)),'th column is repeatedly picked at iter=',num2str(iter)))
    sparsity=k-iter+1;
    display(strcat('sparsity=',(num2str(sparsity))));
    break;
    %indexrep(indrep)=indexpic(iter);%For testing repeated column uses
end
%Calculate/update representation by least square:
Acal=Atran(:,indexpic(1:iter,1));
S=Acal\ycompress;
%Calculate/update residual:
resi=ycompress-Atran(:,indexpic(1:iter,1))*S;
sparsity=k-iter;
display(strcat('sparsity=',(num2str(sparsity))));
reconerr=norm(resi);
if nnz(S)==k-sparsity1
    break;
end
if reconerr<=threshold
    break;
end
end
Xsparse(induniq,1)=S(1:length(induniq),1);
