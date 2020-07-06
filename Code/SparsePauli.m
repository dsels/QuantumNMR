function [Sxloc,Syloc,Szloc]=SparsePauli(L)

N=2^L; %Hilbert space dimension
Mbin=dec2bin(0:(2^L)-1);
Mdec=bin2dec(Mbin)+1;

Szloc=cell(1,L);Sxloc=cell(1,L);Syloc=cell(1,L);

    
for i=1:L
   %Generate sz_i
   f=Mbin(:,i)-'0';%Turn it into a double
   Z=2*f-1;
   Szloc{i}=sparse(1:N,1:N,Z);
   %Generate sx_i/sy_i
   
   %Flip spin
   flip=num2str(mod(f+1,2));
   %Insert in configuration
   Mf=Mbin;Mf(:,i)=flip;
   %Convert to dec
   Mfdec=bin2dec(Mf)+1;
   %Make sparse matrix that couples Mdec to Mfdec
   Sxloc{i}=sparse(Mdec,Mfdec,1);
   Syloc{i}=sparse(Mdec,Mfdec,1i*(1-2*f));
end


end