function [w,A]=NMRSpectrum(k)
%Compute the NMR spectrum for molecule number k in the dataset
%Input: k, the label of the molecule 
%Output: -frequency w and spectral weight A

%%%%%%%%%%%%%%%%%%%%%%%
%Get quantum operators%
%%%%%%%%%%%%%%%%%%%%%%%
N=4;DH=2^N;% Number of spins/Hilbert space dim
[Sx,Sy,Sz]=SparsePauli(N); %Pauli matrices

%%%%%%%%%%%%%%%%%%
%Load spin matrix%
%%%%%%%%%%%%%%%%%%
name=strcat('matrix',num2str(k),'.csv'); 
J=load(name); %Load spin matrix number k
h=diag(J);h=h-mean(h); %Get chemical shifts and remove the mean
wext=40;% External magnetic field strength in Mhz
h=h*wext;%Convert from ppm to Mhz
J=J-diag(diag(J)); %Remove diagonal

%%%%%%%%%%%%%%%%%%%%%%
%Generate Hamiltonian%
%%%%%%%%%%%%%%%%%%%%%%
H=zeros(DH,DH);n=H;X=H;
for i=1:N
    if i~=N
    for j=(i+1):N
        H=H+J(i,j)/4*(Sz{i}*Sz{j}+Sx{i}*Sx{j}+Sy{i}*Sy{j});
    end
    end
    n=n+(h(i)+eps)*Sz{i}/2;
    X=X+Sx{i}/2;
end
H=H+n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diagonalize and construct spectrum%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V,D]=eig(H);
xe=V'*X*V;

E=diag(D);
de=E*ones(1,length(D))-ones(length(D),1)*E';
xr=[];dr=[];
for i=1:(2^N-1)
    xr=[xr diag(xe,i)' diag(xe,-i)'];
    dr=[dr diag(de,i)' diag(de,-i)'];
end
BW=max(dr);%Bandwidth

gamma=1; %Decoherence rate in Hz
w=linspace(-1.2*BW,1.2*BW,10000);
f=gamma./(gamma^2+(w'*ones(1,length(dr))-ones(length(w),1)*dr).^2)/pi/DH;
A=f*xr.^2';%Spectral function 

