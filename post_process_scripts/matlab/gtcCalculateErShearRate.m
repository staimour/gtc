function A = gtcCalculateErShearRate(pfile,spdata,inputOption)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if inputOption==1
  p=gtcReadProfile(pfile);
  er=p.Er;
  psip=linspace(0,1,length(p.x));
elseif inputOption==0
    er=1000*pfile;
    psip=linspace(0,1,length(er))';
end
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sp=spdata;%readVmecSpdat(spdata);
r=sp.xcn(1,1,1).*sp.rpsi';
lsp=sp.lsp;
omegaP=9.57883357*(10^7)*sp.bcn(1,1,1);
b=sp.bcn(1,1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi=linspace(0,1,lsp)';
psint=linspace(1/(2*(lsp-1)),1-1/(2*(lsp-1)),lsp-1)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
er=interp1(psip,er,psi);
%er=gauss(er,11,311,393); vmec nonlinear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dpsidr=diff(psi)./diff(r);
dedp=diff(er)./diff(psi);

r=interp1(psi,r,psint);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma = r.*dpsidr.*(dedp./r)./b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(psint,gamma./omegaP)
grid on
title('|\gamma_{E\times{B}}|')
xlabel('\psi_N')
ylabel('|\gamma_{E\times{B}}|/\Omega_{cp}')
set(gca,'fontsize',20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(psint,gamma)
grid on
title('|\gamma_{E\times{B}}|')
xlabel('\psi_N')
ylabel('1/s')
set(gca,'fontsize',20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=gamma./omegaP;
B=dedp;
