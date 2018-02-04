function [A,B]=fftCal(y,dx,iplot)
% plot magnitude of fft with correct horizontal 
%axis. Please use odd number of array length for
% proper nyquist number

    lsy=length(y); % length of array
    
    nyq=floor(lsy/2+1);  %nyquist number

    pow=fft(y); % take fft
    pow=abs(pow); % take modulus

    yp=zeros(lsy,1); % initialize yp
    
    yp(1:nyq)=pow(nyq:end);        %fft result is backwards
    yp(nyq:end)=pow(1:nyq);

    dk=2.*pi./(dx.*lsy);
    k=dk.*(-(nyq-1):(nyq-1));

A=yp;
B=k;

if(iplot==1)
    figure
    plot(k,yp);
end