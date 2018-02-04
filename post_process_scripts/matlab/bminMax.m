mthetamax=1500;
b=sp3.bcn(:,:,1);


minbfield=100;
maxbfield=0;
nb=mthetamax/2;
%deltat=1/(length(b(1,:,1))-1);

for ii=300:10:400
    
pdum=ii;
    
    for j=1:length(b(ii,:))
        tdum=j;
        bdum=b(pdum,tdum);
        if(bdum<minbfield)
            minbfield=bdum;
            thetabmin=tdum;
        end
         if(bdum>maxbfield)
            maxbfield=bdum;
            thetabmax=tdum;
         end
    end
    
    deltab=(maxbfield-minbfield)/(nb-1);
    for i=1:nb
        testb=minbfield+deltab*(i-1);
        
        bdum=minbfield;
        j=1;
        tdum=thetabmin;
        while(bdum<testb)
            j=j+1;
            tdum=floor(mod(thetabmin+j,299));
            bdum=b(pdum,tdum);
        end
        thetadwn(i)=tdum;
        
        bdum=minbfield;
        j=1;
        tdum=thetabmin;
        while(bdum<testb)
            j=j+1;
            tdum=floor(mod(299+thetabmin+j,299));
            bdum=b(pdum,tdum);
        end
        thetaupp(i)=tdum;
    end
    thetabmax=floor(mod(299+thetabmax-thetabmin,299));


for i=1:299
  tdum=i;  
  bdum=b(ii,tdum);
  bind=floor((bdum-minbfield)/deltab)+1;
  tdum=floor(mod(299+tdum-thetabmin,299));
  if(tdum<thetabmax)
      tnew=thetaupp(bind);
  else
      tnew=thetadwn(bind);
  end
  db(i)=100.*abs(b(ii,tnew)-bdum)./bdum;
end

th=linspace(0,2*pi,299);
figure
plot(th,db,'dk','linewidth',2.5,'markerfacecolor','auto')
grid on 
xlim([0 2*pi])
set(gca,'linewidth',2.5)
set(gca,'fontsize',25)
xlabel('\theta')
ylabel('percent difference')
title(['Change in B, ' num2str(ii)])

end
