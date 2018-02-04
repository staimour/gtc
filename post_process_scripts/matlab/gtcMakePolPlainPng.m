function gtcMakePolPlainPng



disp('determining snapshot files...')
snapfiles = dir('snap*.out');
disp(['there are ' num2str(length(snapfiles)) ' snapshot files:'])


for i=1:length(snapfiles)
  
  plotphipol(snapfiles(i).name,60)
  title('\delta\phi','FontSize',16)
  set(gca,'fontsize',16)
  xlabel('R(m)','fontsize',16)
  ylabel('Z(m)','fontsize',16)
  
  print(gcf,char(strcat('snap',regexp(snapfiles(i).name,'\d+','match'),'.png')),'-dpng')
  disp(strcat(regexp(snapfiles(i).name,'\d+','match'),' saved'))
  
end