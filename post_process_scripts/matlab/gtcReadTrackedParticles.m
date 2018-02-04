function p = gtcReadTrackedParticles(trackPDir,nMPI)

trackPDir = 'trackp_dir';

% list all TRACKP files, one per MPI task
tmp = ls([trackPDir '/TRACKP.*']);
% parse string to trim spaces, and list one filename per row
fileNames = sscanf(tmp,'%s',[length([trackPDir '/TRACKP.'])+5,inf])';
% last five characters are the MPI processor number, find largest for nMPI
nMPI = max(str2num(fileNames(:,end-4:end)))+1;
% should be consistent with number of rows if no files are missing
tmp = size(fileNames,1)

if (nMPI ~= tmp)
    error('Number of files not consistent with maximum file number.  Are you missing files?')
end
clear tmp

p.nMPI = nMPI;
p.nTStep = 0;
p.nParticle = [];
p.nData = 1;

%R0=1;
%PE=32;
%ndiag=1;
% 
%for i=1:PE
%    tracking=[path,'trackp_dir/TRACKP.',num2str(i-1,'%5.5d')];
%    fid(i)=fopen(tracking);
%end

for i = 1:nMPI
    fileID = fopen(fileNames(i,:))
    istep = fgetl(fileID);
    

for i=1:PE
    istep=fgetl(fid(i));
    npp=str2num(fgetl(fid(i)));
    
    for j=1:npp
    data=str2num(fgetl(fid(i)));
    line1=str2num(data);
 %   data=fgetl(fid(i));
 %   line2=str2num(data);
 %   data=fgetl(fid(i));
 %   line3=str2num(data);
    
    x=line1(1)*R0;
    y=line1(2)*R0;
    z=line1(3);
 %   energy=line2(1);
    
    X(1)=x;
    Y(1)=y;
    Z(1)=z;
    
    
%    tag(1)=line2(4);
 %   tag(2)=line3(1);
    end
end

for i=1:PE
    istep=fgetl(fid(i));
  while istep~=-1  
    istep=str2num(istep);
    npp=fgetl(fid(i));
    npp=str2num(npp);
    
    for j=1:npp
      data=fgetl(fid(i));
      line1=str2num(data);
   %   data=fgetl(fid(i));
   %   line2=str2num(data);
   %   data=fgetl(fid(i));
   %   line3=str2num(data);
    
      x=line1(1)*R0;
      y=line1(2)*R0;
      z=line1(3);
     
    
      X(istep)=x;
      Y(istep)=y;
      Z(istep)=z;
      
      
    end
   %while  (line2(4)~=tag(1))
   %end
    
    
   
    istep=fgetl(fid(i));
   end  
  %end
end
