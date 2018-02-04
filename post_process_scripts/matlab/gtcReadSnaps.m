function A = gtcReadSnaps



disp('determining snapshot files...')
snapfiles = dir('snap*.out');
disp(['there are ' num2str(length(snapfiles)) ' snapshot files:'])


for i=1:length(snapfiles)
  sn.m(i)=gtcReadSnap(snapfiles(i).name);  
end

A=sn;