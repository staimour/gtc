function A = gtcReadHist(fname)
% Reads GTC history output files.
%
% Syntax
%     A = gtcReadHist(fname)
%
% The single argument, fname, should be a file system directory path to the
% desired GTC history.out file.  If no argument is supplied, the default value
% is 'history.out' in Matlab's current directory.
% 
% The returned variable, A, is a struct containing data from the history file.
%

% default arguments
if nargin < 1
	fname = 'history.out';
end

	ihist = fopen(fname,'r');
	if ihist < 0 % error opening the file
		error(['Could not open the history file ', fname]);
	end

% # of time steps
	A.ndstep = 1;
% # of species: ion, electron, EP, impuries
	A.nspecies = 1;
% # of quantities per species: density,entropy,momentum,energy, fluxes
	A.mpdiag = 1;	
% # of field variables: phi, a_par, fluidne
	A.nfield = 1;
% # of modes per field: (n,m)
	A.modes = 1;
% # of quantities per field: rms, single spatial point, zonal components
	A.mfdiag = 1;
% # time step size
	A.tstep = 0.0;

% read the history data

	[A.ndstep tmp] = fscanf(ihist, '%f',1);
        [A.nspecies tmp] = fscanf(ihist, '%f',1);
	[A.mpdiag tmp] = fscanf(ihist, '%f',1);
	[A.nfield tmp] = fscanf(ihist, '%f',1);
	[A.modes tmp] = fscanf(ihist, '%f',1);
	[A.mfdiag tmp] = fscanf(ihist, '%f',1);
	[A.tstep tmp] = fscanf(ihist, '%f',1);


%        print, ndstep,nspecies,mpdiag,nfield,modes,mfdiag,tstep

%        print,'total time steps=',ndstep
%        print,'read # of last time step'
%        new=''
%        read,new
%        if new ne '' then ndstep=0+new

i=1;
%while ~feof(ihist)
for i=1:A.ndstep
%for i=1:76
	A.parthist(i,1:A.mpdiag,1:A.nspecies) = fscanf(ihist, '%f', [A.mpdiag,A.nspecies]);
	A.fieldhist(i,1:A.mfdiag,1:A.nfield) = fscanf(ihist, '%f', [A.mfdiag,A.nfield]);
	for j=1:A.nfield
		A.modehist(i,1:2,1:A.modes,j) = fscanf(ihist, '%f', [2,A.modes]);
	end
%	i=i+1;
end

%	ntime=ndstep
%	xtime=indgen(ntime)
%	nstart=0
%	nend=ntime-1
%	nfreq=(nend-nstart)/10


fclose(ihist);
