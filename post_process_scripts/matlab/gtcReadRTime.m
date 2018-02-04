function A = gtcReadRTime(fname)
% arguments
% fname: name of data1d.out file to open and read
%

% default arguments
if nargin < 1
	fname = 'data1d.out';
end

	fid = fopen(fname,'r');
	if fid < 0 % error opening the file
		error(['Could not open the rtime file ', fname]);
	end

% # of time steps
	A.ndstep = 1;
% # of radial grids: mpsi + 1
    A.mpsi = 1;
% # of species: ion, electron, EP, impuries
	A.nspecies = 1;
% whether or not electron is loaded
	A.nhybrid = 0;
% # of quantities per species: particle flux, energy flux, ...
	A.mpdata1d = 1;	
% # of field variables: phi, a_par, fluidne
	A.nfield = 1;
% # of quantities per field: phi00, phi_rms, ...
	A.mfdata1d = 1;


% read the rtime data

	[A.ndstep tmp] = fscanf(fid, '%f',1);
    [A.mpsi tmp] = fscanf(fid, '%f',1);
	[A.nspecies tmp] = fscanf(fid, '%f',1);
	[A.nhybrid tmp] = fscanf(fid, '%f',1);
	[A.mpdata1d tmp] = fscanf(fid, '%f',1);
	[A.nfield tmp] = fscanf(fid, '%f',1);
	[A.mfdata1d tmp] = fscanf(fid, '%f',1);
	A

i=1;
%while ~feof(fid)
for i=1:A.ndstep

	A.data1di(i,1:A.mpsi,1:A.mpdata1d) = fscanf(fid, '%f', [A.mpsi,A.mpdata1d]);
	
	if A.nspecies == 2
		if A.nhybrid > 0 
			A.data1de(i,1:A.mpsi,1:A.mpdata1d) = fscanf(fid, '%f', [A.mpsi,A.mpdata1d]);
	    end
	    if A.nhybrid == 0
            A.data1df(i,1:A.mpsi,1:A.mpdata1d) = fscanf(fid, '%f', [A.mpsi,A.mpdata1d]);
		end	
	end

	if A.nspecies == 3
		A.data1df(i,1:A.mpsi,1:A.mpdata1d) = fscanf(fid, '%f', [A.mpsi,A.mpdata1d]);
		A.data1de(i,1:A.mpsi,1:A.mpdata1d) = fscanf(fid, '%f', [A.mpsi,A.mpdata1d]);
	end
	
	A.field00(i,1:A.mpsi,1:A.nfield) = fscanf(fid, '%f', [A.mpsi,A.nfield]);
	A.fieldrms(i,1:A.mpsi,1:A.nfield) = fscanf(fid, '%f', [A.mpsi,A.nfield]);

	i=i+1;
end

fclose(fid);

clear tmp
