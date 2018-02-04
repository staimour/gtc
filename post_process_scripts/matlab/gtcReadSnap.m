function A = gtcReadSnap(snapshotFilename)
% Reads GTC snapshot output files.
%
% Syntax
%     A = gtcReadSnap(snapshotFilename)
%
% The single argument, snapshotFilename, should be a file system directory path to the
% desired GTC snap00000.out file.  If no argument is supplied, the default value
% is 'snap00005.out' in Matlab's current directory.
% 
% The returned variable, A, is a struct containing data from the snapshot file.
%

% default arguments
	if nargin < 1
		snapshotFilename = 'snap00005.out';
	end

% Initialize some default values before reading in
	A.nspecies = 9;		% # of species: ion, electron, EP, impurities
	A.nfield = 1;		% # of field variables: phi, a_para, fluidne
	A.nvgrid = 1;		% # of grids in energy and pitch
	A.mpsi = 1;			% # of radial grids: mpsi+1
	A.mtgrid = 1;		% # of poloidal grids
	A.mtoroidal = 1;	% # of toroidal (parallel) grids
	A.tmax = 1.0;		% upper bound of temperature

% open file for reading
	fid = fopen(snapshotFilename,'r');
	if fid < 0 % error opening the file
		error(['Could not open the snapshot file ', snapshotFilename]);
	end

% read in snapshot data

	% scalar parameters
	A.nspecies = fscanf(fid, '%f',1);
	A.nfield = fscanf(fid, '%f',1);
	A.nvgrid = fscanf(fid, '%f',1);
	A.mpsi = fscanf(fid, '%f',1);
	A.mtgrid = fscanf(fid, '%f',1);
	A.mtoroidal = fscanf(fid, '%f',1);
	A.tmax = fscanf(fid, '%f',1);

	% read in grid data
	for i = 1:A.nspecies
		tmpprofile(:,:,i) = fscanf(fid, '%f', [A.mpsi 6]);
	end
	
	for i = 1:A.nspecies
		tmppdf(:,:,i) = fscanf(fid, '%f', [A.nvgrid 4]);
	end		
	
	for i = 1:A.nfield+2 %extra +1 to read apara field
		tmppoloidata(:,:,i) = fscanf(fid, '%f', [A.mtgrid A.mpsi]);
	end	
	
	for i = 1:A.nfield
		tmpfluxdata(:,:,i) = fscanf(fid, '%f', [A.mtgrid A.mtoroidal]);
	end	

% close snapshot file.
	fclose(fid);

% store grid data in human readable struct fields

	% profiles
		A.EnPDF_ind = 1:A.tmax/(A.nvgrid-1.0):A.tmax;
	
		% ions
		A.iDens_fullf = tmpprofile(:,1,1);
	    A.iDens_delf = tmpprofile(:,2,1);
	    A.iFlow_fullf = tmpprofile(:,3,1);
		A.iFlow_delf = tmpprofile(:,4,1);
		A.iEn_fullf = tmpprofile(:,5,1);
		A.iEn_delf = tmpprofile(:,6,1);

		A.iEnPDF_fullf = tmppdf(:,1,1);
		A.iEnPDF_delf = tmppdf(:,2,1);

		A.iPitchPDF_fullf = tmppdf(:,3,1);
		A.iPitchPDF_delf = tmppdf(:,4,1);



		% electrons
        if A.nspecies > 1
	    	A.eDens_fullf = tmpprofile(:,1,2);
	    	A.eDens_delf = tmpprofile(:,2,2);
	    	A.eFlow_fullf = tmpprofile(:,3,2);
	    	A.eFlow_delf = tmpprofile(:,4,2);
	    	A.eEn_fullf = tmpprofile(:,5,2);
	    	A.eEn_delf = tmpprofile(:,6,2);
    
	    	A.eEnPDF_fullf = tmppdf(:,1,2);
	    	A.eEnPDF_delf = tmppdf(:,2,2);
    
    		A.ePitchPDF_fullf = tmppdf(:,3,2);
		    A.ePitchPDF_delf = tmppdf(:,4,2);
        end

		% energetic pArticles
		if A.nspecies > 2
			A.epDens_fullf = tmpprofile(:,1,3);
			A.epDens_delf = tmpprofile(:,2,3);
			A.epFlow_fullf = tmpprofile(:,3,3);
			A.epFlow_delf = tmpprofile(:,4,3);
			A.epEn_fullf = tmpprofile(:,5,3);
			A.epEn_delf = tmpprofile(:,6,3);

			A.epEnPDF_fullf = tmppdf(:,1,3);
			A.epEnPDF_delf = tmppdf(:,2,3);

			A.epPitchPDF_fullf = tmppdf(:,3,3);
			A.epPitchPDF_delf = tmppdf(:,4,3);
		end

    A.phifluxsurf = tmpfluxdata(:,:,1);
    A.x = tmppoloidata(:,:,A.nfield+1);
    A.y = tmppoloidata(:,:,A.nfield+2);
    A.phipolplane = tmppoloidata(:,:,1);
    A.aparapolplane = tmppoloidata(:,:,2);  
    A.denepolplane = tmppoloidata(:,:,3);  



	%{
% flux surfaces

 15: begin		;%phi on flux surface
	f=fluxdata(*,*,0)
	flux,f
 end

 16: begin		;poloidal and parallel spectra                
	f=fluxdata(*,*,0)
	spectrum,f
 end

 17: begin		;phi on ploidal plane
	x=poloidata(*,*,nfield)
	y=poloidata(*,*,nfield+1)
	f=poloidata(*,*,0)
	poloidal,x,y,f
 end	

 18: begin		;radius profile of field and rms
	f=poloidata(*,*,0)
	cut1d,f,1
 end

 19: begin		;poloidal profile of field and rms
	f=poloidata(*,*,0)
	cut1d,f,2
 end

 20: begin		;a_para on flux surface
	f=fluxdata(*,*,1)
	flux,f
 end

 21: begin		;poloidal and parallel spectra                
	f=fluxdata(*,*,1)
	spectrum,f
 end

 22: begin		;a_para on ploidal plane
	x=poloidata(*,*,nfield)
	y=poloidata(*,*,nfield+1)
	f=poloidata(*,*,1)
	poloidal,x,y,f
 end	

 23: begin		;radius profile of field and rms
	f=poloidata(*,*,1)
	cut1d,f,1
 end

 24: begin		;poloidal profile of field and rms
	f=poloidata(*,*,1)
	cut1d,f,2
 end

 25: begin		;fluidne on flux surface
	f=fluxdata(*,*,2)
	flux,f
 end

 26: begin		;poloidal and parallel spectra                
	f=fluxdata(*,*,2)
	spectrum,f
 end

 27: begin		;fluidne on ploidal plane
	x=poloidata(*,*,nfield)
	y=poloidata(*,*,nfield+1)
	f=poloidata(*,*,2)
	poloidal,x,y,f
 end	

 28: begin		;radius profile of field and rms
	f=poloidata(*,*,2)
	cut1d,f,1
 end

 29: begin		;poloidal profile of field and rms
%	f=poloidata(*,*,2)
%	cut1d,f,2
% end

%}
