;pref_set,'idl_path','+${HOME}/idl/gtc:<idl_default>',/commit
  ;.r ${HOME}/idl/gtc/hplt
  ;.r ${HOME}/idl/gtc/splt
  ;.r ${HOME}/idl/gtc/shear
  ;.r ${HOME}/idl/gtc/history
  ;.r ${HOME}/idl/gtc/snap
  ;.r ${HOME}/idl/gtc/equilibrium
  ;.r ${HOME}/idl/gtc/rtime
  ;.r ${HOME}/idl/gtc/plt

  
.r ${HOME}/idl/hplt
.r ${HOME}/idl/splt 
.r ${HOME}/idl/shear
.r ${HOME}/idl/history
.r ${HOME}/idl/snapshot
.r ${HOME}/idl/equilibrium
.r ${HOME}/idl/rtime
.r ${HOME}/idl/plt
  

;  qsub -W depend=afterok:6211884 job-hopper
