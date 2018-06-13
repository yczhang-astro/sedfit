
pro bestmodel,name,THRESOPTION=thresop,NMODEL=nmodel,NFULLMODEL=nmodel_full,CHISQ_THRES=chisq_thres

;output selected best models.

;INPUTS:
;   name: used in input and output files
;
;OPTIONAL INPUTS:
;   thresoption: 1 to output certain number of models, 2 to output
;                models with chisq within certain range. (default 1)
;   nmodel: number of output models with different (Mc, Sigma, ms)
;           (default 10)
;   nfullmodel: number of output models with different (Mc, Sigma,
;                ms, inc, dis) (default 20)
;   chisq_thres: chisq threshold for output models, relative to the
;                minimum chisq. (default 3)

  print,''
  print,'output best models'
  print,''

;;===============================================================
;; check the inputs

  ;; check the name input
  if n_elements(name) eq 0 or n_elements(name) gt 1 then begin
     message,'input name wrong',/con
     stop
  endif 

  ;; check thresoption input
  if n_elements(thresop) eq 0 then thresop=1
  if n_elements(thresop) gt 1 then begin
     message,'thresop must be 1 or 2',/con
     stop     
  endif
  if thresop ne 1 and thresop ne 2 then begin
     message,'thresop must be 1 or 2',/con
     stop     
  endif

  ;; check nmodel and nmodel_ful inputs
  if n_elements(nmodel) eq 0 or n_elements(nmodel) gt 1 then nmodel=10
  if n_elements(nmodel_full) eq 0 or n_elements(nmodel_full) gt 1 then nmodel_full=20
  if nmodel le 0 then nmodel=10
  if nmodel_full le 0 then nmodel_full=20

  ;; check chisq_thres
  if n_elements(chisq_thres) eq 0 or n_elements(nmodel_full) gt 1 then chisq_thres=3.
  if chisq_thres le 0. then chisq_thres=3.

;;===============================================================




  
;;===============================================================
;; output best models

  dirout=name

  filein=dirout+'/'+name+'.output.allmodel.parameter.dat'
  filein_full=dirout+'/'+name+'.output.allmodel.full.parameter.dat'
  
  fileout=dirout+'/'+name+'.output.bestmodel.parameter.dat'
  fileout1=dirout+'/'+name+'.output.bestmodel.parameter_range.dat'

  fileout_full=dirout+'/'+name+'.output.bestmodel.full.parameter.dat'
  fileout1_full=dirout+'/'+name+'.output.bestmodel.full.parameter_range.dat'
  
  readcol,filein,modelno,sedno,chisq,chisq_nonlimit,mcore,sigma,mstar,inc,dis,av,rcore,menv,theta,rstar,lstar,tstar,mdisk,rdisk,mdotd,ltot,ltot_inc,ltot_av,tnow,$
          format='(I,A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)',skipline=1,/silent
  nmodel_all=n_elements(chisq)
  chisq_min=chisq(0)
  openw,1,fileout
  printf,1,'model No., SED No., chisq, chisq(non-limit), mc (msun), Sigma (g/cm^2), ms (msun), mu (deg), dis (pc), av, rcore (AU), menv (msun), theta_w (deg), rstar (rsun), lstar (lsun), tstar (K), mdisk (msun), rdisk (AU), mdotd (msun/yr), ltot (Lsun), ltot_inc (Lsun), ltot_av (Lsun), tnow (yr)'
  imodel=1
  for imodel=1l,nmodel_all do begin
     if thresop eq 1 and imodel gt nmodel then break
     if thresop eq 2 and chisq(imodel-1) gt chisq_min then break
     printf,1,strtrim(string(modelno(imodel-1)),2)+' '+sedno(imodel-1),chisq(imodel-1),chisq_nonlimit(imodel-1),mcore(imodel-1),sigma(imodel-1),$
            mstar(imodel-1),inc(imodel-1),dis(imodel-1),av(imodel-1),rcore(imodel-1),menv(imodel-1),theta(imodel-1),rstar(imodel-1),$
            lstar(imodel-1),tstar(imodel-1),mdisk(imodel-1),rdisk(imodel-1),mdotd(imodel-1),ltot(imodel-1),ltot_inc(imodel-1),$
            ltot_av(imodel-1),tnow(imodel-1),FORMAT='(A,21e14.6)'
  endfor
  close,1

  readcol,filein_full,modelno,sedno,chisq,chisq_nonlimit,mcore,sigma,mstar,inc,dis,av,rcore,menv,theta,rstar,lstar,tstar,mdisk,rdisk,mdotd,ltot,ltot_inc,ltot_av,tnow,$
          format='(I,A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)',skipline=1,/silent
  nmodel_all=n_elements(chisq)
  chisq_min=chisq(0)
  openw,1,fileout_full
  printf,1,'model No., SED No., chisq, chisq(non-limit), mc (msun), Sigma (g/cm^2), ms (msun), mu (deg), dis (pc), av, rcore (AU), menv (msun), theta_w (deg), rstar (rsun), lstar (lsun), tstar (K), mdisk (msun), rdisk (AU), mdotd (msun/yr), ltot (Lsun), ltot_inc (Lsun), ltot_av (Lsun), tnow (yr)'
  imodel=1
  for imodel=1l,nmodel_all do begin
     if thresop eq 1 and imodel gt nmodel_full then break
     if thresop eq 2 and chisq(imodel-1) gt chisq_min then break
     printf,1,strtrim(string(modelno(imodel-1)),2)+' '+sedno(imodel-1),chisq(imodel-1),chisq_nonlimit(imodel-1),mcore(imodel-1),sigma(imodel-1),$
            mstar(imodel-1),inc(imodel-1),dis(imodel-1),av(imodel-1),rcore(imodel-1),menv(imodel-1),theta(imodel-1),rstar(imodel-1),$
            lstar(imodel-1),tstar(imodel-1),mdisk(imodel-1),rdisk(imodel-1),mdotd(imodel-1),ltot(imodel-1),ltot_inc(imodel-1),$
            ltot_av(imodel-1),tnow(imodel-1),FORMAT='(A,21e14.6)'
  endfor
  close,1

;;===============================================================




  
;;===============================================================
;; list the range of parameters

  for ifile=1,2 do begin
     if ifile eq 1 then begin
        filei=fileout
        fileo=fileout1
     endif else begin
        filei=fileout_full
        fileo=fileout1_full
     endelse
     
     readcol,filei,modelno,sedno,chisq,chisq_nonlimit,mcore,sigma,mstar,inc,dis,av,rcore,menv,theta,rstar,lstar,tstar,mdisk,rdisk,mdotd,ltot,ltot_inc,ltot_av,tnow,$
             format='(I,A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)',skipline=1,/silent
     
     chisq_best=chisq(0)
     chisq_nonlimit_best=chisq_nonlimit(0)
     mcore_best=mcore(0)
     sigma_best=sigma(0)
     mstar_best=mstar(0)
     inc_best=inc(0)
     dis_best=dis(0)
     av_best=av(0)
     rcore_best=rcore(0)
     menv_best=menv(0)
     theta_best=theta(0)
     rstar_best=rstar(0)
     lstar_best=lstar(0)
     tstar_best=tstar(0)
     mdisk_best=mdisk(0)
     rdisk_best=rdisk(0)
     mdotd_best=mdotd(0)
     ltot_best=ltot(0)
     ltot_inc_best=ltot_inc(0)
     ltot_av_best=ltot_av(0)
     tnow_best=tnow(0)
     
     chisq_min=min(chisq)
     chisq_nonlimit_min=min(chisq_nonlimit)
     mcore_min=min(mcore)
     sigma_min=min(sigma)
     mstar_min=min(mstar)
     inc_min=min(inc)
     dis_min=min(dis)
     av_min=min(av)
     rcore_min=min(rcore)
     menv_min=min(menv)
     theta_min=min(theta)
     rstar_min=min(rstar)
     lstar_min=min(lstar)
     tstar_min=min(tstar)
     mdisk_min=min(mdisk)
     rdisk_min=min(rdisk)
     mdotd_min=min(mdotd)
     ltot_min=min(ltot)
     ltot_inc_min=min(ltot_inc)
     ltot_av_min=min(ltot_av)
     tnow_min=min(tnow)
     
     chisq_max=max(chisq)
     chisq_nonlimit_max=max(chisq_nonlimit)
     mcore_max=max(mcore)
     sigma_max=max(sigma)
     mstar_max=max(mstar)
     inc_max=max(inc)
     dis_max=max(dis)
     av_max=max(av)
     rcore_max=max(rcore)
     menv_max=max(menv)
     theta_max=max(theta)
     rstar_max=max(rstar)
     lstar_max=max(lstar)
     tstar_max=max(tstar)
     mdisk_max=max(mdisk)
     rdisk_max=max(rdisk)
     mdotd_max=max(mdotd)
     ltot_max=max(ltot)
     ltot_inc_max=max(ltot_inc)
     ltot_av_max=max(ltot_av)
     tnow_max=max(tnow)
     
     close,1
     openw,1,fileo
     printf,1,'chisq, chisq_nonlimit, mc (msun), Sigma (g/cm^2), ms (msun), mu (deg), dis (pc), av, rcore (AU), menv (msun), theta_w (deg), rstar (rsun), lstar (lsun), tstar (K), mdisk (msun), rdisk (AU), mdotd (msun/yr), ltot (Lsun), ltot_inc (Lsun), ltot_av (Lsun), tnow (yr)'
     printf,1,chisq_best,chisq_nonlimit_best,mcore_best,sigma_best,mstar_best,inc_best,dis_best,av_best,rcore_best,menv_best,theta_best,rstar_best,lstar_best,tstar_best,mdisk_best,rdisk_best,mdotd_best,ltot_best,ltot_inc_best,ltot_av_best,tnow_best,FORMAT='(21e14.6)'
     printf,1,chisq_min,chisq_nonlimit_min,mcore_min,sigma_min,mstar_min,inc_min,dis_min,av_min,rcore_min,menv_min,theta_min,rstar_min,lstar_min,tstar_min,mdisk_min,rdisk_min,mdotd_min,ltot_min,ltot_inc_min,ltot_av_min,tnow_min,FORMAT='(21e14.6)'
     printf,1,chisq_max,chisq_nonlimit_max,mcore_max,sigma_max,mstar_max,inc_max,dis_max,av_max,rcore_max,menv_max,theta_max,rstar_max,lstar_max,tstar_max,mdisk_max,rdisk_max,mdotd_max,ltot_max,ltot_inc_max,ltot_av_max,tnow_max,FORMAT='(21e14.6)'
     close,1

  endfor
  
;;===============================================================

  return
end
