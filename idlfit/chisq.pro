
pro chisq,name,filt_arr,filt_wav_arr,flux_arr,errup,$
          errlo,ERROPTION=errop,AVOPTION=avop,$
          LIMIT=limit,NDIS=ndis,DISTANCE=dis,DISERR=dis_err

;calculate chisq for each set of (Mc, Sigma, ms, inc, dis) by changing
;Av, and output models based on the chisq
;
;INPUTS:
;   NAME: used in output files
;   FILT_ARR: array of the filter names
;   FILT_WAV_ARR: array of the filter wavelength
;   FLUX_ARR: array of the fluxes in Jy
;   ERRUP: fractional upper error of the fluxes (array or scalar)
;   ERRLO: fractional lower error of the fluxes (array or scalar)
;   IF ERRLO is not set, or ERROPTION is set to be 1, ONLY ERRUP is
;   used for both upper and lower errors
;  
;OUTPUT:
;   files with "model": best models with different (Mc, Sigma, ms)
;   files with "fullmodel": best models with different (Mc, Sigma, ms,
;                           inc, dis)
;   parameter.dat files: lists parameters for the best models
;   parameter_range.dat files: lists parameter range for the best models
;   eps files: figures showing the fitting SEDs
;
;OPTIONAL INPUTS:
;   erroption: (array or scalar, default 1) 1 to convert flux and
;              error to log, 2 to convert the upper and lower errors
;              differently, 3 similar to 2 but has a special treatment
;              for the large (>=100%) errors
;   avoption: set to 1 to search Av in the range corresponding to
;             [0,5Sigma], set to 0 to search Av in the range [0,1000].
;             Default is 0.
;   limit: (array or scalar) 0 for exact values, 1 for upper
;          limits, -1 for lower limits. Default is 0
;   ndis: number of distances to search. 1 to use the exact provided
;         dis, 2 to use dis+-dis_err, 3 or more to search evenly within
;         [dis-dis_err, dis+dis_err].
;   dis: the central value of the distance range to search
;   diserr: the half width of the distance range to search
;   If dis is not provided, [0.1,10] kpc
;   range will be searched with 0.1 kpc interval.
;


  print,''
  print,'calculating chi square'
  print,''
  
;;===============================================================
;; Update following parameters if the model grid is updated.
  
;; the current model grid has 9 core masses (mc), 4 ambient surface
;; densities (sigma), 12 stellar masses (ms), and 40 inclinations (mu).
;; the inclinations are even in cos space from 1 to -1 
;; mu=1 and 40 (imu=2 and 39, etc.) give same SED and can be added to 
;; give better quality.

  dir0='../Model_SEDs/'


  nmc=15
  mc_str_arr=['10','20','30','40','50','60','80','100','120','160','200','240','320','400','480']
  mc_arr=float(mc_str_arr)

  nsigma=4
  sigma_str_arr=['0.1','0.316','1','3.16']
  sigma_arr=float(sigma_str_arr)

  nms=14
  ms_str_arr=['0.5','1.0','2.0','4.0','8.0','12.0','16.0','24.0','32.0','48.0','64.0','96.0','128.0','160.0']
  ms_arr=float(ms_str_arr)
  
  nmu=20
  mu_arr=reverse(findgen(nmu)/float(nmu)+1./float(nmu)/2.)
  theta_arr=acos(mu_arr)/!pi*180.

  nlambda=500 ;; number of wavelength for the SED. This is fixed.

  nav=100
  
  pc=3.0857d18
  lsun=3.845d33
  mH=1.6733d-24
  clight=2.9979e14
;;===============================================================



  

;;===============================================================
;; check the inputs

  ;; check the filt_arr and flux_arr inputs
  if n_elements(filt_arr) ne n_elements(flux_arr) then begin
     print,filt_arr
     print,flux_arr
     message,'filt_arr and flux_arr do not match',/con
     stop
  endif
  nfilt=n_elements(filt_arr)
  for ifilt=1,nfilt do begin
     if (file_test('../Model_SEDs/flux_filt/'+filt_arr(ifilt-1)+'.fits') eq 0) then begin
        message,'filt_arr contains filter not included in Model_SEDs/flux_filt/',/con
        stop
     endif
     if flux_arr(ifilt-1) le 0. then begin
        message,'negative flux_arr',/con
        stop        
     endif
  endfor

  ;; check the filt_wav_arr
  if n_elements(filt_wav_arr) ne nfilt then begin
     message,'filt_wav_arr and filt_arr do not match',/con
     stop
  endif
  a=where(filt_wav_arr le 0.)
  if a(0) ge 0 then begin
     message,'negative flux_wav_arr',/con
     stop        
  endif
     
  
  ;; check the err input
  if n_elements(errup) ne 1 and n_elements(errup) ne nfilt then begin
     message,'errup does not match filt_arr',/con
     stop
  endif
  errup_arr=fltarr(nfilt)+errup
  a=where(errup_arr lt 0.)
  if a(0) ge 0 then begin
     message,'negative err',/con
     stop             
  endif
  if n_elements(errlo) ne 0 then begin
     if n_elements(errlo) ne 1 and n_elements(errlo) ne nfilt then begin
        message,'errlo does not match filt_arr',/con
        stop
     endif
     errlo_arr=fltarr(nfilt)+errlo
     a=where(errlo_arr lt 0.)
     if a(0) ge 0 then begin
        message,'negative err',/con
        stop             
     endif
  endif else begin
     errlo_arr=errup_arr
  endelse
  a=where(errlo_arr ge 1. or errup_arr ge 1.)
  if a(0) ge 0 then begin
     print,''
     print,'WARN: >100% err'
     print,''
  endif

  ;; check the erroption input
  if n_elements(errop) eq 0 then errop=1
  if n_elements(errop) ne 1 and n_elements(errop) ne nfilt then begin
     message,'erroption does not match filt_arr',/con
     stop
  endif
  errop_arr=intarr(nfilt)+errop
  a=where(errop_arr ne 1 and errop_arr ne 2 and errop_arr ne 3)
  if a(0) ge 0 then begin
     message,'errop must be 1, 2 or 3',/con
     stop     
  endif

  ;; check the avoption input
  if n_elements(avop) eq 0 then avop=0
  if n_elements(avop) ne 1 then begin
     message,'avop must a scalar',/con
     stop
  endif
  if avop ne 0 and avop ne 1 then begin
     message,'avop must be 0 or 1',/con
     stop     
  endif

  ;; check the name input
  if n_elements(name) eq 0 or n_elements(name) gt 1 then begin
     message,'input name wrong',/con
     stop
  endif 

  ;; check the limit input
  if n_elements(limit) eq 0 then limit=0
  if n_elements(limit) ne 1 and n_elements(limit) ne nfilt then begin
     message,'limit does not match filt_arr',/con
     stop
  endif
  limit_arr=intarr(nfilt)+limit
  a=where(limit_arr ne 1 and limit_arr ne -1 and limit_arr ne 0)
  if a(0) ge 0 then begin
     message,'errop must be 0, 1 or -1',/con
     stop     
  endif

  ;; check the distance input and set up the dis array
  if n_elements(dis) eq 0 then begin
     ndis=100
     dis_arr=(findgen(ndis)+1)*100.
  endif else if n_elements(dis) gt 1 then begin
     message,'dis must a scalar',/con
     stop
  endif else if dis le 0. then begin
     message,'invalid dis',/con
     stop
  endif else begin
     if n_elements(dis_err) eq 0 or n_elements(dis_err) gt 1 then begin
        print,''
        print,'dis_err error, use dis as exact value'
        print,''
        ndis=1
        dis_arr=[dis]
     endif else if dis_err le 0. then begin
        print,''
        print,'dis_err invalid, use dis as exact value'
        print,''
        ndis=1
        dis_arr=[dis]
     endif else begin
        dis_max=dis+dis_err
        dis_min=dis-dis_err
        if dis_min le 0. then begin
           print,''
           print,'dis_min less than 0, resetting it to dis/10'
           print,''
           dis_min=dis/10.
        endif
        if n_elements(ndis) eq 0 or n_elements(ndis) gt 1 then begin
           print,''
           print,'ndis error, resetting ndis to 1'
           print,''
           ndis=1
        endif else if ndis le 0 then begin
           print,''
           print,'ndis invalid, resetting ndis to 1'
           print,''
           ndis=1
        endif
        if ndis eq 1 then begin
           dis_arr=[dis]
        endif else if ndis eq 2 then begin
           dis_arr=[dis_min,dis_max]
        endif else if ndis eq 3 then begin
           dis_arr=[dis_min,dis,dis_max]
        endif else begin
           dis_arr=dis_min+(findgen(ndis)/float(ndis-1))*(dis_max-dis_min)
        endelse
     endelse
  endelse
  
  ;; output the input parameters
  dirout=name
  if (file_test(name,/directory) eq 0) then file_mkdir,dirout

  fileintxt=dirout+'/'+name+'.input.dat'
  openw,1,fileintxt
  printf,1,'number of wavelengths,   number of distance,   Av option'
  printf,1,nfilt,ndis,avop
  printf,1,''
  printf,1,'filter name,   wavelength (micron),   flux (Jy),   nuF_nu (erg/s/cm2),   upper percentage error,   lower percentage error,   error option,   limit'
  for ifilt=1,nfilt do begin
     printf,1,filt_arr(ifilt-1),filt_wav_arr(ifilt-1),flux_arr(ifilt-1),flux_arr(ifilt-1)/1.e23*3.e14/filt_wav_arr(ifilt-1),$
            errup_arr(ifilt-1),errlo_arr(ifilt-1),errop_arr(ifilt-1),limit_arr(ifilt-1),format='(A,F,e14.6,e14.6,F,F,I,I)'
  endfor
  printf,1,''
  printf,1,'distance array (pc)'
  for idis=1,ndis do begin
     printf,1,dis_arr(idis-1)
  endfor
  close,1

  a=where(finite(flux_arr) eq 0 or finite(errup_arr) eq 0 or finite(errlo_arr) eq 0)
  if a(0) ge 0 then begin
     print,'Input flux/error contains INF or NAN, please check!'
     stop
  endif

;;===============================================================





;;===============================================================

;; read in the fluxes convolved with filters  
  filtdata=fltarr(nmc,nsigma,nms,nmu,nfilt)
  for ifilt=1,nfilt do begin
     fitsname=dir0+'/flux_filt/'+filt_arr(ifilt-1)+'.fits'
     filtdata_temp=readfits(fitsname,/silent)
     filtdata(*,*,*,*,ifilt-1)=filtdata_temp
  endfor
   
;; This is opacity file for the foreground extinction
  readcol,'../Model_SEDs/parfiles/kmh.par',klam,junk1,junk2,kkap,/silent
  kapv=interpol(kkap,klam,0.55)                    

;; initializing some arrays
  chisq_arr=fltarr(nmc,nsigma,nms,nmu,ndis)
  av_chisq_arr=fltarr(nmc,nsigma,nms,nmu,ndis)
  chisq_nonlimit_arr=fltarr(nmc,nsigma,nms,nmu,ndis)
  I_chisq_arr=fltarr(nmc,nsigma,nms,nmu,ndis,nfilt)
  
  for imc=1,nmc do begin
     mc_str=mc_str_arr(imc-1)
     mc=mc_arr(imc-1)
     
     for isigma=1,nsigma do begin
        sigma_str=sigma_str_arr(isigma-1)
        sigma=sigma_arr(isigma-1)

        if avop eq 0 then begin ; if avop=0, then Av is sampled evenly from 0 to 100 (This can be adjusted)
           av_min=0.
           av_max=100.
           av_arr=av_min+findgen(nav)/float(nav-1)*(av_max-av_min)
        endif else begin ; if avop=1, then Av is sampled within a range which is related to the Sigma parameter (ambient surface density)
           av_clump=sigma/mH/1.8d21/2.
           av_min=0.
           av_max=av_clump*5.
           av_arr=av_min+findgen(nav)/float(nav-1)*(av_max-av_min)
        endelse

        for ims=1,nms do begin
           ms_str=ms_str_arr(ims-1)
           ms=ms_arr(ims-1)
                         
           for imu=1,nmu do begin
              
              filtflux=filtdata(imc-1,isigma-1,ims-1,imu-1,*)
              a=where(filtflux gt 0.) ;; in case there is bad model data.
              if a(0) ge 0 then begin

                 index_wav=a
                 
                 filtflux1=filtflux(a)
                 filt_wav_arr1=filt_wav_arr(a)
                 flux_arr1=flux_arr(a)
                 errup_arr1=errup_arr(a)
                 errlo_arr1=errlo_arr(a)
                 errop_arr1=errop_arr(a)
                 limit_arr1=limit_arr(a)

                 nfit=n_elements(flux_arr1); number of data points to fit, including upper/lower limits
                 filt_wav_fit_arr=filt_wav_arr1

                 kap_fit_arr=interpol(kkap,klam,filt_wav_fit_arr) ;; foreground opacities at those wavelengths

                 flux_fit_log_arr=fltarr(nfit); observed fluxes in log used to fit
                 errup_fit_log_arr=fltarr(nfit); observed upper error in log used to fit
                 errlo_fit_log_arr=fltarr(nfit) ; observed lower error in log used to fit
                 
                 a=where(errop_arr1 eq 1); first way to convert fluxes and errors into log
                 if a(0) ge 0 then begin
                    flux_fit_log_arr(a)=alog10(flux_arr1(a))-0.5*errup_arr1(a)^2/alog(10.); errup_arr contains fractional errors
                    errup_fit_log_arr(a)=errup_arr1(a)/alog(10.); this is absolute error in log space
                    errlo_fit_log_arr(a)=errup_arr1(a)/alog(10.)
                 endif
                 a=where(errop_arr1 ge 2) ; second way to convert fluxes and errors into log
                 if a(0) ge 0 then begin
                    flux_fit_log_arr(a)=alog10(flux_arr1(a))
                    errup_fit_log_arr(a)=alog10(1.+errup_arr1(a)); this is absolute error in log space                    
                    errlo_fit_log_arr(a)=-alog10(1.-errup_arr1(a)); note 100% error means a infinite lower error in log
                 endif
                 a=where(finite(errlo_fit_log_arr) eq 0)
                 if a(0) ge 0 then errlo_fit_log_arr(a)=1.e33; set to a very large number

                 a=where(limit_arr1 eq 1) ; for upper limits
                 if a(0) ge 0 then begin
                    errlo_fit_log_arr(a)=1.e33; set to a very large number, no constraint on the lower side
                 endif
                 a=where(limit_arr1 eq -1) ; for lower limits
                 if a(0) ge 0 then begin
                    errup_fit_log_arr(a)=1.e33 ; set to a very large number, no constraint on the upper side
                 endif

                 a=where(errlo_fit_log_arr eq 0.)
                 if a(0) ge 0 then errlo_fit_log_arr(a)=1.e-33 ; set to a very small number to avoid singularity
                 a=where(errup_fit_log_arr eq 0.)
                 if a(0) ge 0 then errup_fit_log_arr(a)=1.e-33 ; set to a very small number to avoid singularity
                 
                 for idis=1,ndis do begin
                    
                    dis=dis_arr(idis-1)
                    fnorm=1./(4.d0*!pi)/dis^2/pc/pc*lsun
                    modelflux_arr1=filtflux1*fnorm                          ; now in erg/s/cm^2. it was in lsun
                    modelflux_arr1=modelflux_arr1/3.e14*filt_wav_arr1*1.e23 ; now in Jy
                    modelflux_fit_log_arr=alog10(modelflux_arr1)
                                                          
                    for iav=1,nav do begin ; searching the best Av for particular set of (mc, sigma, ms, mu, dis)
                       
                       av=av_arr(iav-1)
                       modelflux_fit_log_av_arr=modelflux_fit_log_arr-0.4*av*kap_fit_arr/kapv ;extincted

                       chisq_fit_arr=fltarr(nfit)
                       iflimit_arr=fltarr(nfit)
                       for ifit=1,nfit do begin
                          flux_t=flux_fit_log_arr(ifit-1)
                          modelflux_t=modelflux_fit_log_av_arr(ifit-1)
                          errup_t=errup_fit_log_arr(ifit-1)
                          errlo_t=errlo_fit_log_arr(ifit-1)
                          limit_t=limit_arr1(ifit-1)
                          errop_t=errop_arr1(ifit-1)

                          iflimit_t=0
                          if modelflux_t ge flux_t then begin; use upper error
                             chisq_t=(modelflux_t-flux_t)^2/errup_t^2
                             if limit_t eq -1 then iflimit_t=1
                          endif else begin
                             if errop_t eq 3 and limit_t eq 0 and errlo_t gt 1.e30 then begin; only when errop set to be 3, for not upper/lower limit, and with very large errors   
                                sigma0=2.
                                x_temp=(flux_t-modelflux_t)^2/sigma0^2
                                chisq_t=alog(x_temp+sqrt(x_temp^2+1.))
                             endif else begin
                                chisq_t=(modelflux_t-flux_t)^2/errlo_t^2
                             endelse
                             if limit_t eq 1 then iflimit_t=1
                          endelse

                          chisq_fit_arr(ifit-1)=chisq_t
                          iflimit_arr(ifit-1)=iflimit_t
                       endfor

                       chisq=total(chisq_fit_arr)/float(nfit) ; reduced chisq, used for ranking
                       a=where(iflimit_arr eq 0)
                       nfit_nonlimit=n_elements(a)
                       chisq_nonlimit=total(chisq_fit_arr)/float(nfit_nonlimit) ; reduced chisq, used for estimate the average chisq contribution from each data points (non counting limits which have 0 contribution)

                       if iav eq 1 then begin
                          chisq_min=chisq; search for minimum chisq for different Av
                          av_chisq_min=av
                          chisq_nonlimit_min=chisq_nonlimit
                       endif else begin
                          if chisq lt chisq_min then begin
                             chisq_min=chisq ; search for minimum chisq for different Av
                             av_chisq_min=av
                             chisq_nonlimit_min=chisq_nonlimit
                          endif
                       endelse
                       
                    endfor                       
                                           
                    chisq_arr(imc-1,isigma-1,ims-1,imu-1,idis-1)=chisq_min ; minimum chisq for this set of (mc, sigma, ms, mu, dis)
                    av_chisq_arr(imc-1,isigma-1,ims-1,imu-1,idis-1)=av_chisq_min ; the Av to reach this chisq for this set of (mc, sigma, ms, mu, dis)
                    chisq_nonlimit_arr(imc-1,isigma-1,ims-1,imu-1,idis-1)=chisq_nonlimit_min;
                    I_chisq_arr(imc-1,isigma-1,ims-1,imu-1,idis-1,index_wav)=modelflux_fit_log_arr-0.4*av_chisq_min*kap_fit_arr/kapv ; the model flux (log of Jy) for minimum chsq for this set of (mc, sigma, ms, mu, dis)
                    
                 endfor

              endif else begin
                 chisq_arr(imc-1,isigma-1,ims-1,*,*)=-1.
                 av_chisq_arr(imc-1,isigma-1,ims-1,*,*)=-1.
                 chisq_nonlimit_arr(imc-1,isigma-1,ims-1,*,*)=-1.
                 I_chisq_arr(imc-1,isigma-1,ims-1,*,*,*)=0.
              endelse
                
           endfor
        endfor
     endfor
  endfor

;;===============================================================




  
;;===============================================================
;; Output the results

  mkhdr,head,chisq_arr
  writefits,dirout+'/'+name+'.output.chisqarr.fits',chisq_arr,head 
  writefits,dirout+'/'+name+'.output.avchisqarr.fits',av_chisq_arr,head
  
  mkhdr,head,I_chisq_arr
  writefits,dirout+'/'+name+'.output.Ichisqarr.fits',I_chisq_arr,head 
  
  chisq_max=max(chisq_arr)
  a=where(chisq_arr lt 0.)
  if a(0) ge 0 then chisq_arr(a)=chisq_max*2.

  fileout=dirout+'/'+name+'.output.allmodel.parameter.dat'
  fileout_full=dirout+'/'+name+'.output.allmodel.full.parameter.dat'
  
  index_sort=sort(chisq_arr)    ; sort the chisq array.

  openw,1,fileout
  openw,2,fileout_full

  imodel=1l
  imodel_full=1l
  repeatcheck_arr=intarr(imc,isigma,ims)
  chisq_min=chisq_arr(index_sort(0))
  n_index_sort=n_elements(index_sort)
  for ii=0l,n_index_sort-1 do begin

     index=index_sort(ii)
     chisq=chisq_arr(index)
     av=av_chisq_arr(index)
     chisq_nonlimit=chisq_nonlimit_arr(index)
     index_cube=array_indices(chisq_arr,index)
     imc=index_cube(0)
     isigma=index_cube(1)
     ims=index_cube(2)
     imu=index_cube(3)
     if ndis eq 1 then begin
        idis=0
     endif else begin
        idis=index_cube(4)
     endelse
     dis=dis_arr(idis)
     repeatcheck=repeatcheck_arr(imc,isigma,ims)

     if chisq gt chisq_max then break
          
     suffix=strtrim(string(imc+1,'(I02)'),2)+'_'+$
            strtrim(string(isigma+1,'(I02)'),2)+'_'+$
            strtrim(string(ims+1,'(I02)'),2)
     suffix1=suffix+'_'+strtrim(string(imu+1,'(I02)'),2)
     
     imodel_str=strtrim(string(imodel),2)
     imodel_full_str=strtrim(string(imodel_full),2)
     mc_str=mc_str_arr(imc)
     sigma_str=sigma_str_arr(isigma)
     ms_str=ms_str_arr(ims)
     theta_str=string(theta_arr(imu))
     dis_str=string(dis)
     av_str=string(av)
     chisq_str=string(chisq)

     filename=dir0+'/model_info/'+suffix+'.dat'
     readcol,filename,mcore,sigma,mstar,rcore,menv,theta,rstar,lstar,tstar,mdisk,rdisk,mdotd,ltot,tnow,skipline=1,/silent ; model informations
     rcore=rcore*206265.

     sedname1=dir0+'/sed/'+suffix1+'.dat'
     readcol,sedname1,lambda_arr,I_arr,/silent 
     nu_arr=clight/lambda_arr
     ltot_inc=tsum(nu_arr,I_arr/nu_arr)
     
     kap_arr=interpol(kkap,klam,lambda_arr)
     I_arr=I_arr*10.^(-0.4*av*kap_arr/kapv)
     ltot_av=tsum(nu_arr,I_arr/nu_arr)
     
     if ii eq 0 then begin
        printf,1,'model No., SED No., chisq, chisq(non-limit), mc (msun), Sigma (g/cm^2), ms (msun), mu (deg), dis (pc), av, rcore (AU), menv (msun), theta_w (deg), rstar (rsun), lstar (lsun), tstar (K), mdisk (msun), rdisk (AU), mdotd (msun/yr), ltot (Lsun), ltot_inc (Lsun), ltot_av (Lsun), tnow (yr)'
        printf,2,'model No., SED No., chisq, chisq(non-limit), mc (msun), Sigma (g/cm^2), ms (msun), mu (deg), dis (pc), av, rcore (AU), menv (msun), theta_w (deg), rstar (rsun), lstar (lsun), tstar (K), mdisk (msun), rdisk (AU), mdotd (msun/yr), ltot (Lsun), ltot_inc (Lsun), ltot_av (Lsun), tnow (yr)'
     endif
     if repeatcheck eq 0 then begin
        printf,1,imodel_str+' '+suffix1+' ',chisq,chisq_nonlimit,mcore,sigma,mstar,theta_arr(imu),dis,av,$
               rcore,menv,theta,rstar,lstar,tstar,mdisk,rdisk,mdotd,ltot,ltot_inc,ltot_av,tnow,FORMAT='(A,21e14.6)'
        imodel=imodel+1
        repeatcheck_arr(imc,isigma,ims)=1
     endif
     printf,2,imodel_full_str+' '+suffix1+' ',chisq,chisq_nonlimit,mcore,sigma,mstar,theta_arr(imu),dis,av,$
            rcore,menv,theta,rstar,lstar,tstar,mdisk,rdisk,mdotd,ltot,ltot_inc,ltot_av,tnow,FORMAT='(A,21e14.6)'
     imodel_full=imodel_full+1

  endfor
  close,1
  close,2
;;===============================================================

  return
end
