
pro plotsed,name

;plot seds of the selected best models.

;INPUTS:
;   name: used in input and output files
;

;;===============================================================
;; check the inputs

  ;; check the name input
  if n_elements(name) eq 0 or n_elements(name) gt 1 then begin
     message,'input name wrong',/con
     stop
  endif 

;;===============================================================



  

;;===============================================================
;; plot the SEDs of best models

  pc=3.0857d18
  lsun=3.845d33

  dirout=name
  fileinput=dirout+'/'+name+'.input.dat'

  filein=dirout+'/'+name+'.output.bestmodel.parameter.dat'
  fileout=dirout+'/'+name+'.plotsed.bestmodel.eps'

  filein_full=dirout+'/'+name+'.output.bestmodel.full.parameter.dat'
  fileout_full=dirout+'/'+name+'.plotsed.bestmodel.full.eps'
  
  readcol,fileinput,nfilt,ndis,dump,format='(I,I,I)',skipline=1,numline=1,/silent
  nfilt=nfilt(0)
  ndis=ndis(0)
  readcol,fileinput,filt_arr,filt_wav_arr,flux_arr,flux_nufnu_arr,errup_arr,errlo_arr,errop_arr,limit_arr,format='(A,F,F,F,F,F,I,I)',skipline=4,numline=nfilt,/silent
  
  readcol,'../Model_SEDs/parfiles/kmh.par',klam,junk1,junk2,kkap,/silent

  flux_max=max(flux_nufnu_arr)
  flux_min=min(flux_nufnu_arr)
  flux_max_log=alog10(flux_max)
  flux_min_log=alog10(flux_min)

  ymax=ceil(flux_max_log)+1
  ymin=floor(flux_min_log)-2
  ymax=10.^ymax
  ymin=10.^ymin
  xtit='!4k!3 (!4l!3m)'
  ytit='!4m!3F!d!4m!3!n (ergs s!u-1!n cm!u-2!n)'

  xs_img=6.
  ys_img=4.
  
  l_margin=0.9
  r_margin=0.3
  b_margin=0.7
  t_margin=0.5

  l_img=l_margin
  r_img=l_img+xs_img
  b_img=b_margin
  t_img=b_img+ys_img
  
  xs_ps=xs_img+l_margin+r_margin
  ys_ps=ys_img+b_margin+t_margin

  pos=[l_img/xs_ps,b_img/ys_ps,r_img/xs_ps,t_img/ys_ps]
  
  for ifile=1,2 do begin
     
     if ifile eq 1 then begin
        filei=filein
        fileo=fileout
     endif else begin
        filei=filein_full
        fileo=fileout_full
     endelse
        
     readcol,filei,modelno,sedno,chisq,chisq_nonlimit,mcore,sigma,mstar,inc,dis,av,rcore,menv,theta,rstar,lstar,tstar,mdisk,rdisk,mdotd,ltot,ltot_inc,ltot_av,tnow,$
             format='(I,A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)',skipline=1,/silent
     ntop=n_elements(modelno)
     
     set_plot,'ps'
     device,xsize=xs_ps,ysize=ys_ps,/inches,/color,xoffset=0.,yoffset=0.,$
            bits=24,filename=fileo,/encapsulate
     
     plot,[0.,0.],[1.,1.],/nodata,/xlog,/ylog,position=pos,yrange=[ymin,ymax],xrange=[1.,1000.],xstyle=1+4,ystyle=1+4

     ;; plot the model SEDs
     for ii=ntop-1,0,-1 do begin

        dis_t=dis(ii)
        suffix1=sedno(ii)
        filename1='../Model_SEDs/sed/'+suffix1+'.dat'
        readcol,filename1,lambda_arr,I_arr,/silent 
        
        fnorm=1./(4.d0*!pi)/dis_t^2/pc/pc*lsun
        I_arr=I_arr*fnorm       ; nfu
        
        kap_arr=interpol(kkap,klam,lambda_arr)
        kapv=interpol(kkap,klam,0.55)
        I_arr=I_arr*10.^(-0.4*av(ii)*kap_arr/kapv)
         
        loadct,0,/silent
        if ii eq 0 then begin
           color=0
           thick=3
        endif else begin
           color=120
           thick=1
        endelse
        oplot,lambda_arr,I_arr,color=color,linestyle=0,thick=thick

     endfor
     
     ;; plot the original data points
     flux_nufnu_arr=flux_arr/1.e23*3.e14/filt_wav_arr ; in erg /s /cm^2
     errup_nufnu_arr=errup_arr
     errlo_nufnu_arr=errlo_arr

     ;loadcolors
     ;plotsym,8,0.5,/FILL
     ;plots,filt_wav_arr,flux_nufnu_arr,psym=8,color=6

     ;; plot the data points used to fit
     flux_nufnu_log_arr=fltarr(nfilt)
     errup_nufnu_log_arr=fltarr(nfilt)
     errlo_nufnu_log_arr=fltarr(nfilt)
     
     a=where(errop_arr eq 1)
     if a(0) ge 0 then begin
        flux_nufnu_log_arr(a)=alog10(flux_nufnu_arr(a))-0.5*errup_nufnu_arr(a)^2/alog(10.)
        errup_nufnu_log_arr(a)=errup_nufnu_arr(a)/alog(10.)
        errlo_nufnu_log_arr(a)=errup_nufnu_arr(a)/alog(10.)
     endif
     a=where(errop_arr ge 2)
     if a(0) ge 0 then begin
        flux_nufnu_log_arr(a)=alog10(flux_nufnu_arr(a))
        errup_nufnu_log_arr(a)=alog10(1.+errup_nufnu_arr(a))
        errlo_nufnu_log_arr(a)=-alog10(1.-errlo_nufnu_arr(a))
     endif
     a=where(finite(errlo_nufnu_log_arr) eq 0)
     if a(0) ge 0 then errlo_nufnu_log_arr(a)=1.e33
     
     a=where(limit_arr eq 1)
     if a(0) ge 0 then begin
        errlo_nufnu_log_arr(a)=1.e33
     endif
     a=where(limit_arr eq -1)
     if a(0) ge 0 then begin
        errup_nufnu_log_arr(a)=1.e33
     endif

     flux_nufnu_lo_log_arr=flux_nufnu_log_arr-errlo_nufnu_log_arr
     flux_nufnu_up_log_arr=flux_nufnu_log_arr+errup_nufnu_log_arr

     flux_plot_arr=10.^flux_nufnu_log_arr
     errup_plot_arr=abs(10.^flux_nufnu_up_log_arr-10.^flux_nufnu_log_arr)
     errlo_plot_arr=abs(10.^flux_nufnu_lo_log_arr-10.^flux_nufnu_log_arr)
     
     loadcolors
     a=where(limit_arr eq 0)
     plotsym,8,0.5,/FILL
     plots,filt_wav_arr,flux_plot_arr,psym=8,color=5
     oploterror,filt_wav_arr(a),flux_plot_arr(a),errup_plot_arr(a),errthick=2,psym=3,errcolor=5,/hibar
     oploterror,filt_wav_arr(a),flux_plot_arr(a),errlo_plot_arr(a),errthick=2,psym=3,errcolor=5,/lobar
     
     ;; a=where(errup_arr gt 0.1 or errlo_arr gt 0.1 or limit_arr eq 0)
     ;; if a(0) ge 0 then begin
     ;;    plotsym,8,0.5,/FILL
     ;;    plots,filt_wav_arr(a),flux_plot_arr(a),psym=8,color=5
     ;; endif
  
     ;; plot the upper and lower limits
     a=where(limit_arr eq 1)
     if a(0) ge 0 then begin
        plotsym,1,2,thick=5
        plots,filt_wav_arr(a),flux_plot_arr(a),psym=8,color=5
        oploterror,filt_wav_arr(a),flux_plot_arr(a),errup_plot_arr(a),errthick=2,psym=3,errcolor=5,/hibar
     endif
     a=where(limit_arr eq -1)
     if a(0) ge 0 then begin
        plotsym,2,2,thick=5
        plots,filt_wav_arr(a),flux_plot_arr(a),psym=8,color=5
        oploterror,filt_wav_arr(a),flux_plot_arr(a),errlo_plot_arr(a),errthick=2,psym=3,errcolor=5,/lobar
     endif
     
     plot,[0.,0.],[1.,1.],/nodata,/noerase,/xlog,/ylog,position=pos,yrange=[ymin,ymax],xrange=[1.,1000.],$
          charthick=3,xstyle=1,ystyle=1,ytitle=ytit,xtitle=xtit,$
          charsize=1.,xthick=3,ythick=3,title=name,yminor=9

     device,/close
  endfor
;;===============================================================

  
  return
end
