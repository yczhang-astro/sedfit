
pro plotsed1,name,chisq_max=chisq_show_max

;plot seds of models with wide range of chisq

;INPUTS:
;   name: used in input and output files
;
;OPTIONAL INPUTS:
;   chisq_max: maximum chisq to show, default 50

  print,''
  print,'plot SEDs of more models'
  print,''

;;===============================================================
;; check the inputs
  
  ;; check the name input
  if n_elements(name) eq 0 or n_elements(name) gt 1 then begin
     message,'input name wrong',/con
     stop
  endif

  if n_elements(chisq_show_max) eq 0 then begin
     chisq_show_max=50.
  endif else if n_elements(chisq_show_max) gt 1 then begin
     message,'chisq_show_max wrong',/con
     stop
  endif

;;===============================================================



  

;;===============================================================
;; plot the SEDs of best models

  pc=3.0857d18
  lsun=3.845d33

  dirout=name
  fileinput=dirout+'/'+name+'.input.dat'

  filein=dirout+'/'+name+'.output.allmodel.full.parameter.dat'
  fileout=dirout+'/'+name+'.plotsed.modelrange.eps'
  
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
        
  readcol,filein,modelno,sedno,chisq,chisq_nonlimit,mcore,sigma,mstar,inc,dis,av,rcore,menv,theta,rstar,lstar,tstar,mdisk,rdisk,mdotd,ltot,ltot_inc,ltot_av,tnow,$
          format='(I,A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)',skipline=1,/silent
  chisq_min=chisq(0)
  a=where(chisq lt chisq_show_max)
  if a(0) lt 0 then a=[0]
  ntop=n_elements(a)
  sedno=sedno(a)
  chisq=chisq(a)
  dis=dis(a)
  av=av(a)

  xs_img=6.
  ys_img=4.
  
  l_margin=0.9
  r_margin=1.4
  b_margin=0.7
  t_margin=0.5

  l_img=l_margin
  r_img=l_img+xs_img
  b_img=b_margin
  t_img=b_img+ys_img
  
  xs_ps=xs_img+l_margin+r_margin
  ys_ps=ys_img+b_margin+t_margin

  pos=[l_img/xs_ps,b_img/ys_ps,r_img/xs_ps,t_img/ys_ps]
  
  xs_bar=0.15
  ys_bar=ys_img
  l_bar=r_img+0.2
  r_bar=l_bar+xs_bar
  b_bar=b_img
  t_bar=b_bar+ys_bar

  pos_bar=[l_bar/xs_ps,b_bar/ys_ps,r_bar/xs_ps,t_bar/ys_ps]
  
  set_plot,'ps'
  device,xsize=xs_ps,ysize=ys_ps,/inches,/color,xoffset=0.,yoffset=0.,$
         bits=24,filename=fileout,/encapsulate
     
  plot,[0.,0.],[1.,1.],/nodata,/xlog,/ylog,position=pos,yrange=[ymin,ymax],xrange=[1.,1000.],xstyle=1+4,ystyle=1+4
  
  ;; plot the model SEDs
  for ii=ntop-1,0,-1 do begin
     
     dis_t=dis(ii)
     suffix1=sedno(ii)
     filename1='../Model_SEDs/sed/'+suffix1+'.dat'
     readcol,filename1,lambda_arr,I_arr,/silent 
        
     fnorm=1./(4.d0*!pi)/dis_t^2/pc/pc*lsun
     I_arr=I_arr*fnorm          ; nfu
     
     kap_arr=interpol(kkap,klam,lambda_arr)
     kapv=interpol(kkap,klam,0.55)
     I_arr=I_arr*10.^(-0.4*av(ii)*kap_arr/kapv)
     
     loadct,34,/silent
     if ii eq 0 then begin
        loadcolors
        color=0
        thick=3
        oplot,lambda_arr,I_arr,color=color,linestyle=0,thick=thick
     endif else begin
        oplot,lambda_arr,I_arr,color=255-bytscl(alog10(chisq(ii)),min=alog10(chisq_min),max=alog10(chisq_show_max)),linestyle=0,thick=1
     endelse
        
  endfor
     
  ;; plot the original data points
  flux_nufnu_arr=flux_arr/1.e23*3.e14/filt_wav_arr ; in erg /s /cm^2
  errup_nufnu_arr=errup_arr
  errlo_nufnu_arr=errlo_arr
  
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

  symbolcolor=0
  
  loadcolors
  a=where(limit_arr eq 0)
  plotsym,8,0.5
  plots,filt_wav_arr,flux_plot_arr,psym=8,color=symbolcolor
  oploterror,filt_wav_arr(a),flux_plot_arr(a),errup_plot_arr(a),errthick=2,psym=3,errcolor=symbolcolor,/hibar
  oploterror,filt_wav_arr(a),flux_plot_arr(a),errlo_plot_arr(a),errthick=2,psym=3,errcolor=symbolcolor,/lobar
  
  a=where(errup_arr gt 0.1 or errlo_arr gt 0.1 or limit_arr eq 0)
  if a(0) ge 0 then begin
     plotsym,8,0.5,/FILL
     plots,filt_wav_arr(a),flux_plot_arr(a),psym=8,color=symbolcolor
  endif
  
  ;; plot the upper and lower limits
  a=where(limit_arr eq 1)
  if a(0) ge 0 then begin
     plotsym,1,2,thick=5
     plots,filt_wav_arr(a),flux_plot_arr(a),psym=8,color=symbolcolor
     oploterror,filt_wav_arr(a),flux_plot_arr(a),errup_plot_arr(a),errthick=2,psym=3,errcolor=symbolcolor,/hibar
  endif
  a=where(limit_arr eq -1)
  if a(0) ge 0 then begin
     plotsym,2,2,thick=5
     plots,filt_wav_arr(a),flux_plot_arr(a),psym=8,color=symbolcolor
     oploterror,filt_wav_arr(a),flux_plot_arr(a),errlo_plot_arr(a),errthick=2,psym=3,errcolor=symbolcolor,/lobar
  endif

  plot,[0.,0.],[1.,1.],/nodata,/noerase,position=pos,/xlog,/ylog,yrange=[ymin,ymax],xrange=[1.,1000.],$
       charthick=3,xstyle=1,ystyle=1,ytitle=ytit,xtitle=xtit,$
       charsize=1.,xthick=3,ythick=3,title=name,yminor=9

  chisq_tickv=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,$
               1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000]
  chisq_tickname=['0.01',' ',' ',' ','0.05',' ',' ',' ',' ','0.1',' ',' ',' ','0.5',' ',' ',' ',' ',$
                  '1',' ',' ',' ','5',' ',' ',' ',' ','10',' ',' ',' ','50',' ',' ',' ',' ','100',' ',' ',' ','500',' ',' ',' ',' ','1000']
  a=where(chisq_tickv ge chisq_min and chisq_tickv le chisq_show_max)
  y_tickv=chisq_tickv(a)
  y_tickname=chisq_tickname(a)

  if a(0) ge 0 then begin
     loadct,34,/silent
     bar=reverse(bindgen(256))#replicate(1b,10)
     bar=transpose(bar)
     tv,bar,l_bar,b_bar,xsize=xs_bar,ysize=ys_bar,/inches
     
     loadcolors
     contour,bar,position=pos_bar,color=0,/noerase,/nodata,xstyle=4,ystyle=4
     axis,xstyle=4,ystyle=1,yaxis=1,yrange=[chisq_min,chisq_show_max],ytitle=textoidl('\chi^2'),/ylog,$
          ythick=3,yticklen=0.5,charthick=3,charsize=1.2,ytickv=y_tickv,ytickname=y_tickname,yticks=n_elements(y_tickv)-1
  endif
     
  device,/close

;;===============================================================


  return
end
