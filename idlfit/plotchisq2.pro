
pro plotchisq2,name,xwav_arr,ywav_arr,chisq_max=chisq_show_max,nxp=nxp,nyp=nyp

;plot distribution of chisq with color-luminosity diagram.

;INPUTS:
;   name: used in input and output files
;   xwav_arr, ywav_arr: indicate the wavelength to draw the color-flux
;   diagram, x axis will be log(F(xwav)-F(ywav)), y axis will by log(F(ywav))
;  
;OPTIONAL INPUTS:
;   chisq_max: maximum chisq to show, default 50
;   nxp: number of panels in x
;   nyp: number of panels in y

  print,''
  print,'plotchisq2'
  print,''
  
;;===============================================================
;; check the inputs

  ;; check the name input
  if n_elements(name) eq 0 or n_elements(name) gt 1 then begin
     message,'input name wrong',/con
     stop
  endif

  if n_elements(xwav_arr) ne n_elements(ywav_arr) then begin
     message,'xwav_arr and ywav_arr do not match',/con
     stop
  endif
  npanel=n_elements(xwav_arr)
  if npanel eq 0 then begin
     message,'xwav_arr and ywav_arr missing',/con
     stop
  endif

  if n_elements(chisq_show_max) eq 0 then begin
     chisq_show_max=50.
  endif else if n_elements(chisq_show_max) gt 1 then begin
     message,'chisq_show_max wrong',/con
     stop
  endif

  if n_elements(nxp) gt 1 or n_elements(nyp) gt 1 then begin
     message,'nxp, nyp wrong',/con
     stop
  endif
  if n_elements(nxp) eq 0 then begin
     if n_elements(nyp) eq 0 then begin
        nxp=3
        nyp=ceil(float(npanel)/float(nxp))
     endif
     if nyp eq 0 then begin
        nxp=3
        nyp=ceil(float(npanel)/float(nxp))
     endif else begin
        nxp=ceil(float(npanel)/float(nyp))
     endelse
  endif
  if nxp eq 0 then begin
     if n_elements(nyp) eq 0 then begin
        nxp=3
        nyp=ceil(float(npanel)/float(nxp))
     endif
     if nyp eq 0 then begin
        nxp=3
        nyp=ceil(float(npanel)/float(nxp))
     endif else begin
        nxp=ceil(float(npanel)/float(nyp))
     endelse
  endif else begin
     if n_elements(nyp) eq 0 then begin
        nyp=ceil(float(npanel)/float(nxp))
     endif
     if nyp eq 0 then begin
        nyp=ceil(float(npanel)/float(nxp))
     endif
  endelse
  
;;===============================================================




  

;;===============================================================
;; calcualte the fluxes and colors of the models
  
  nbest=5
  ifctreverse=1
  ctable=34
  
  dirout=name
  filein1=dirout+'/'+name+'.output.chisqarr.fits'
  filein2=dirout+'/'+name+'.output.Ichisqarr.fits'
  filein3=dirout+'/'+name+'.input.dat'
  fileout=dirout+'/'+name+'.plotchisq.color.eps'
  
  chisq_arr=readfits(filein1,head,/silent)
  I_chisq_arr=readfits(filein2,head,/silent)
  
  nmc=sxpar(head,'NAXIS1')
  nsigma=sxpar(head,'NAXIS2')
  nms=sxpar(head,'NAXIS3')
  nmu=sxpar(head,'NAXIS4')
  ndis=sxpar(head,'NAXIS5')
  nfilt=sxpar(head,'NAXIS6')
  
  readcol,filein3,nfilt1,dump,dump,format='(I,I,I)',skipline=1,numline=1,/silent
  nfilt1=nfilt1(0)
  readcol,filein3,filt_arr,filt_wav_arr,flux_obs_arr,flux_nufnu_obs_arr,errup_arr,errlo_arr,errop_arr,limit_arr,format='(A,F,F,F,F,F,I,I)',skipline=4,numline=nfilt1,/silent
  if nfilt1 ne nfilt then begin
     print,'nfilt wrong'
     stop
  endif
  
  chisq_arr1=fltarr(nmc,nsigma,nms)
  I_chisq_arr1=fltarr(nmc,nsigma,nms,nfilt)
  for imc=1,nmc do begin
     for isigma=1,nsigma do begin
        for ims=1,nms do begin
           chisq_min_temp=min(chisq_arr(imc-1,isigma-1,ims-1,*,*))
           chisq_arr1(imc-1,isigma-1,ims-1)=chisq_min_temp
           a=where(chisq_arr(imc-1,isigma-1,ims-1,*,*) eq chisq_min_temp)
           for ifilt=1,nfilt do begin
              I_chisq_arr_temp=I_chisq_arr(imc-1,isigma-1,ims-1,*,*,ifilt-1)
              I_chisq_arr1(imc-1,isigma-1,ims-1,ifilt-1)=I_chisq_arr_temp(a(0))
           endfor
        endfor
     endfor
  endfor
  chisq_arr=chisq_arr1
  I_chisq_arr=I_chisq_arr1

  chisq_1d_arr=reform(chisq_arr,nmc*nsigma*nms,1)
  I_chisq_1d_arr=fltarr(nmc*nsigma*nms,nfilt)
  for ifilt=1,nfilt do begin
     I_chisq_1d_arr(*,ifilt-1)=reform(I_chisq_arr(*,*,*,ifilt-1),nmc*nsigma*nms,1)
  endfor
  a=where(chisq_1d_arr gt 0. and chisq_1d_arr lt chisq_show_max)
  chisq_1d_arr=chisq_1d_arr(a)
  I_chisq_1d_arr=I_chisq_1d_arr(a,*)
  nmodel=n_elements(chisq_1d_arr)
  index_sort=sort(chisq_1d_arr)
  
  x_all_arr=fltarr(nmodel,npanel)
  y_all_arr=fltarr(nmodel,npanel)
  xtitle_arr=fltarr(npanel)
  ytitle_arr=fltarr(npanel)
  for ipanel=1,npanel do begin
     a=where(filt_arr eq xwav_arr(ipanel-1))
     x_all_arr(*,ipanel-1)=I_chisq_1d_arr(*,a(0))
     xtitle_arr(ipanel-1)=filt_wav_arr(a(0))
     a=where(filt_arr eq ywav_arr(ipanel-1))
     y_all_arr(*,ipanel-1)=I_chisq_1d_arr(*,a(0))
     ytitle_arr(ipanel-1)=filt_wav_arr(a(0))
  endfor
  x_all_arr=x_all_arr-y_all_arr
  xtitle_arr=strtrim(string(xtitle_arr,'(F6.1)'),2)
  ytitle_arr=strtrim(string(ytitle_arr,'(F6.1)'),2)
  
  x_obs_arr=fltarr(npanel)
  y_obs_arr=fltarr(npanel)
  x_errup_arr=fltarr(npanel)
  y_errup_arr=fltarr(npanel)
  x_errlo_arr=fltarr(npanel)
  y_errlo_arr=fltarr(npanel)
  for ipanel=1,npanel do begin
     a=where(filt_arr eq xwav_arr(ipanel-1))
     ifilt=a(0)+1
     flux_obs=flux_obs_arr(ifilt-1)
     errup=errup_arr(ifilt-1)
     errlo=errlo_arr(ifilt-1)
     errop=errop_arr(ifilt-1)
     limit=limit_arr(ifilt-1)
     if errop eq 1 then begin
        flux_fit_log=alog10(flux_obs)-0.5*errup^2/alog(10.)
        errup_fit_log=errup/alog(10.)
        errlo_fit_log=errlo/alog(10.)
     endif else begin
        flux_fit_log=alog10(flux_obs)
        errup_fit_log=alog10(1.+errup)
        errlo_fit_log=-alog10(1.-errup)
     endelse
     if finite(errlo_fit_log) eq 0 then errlo_fit_log=1.e10
     if limit eq 1 then errlo_fit_log=1.e10
     if limit eq -1 then errup_fit_log=1.e10
     x_obs_arr(ipanel-1)=flux_fit_log
     x_errlo_arr(ipanel-1)=errlo_fit_log
     x_errup_arr(ipanel-1)=errup_fit_log
  endfor
  for ipanel=1,npanel do begin
     a=where(filt_arr eq ywav_arr(ipanel-1))
     ifilt=a(0)+1
     flux_obs=flux_obs_arr(ifilt-1)
     errup=errup_arr(ifilt-1)
     errlo=errlo_arr(ifilt-1)
     errop=errop_arr(ifilt-1)
     limit=limit_arr(ifilt-1)
     if errop eq 1 then begin
        flux_fit_log=alog10(flux_obs)-0.5*errup^2/alog(10.)
        errup_fit_log=errup/alog(10.)
        errlo_fit_log=errlo/alog(10.)
     endif else begin
        flux_fit_log=alog10(flux_obs)
        errup_fit_log=alog10(1.+errup)
        errlo_fit_log=-alog10(1.-errup)
     endelse
     if finite(errlo_fit_log) eq 0 then errlo_fit_log=1.e10
     if limit eq 1 then errlo_fit_log=1.e10
     if limit eq -1 then errup_fit_log=1.e10
     y_obs_arr(ipanel-1)=flux_fit_log
     y_errlo_arr(ipanel-1)=errlo_fit_log
     y_errup_arr(ipanel-1)=errup_fit_log
  endfor
  x_obs_arr=x_obs_arr-y_obs_arr
  x_errup_arr=sqrt((x_errup_arr)^2+(y_errlo_arr)^2)
  x_errlo_arr=sqrt((x_errlo_arr)^2+(y_errup_arr)^2)
  
  chisq_rag=[min(chisq_1d_arr),chisq_show_max]
  if ifctreverse eq 1 then begin
     mag_arr=255-bytscl(alog10(chisq_1d_arr),min=alog10(chisq_rag(0)),max=alog10(chisq_rag(1)))
  endif else begin
     mag_arr=bytscl(alog10(chisq_1d_arr),min=alog10(chisq_rag(0)),max=alog10(chisq_rag(1)))
  endelse
  
;;===============================================================





;;===============================================================
;; plot the figures

  xs_img=3.
  ys_img=3.
  
  l_margin=0.9
  r_margin=1.4
  b_margin=0.8
  t_margin=0.4
  
  xs_ps=(xs_img+l_margin)*float(nxp)+r_margin
  ys_ps=(ys_img+b_margin)*float(nyp)+t_margin
  
  l_img_arr=l_margin+findgen(nxp)*(xs_img+l_margin)
  r_img_arr=l_img_arr+xs_img
  b_img_arr=b_margin+findgen(nyp)*(ys_img+b_margin)
  t_img_arr=b_img_arr+ys_img
  
  xs_bar=0.15
  ys_bar=t_img_arr(nyp-1)-b_img_arr(0)
  l_bar=r_img_arr(nxp-1)+0.2
  r_bar=l_bar+xs_bar
  b_bar=b_img_arr(0)
  t_bar=b_bar+ys_bar
  
  set_plot,'PS'
  device,xsize=xs_ps,ysize=ys_ps,/inches,/encapsulate,$
         /portrait,xoffset=0.,yoffset=0.,bits=24,filename=fileout,/color
  
  for iyp=nyp,1,-1 do begin
     for ixp=1,nxp do begin
        
        l_img=l_img_arr(ixp-1)
        r_img=r_img_arr(ixp-1)
        b_img=b_img_arr(iyp-1)
        t_img=t_img_arr(iyp-1)
        pos=[l_img/xs_ps,b_img/ys_ps,r_img/xs_ps,t_img/ys_ps]
        
        ipanel=(nyp-iyp)*nxp+ixp
        
        if ipanel le npanel then begin
           
           x_arr=x_all_arr(*,ipanel-1)
           y_arr=y_all_arr(*,ipanel-1)
           xtitle=xtitle_arr(ipanel-1)
           ytitle=ytitle_arr(ipanel-1)
           xtit=textoidl('log(F_{'+xtitle+'\mum}/F_{'+ytitle+'\mum})')
           ytit=textoidl('log(F_{'+ytitle+'\mum})')
           x_obs=x_obs_arr(ipanel-1)
           y_obs=y_obs_arr(ipanel-1)
           x_errup=x_errup_arr(ipanel-1)
           y_errup=y_errup_arr(ipanel-1)
           x_errlo=x_errlo_arr(ipanel-1)
           y_errlo=y_errlo_arr(ipanel-1)
           
           xrag=x_obs+[-1.5,1.5]
           yrag=y_obs+[-1.5,1.5]
           
           loadcolors
           plot,[0,0],[1,1],/nodata,/noerase,position=pos,xstyle=4+1,ystyle=4+1,$
                xrange=xrag,yrange=yrag
           
           loadct,0,/silent
           x_errup1=min([3.*x_errup,xrag(1)-x_obs])/3.
           x_errlo1=min([3.*x_errlo,abs(xrag(0)-x_obs)])/3.
           y_errup1=min([3.*y_errup,yrag(1)-y_obs])/3.
           y_errlo1=min([3.*y_errlo,abs(yrag(0)-y_obs)])/3.        
           polyfill,[x_obs-3.*x_errlo1, x_obs+3.*x_errup1, x_obs+3.*x_errup1, x_obs-3.*x_errlo1, x_obs-3.*x_errlo1],$
                    [y_obs-3.*y_errlo1, y_obs-3.*y_errlo1, y_obs+3.*y_errup1, y_obs+3.*y_errup1, y_obs-3.*y_errlo1],$
                    color=200
           
           x_errup2=min([x_errup,xrag(1)-x_obs])
           x_errlo2=min([x_errlo,abs(xrag(0)-x_obs)])
           y_errup2=min([y_errup,yrag(1)-y_obs])
           y_errlo2=min([y_errlo,abs(yrag(0)-y_obs)])        
           polyfill,[x_obs-x_errlo2, x_obs+x_errup2, x_obs+x_errup2, x_obs-x_errlo2, x_obs-x_errlo2],$
                    [y_obs-y_errlo2, y_obs-y_errlo2, y_obs+y_errup2, y_obs+y_errup2, y_obs-y_errlo2],$
                    color=160
           oplot,x_obs*[1,1],yrag,linestyle=2,thick=1,color=0
           oplot,xrag,y_obs*[1,1],linestyle=2,thick=1,color=0
           
           plot,[0,0],[1,1],/nodata,/noerase,position=pos,xstyle=1,ystyle=1,$
                xrange=xrag,yrange=yrag,ytitle=ytit,xtitle=xtit,$
                charthick=3,charsize=1.2,xthick=3,ythick=3
           
           plotsym,0,0.6,/fill
           loadct,ctable,/silent
           plots,reverse(x_arr),reverse(y_arr),psym=8,color=reverse(mag_arr),thick=3
           
           loadcolors
           plots,x_arr(index_sort(0)),y_arr(index_sort(0)),psym=1,symsize=3,thick=5,color=0
           plots,x_arr(index_sort(1:4)),y_arr(index_sort(1:4)),psym=1,symsize=2,thick=5,color=0
           
        endif
     endfor
  endfor
  
  pos_bar=[l_bar/xs_ps,b_bar/ys_ps,r_bar/xs_ps,t_bar/ys_ps]

  chisq_tickv=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,$
               1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000]
  chisq_tickname=['0.01',' ',' ',' ','0.05',' ',' ',' ',' ','0.1',' ',' ',' ','0.5',' ',' ',' ',' ',$
                  '1',' ',' ',' ','5',' ',' ',' ',' ','10',' ',' ',' ','50',' ',' ',' ',' ','100',' ',' ',' ','500',' ',' ',' ',' ','1000']
  a=where(chisq_tickv ge chisq_rag(0) and chisq_tickv le chisq_rag(1))
  y_tickname=chisq_tickname(a)
  y_tickv=chisq_tickv(a)

  loadct,ctable,/silent
  if ifctreverse eq 1 then begin
     bar=reverse(bindgen(256))#replicate(1b,10)
  endif else begin
     bar=bindgen(256)#replicate(1b,10)
  endelse
  bar=transpose(bar)
  tv,bar,l_bar,b_bar,xsize=xs_bar,ysize=ys_bar,/inches
  
  loadcolors
  contour,bar,position=pos_bar,color=0,/noerase,/nodata,xstyle=4,ystyle=4
  axis,xstyle=4,ystyle=1,yaxis=1,yrange=chisq_rag,ytitle=textoidl('\chi^2'),/ylog,$
       ythick=3,yticklen=0.5,charthick=3,charsize=1.2,ytickname=y_tickname,ytickv=y_tickv,yticks=n_elements(y_tickv)-1

  xyouts,(l_margin+xs_ps-r_margin)/2./xs_ps,0.97,name,alignment=0.5,charthick=3,charsize=1.2,/normal

  device,/close
  
;;===============================================================

  return

end
