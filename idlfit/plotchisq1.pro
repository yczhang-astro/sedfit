
pro plotchisq1,name,chisq_max=chisq_show_max

;plot distribution of chisq with secondary parameters.

;INPUTS:
;   name: used in input and output files
;
;OPTIONAL INPUTS:
;   chisq_max: maximum chisq to show, default 50

  print,''
  print,'plotchisq1'
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

  dirout=name
  filein=dirout+'/'+name+'.output.allmodel.parameter.dat'
  fileout=dirout+'/'+name+'.plotchisq.secondary.eps'

  nbest=5
  ifctreverse=1
  ctable=34
  
  ;; nxp=4
  ;; nyp=3
  ;; npanel=12
  ;; xsuffix_arr=['menv', 'menv', 'mstar', 'mstar', 'mstar', 'inc', 'inc', 'sigma', 'inc', 'inc', 'av', 'av']
  ;; ysuffix_arr=['ltot', 'mstar', 'ltot', 'mdotd', 'tnow', 'thetaw', 'av', 'av_sigma', 'ltot', 'flash2', 'ltot', 'flash1']
  ;; xlog_arr=[1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0]
  ;; ylog_arr=[1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1]
  ;; xtit_arr=textoidl(['M_{env} (M_\odot)', 'M_{env} (M_\odot)', 'm_* (M_\odot)', 'm_* (M_\odot)', 'm_* (M_\odot)', '\theta_{view} (deg)', '\theta_{view} (deg)', '\Sigma_{cl} (g cm^{-2})', '\theta_{view} (deg)', '\theta_{view} (deg)', 'A_V (mag)', 'A_V (mag)'])
  ;; ytit_arr=textoidl(['L_{tot} (L_\odot)', 'm_* (M_\odot)', 'L_{tot} (L_\odot)', 'dm_*/dt (M_\odot yr^{-1})', 'Age (yr)', '\theta_{w,esc} (deg)', 'A_V (mag)', 'A_V/A_V(\Sigma_{cl}/2)', 'L_{tot} (L_\odot)', 'L_{tot,inc}/L_{tot}', 'L_{tot} (L_\odot)', 'L_{tot,av}/L_{tot,inc}'])
  ;; xragl_arr=[1., 1., 1.0, 1.0, 1.0, 0., 0., 0.04, 0., 0., 0., 0.]
  ;; xragh_arr=[1000., 1000., 200., 200., 200., 90., 90., 4., 90., 90., 300., 300.]
  ;; yragl_arr=[1.e3, 1.0, 1.e3, 1.e-5, 1.e4, 0., 0., 0., 1.e3, 0.01, 1.e3, 0.1]
  ;; yragh_arr=[1.e7, 200., 1.e7, 1.e-2, 1.e6, 90., 300., 5., 1.e7, 10., 1.e7, 2.]
  
  nxp=5
  nyp=2
  npanel=12
  xsuffix_arr=['menv', 'menv', 'mstar', 'mstar', 'mstar', 'inc', 'inc', 'sigma', 'inc', 'av']
  ysuffix_arr=['ltot', 'mstar', 'ltot', 'mdotd', 'tnow', 'thetaw', 'av', 'av_sigma', 'flash', 'flash']
  xlog_arr=[1, 1, 1, 1, 1, 0, 0, 1, 0, 0]
  ylog_arr=[1, 1, 1, 1, 1, 0, 0, 0, 1, 1]
  xtit_arr=textoidl(['M_{env} (M_\odot)', 'M_{env} (M_\odot)', 'm_* (M_\odot)', 'm_* (M_\odot)', 'm_* (M_\odot)', '\theta_{view} (deg)', '\theta_{view} (deg)', '\Sigma_{cl} (g cm^{-2})', '\theta_{view} (deg)', 'A_V (mag)'])
  ytit_arr=textoidl(['L_{tot} (L_\odot)', 'm_* (M_\odot)', 'L_{tot} (L_\odot)', 'dm_*/dt (M_\odot yr^{-1})', 'Age (yr)', '\theta_{w,esc} (deg)', 'A_V (mag)', 'A_V/A_V(\Sigma_{cl}/2)', 'L_{inc,av}/L_{tot}', 'L_{inc,av}/L_{tot}'])
  xragl_arr=[1., 1., 1.0, 1.0, 1.0, 0., 0., 0.04, 0., 0.]
  xragh_arr=[1000., 1000., 200., 200., 200., 90., 90., 4., 90., 300.]
  yragl_arr=[1.e3, 1.0, 1.e3, 1.e-5, 1.e4, 0., 0., 0., 0.01, 0.01]
  yragh_arr=[1.e7, 200., 1.e7, 1.e-2, 1.e6, 90., 300., 5., 10., 10.]
;;===============================================================




  
;;===============================================================
;; readin parameter files
  
  mH=1.6733d-24
  
  readcol,filein,modelno_arr,sedno_arr,chisq_arr,chisq_nonlimit_arr,mcore_arr,sigma_arr,mstar_arr,$
          inc_arr,dis_arr,av_arr,rcore_arr,menv_arr,thetaw_arr,rstar_arr,lstar_arr,tstar_arr,mdisk_arr,$
          rdisk_arr,mdotd_arr,ltot_arr,ltot_inc_arr,ltot_av_arr,tnow_arr,$
          format='(I,A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)',skipline=1,/silent
  
  a=where(chisq_arr lt chisq_show_max)
  modelno_arr=modelno_arr(a)
  sedno_arr=sedno_arr(a)
  chisq_arr=chisq_arr(a)
  chisq_nonlimit_arr=chisq_nonlimit_arr(a)
  mcore_arr=mcore_arr(a)
  sigma_arr=sigma_arr(a)
  mstar_arr=mstar_arr(a)
  inc_arr=inc_arr(a)
  dis_arr=dis_arr(a)
  av_arr=av_arr(a)
  rcore_arr=rcore_arr(a)
  menv_arr=menv_arr(a)
  thetaw_arr=thetaw_arr(a)
  rstar_arr=rstar_arr(a)
  lstar_arr=lstar_arr(a)
  tstar_arr=tstar_arr(a)
  mdisk_arr=mdisk_arr(a)
  rdisk_arr=rdisk_arr(a)
  mdotd_arr=mdotd_arr(a)
  ltot_arr=ltot_arr(a)
  ltot_inc_arr=ltot_inc_arr(a)
  ltot_av_arr=ltot_av_arr(a)
  tnow_arr=tnow_arr(a)
  
  lacc_arr=ltot_arr-lstar_arr
  lratio_arr=lstar_arr/ltot_arr
  theta_diff_arr=inc_arr-thetaw_arr
  flash_arr=ltot_av_arr/ltot_arr
  flash1_arr=ltot_av_arr/ltot_inc_arr
  flash2_arr=ltot_inc_arr/ltot_arr
  av_clump_arr=sigma_arr/mH/1.8d21/2.
  av_sigma_arr=av_arr/av_clump_arr

  chisq_rag=[chisq_arr(0),chisq_show_max]
  if ifctreverse eq 1 then begin
     mag_arr=255-bytscl(alog10(chisq_arr),min=alog10(chisq_rag(0)),max=alog10(chisq_rag(1)))
  endif else begin
     mag_arr=bytscl(alog10(chisq_arr),min=alog10(chisq_rag(0)),max=alog10(chisq_rag(1)))
  endelse
  
;;===============================================================





;;===============================================================
;; plot figures

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
           
           xlog=xlog_arr(ipanel-1)
           ylog=ylog_arr(ipanel-1)
           xrag=[xragl_arr(ipanel-1),xragh_arr(ipanel-1)]
           yrag=[yragl_arr(ipanel-1),yragh_arr(ipanel-1)]
           xtit=xtit_arr(ipanel-1)
           ytit=ytit_arr(ipanel-1)
           xsuffix=xsuffix_arr(ipanel-1)
           ysuffix=ysuffix_arr(ipanel-1)
           void=execute('x_arr='+xsuffix+'_arr')
           void=execute('y_arr='+ysuffix+'_arr')
           
           loadcolors
           plot,[0,0],[1,1],/nodata,/noerase,position=pos,xstyle=1,ystyle=1,$
                xlog=xlog,ylog=ylog,xrange=xrag,yrange=yrag,ytitle=ytit,xtitle=xtit,$
                charthick=3,charsize=1.2,xthick=3,ythick=3

           if ipanel eq 1 then begin
              xx=xrag
              yy1=10.^(2.5+1.8*alog10(xrag))
              yy2=10.^(0.54+1.9*alog10(xrag))
              loadct,0,/silent
              oplot,xx,yy1,linestyle=2,thick=3,color=120
              oplot,xx,yy2,linestyle=2,thick=3,color=120
           endif
           if ipanel eq 2 then begin
              xx=xrag
              yy1=xrag*0.01
              yy2=xrag*0.1
              yy3=xrag*1.
              loadct,0,/silent
              oplot,xx,yy1,linestyle=2,thick=3,color=120
              oplot,xx,yy2,linestyle=2,thick=3,color=120
              oplot,xx,yy3,linestyle=2,thick=3,color=120
           endif

           if ipanel eq 6 then begin
              xx=xrag
              yy=xrag
              loadct,0,/silent
              oplot,xx,yy,linestyle=2,thick=3,color=120
           endif

           plotsym,0,0.8,/fill
           loadct,ctable,/silent
           plots,reverse(x_arr),reverse(y_arr),psym=8,color=reverse(mag_arr),thick=3
           
           loadcolors
           plots,x_arr(0),y_arr(0),psym=1,symsize=3,thick=5,color=0
           plots,x_arr(1:4),y_arr(1:4),psym=1,symsize=2,thick=5,color=0
                     
           if xlog eq 0 then begin
              xpos=xrag(0)+(xrag(1)-xrag(0))*0.05
           endif else begin
              xpos=10.^(alog10(xrag(0))+(alog10(xrag(1))-alog10(xrag(0)))*0.05)
           endelse
           if ylog eq 0 then begin
              ypos=yrag(1)-(yrag(1)-yrag(0))*0.1
           endif else begin
              ypos=10.^(alog10(yrag(1))-(alog10(yrag(1))-alog10(yrag(0)))*0.1)
           endelse
           xyouts,xpos,ypos,'('+string(96B+byte(ipanel))+')',charthick=3,charsize=1.2,color=0,alignment=0
           
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

  xyouts,(l_margin+xs_ps-r_margin)/2./xs_ps,0.98,name,alignment=0.5,charthick=3,charsize=1.2,/normal

  device,/close

;;===============================================================

  return
end
