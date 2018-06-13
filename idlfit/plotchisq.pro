
pro plotchisq,name,chisq_max=chisq_show_max

;plot distribution of chisq with mc, sigma and ms.
;The distribution of chisq with mc and sigma shows the best chisq for
;each set of (mc,sigma), and the values of (ms,mu,av) to achieve them.
;Same for the distribution of chisq with sigma and ms, or with mc and ms.

;INPUTS:
;   name: used in input and output files
;
;OPTIONAL INPUTS:
;   chisq_max: maximum chisq to show, default 50

  print,''
  print,'plotchisq'
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
  
  av_show_max=200.
  nbest=5
  av_rag=[0.,av_show_max]
  
;;===============================================================




  
;;===============================================================

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

;;===============================================================




  
;;===============================================================
;; prepare the chisq cube and find the best models

  dirout=name
  filein1=dirout+'/'+name+'.output.chisqarr.fits'
  filein2=dirout+'/'+name+'.output.avchisqarr.fits'
  fileout=dirout+'/'+name+'.plotchisq.primary.eps'
  fileout2=dirout+'/'+name+'.plotchisq.primary.simplified.eps'
  
  chisq_arr=readfits(filein1,head,/silent)
  avchisq_arr=readfits(filein2,head,/silent)

  nmc_temp=sxpar(head,'NAXIS1')
  nsigma_temp=sxpar(head,'NAXIS2')
  nms_temp=sxpar(head,'NAXIS3')
  nmu_temp=sxpar(head,'NAXIS4')
  
  if nmc ne nmc_temp then begin
     print,'check nmc'
     stop
  endif
  if nsigma ne nsigma_temp then begin
     print,'check nsigma'
     stop
  endif
  if nms ne nms_temp then begin
     print,'check nms'
     stop
  endif
  if nmu ne nmu_temp then begin
     print,'check nmu'
     stop
  endif
  
  chisq_arr1=fltarr(nmc,nsigma,nms)
  avchisq_arr1=fltarr(nmc,nsigma,nms)
  chisqmu_arr=fltarr(nmc,nsigma,nms)
  for imc=1,nmc do begin
     for isigma=1,nsigma do begin
        for ims=1,nms do begin
           chisq_min_temp=min(chisq_arr(imc-1,isigma-1,ims-1,*,*))
           chisq_arr1(imc-1,isigma-1,ims-1)=chisq_min_temp
           a=where(chisq_arr(imc-1,isigma-1,ims-1,*,*) eq chisq_min_temp)
           index_cube=array_indices(chisq_arr(imc-1,isigma-1,ims-1,*,*),a(0))
           chisqmu_arr(imc-1,isigma-1,ims-1)=mu_arr(index_cube(3))
           if n_elements(index_cube) eq 4 then begin
              avchisq_arr1(imc-1,isigma-1,ims-1)=avchisq_arr(imc-1,isigma-1,ims-1,index_cube(3))
           endif else begin
              avchisq_arr1(imc-1,isigma-1,ims-1)=avchisq_arr(imc-1,isigma-1,ims-1,index_cube(3),index_cube(4))
           endelse
           if chisq_min_temp eq -1. then begin
              chisqmu_arr(imc-1,isigma-1,ims-1)=-1.
           endif
        endfor
     endfor
  endfor
  chisq_arr=chisq_arr1
  avchisq_arr=avchisq_arr1
  
  chisq_min_arr=fltarr(nbest)
  mc_min_arr=fltarr(nbest)
  sigma_min_arr=fltarr(nbest)
  ms_min_arr=fltarr(nbest)
  
  chisq_arr1=chisq_arr
  chisq_max=max(chisq_arr1)
  a=where(chisq_arr1 lt 0.)
  if a(0) ge 0 then chisq_arr1(a)=chisq_max+10.
  index_sort=sort(chisq_arr1)
  for ibest=1,nbest do begin
     index=index_sort(ibest-1)
     index_cube=array_indices(chisq_arr,index)
     imc_min=index_cube(0)+1
     isigma_min=index_cube(1)+1
     ims_min=index_cube(2)+1
     chisq_min_arr(ibest-1)=chisq_arr1(index)
     mc_min_arr(ibest-1)=alog10(mc_arr(imc_min-1))
     sigma_min_arr(ibest-1)=alog10(sigma_arr(isigma_min-1))
     ms_min_arr(ibest-1)=alog10(ms_arr(ims_min-1))
  endfor

  chisq_rag=[chisq_min_arr(0),chisq_show_max]
  
;;===============================================================




  
;;===============================================================
;; make 2-parameter maps from the chisq cube

  chisq_mc_sigma_arr=fltarr(nmc,nsigma)
  chisqms_mc_sigma_arr=fltarr(nmc,nsigma)
  avchisq_mc_sigma_arr=fltarr(nmc,nsigma)
  chisqmu_mc_sigma_arr=fltarr(nmc,nsigma)
  
  chisq_sigma_ms_arr=fltarr(nsigma,nms)
  chisqmc_sigma_ms_arr=fltarr(nsigma,nms)
  avchisq_sigma_ms_arr=fltarr(nsigma,nms)
  chisqmu_sigma_ms_arr=fltarr(nsigma,nms)
  
  chisq_mc_ms_arr=fltarr(nmc,nms)
  chisqsigma_mc_ms_arr=fltarr(nmc,nms)
  avchisq_mc_ms_arr=fltarr(nmc,nms)
  chisqmu_mc_ms_arr=fltarr(nmc,nms)
  
  for imc=1,nmc do begin
     for isigma=1,nsigma do begin
        chisq_arr_temp=chisq_arr(imc-1,isigma-1,*)
        a=where(chisq_arr_temp gt 0.)
        if a(0) ge 0 then begin
           chisq_min_temp=min(chisq_arr_temp(a))
           b=where(chisq_arr_temp eq chisq_min_temp)
           chisq_mc_sigma_arr(imc-1,isigma-1)=chisq_min_temp
           chisqms_mc_sigma_arr(imc-1,isigma-1)=ms_arr(b(0))
           avchisq_mc_sigma_arr(imc-1,isigma-1)=avchisq_arr(imc-1,isigma-1,b(0))
           chisqmu_mc_sigma_arr(imc-1,isigma-1)=chisqmu_arr(imc-1,isigma-1,b(0))
        endif else begin
           chisq_mc_sigma_arr(imc-1,isigma-1)=-1.
           chisqms_mc_sigma_arr(imc-1,isigma-1)=-1.
           avchisq_mc_sigma_arr(imc-1,isigma-1)=-1.
           chisqmu_mc_sigma_arr(imc-1,isigma-1)=-1.
        endelse
     endfor
  endfor
  
  for isigma=1,nsigma do begin
     for ims=1,nms do begin
        chisq_arr_temp=chisq_arr(*,isigma-1,ims-1)
        a=where(chisq_arr_temp gt 0.)
        if a(0) ge 0 then begin
           chisq_min_temp=min(chisq_arr_temp(a))
           b=where(chisq_arr_temp eq chisq_min_temp)
           chisq_sigma_ms_arr(isigma-1,ims-1)=chisq_min_temp
           chisqmc_sigma_ms_arr(isigma-1,ims-1)=mc_arr(b(0))
           avchisq_sigma_ms_arr(isigma-1,ims-1)=avchisq_arr(b(0),isigma-1,ims-1)
           chisqmu_sigma_ms_arr(isigma-1,ims-1)=chisqmu_arr(b(0),isigma-1,ims-1)
        endif else begin
           chisq_sigma_ms_arr(isigma-1,ims-1)=-1.
           chisqmc_sigma_ms_arr(isigma-1,ims-1)=-1.
           avchisq_sigma_ms_arr(isigma-1,ims-1)=-1.
           chisqmu_sigma_ms_arr(isigma-1,ims-1)=-1.
        endelse
     endfor
  endfor
  
  for imc=1,nmc do begin
     for ims=1,nms do begin
        chisq_arr_temp=chisq_arr(imc-1,*,ims-1)
        a=where(chisq_arr_temp gt 0.)
        if a(0) ge 0 then begin
           chisq_min_temp=min(chisq_arr_temp(a))
           b=where(chisq_arr_temp eq chisq_min_temp)
           chisq_mc_ms_arr(imc-1,ims-1)=chisq_min_temp
           chisqsigma_mc_ms_arr(imc-1,ims-1)=sigma_arr(b(0))
           avchisq_mc_ms_arr(imc-1,ims-1)=avchisq_arr(imc-1,b(0),ims-1)
           chisqmu_mc_ms_arr(imc-1,ims-1)=chisqmu_arr(imc-1,b(0),ims-1)
        endif else begin
           chisq_mc_ms_arr(imc-1,ims-1)=-1.
           chisqsigma_mc_ms_arr(imc-1,ims-1)=-1.
           avchisq_mc_ms_arr(imc-1,ims-1)=-1.
           chisqmu_mc_ms_arr(imc-1,ims-1)=-1.
        endelse
     endfor
  endfor
  
  if chisq_show_max gt 0. then begin
     a=where(chisq_mc_sigma_arr gt chisq_show_max)
     chisq_mc_sigma_arr(a)=-.5
     chisqms_mc_sigma_arr(a)=-.5
     avchisq_mc_sigma_arr(a)=-.5
     chisqmu_mc_sigma_arr(a)=-.5
  endif
  if chisq_show_max gt 0. then begin
     a=where(chisq_sigma_ms_arr gt chisq_show_max)
     chisq_sigma_ms_arr(a)=-.5
     chisqmc_sigma_ms_arr(a)=-.5
     avchisq_sigma_ms_arr(a)=-.5
     chisqmu_sigma_ms_arr(a)=-.5
  endif
  if chisq_show_max gt 0. then begin
     a=where(chisq_mc_ms_arr gt chisq_show_max)
     chisq_mc_ms_arr(a)=-.5
     chisqsigma_mc_ms_arr(a)=-.5
     avchisq_mc_ms_arr(a)=-.5
     chisqmu_mc_ms_arr(a)=-.5
  endif

;;===============================================================




  
;;===============================================================
;; new grid and interpolate

  mc_grid_max=500.
  mc_grid_min=10.
  sigma_grid_max=4.
  sigma_grid_min=0.07
  ms_grid_max=100.
  ms_grid_min=0.5
  mu_grid_max=1.
  mu_grid_min=0.
  
  nmcgrid=500
  nsigmagrid=500
  nmsgrid=500
  nmugrid=500
  
  mc_grid_arr=findgen(nmcgrid)/float(nmcgrid-1)*(alog10(mc_grid_max)-alog10(mc_grid_min))+alog10(mc_grid_min)
  sigma_grid_arr=findgen(nsigmagrid)/float(nsigmagrid-1)*(alog10(sigma_grid_max)-alog10(sigma_grid_min))+alog10(sigma_grid_min)
  ms_grid_arr=findgen(nmsgrid)/float(nmsgrid-1)*(alog10(ms_grid_max)-alog10(ms_grid_min))+alog10(ms_grid_min)
  mu_grid_arr=findgen(nmugrid)/float(nmugrid-1)*(mu_grid_max-mu_grid_min)+mu_grid_min
  
  delt_mc_grid=(alog10(mc_grid_max)-alog10(mc_grid_min))/float(nmcgrid-1)
  delt_sigma_grid=(alog10(sigma_grid_max)-alog10(sigma_grid_min))/float(nsigmagrid-1)
  delt_ms_grid=(alog10(ms_grid_max)-alog10(ms_grid_min))/float(nmsgrid-1)
  
  mc_rag=[alog10(mc_grid_min)-delt_mc_grid/2.,alog10(mc_grid_max)+delt_mc_grid/2.]
  sigma_rag=[alog10(sigma_grid_min)-delt_sigma_grid/2.,alog10(sigma_grid_max)+delt_sigma_grid/2.]
  ms_rag=[alog10(ms_grid_min)-delt_ms_grid/2.,alog10(ms_grid_max)+delt_ms_grid/2.]
  
  mc_tickv=alog10([10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,200.,300.,400.,500.])
  mc_tickname=['10',' ',' ',' ',' ',' ',' ',' ',' ','100',' ',' ',' ',' ']
  sigma_tickv=alog10([0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,2.,3.,4.])
  sigma_tickname=[' ',' ',' ','0.1',' ',' ',' ',' ',' ',' ',' ',' ','1',' ',' ','4']
  ms_tickv=alog10([0.5, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 20., 30., 40., 50., 60., 70., 80., 90.])
  ms_tickname=[' ','1',' ',' ',' ',' ',' ',' ',' ',' ','10',' ',' ',' ',' ',' ',' ',' ',' ']
  
  mc_arr=alog10(mc_arr)
  sigma_arr=alog10(sigma_arr)
  ms_arr=alog10(ms_arr)
  
  chisq_mc_sigma_grid_arr=newgrid2d(mc_arr,sigma_arr,chisq_mc_sigma_arr,mc_grid_arr,sigma_grid_arr,missing=-1.)
  chisq_mc_ms_grid_arr=newgrid2d(mc_arr,ms_arr,chisq_mc_ms_arr,mc_grid_arr,ms_grid_arr,missing=-1.)
  chisq_sigma_ms_grid_arr=newgrid2d(sigma_arr,ms_arr,chisq_sigma_ms_arr,sigma_grid_arr,ms_grid_arr,missing=-1.)
  
  chisqms_mc_sigma_grid_arr=newgrid2d(mc_arr,sigma_arr,chisqms_mc_sigma_arr,mc_grid_arr,sigma_grid_arr,missing=-1.)
  chisqsigma_mc_ms_grid_arr=newgrid2d(mc_arr,ms_arr,chisqsigma_mc_ms_arr,mc_grid_arr,ms_grid_arr,missing=-1.)
  chisqmc_sigma_ms_grid_arr=newgrid2d(sigma_arr,ms_arr,chisqmc_sigma_ms_arr,sigma_grid_arr,ms_grid_arr,missing=-1.)
  
  avchisq_mc_sigma_grid_arr=newgrid2d(mc_arr,sigma_arr,avchisq_mc_sigma_arr,mc_grid_arr,sigma_grid_arr,missing=-1.)
  avchisq_mc_ms_grid_arr=newgrid2d(mc_arr,ms_arr,avchisq_mc_ms_arr,mc_grid_arr,ms_grid_arr,missing=-1.)
  avchisq_sigma_ms_grid_arr=newgrid2d(sigma_arr,ms_arr,avchisq_sigma_ms_arr,sigma_grid_arr,ms_grid_arr,missing=-1.)
  
  chisqmu_mc_sigma_grid_arr=newgrid2d(mc_arr,sigma_arr,chisqmu_mc_sigma_arr,mc_grid_arr,sigma_grid_arr,missing=-1.)
  chisqmu_mc_ms_grid_arr=newgrid2d(mc_arr,ms_arr,chisqmu_mc_ms_arr,mc_grid_arr,ms_grid_arr,missing=-1.)
  chisqmu_sigma_ms_grid_arr=newgrid2d(sigma_arr,ms_arr,chisqmu_sigma_ms_arr,sigma_grid_arr,ms_grid_arr,missing=-1.)
  
  valmin=chisq_rag(0)
  valmax=chisq_rag(1)
  chisq_mc_sigma_grid_arr_new=bytscl(alog10(chisq_mc_sigma_grid_arr),min=alog10(valmin),max=alog10(valmax))
  chisq_mc_ms_grid_arr_new=bytscl(alog10(chisq_mc_ms_grid_arr),min=alog10(valmin),max=alog10(valmax))
  chisq_sigma_ms_grid_arr_new=bytscl(alog10(chisq_sigma_ms_grid_arr),min=alog10(valmin),max=alog10(valmax))
  
  chisqms_mc_sigma_grid_arr_new=bytscl(alog10(chisqms_mc_sigma_grid_arr),min=alog10(ms_grid_min),max=alog10(ms_grid_max))
  chisqsigma_mc_ms_grid_arr_new=bytscl(alog10(chisqsigma_mc_ms_grid_arr),min=alog10(sigma_grid_min),max=alog10(sigma_grid_max))
  chisqmc_sigma_ms_grid_arr_new=bytscl(alog10(chisqmc_sigma_ms_grid_arr),min=alog10(mc_grid_min),max=alog10(mc_grid_max))
  
  valmin=av_rag(0)
  valmax=av_rag(1)
  avchisq_mc_sigma_grid_arr_new=bytscl(avchisq_mc_sigma_grid_arr,min=valmin,max=valmax)
  avchisq_mc_ms_grid_arr_new=bytscl(avchisq_mc_ms_grid_arr,min=valmin,max=valmax)
  avchisq_sigma_ms_grid_arr_new=bytscl(avchisq_sigma_ms_grid_arr,min=valmin,max=valmax)
  
  valmin=mu_grid_min
  valmax=mu_grid_max
  chisqmu_mc_sigma_grid_arr_new=bytscl(chisqmu_mc_sigma_grid_arr,min=valmin,max=valmax)
  chisqmu_mc_ms_grid_arr_new=bytscl(chisqmu_mc_ms_grid_arr,min=valmin,max=valmax)
  chisqmu_sigma_ms_grid_arr_new=bytscl(chisqmu_sigma_ms_grid_arr,min=valmin,max=valmax)

;;===============================================================




  
;;===============================================================
;; plot figure

  nxp=4
  nyp=3
  
  xs_img=3.
  ys_img=3.
  
  l_margin=0.8
  r_margin=0.3
  b_margin=1.4
  t_margin=0.4
  
  xs_ps=(xs_img+r_margin)*float(nxp)+l_margin
  ys_ps=(ys_img+b_margin+t_margin)*float(nyp)
  
  l_img_arr=l_margin+findgen(nxp)*(xs_img+r_margin)
  r_img_arr=l_img_arr+xs_img
  b_img_arr=b_margin+findgen(nyp)*(ys_img+b_margin+t_margin)
  t_img_arr=b_img_arr+ys_img
  
  xs_bar=xs_img-0.2
  ys_bar=0.15
  l_bar_arr=l_img_arr+0.1
  r_bar_arr=l_bar_arr+xs_bar
  t_bar_arr=b_img_arr-0.6
  b_bar_arr=t_bar_arr-ys_bar
  
  chisq_tickv=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,$
               1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000]
  chisq_tickname=['0.01',' ',' ',' ','0.05',' ',' ',' ',' ','0.1',' ',' ',' ','0.5',' ',' ',' ',' ',$
                  '1',' ',' ',' ','5',' ',' ',' ',' ','10',' ',' ',' ','50',' ',' ',' ',' ','100',' ',' ',' ','500',' ',' ',' ',' ','1000']
  
  set_plot,'PS'
  device,xsize=xs_ps,ysize=ys_ps,/inches,/encapsulate,$
         /portrait,xoffset=0.,yoffset=0.,bits=24,filename=fileout,/color
  
  for iyp=nyp,1,-1 do begin
     
     if iyp eq 1 then begin
        xsuffix='sigma'
        ysuffix='ms'
        xtitle=textoidl('\Sigma_{cl} (g cm^{-2})')
        ytitle=textoidl('m_* (M_\odot)')
     endif else if iyp eq 2 then begin
        xsuffix='mc'
        ysuffix='ms'
        xtitle=textoidl('M_c (M_\odot)')
        ytitle=textoidl('m_* (M_\odot)')
     endif else if iyp eq 3 then begin
        xsuffix='mc'
        ysuffix='sigma'
        xtitle=textoidl('M_c (M_\odot)')
        ytitle=textoidl('\Sigma_{cl} (g cm^{-2})')
     endif
     
     for ixp=1,nxp do begin
        
        l_img=l_img_arr(ixp-1)
        r_img=r_img_arr(ixp-1)
        b_img=b_img_arr(iyp-1)
        t_img=t_img_arr(iyp-1)
        pos_img=[l_img/xs_ps,b_img/ys_ps,r_img/xs_ps,t_img/ys_ps]
        
        l_bar=l_bar_arr(ixp-1)
        r_bar=r_bar_arr(ixp-1)
        b_bar=b_bar_arr(iyp-1)
        t_bar=t_bar_arr(iyp-1)
        pos_bar=[l_bar/xs_ps,b_bar/ys_ps,r_bar/xs_ps,t_bar/ys_ps]

        if ixp eq 1 then begin
           zsuffix='chisq'
           ztitle=textoidl('\chi^2')
           z_rag=chisq_rag
           ifctreverse=1
           iflog=1
           bararr=reverse(bindgen(256))
           a=where(chisq_tickv ge chisq_rag(0) and chisq_tickv le chisq_rag(1))
           z_tickname=chisq_tickname(a)
           z_tickv=chisq_tickv(a)
        endif else if ixp eq 2 then begin
           if iyp eq 1 then begin
              zsuffix='chisqmc'
              ztitle=textoidl('M_c (M_\odot)')
              z_rag=[mc_grid_min,mc_grid_max]
              ifctreverse=0
              iflog=1
              bararr=bytscl(newgrid(mc_arr,mc_arr,mc_grid_arr),min=alog10(mc_grid_min),max=alog10(mc_grid_max))
              z_tickname=mc_tickname
              z_tickv=10.^mc_tickv
           endif else if iyp eq 2 then begin
              zsuffix='chisqsigma'
              ztitle=textoidl('\Sigma_{cl} (g cm^{-1})')
              z_rag=[sigma_grid_min,sigma_grid_max]
              ifctreverse=0
              iflog=1
              bararr=bytscl(newgrid(sigma_arr,sigma_arr,sigma_grid_arr),min=alog10(sigma_grid_min),max=alog10(sigma_grid_max))
              z_tickname=sigma_tickname
              z_tickv=10.^sigma_tickv
           endif else if iyp eq 3 then begin
              zsuffix='chisqms'
              ztitle=textoidl('m_* (M_\odot)')
              z_rag=[ms_grid_min,ms_grid_max]
              ifctreverse=0
              iflog=1
              bararr=bytscl(newgrid(ms_arr,ms_arr,ms_grid_arr),min=alog10(ms_grid_min),max=alog10(ms_grid_max))
              z_tickname=ms_tickname
              z_tickv=10.^ms_tickv
           endif
        endif else if ixp eq 3 then begin
           zsuffix='chisqmu'
           ztitle=textoidl('cos\theta_{view}')
           z_rag=[mu_grid_min,mu_grid_max]
           ifctreverse=1
           iflog=0
           bararr=reverse(bytscl(newgrid(mu_arr,mu_arr,mu_grid_arr),min=mu_grid_min,max=mu_grid_max))
           z_tickname=textoidl(['0','0.2','0.4','0.6','0.8','1'])
           z_tickv=[0.,0.2,0.4,0.6,0.8,1.]
        endif else if ixp eq 4 then begin
           zsuffix='avchisq'
           ztitle=textoidl('A_V')
           z_rag=av_rag
           ifctreverse=0
           iflog=0
           bararr=bindgen(256)
           z_tickname=textoidl(['0'])
           z_tickv=[0.]
        endif
        ctable=33
        
        void=execute('z_grid_arr='+zsuffix+'_'+xsuffix+'_'+ysuffix+'_grid_arr')
        void=execute('z_grid_arr_new='+zsuffix+'_'+xsuffix+'_'+ysuffix+'_grid_arr_new')
        void=execute('x_grid_arr='+xsuffix+'_grid_arr')
        void=execute('y_grid_arr='+ysuffix+'_grid_arr')
        void=execute('xrag='+xsuffix+'_rag')
        void=execute('yrag='+ysuffix+'_rag')
        void=execute('xtickv='+xsuffix+'_tickv')
        void=execute('ytickv='+ysuffix+'_tickv')
        void=execute('xtickname='+xsuffix+'_tickname')
        void=execute('ytickname='+ysuffix+'_tickname')
        void=execute('x_min_arr='+xsuffix+'_min_arr')
        void=execute('y_min_arr='+ysuffix+'_min_arr')
        
        if ixp eq 1 then begin
           z_grid_arr1=z_grid_arr
           z_max=max(z_grid_arr1)
           a=where(z_grid_arr1 lt -0.9)
           if a(0) ge 0 then begin
              z_grid_arr1(a)=z_max+100.
           endif
           a=where(z_grid_arr1 lt -0.4)
           if a(0) ge 0 then begin
              z_grid_arr1(a)=z_max+10.
           endif
        endif
        
        nx=n_elements(x_grid_arr)
        ny=n_elements(y_grid_arr)
        
        if ifctreverse eq 1 then begin
           z_grid_arr_new=255-z_grid_arr_new
        endif
        
        loadct,ctable,/silent
        tvlct,red,green,blue,/get
        redmap=bytarr(nx,ny)
        greenmap=bytarr(nx,ny)
        bluemap=bytarr(nx,ny)
        for ix=1,nx do begin
           for iy=1,ny do begin
              redmap(ix-1,iy-1)=red(z_grid_arr_new(ix-1,iy-1))
              greenmap(ix-1,iy-1)=green(z_grid_arr_new(ix-1,iy-1))
              bluemap(ix-1,iy-1)=blue(z_grid_arr_new(ix-1,iy-1))
           endfor
        endfor
        a=where(z_grid_arr lt -0.9)
        if a(0) ge 0 then begin
           redmap(a)=120
           greenmap(a)=120
           bluemap(a)=120
        endif
        a=where(z_grid_arr lt -0.4 and z_grid_arr gt -0.9)
        if a(0) ge 0 then begin
           redmap(a)=255
           greenmap(a)=255
           bluemap(a)=255
        endif
        loadct,0,/silent
        tv,[[[redmap]],[[greenmap]],[[bluemap]]],l_img,b_img,xsize=xs_img,ysize=ys_img,/inches,true=3
        
        xtit=xtitle
        xcharsize=1.2
        if ixp eq 1 then begin
           ytit=ytitle
           ycharsize=1.2
        endif else begin
           ytit=''
           ycharsize=0.0001
        endelse
        
        loadcolors
        contour,z_grid_arr,x_grid_arr,y_grid_arr,position=pos_img,$
                /nodata,/noerase,xstyle=1,ystyle=1,$
                xrange=xrag,yrange=yrag,xtitle=xtit,ytitle=ytit,$
                xtickname=xtickname,xtickv=xtickv,xticks=n_elements(xtickv)-1,$
                ytickname=ytickname,ytickv=ytickv,yticks=n_elements(ytickv)-1,$
                charthick=3,xthick=3,ythick=3,xcharsize=xcharsize,ycharsize=ycharsize
        
        contour,z_grid_arr1,x_grid_arr,y_grid_arr,/overplot,levels=[chisq_min_arr(0)+5.],color=7,thick=10
        contour,z_grid_arr1,x_grid_arr,y_grid_arr,/overplot,levels=[chisq_min_arr(0)+5.],color=5,thick=4
        
        contour,z_grid_arr1,x_grid_arr,y_grid_arr,/overplot,levels=[z_max+50.],color=0,thick=5
        
        for ibest=1,nbest do begin
           if ibest eq 1 then begin
              symsize=3
           endif else begin
              symsize=2
           endelse
           plots,x_min_arr(ibest-1),y_min_arr(ibest-1),psym=1,symsize=symsize,thick=5,color=7,/data
        endfor
        
        bar=bararr#replicate(1b,10)
        loadct,ctable,/silent
        tv,bar,l_bar,b_bar,xsize=xs_bar,ysize=ys_bar,/inches
        
        loadcolors
        contour,bar,position=pos_bar,color=0,/noerase,/nodata,xstyle=4,ystyle=4
        axis,xstyle=1,ystyle=4,xaxis=0,xrange=z_rag,xtitle=ztitle,xlog=iflog,$
             xthick=3,xticklen=0.3,charthick=3,charsize=1.2,$
             xtickname=z_tickname,xtickv=z_tickv,xticks=n_elements(z_tickv)-1
        axis,xstyle=1,ystyle=4,xaxis=1,xrange=z_rag,xtitle=ztitle,xlog=iflog,$
             xthick=3,xticklen=0.3,charthick=3,charsize=0.0001,$
             xtickname=z_tickname,xtickv=z_tickv,xticks=n_elements(z_tickv)-1         
        
     endfor
  endfor

  xyouts,(l_margin+xs_ps-r_margin)/2./xs_ps,0.98,name,alignment=0.5,charthick=3,charsize=1.2,/normal
  
  device,/close

;;===============================================================




;;===============================================================
;; plot a simple version of the figure

  np=3
  
  xs_img=3.
  ys_img=3.
  
  l_margin=0.8
  r_margin=1.1
  b_margin=0.6
  t_margin=0.4
  
  xs_ps=(xs_img+l_margin)*float(np)+r_margin
  ys_ps=ys_img+b_margin+t_margin
  
  l_img_arr=l_margin+findgen(np)*(xs_img+l_margin)
  r_img_arr=l_img_arr+xs_img
  b_img=b_margin
  t_img=b_img+ys_img
    
  ys_bar=ys_img
  xs_bar=0.15
  l_bar=max(r_img_arr)+0.2
  r_bar=l_bar+xs_bar
  b_bar=b_img
  t_bar=t_img
  pos_bar=[l_bar/xs_ps,b_bar/ys_ps,r_bar/xs_ps,t_bar/ys_ps]

  chisq_tickv=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,$
               1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000]
  chisq_tickname=['0.01',' ',' ',' ','0.05',' ',' ',' ',' ','0.1',' ',' ',' ','0.5',' ',' ',' ',' ',$
                  '1',' ',' ',' ','5',' ',' ',' ',' ','10',' ',' ',' ','50',' ',' ',' ',' ','100',' ',' ',' ','500',' ',' ',' ',' ','1000']
  
  set_plot,'PS'
  device,xsize=xs_ps,ysize=ys_ps,/inches,/encapsulate,$
         /portrait,xoffset=0.,yoffset=0.,bits=24,filename=fileout2,/color
  
  for ip=1,np do begin
     
     if ip eq 1 then begin
        xsuffix='mc'
        ysuffix='sigma'
        xtitle=textoidl('M_c (M_\odot)')
        ytitle=textoidl('\Sigma_{cl} (g cm^{-2})')
     endif else if ip eq 2 then begin
        xsuffix='mc'
        ysuffix='ms'
        xtitle=textoidl('M_c (M_\odot)')
        ytitle=textoidl('m_* (M_\odot)')
     endif else if ip eq 3 then begin
        xsuffix='sigma'
        ysuffix='ms'
        xtitle=textoidl('\Sigma_{cl} (g cm^{-2})')
        ytitle=textoidl('m_* (M_\odot)')
     endif
     
     l_img=l_img_arr(ip-1)
     r_img=r_img_arr(ip-1)
     pos_img=[l_img/xs_ps,b_img/ys_ps,r_img/xs_ps,t_img/ys_ps]

     zsuffix='chisq'
     ztitle=textoidl('\chi^2')
     z_rag=chisq_rag
     ifctreverse=1
     iflog=1
     bararr=reverse(bindgen(256))
     a=where(chisq_tickv ge chisq_rag(0) and chisq_tickv le chisq_rag(1))
     z_tickname=chisq_tickname(a)
     z_tickv=chisq_tickv(a)
     ctable=33
        
     void=execute('z_grid_arr='+zsuffix+'_'+xsuffix+'_'+ysuffix+'_grid_arr')
     void=execute('z_grid_arr_new='+zsuffix+'_'+xsuffix+'_'+ysuffix+'_grid_arr_new')
     void=execute('x_grid_arr='+xsuffix+'_grid_arr')
     void=execute('y_grid_arr='+ysuffix+'_grid_arr')
     void=execute('xrag='+xsuffix+'_rag')
     void=execute('yrag='+ysuffix+'_rag')
     void=execute('xtickv='+xsuffix+'_tickv')
     void=execute('ytickv='+ysuffix+'_tickv')
     void=execute('xtickname='+xsuffix+'_tickname')
     void=execute('ytickname='+ysuffix+'_tickname')
     void=execute('x_min_arr='+xsuffix+'_min_arr')
     void=execute('y_min_arr='+ysuffix+'_min_arr')
        
     z_grid_arr1=z_grid_arr
     z_max=max(z_grid_arr1)
     a=where(z_grid_arr1 lt -0.9)
     if a(0) ge 0 then begin
        z_grid_arr1(a)=z_max+100.
     endif
     a=where(z_grid_arr1 lt -0.4)
     if a(0) ge 0 then begin
        z_grid_arr1(a)=z_max+10.
     endif

     nx=n_elements(x_grid_arr)
     ny=n_elements(y_grid_arr)
        
     if ifctreverse eq 1 then begin
        z_grid_arr_new=255-z_grid_arr_new
     endif
        
     loadct,ctable,/silent
     tvlct,red,green,blue,/get
     redmap=bytarr(nx,ny)
     greenmap=bytarr(nx,ny)
     bluemap=bytarr(nx,ny)
     for ix=1,nx do begin
        for iy=1,ny do begin
           redmap(ix-1,iy-1)=red(z_grid_arr_new(ix-1,iy-1))
           greenmap(ix-1,iy-1)=green(z_grid_arr_new(ix-1,iy-1))
           bluemap(ix-1,iy-1)=blue(z_grid_arr_new(ix-1,iy-1))
        endfor
     endfor
     a=where(z_grid_arr lt -0.9)
     if a(0) ge 0 then begin
        redmap(a)=120
        greenmap(a)=120
        bluemap(a)=120
     endif
     a=where(z_grid_arr lt -0.4 and z_grid_arr gt -0.9)
     if a(0) ge 0 then begin
        redmap(a)=255
        greenmap(a)=255
        bluemap(a)=255
     endif
     loadct,0,/silent
     tv,[[[redmap]],[[greenmap]],[[bluemap]]],l_img,b_img,xsize=xs_img,ysize=ys_img,/inches,true=3
     
     xtit=xtitle
     xcharsize=1.2
     ytit=ytitle
     ycharsize=1.2

     loadcolors
     contour,z_grid_arr,x_grid_arr,y_grid_arr,position=pos_img,$
             /nodata,/noerase,xstyle=1,ystyle=1,$
             xrange=xrag,yrange=yrag,xtitle=xtit,ytitle=ytit,$
             xtickname=xtickname,xtickv=xtickv,xticks=n_elements(xtickv)-1,$
             ytickname=ytickname,ytickv=ytickv,yticks=n_elements(ytickv)-1,$
             charthick=3,xthick=3,ythick=3,xcharsize=xcharsize,ycharsize=ycharsize
     
     contour,z_grid_arr1,x_grid_arr,y_grid_arr,/overplot,levels=[chisq_min_arr(0)+5.],color=7,thick=10
     contour,z_grid_arr1,x_grid_arr,y_grid_arr,/overplot,levels=[chisq_min_arr(0)+5.],color=5,thick=4
     
     contour,z_grid_arr1,x_grid_arr,y_grid_arr,/overplot,levels=[z_max+50.],color=0,thick=5
     
     for ibest=1,nbest do begin
        if ibest eq 1 then begin
           symsize=3
        endif else begin
           symsize=2
        endelse
        plots,x_min_arr(ibest-1),y_min_arr(ibest-1),psym=1,symsize=symsize,thick=5,color=7,/data
     endfor
     
  endfor
  
  loadct,33,/silent
  bar=reverse(bytscl(newgrid(mu_arr,mu_arr,mu_grid_arr),min=mu_grid_min,max=mu_grid_max))#replicate(1b,10)
  bar=transpose(bar)
  tv,bar,l_bar,b_bar,xsize=xs_bar,ysize=ys_bar,/inches
  
  loadcolors
  contour,bar,position=pos_bar,color=0,/noerase,/nodata,xstyle=4,ystyle=4
  axis,xstyle=4,ystyle=1,yaxis=1,yrange=z_rag,ytitle=ztitle,ylog=iflog,$
       ythick=3,yticklen=0.3,charthick=3,charsize=1.2,$
       ytickname=z_tickname,ytickv=z_tickv,yticks=n_elements(z_tickv)-1
  axis,xstyle=4,ystyle=1,yaxis=0,yrange=z_rag,ytitle=ztitle,ylog=iflog,$
       ythick=3,yticklen=0.3,charthick=3,charsize=0.001,$
       ytickname=z_tickname,ytickv=z_tickv,yticks=n_elements(z_tickv)-1
  
  xyouts,(l_margin+xs_ps-r_margin)/2./xs_ps,0.92,name,alignment=0.5,charthick=3,charsize=1.2,/normal
  
  device,/close

;;===============================================================


  return
  
end
