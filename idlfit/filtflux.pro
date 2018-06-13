
;; from model SEDs to filter fluxes

;;===============================================================
;; Input

filt_name='X2'; name of the filter
filt_wav=7.5; the wavelength of the filter, in micron
filt_width=1.; the width of the filter, in micron. If ifiltfile=1, this is not used.
iffiltfile=0 ; 1 to use filter files in the parfiles directory, 0 to use filt_wav and filt_width to define a square filter.
; If iffiltfile=0, a filter file named 'filt_name.txt" will
; also be added to parfiles directory.

;;===============================================================

dirin='../Model_SEDs/sed/'
dir_filt='../Model_SEDs/parfiles/'
dirout='../Model_SEDs/flux_filt/'

nmc=15
mc_str_arr=['10','20','30','40','50','60','80','100','120','160','200','240','320','400','480']
mc_arr=[10.,20.,30.,40.,50.,60.,80.,100.,120.,160.,200.,240.,320.,400.,480.]

nsigma=4
sigma_str_arr=['0.1','0.316','1','3.16']
sigma_arr=[0.1,0.316,1.,3.16]

nms=14
ms_str_arr=['0.5','1.0','2.0','4.0','8.0','12.0','16.0','24.0','32.0','48.0','64.0','96.0','128.0','160.0']
ms_arr=[0.5,1.,2.,4.,8.,12.,16.,24.,32.,48.,64.,96.,128.,160.]

nmu=20

if iffiltfile eq 0 then begin
   if filt_width le 0. then begin
      print,'filt_width wrong'
      stop
   endif
   dwave=filt_width/100.
   fwave=reverse((findgen(110)-55.)*dwave+filt_wav)
   fresponse=fltarr(110)
   a=where(abs(fwave-filt_wav) lt filt_width/2.)
   fresponse(a)=1.
   if (file_test(dir_filt+filt_name+'.txt') eq 1) then begin
      print,'filter file named as '+filt_name+'.txt already exists.'
      stop
   endif else begin
      openw,1,dir_filt+filt_name+'.txt'
      printf,1,'#110 '+filt_name+' '+strtrim(string(filt_wav),2)
      for ii=1,110 do begin
         printf,1,fwave(ii-1),fresponse(ii-1)
      endfor
      close,1
   endelse
endif

readcol,dir_filt+filt_name+'.txt',fwave,fresponse,skipline=1,/silent
fnu=3.e14/fwave
nf=n_elements(fnu)

if fnu(0) gt fnu(1) then begin
   fnu=reverse(fnu)
   fresponse=reverse(fresponse)
endif

dfnu=fnu(1:nf-1)-fnu(0:nf-2)
fint=total(0.5*(fresponse(0:nf-2)+fresponse(1:nf-1))*dfnu)
fresponse=fresponse(1:nf-1)/fint
fwave=fwave(1:nf-1)

if (file_test(dirout,/directory) eq 0) then file_mkdir,dirout

filtdata=fltarr(nmc,nsigma,nms,nmu)
   
for imc=1,nmc do begin
   for isigma=1,nsigma do begin   
      for ims=1,nms do begin
         for imu=1,nmu do begin
            
            suffix=strtrim(string(imc,'(I02)'),2)+'_'+$
                   strtrim(string(isigma,'(I02)'),2)+'_'+$
                   strtrim(string(ims,'(I02)'),2)+'_'+$
                   strtrim(string(imu,'(I02)'),2)
            filename=suffix+'.dat'
            if (file_test(dirin+filename) eq 1) then begin
               
               readcol,dirin+filename,lambda_arr,I_arr,/silent

               I_arr1=interpol(I_arr,lambda_arr,fwave)
               I_arr1=I_arr1/3.e14*fwave
               I_filt=total(I_arr1*dfnu*fresponse)
               I_filt=I_filt*3.e14/filt_wav
                  
               filtdata(imc-1,isigma-1,ims-1,imu-1)=I_filt
                 
            endif
         endfor
      endfor
   endfor
endfor

if file_test(dirout+filt_name+'.fits') eq 0 then begin
   mkhdr,header,filtdata
   sxaddpar,header,'FWAVE',filt_wav
   writefits,dirout+filt_name+'.fits',filtdata,header
endif else begin
   print,filt_name+'.fits already exists.'
   stop
endelse

end
