
function newgrid,x_arr,data_arr,x_grid_arr,missing=missing

  data_arr1=data_arr
  
  nx=n_elements(x_arr)
  nxgrid=n_elements(x_grid_arr)

  x0=x_arr(0)
  x1=x_arr(nx-1)
  if x1 lt x0 then begin
     x_arr=reverse(x_arr)
     data_arr1=reverse(data_arr1)
  endif
  data_arr=data_arr1
  
  data_grid_arr=fltarr(nxgrid)
  
  for ixgrid=1,nxgrid do begin

     x_grid=x_grid_arr(ixgrid-1)
     a=where(x_arr gt x_grid)
     if a(0) eq 0 then begin
        if x_arr(0)-x_grid le (x_arr(1)-x_arr(0))/2.*1.00001 then begin
           indx=0
        endif else begin
           indx=-1
        endelse
     endif else if a(0) eq -1 then begin
        if x_grid-x_arr(nx-1) le (x_arr(nx-1)-x_arr(nx-2))/2.*1.00001 then begin
           indx=nx-1
        endif else begin
           indx=-1
        endelse        
     endif else begin
        indl=a(0)-1
        indh=a(0)
        xl=x_arr(indl)
        xh=x_arr(indh)
        if x_grid-xl le xh-x_grid then begin
           indx=indl
        endif else begin
           indx=indh
        endelse
     endelse
     
     if indx ge 0 then begin
        data_grid_arr(ixgrid-1)=data_arr(indx)
     endif else begin
        data_grid_arr(ixgrid-1)=missing
     endelse
        
  endfor
  
  return,data_grid_arr
  
end
