
function newgrid2d,x_arr,y_arr,data_arr,x_grid_arr,y_grid_arr,missing=missing

  nx=n_elements(x_arr)
  ny=n_elements(y_arr)

  nxgrid=n_elements(x_grid_arr)
  nygrid=n_elements(y_grid_arr)

  x0=x_arr(0)
  x1=x_arr(nx-1)
  if x1 lt x0 then begin
     x_arr=reverse(x_arr)
     data_arr=reverse(data_arr,1)
  endif
  y0=y_arr(0)
  y1=y_arr(ny-1)
  if y1 lt y0 then begin
     y_arr=reverse(y_arr)
     data_arr=reverse(data_arr,2)
  endif
  
  data_grid_arr=fltarr(nxgrid,nygrid)
  
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
        
     for iygrid=1,nygrid do begin

        y_grid=y_grid_arr(iygrid-1)
        a=where(y_arr gt y_grid)
        if a(0) eq 0 then begin
           if y_arr(0)-y_grid le (y_arr(1)-y_arr(0))/2.*1.00001 then begin
              indy=0
           endif else begin
              indy=-1
           endelse
        endif else if a(0) eq -1 then begin
           if y_grid-y_arr(ny-1) le (y_arr(ny-1)-y_arr(ny-2))/2.*1.00001 then begin
              indy=ny-1
           endif else begin
              indy=-1
           endelse        
        endif else begin
           indl=a(0)-1
           indh=a(0)
           yl=y_arr(indl)
           yh=y_arr(indh)
           if y_grid-yl le yh-y_grid then begin
              indy=indl
           endif else begin
              indy=indh
           endelse
        endelse

        if indx ge 0 and indy ge 0 then begin
           data_grid_arr(ixgrid-1,iygrid-1)=data_arr(indx,indy)
        endif else begin
           data_grid_arr(ixgrid-1,iygrid-1)=missing
        endelse
        
     endfor
  endfor
  
  return,data_grid_arr

end
