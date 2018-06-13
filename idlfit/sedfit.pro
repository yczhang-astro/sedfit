
;;===============================================================
;; Edit the input here

sourcename='example'
filt_arr=['I1','I2','I3','I4','F8','L1','L4','P1','P3','P5','P6']
filt_wav_arr=[3.6, 4.5, 5.8, 8., 19.7, 31.5, 37.1, 70, 160, 350, 500]
flux_arr=[0.497, 1.24, 1.84, 3.22, 68.18, 552.8, 1192.5, 2627.6, 2385.6, 428.6, 127.4]
errup=[0.216, 0.1, 0.434, 0.522, 0.18, 0.1, 0.1, 0.1, 0.176, 0.387, 0.536]; relative errors
errlo=errup
limit=[1,1,0,0,0,0,0,0,0,1,1]

errop=2 ; 1 to convert flux and error to log, 2 to convert with different upper and lower errors, 3 same as 2 but treat >=100% lower error specially. Single number or an array.
avop=1 ; 1 if Av is correlated to the sigma parater, 0 if not
ndis=1; number of distance grid, set to 1 to use exact distance
dis=2200.                      
dis_err=0.; error on distance if disop=1
thresop=1 ; output specified number of models
nmodel=10
nmodel_full=20
chisq_thres=3.

xfilt_arr=['L1', 'P1', 'F8', 'L4', 'I4', 'P1']
yfilt_arr=['L4', 'P3', 'L4', 'P1', 'L4', 'P5']

;; End of input section
;;===============================================================








;;===============================================================
;; prepare the subroutines

!PATH=EXPAND_PATH('+lib')+':'+!PATH
RESOLVE_ROUTINE,'chisq'
RESOLVE_ROUTINE,'bestmodel'
RESOLVE_ROUTINE,'plotsed'
RESOLVE_ROUTINE,'plotsed1'
RESOLVE_ROUTINE,'plotchisq'
RESOLVE_ROUTINE,'plotchisq1'
RESOLVE_ROUTINE,'plotchisq2'

;;===============================================================
;; calculate the chisq for each model and output the results
chisq,sourcename,filt_arr,filt_wav_arr,flux_arr,errup,errlo,$
      ERROPTION=errop,AVOPTION=avop,$
      LIMIT=limit,NDIS=ndis,DISTANCE=dis,DISERR=dis_err

;; From the output of chisq program, select the best models and output
;; the parameters
bestmodel,sourcename,THRESOPTION=thresop,NMODEL=nmodel,NFULLMODEL=nmodel_full,CHISQ_THRES=chisq_thres

;; Plot the SEDs of the best models
plotsed,sourcename

;;===============================================================
;; Additional plots for analyzing the results

;; Plot the SEDs with a wider range of chisq
plotsed1,sourcename,CHISQ_MAX=50

;; Plot the chisq distribution in the parameter space
;; For the primary parameters (Mc, Sigma, ms, and mu)
plotchisq,sourcename,CHISQ_MAX=50
;; For the secondary parameters
plotchisq1,sourcename,CHISQ_MAX=50

;; Plot the chisq distribution in the color-flux diagrams
plotchisq2,sourcename,xfilt_arr,yfilt_arr,CHISQ_MAX=50,NXP=3,NYP=2

;;===============================================================

end
