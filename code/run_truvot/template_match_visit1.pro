pro template_match_visit1

  ;optimal_shift_time,(operate on 1st image extension=1),(sky background
  ;estimate=50),(filename of template image),(output image),(filename
  ;of template image again),(x position of zeroth order feature in 1st
  ;image extension),(y position of zeroth order feature in 1st image extension),(x position of zeroth order feature in 2nd
  ;image extension),(y position of zeroth order feature in 2nd image
  ;extension),(/coadd flag indicates that you wish to coadd image
  ;extensions 1 and 2),(ext2=2 indicates that you wish to coadd the 2nd
  ;image extension onto the first image extension)


  ;optimal_shift_time,(operate on the third image extension),(sky
  ;background estimate=40),(template image),(output image),(data
  ;image),(x position of zeroth order feature in data image),(y position
  ;of zeroth order feature in data image),(x position of zeroth order
  ;feature in template image),(y position of zeroth order feature in
  ;template image),(output a residual image of data image minus template image)
  

;-----------------------
;PA 248
;-----------------------
;COADD THE 2 TEMPLATE SNAPSHOTS INTO A SINGLE, HIGHER SIGNAL-TO-NOISE
;TEMPLATE
base_template_str = '/Users/bostroem/Desktop/research/asassn15oz/data/swiftuvot/test_truvot_extraction/truvot/00034040013/uvot/image/sw00034040013ugu_dt'
;COADD first and second extensions
ext =1
skycounts=15
backfile=base_template_str+'.img'
backfileoutput=base_template_str+'_coadd1.fits'
datafile =base_template_str+'.img'
x1=1043.
y1=554.
x2=1030.
y2=594
ext2=2
optimal_shift_time,ext,skycounts,backfile,backfileoutput,datafile,x1,y1,x2,y2,/coadd,ext2=ext2

;COADD third extension to first and second
ext =0
skycounts=18
backfile=base_template_str+'.img'
backfileoutput=base_template_str+'_coadd2.fits'
datafile =base_template_str+'_coadd1.fits'
x1=1043.
y1=554.
x2=1033.
y2=583.
ext2=3
optimal_shift_time,ext,skycounts,backfile,backfileoutput,datafile,x1,y1,x2,y2,/coadd,ext2=ext2

;----------------------
; Visit 1
;----------------------
base_obs_str1 = '/Users/bostroem/Desktop/research/asassn15oz/data/swiftuvot/test_truvot_extraction/truvot/00034040001/uvot/image/sw00034040001ugu_dt'

;TRANSLATE THE TEMPLATE IMAGE TO OVERLAY THE DATA IMAGE AND FLUX SCALE
backfile=base_template_str+'_coadd2.fits'
datafile =base_obs_str1+'.img'
x2=1476.
y2=602.
;Scale the first extension
ext =1
skycounts=3.
backfileoutput=base_obs_str1+'_template_1.fits'
subfile = base_obs_str1+'_subtracted_1.fits'
x1=1285.
y1=698.
optimal_shift_time,ext,skycounts,backfile,backfileoutput,datafile,x1,y1, x2,y2,outputsubtract=subfile

;Scale the second extension
ext =2
skycounts=25.
backfileoutput=base_obs_str1+'_template_2.fits'
subfile = base_obs_str1+'_subtracted_2.fits'
x1=1281.
y1=681.
optimal_shift_time,ext,skycounts,backfile,backfileoutput,datafile,x1,y1, x2,y2,outputsubtract=subfile

;Scale the third extension
ext=3
skycounts=29.
backfileoutput=base_obs_str1+'_template_3.fits'
subfile=base_obs_str1+'_subtracted_3.fits'
x1=1294.
y1=752.
optimal_shift_time,ext,skycounts,backfile, backfileoutput,datafile,x1,y1, x2,y2,outputsubtract=subfile



end



