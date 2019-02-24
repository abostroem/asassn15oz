pro template_match_visit2

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
  ;PA 260
  ;COADD THE 2 TEMPLATE SNAPSHOTS INTO A SINGLE, HIGHER SIGNAL-TO-NOISE
  ;TEMPLATE
  ;-----------------------
  base_template_str2 = '/Users/bostroem/Desktop/research/asassn15oz/data/swiftuvot/test_truvot_extraction/truvot/00034040015/uvot/image/sw00034040015ugu_dt'
  base_obs_str1 = '/Users/bostroem/Desktop/research/asassn15oz/data/swiftuvot/test_truvot_extraction/truvot/00034040002/uvot/image/sw00034040002ugu_dt'

  ;Visit 15
  x1=1215.
  y1=689.
  
  ;COADD first and second extensions
  print, 'COADD first and second extensions visit 15'
  ext =1
  skycounts=11.
  backfile=base_template_str2+'.img'
  backfileoutput=base_template_str2+'_coadd1.fits'
  datafile =base_template_str2+'.img'
  x2=1207.
  y2=687.
  ext2=2
  optimal_shift_time,ext,skycounts,backfile,backfileoutput,datafile,x1,y1,x2,y2,/coadd,ext2=ext2

  ;COADD third extension to first and second
  print, 'COADD third extension to first and second visit 15'
  ext =0
  skycounts=11.
  backfile=base_template_str2+'.img'
  backfileoutput=base_template_str2+'_coadd2.fits'
  datafile =base_template_str2+'_coadd1.fits'
  x2=1224.
  y2=690.
  ext2=3
  optimal_shift_time,ext,skycounts,backfile,backfileoutput,datafile,x1,y1,x2,y2,/coadd,ext2=ext2  
  
  ;COADD fourth extension to first three
  print, 'COADD fourth extension to first three visit 15'
  ext =0
  skycounts=12.
  backfile=base_template_str2+'.img'
  backfileoutput=base_template_str2+'_coadd3.fits'
  datafile =base_template_str2+'_coadd2.fits'
  x2=1212.
  y2=747.
  ext2=4
  optimal_shift_time,ext,skycounts,backfile,backfileoutput,datafile,x1,y1,x2,y2,/coadd,ext2=ext2
  
  ;COADD fifth extension to first four
  print, 'COADD fifth extension to first four visit 15'
  ext =0
  skycounts=12.
  backfile=base_template_str2+'.img'
  backfileoutput=base_template_str2+'_coadd4.fits'
  datafile =base_template_str2+'_coadd3.fits'
  x2=1221.
  y2=716.
  ext2=5
  optimal_shift_time,ext,skycounts,backfile,backfileoutput,datafile,x1,y1,x2,y2,/coadd,ext2=ext2
  
  ;----------------------
  ; Visit 2
  ;----------------------
 

  ;TRANSLATE THE TEMPLATE IMAGE TO OVERLAY THE DATA IMAGE AND FLUX SCALE
  backfile=base_template_str2+'_coadd4.fits'
  datafile =base_obs_str1+'.img'
  x2=1440.
  y2=602.
  ;Scale the first extension
  ext =1
  skycounts=23.
  backfileoutput=base_obs_str1+'_template_1.fits'
  subfile = base_obs_str1+'_subtracted_1.fits'
  x1=1455.
  y1=795.
  optimal_shift_time,ext,skycounts,backfile,backfileoutput ,datafile,x1,y1, x2,y2,outputsubtract=subfile

  ;Scale the second extension
  ext =2
  skycounts=25.
  backfileoutput=base_obs_str1+'_template_2.fits'
  subfile = base_obs_str1+'_subtracted_2.fits'
  x1=1262.
  y1=565.
  optimal_shift_time,ext,skycounts,backfile,backfileoutput ,datafile,x1,y1, x2,y2,outputsubtract=subfile

  ;Scale the third extension
  ext=3
  skycounts=2.
  backfileoutput=base_obs_str1+'_template_3.fits'
  subfile=base_obs_str1+'_subtracted_3.fits'
  x1=1259.
  y1=573.
  optimal_shift_time,ext,skycounts,backfile, backfileoutput,datafile,x1,y1, x2,y2,outputsubtract=subfile



end



