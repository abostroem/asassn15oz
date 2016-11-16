PESSTOEFOSC2dSPEC -i -C
cp tASASSN-15oz_20160918_Gr13_Free_slit1.5_57707_2.fits host_tASASSN-15oz_20160918_Gr13_Free_slit1.5_57707_2.fits
cp tASASSN-15oz_20160918_Gr13_Free_slit1.5_57707_1.fits host_tASASSN-15oz_20160918_Gr13_Free_slit1.5_57707_1.fits
ls host*.fits > list_hosts.txt #Make a list of hosts to run first to get the trace
ls t*.fits > list2d.txt #Make a list of all 2D observations
PESSTOEFOSC1dSPEC -i list_hosts.txt -d
PESSTOEFOSC1dSPEC -i list2d.txt -d -t
