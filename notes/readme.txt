print (icl.hdus(file='ALPAQL_Light_4,000secs_00000195.fit'))

print ()
for hdu in icl.hdus(file='ALPAQL_Light_4,000secs_00000195.fit'):
    print (repr(hdu.header))

print ()
for hdu in icl.hdus(*{'date-obs':'2019-09-20T21:06:35.659'}):
    print (repr(hdu.header))
    hdu.header.set('TEST','This is a test')
    print (repr(hdu.header))

print ()
for hdu in icl.hdus(*{'date-obs':'2019-09-20T21:06:35.659'}):
    print (repr(hdu.header))

print (icl.summary.columns)
<TableColumns names=('file','simple','bitpix','naxis','naxis1','naxis2','extend','bzero','bscale','observer','origin','telescop','focallen','aptarea','aptdia','sbuuid','exptime','swcreate','colorccd','dispincr','picttype','imagetyp','xorgsubf','yorgsubf','xbinning','ybinning','ccd-temp','set-temp','instrume','xpixsz','ypixsz','date-obs','localtim','comment','objctra','objctdec','object')>