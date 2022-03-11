from .globals import *

def load_kernels(start_utc, KERNEL_DATAFOLDER):
    ## find and load the kernels for a specific date 
    iks   = sorted(glob.glob(KERNEL_DATAFOLDER+"ik/juno_jiram_v*.ti"))
    cks   = sorted(glob.glob(KERNEL_DATAFOLDER+"ck/juno_sc_rec_*.bc"))
    spks1 = sorted(glob.glob(KERNEL_DATAFOLDER+"spk/spk_rec_*.bsp"))
    spks2 = sorted(glob.glob(KERNEL_DATAFOLDER+"spk/jup*.bsp"))
    spks3 = sorted(glob.glob(KERNEL_DATAFOLDER+"spk/de*.bsp"))
    pcks  = sorted(glob.glob(KERNEL_DATAFOLDER+"pck/pck*.tpc"))
    fks   = sorted(glob.glob(KERNEL_DATAFOLDER+"fk/juno_v*.tf"))
    sclks = sorted(glob.glob(KERNEL_DATAFOLDER+"sclk/JNO_SCLKSCET.*.tsc"))
    lsks  = sorted(glob.glob(KERNEL_DATAFOLDER+"lsk/naif*.tls"))

    year, month, day = start_utc.split('-')
    yy = year[2:]
    mm = month
    dd = day[:2]

    intdate = int("%s%s%s"%(yy,mm,dd))
    kernels = []

    ## find the ck and spk kernels for the given date 
    ckpattern = r'juno_sc_rec_([0-9]{6})_([0-9]{6})\S*'
    nck = 0
    for ck in cks:
        fname = os.path.basename(ck)
        groups = re.findall(ckpattern, fname)
        if(len(groups) == 0):
            continue
        datestart, dateend = groups[0]

        if( (int(datestart) <= intdate) & (int(dateend) >= intdate) ):
            kernels.append(ck)
            nck += 1
    
    ''' use the predicted kernels if there are no rec '''
    if(nck == 0):
        ckpattern = r'juno_sc_pre_([0-9]{6})_([0-9]{6})\S*'
        for ck in cks:
            fname = os.path.basename(ck)
            groups = re.findall(ckpattern, fname)
            if(len(groups) == 0):
                continue
            datestart, dateend = groups[0]

            if( (int(datestart) <= intdate) & (int(dateend) >= intdate) ):
                kernels.append(ck)
                nck += 1

    spkpattern = r'spk_rec_([0-9]{6})_([0-9]{6})\S*'
    nspk = 0
    for spk in spks1:
        fname = os.path.basename(spk)
        groups = re.findall(spkpattern, fname)
        if(len(groups) == 0):
            continue
        datestart, dateend = groups[0]

        if( (int(datestart) <= intdate) & (int(dateend) >= intdate) ):
            kernels.append(spk)
            nspk += 1

    ''' use the predicted kernels if there are no rec '''
    if(nspk == 0):
        spkpattern = r'spk_pre_([0-9]{6})_([0-9]{6})\S*'
        for spk in spks1:
            fname = os.path.basename(spk)
            groups = re.findall(spkpattern, fname)
            if(len(groups) == 0):
                continue
            datestart, dateend = groups[0]

            if( (int(datestart) <= intdate) & (int(dateend) >= intdate) ):
                kernels.append(spk)
                nspk += 1

    if(nck*nspk == 0):
        print("ERROR: Kernels not found for the date range!")

    ## load the latest updates for these 
    kernels.append(iks[-1])
    kernels.append(spks2[-1])
    kernels.append(spks3[-1])
    kernels.append(pcks[-1])
    kernels.append(fks[-1])
    kernels.append(sclks[-1])
    kernels.append(lsks[-1])

    nloaded = spice.ktotal('ALL')
    loaded_kernels = [spice.kdata(i,'ALL')[0] for i in range(nloaded)]

    for kernel in kernels:
        ''' load the kernel if it has not been loaded already '''
        if(kernel not in loaded_kernels):
            furnish_c(kernel.encode('ascii'))
            spice.furnsh(kernel)
    return kernels

class Projector():
    def __init__(self, IMG, LBL, KERNEL_DATAFOLDER):
        metafile      = open(LBL, 'r')
        self.metadata = self.load_LBL(LBL)

        self.KERNEL_DATAFOLDER = KERNEL_DATAFOLDER

        ''' load the image (which is either 256x432 or 128x432 at 32-bit depth) '''
        if 'S1' in self.mode:
            self.kernels = load_kernels(self.start_utc, self.KERNEL_DATAFOLDER)

            ''' read the IMG file '''
            IMGfile = open(IMG, "rb")

            # I1 mode is both L band and M band
            if 'I1' in self.mode:
                nbytes = IMG_HEIGHT*FRAME_WIDTH*4*2
                buffer  = IMGfile.read(nbytes)
                fullimg = np.frombuffer(buffer, dtype=np.uint8, count=nbytes)
                fullimg = np.frombuffer(fullimg, dtype=np.float32, count=IMG_HEIGHT*FRAME_WIDTH*2).reshape((IMG_HEIGHT*2, FRAME_WIDTH))
                self.Lband = fullimg[:128,:]
                self.Mband = fullimg[128:,:]
            if 'I2' in self.mode: # M band only
                nbytes = IMG_HEIGHT*FRAME_WIDTH*4
                buffer  = IMGfile.read(nbytes)
                fullimg = np.frombuffer(buffer, dtype=np.uint8, count=nbytes)
                self.Mband = np.frombuffer(fullimg, dtype=np.float32, count=IMG_HEIGHT*FRAME_WIDTH).reshape((IMG_HEIGHT, FRAME_WIDTH))
                self.Lband = np.zeros_like(self.Mband)
            if 'I3' in self.mode: # L band only
                nbytes = IMG_HEIGHT*FRAME_WIDTH*4
                buffer  = IMGfile.read(nbytes)
                fullimg = np.frombuffer(buffer, dtype=np.uint8, count=nbytes)
                self.Lband = np.frombuffer(fullimg, dtype=np.float32, count=IMG_HEIGHT*FRAME_WIDTH).reshape((IMG_HEIGHT, FRAME_WIDTH))
                self.Mband = np.zeros_like(self.Lband)

            self.re, _, self.rp = spice.bodvar(spice.bodn2c('JUPITER'), 'RADII', 3)
            self.flattening = (self.re - self.rp)/self.re

            ## calculate the start time 
            self.start_et    = spice.str2et(self.start_utc)
            self.end_et    = spice.str2et(self.end_utc)

            self.init = True
        else:
            self.init = False

        metafile.close()

    def load_LBL(self, LBLfile):
        #time_pattern    = r"([0-9]{4})-([0-9]{2})-([0-9]{2})T([0-9]{2}):([0-9]{2}):([0-9.]{6})"
        time_pattern    = r"= ([0-9.:T-]+)"
        id_pattern     = r"\s*= ([A-Z0-9_]+)"
        target_pattern  = r"= \"([0-9]+)\""

        LBLlines = open(LBLfile, "r").readlines()

        for line in LBLlines:
            if(re.search(r"^START_TIME",line)):
                match = re.findall(time_pattern, line)
                self.start_utc = match[0]

            if(re.search(r"^STOP_TIME",line)):
                match = re.findall(time_pattern, line)
                self.end_utc = match[0]

            if(re.search(r"^PRODUCT_ID",line)):
                match = re.findall(id_pattern, line)
                self.fname = match[0]
            
            if(re.search(r"^INSTRUMENT_MODE_ID",line)):
                match = re.findall(id_pattern, line)
                self.mode = match[0]

            if(re.search(r"^JNO:TARGET_PRESENCE_FLAG", line)):
                match = re.findall(target_pattern, line)
                try:
                    target_flag = match[0]
                    self.target_present   = int(target_flag[0]) == 1
                    self.target_present_L = int(target_flag[5]) == 1
                    self.target_present_M = int(target_flag[10]) == 1
                except:
                    self.target_present   = False
                    self.target_present_L = False
                    self.target_present_M = False

        

    def process(self, band='M', ncfolder="nc/", pngfolder="png/"):
        if not os.path.exists(ncfolder):
            os.mkdir(ncfolder)
        if not os.path.exists(pngfolder):
            os.mkdir(pngfolder)

        if(self.init == False):
            raise RuntimeError("This target is not in high image resolution mode!")


        if band=='M':
            image = self.Mband
            band_number = 0
            if not self.target_present_M:
                raise RuntimeError(f"Target is not present in {band} band")
        elif band=='L':
            image = self.Lband
            band_number = 1
            if not self.target_present_M:
                raise RuntimeError(f"Target is not present in {band} band")
        else:
            raise KeyError('band must be either M or L')

        eti   = self.start_et#0.5*(self.start_et + self.end_et)
        
        # initialize arrays to hold the data and pass into the 
        # C function
        lats     = -1000.*np.ones((FRAME_HEIGHT, FRAME_WIDTH))
        lons     = -1000.*np.ones((FRAME_HEIGHT, FRAME_WIDTH))
        phase    = np.ones((FRAME_HEIGHT, FRAME_WIDTH))
        inc      = np.ones((FRAME_HEIGHT, FRAME_WIDTH))
        emission = np.ones((FRAME_HEIGHT, FRAME_WIDTH))
        scloc = np.zeros(3)

        # call the C function to process each pixel into a 
        # map of lat/lon. also get the illumination angles
        process_c(eti, band_number, scloc, lats, lons, phase, inc, emission)

        '''
        fig, ax = plt.subplots(1, 1, figsize=(10,3), dpi=150)

        ax.imshow(image, cmap='hot')
        ax.contour(lons, np.arange(-180, 180, 30), colors='white', linewidths=0.1) 
        ax.contour(lats, np.arange(-90, 90, 30), colors='white', linewidths=0.1) 
        ax.axis('off')
        fig.tight_layout()
        plt.savefig("%s/%s.png"%(pngfolder,self.fname))
        #plt.close()
        #plt.show()
        '''
        return lons, lats, image, scloc, phase, inc, emission

    def project(self, x, y, cam2jup, eti):
        '''
            projects a single pixel in the filter given by 
            the cam object for a given spacecraft location
            and et
        '''
        
        xyvec  = np.array([CENTER_Y-y, CENTER_X-x, F1])
        xyvec  = xyvec/np.linalg.norm(xyvec)


        ## get the vector in the Jupiter frame
        pos_jup = np.matmul(cam2jup, xyvec)
        if y==0:
            print(x, y, xyvec, pos_jup)
        
        try:           
            point, _, srfvec = spice.sincpt("Ellipsoid", "JUPITER", eti, "IAU_JUPITER", "CN+S", "JUNO", "IAU_JUPITER", pos_jup)
            #point = spice.surfpt(scloc, pos_jup, self.re, self.re, self.rp)
            #print(xyvec, xyvec2)
        except: 
            return (-1000., -1000.)
        dist, loni, lati = spice.reclat(point)

        
        return (np.degrees(loni), np.degrees(lati))
    
def vec2pix(spoint, et, scloc, jup2cam):
    
    ''' 
        find the difference vector = point - sclocation
        in this case, scloc is actually position of jupiter 
        from camera so sclocation = -scloc
    '''
    point_vec  = spoint - scloc

    cam_vec = np.matmul(jup2cam, point_vec)
    
    cam_vec = cam_vec/np.linalg.norm(cam_vec)

    alpha = cam_vec[2]/F1
    xx    = CENTER_X - cam_vec[1]/alpha
    yy    = CENTER_Y - cam_vec[0]/alpha

    return (xx, yy)


def map_project_multi(file, KERNEL_DATAFOLDER, pixres=1./25.):


    with nc.Dataset(file, 'r') as dataset:
        lats = dataset.variables['lat'][:]
        lons = dataset.variables['lon'][:]
        imgs = dataset.variables['img'][:]
        
        scloc = dataset.variables['scloc'][:]
        et    = dataset.variables['et'][:]

        nfiles = len(et)

    jup2cam  = np.zeros((nfiles, 9))
    for i in range(nfiles):
        jup2cam[i,:] = spice.pxform('IAU_JUPITER', 'JUNO_JIRAM_I_MBAND', et[i]).flatten()


    latmin = np.min([lati[lati!=-1000].min() for lati in lats])
    latmax = np.max([lati[lati!=-1000].max() for lati in lats])
    lonmin = np.min([loni[loni!=-1000].min() for loni in lons])
    lonmax = np.max([loni[loni!=-1000].max() for loni in lons])

    print("Extents - lon: %.3f %.3f  lat: %.3f %.3f"%(lonmin, lonmax, latmin, latmax))

    newlon = np.arange(lonmin, lonmax, pixres)
    newlat = np.arange(latmin, latmax, pixres)

    LAT, LON = np.meshgrid(newlat, newlon)

    nlat = newlat.size
    nlon = newlon.size 
    print("Mosaic shape: %d x %d"%(nlon, nlat))
    

    print("Getting mask...")
    output = image_mask_c(np.radians(newlat), np.radians(newlon), nlat, nlon,\
                    scloc, jup2cam, et, nfiles)
    mask = ctypes.cast(output, ctypes.POINTER(ctypes.c_int*(nlat*nlon))).contents
    mask = np.asarray(mask, dtype=np.int).reshape((nlat, nlon))

    #lats = lats.flatten()
    #lons = lons.flatten()
    #imgs = imgs.flatten()

    invmask = np.where((lats==-1000.)|(lons==-1000.))[0]
    ## remove pixels that were not projected
    lat = np.delete(lats, invmask)
    lon = np.delete(lons, invmask)
    img = np.delete(imgs, invmask)
    
    fig, ax = plt.subplots(figsize=(10,5))

    ax.imshow(mask, cmap='gray', origin='lower', extent=(LON.min(), LON.max(), LAT.min(), LAT.max()))
    ax.plot(lon, lat, 'r.')

    plt.tight_layout()
    fig.savefig('mask.png')
    
    IMG  = np.zeros((nlat, nlon))
    NPIX = np.zeros((nlat, nlon))

    #IMG = griddata((lon, lat), img, (LON, LAT), method='linear').T
    print("Mosaicing...")
    for n in range(nfiles):
        invmask = np.where((lats[n,:]==-1000.)|(lons[n,:]==-1000.))[0]
        ## remove pixels that were not projected
        lat = np.delete(lats[n,:], invmask)
        lon = np.delete(lons[n,:], invmask)
        img = np.delete(imgs[n,:], invmask)

        print("\r %d/%d"%(n, nfiles), end='')
        #IMGI = map_project(LON, LAT, lons[n,:,:], lats[n,:,:], imgs[n,:,:])
        IMGI = griddata((lon, lat), img, (LON, LAT), method='cubic').T
        IMG[IMGI>0.01]  += IMGI[IMGI>0.01]
        NPIX[IMGI>0.01] += 1

    '''
    IMGI = griddata((lon, lat), img, (LON, LAT), method='linear').T
    IMGI[np.isnan(IMGI)]  = 0.
    IMGI[IMGI<0.] = 0.
    ## do a simple averaging to construct the mosaic
    for jj, latj in enumerate(newlat):
        print("\r %4d/%d"%(jj, newlat.size), end='')

        ## find the pixels that correspond to this latitude 
        masklat = np.abs(latj-lat)<=(pixres/2.)

        ## and subset the image/lon data 
        lonsub  = lon[masklat]
        IMGsub  = img[masklat]
        for ii, loni in enumerate(newlon):
            ## get the points from the image subset that correspond to this lon 
            masklon = np.abs(lonsub-loni)<=(pixres/2.)

            ## get the points from the images which correspond to both this lon and lat
            mx    =  IMGsub[masklon]
            if(len(mx) > 0):
                ## save the points into the image by averaging over all valid pixels 
                IMG[jj,ii] = np.average(mx)
    '''

    IMG[NPIX>0] = IMG[NPIX>0]/NPIX[NPIX>0]
    IMG = IMG*mask

    ## save these parameters to a NetCDF file so that we can plot it later 
    f = nc.Dataset('multi_proj_raw.nc', 'w')

    xdim     = f.createDimension('x',nlon)
    ydim     = f.createDimension('y',nlat)

    ##  create the NetCDF variables 
    latVar  = f.createVariable('lat', 'float64', ('y'))
    lonVar  = f.createVariable('lon', 'float64', ('x'))
    imgVar  = f.createVariable('img', 'float64', ('y','x'))

    latVar[:]  = newlat[:]
    lonVar[:]  = newlon[:]
    imgVar[:]  = IMG[:]

    f.close()
    
    ## normalize the image by the 95% percentile 
    IMG = IMG/(np.percentile(IMG[IMG>0.], 99.))

    plt.imsave('mosaic.png', IMG, vmin=0., vmax=1., origin='lower', cmap='hot')
    return (newlon, newlat, IMG)

def map_project(LON, LAT, lon, lat, img):
    lons = lon.flatten()
    lats = lat.flatten()
    imgs = img.flatten()

    invmask = np.where((lats==-1000.)|(lons==-1000.))[0]
    ## remove pixels that were not projected
    lat = np.delete(lats, invmask)
    lon = np.delete(lons, invmask)
    img = np.delete(imgs, invmask)

    IMGI = griddata((lon, lat), img, (LON, LAT), method='cubic').T
    
    IMGI[np.isnan(IMGI)]  = 0.
    IMGI[IMGI<0.] = 0.

    #plt.imsave("%s.png"%fname, IMGI/np.percentile(IMGI[IMGI>0.], 95.), cmap='hot')

    return IMGI

