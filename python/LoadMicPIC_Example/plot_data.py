from scipy.io import FortranFile
import numpy

def readFile(file):
    # This is where the reading of the binary file is taken care of

    #Contains a list of all the fields-files
    nss = 6
    with open(file, 'r') as f:
        header, time  = np.fromfile(f, dtype='<i4', count=2)
        dt, teti, xmax, zmax,  = np.fromfile(f, dtype='<f4', count=4)
        nnx, nnz = np.fromfile(f, dtype='<i4', count=2)

        vxs = np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz*nss),(nss,nnz,nnx))#(nnx,nnz,nss))
        vys = np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz*nss),(nss,nnz,nnx))#(nnx,nnz,nss))
        vzs = np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz*nss),(nss,nnz,nnx))#(nnx,nnz,nss))

        bx = np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz),(nnz,nnx))
        by = np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz),(nnz,nnx))
        bz = np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz),(nnz,nnx))

        ex = np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz),(nnz,nnx))
        ey = np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz),(nnz,nnx))
        ez = np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz),(nnz,nnx))

        dns = np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz*nss),(nss,nnz,nnx))

        xe = np.fromfile(f, dtype = '<f4', count = nnx)
        ze = np.fromfile(f, dtype = '<f4', count = nnz)

        mass = np.fromfile(f, dtype = '<f4', count = nss)
        q = np.fromfile(f, dtype = '<f4', count = nss)

        time = np.fromfile(f, dtype = '<f8', count = 1)
        wpewce =  np.fromfile(f, dtype = '<f4', count = 1)
        dfac =  np.fromfile(f, dtype = '<f4', count = nss)

        vxx =  np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz*nss),(nss,nnz,nnx))
        vyy =  np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz*nss),(nss,nnz,nnx))
        vzz =  np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz*nss),(nss,nnz,nnx))
        vxy =  np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz*nss),(nss,nnz,nnx))
        vxz =  np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz*nss),(nss,nnz,nnx))
        vyz =  np.reshape(np.fromfile(f, dtype = '<f4', count = nnx*nnz*nss),(nss,nnz,nnx))

    print("Done reading {}".format(file))

    xe = xe/np.sqrt(mass[0])
    ze = ze/np.sqrt(mass[0])

    coords = [xe,ze]
    bx = bx*wpewce[0]
    by = by*wpewce[0]
    bz = bz*wpewce[0]

    ex = ex*np.sqrt(mass[0])*wpewce[0]**2
    ey = ey*np.sqrt(mass[0])*wpewce[0]**2
    ez = ez*np.sqrt(mass[0])*wpewce[0]**2

##### Plasma species 1
    ind = 1;
    n = dns[ind,:,:].squeeze()*dfac[ind]

    jx = vxs[ind,:,:].squeeze()*dfac[ind]
    jy = vys[ind,:,:].squeeze()*dfac[ind]
    jz = vzs[ind,:,:].squeeze()*dfac[ind]

    vx = jx/n
    vy = jy/n
    vz = jz/n
    vx[n < 0.005] = 0
    vy[n < 0.005] = 0
    vz[n < 0.005] = 0

    pxx = vxx[ind,:,:].squeeze()*dfac[ind]
    pyy = vyy[ind,:,:].squeeze()*dfac[ind]
    pzz = vzz[ind,:,:].squeeze()*dfac[ind]
    pxy = vxy[ind,:,:].squeeze()*dfac[ind]
    pxz = vxz[ind,:,:].squeeze()*dfac[ind]
    pyz = vyz[ind,:,:].squeeze()*dfac[ind]

    pxz = mass[ind]*(vxz-vx*jz)
    pxx = mass[ind]*(vxx-vx*jx)
    pxy = mass[ind]*(vxy-vx*jy)
    pyy = mass[ind]*(vyy-vy*jy)
    pyz = mass[ind]*(vyz-vy*jz)
    pzz = mass[ind]*(vzz-vz*jz)


    # Scaled for 2nd time

    jx = jx*np.sqrt(mass[0])*wpewce[0]
    jy = jy*np.sqrt(mass[0])*wpewce[0]
    jz = jz*np.sqrt(mass[0])*wpewce[0]

    vOx = jOx/dnO
    vOy = jOy/dnO
    vOz = jOz/dnO
    vOx[dnO < 0.005] = 0
    vOy[dnO < 0.005] = 0
    vOz[dnO < 0.005] = 0

    jeOx = jeOx*np.sqrt(mass[0])*wpewce[0]
    jeOy = jeOy*np.sqrt(mass[0])*wpewce[0]
    jeOz = jeOz*np.sqrt(mass[0])*wpewce[0]


    veOx = jeOx/dneO
    veOy = jeOy/dneO
    veOz = jeOz/dneO

    veOx[dneO < 0.005] = 0
    veOy[dneO < 0.005] = 0
    veOz[dneO < 0.005] = 0

    pxxO = pxxO*wpewce[0]**2
    pyyO = pyyO*wpewce[0]**2
    pzzO = pzzO*wpewce[0]**2
    pxyO = pxyO*wpewce[0]**2
    pxzO = pxzO*wpewce[0]**2
    pyzO = pyzO*wpewce[0]**2

    pxxeO = pxxeO*wpewce[0]**2
    pyyeO = pyyeO*wpewce[0]**2
    pzzeO = pzzeO*wpewce[0]**2
    pxyeO = pxyeO*wpewce[0]**2
    pxzeO = pxzeO*wpewce[0]**2
    pyzeO = pyzeO*wpewce[0]**2

    pO=(pxxO+pyyO+pzzO)/3
    peO=(pxxeO+pyyeO+pzzeO)/3


    #### Total electron pressure

    jex = vxs[3,:,:].squeeze()*dfac[3] + vxs[1,:,:].squeeze()*dfac[1]
    jey = vys[3,:,:].squeeze()*dfac[3] + vys[1,:,:].squeeze()*dfac[1]
    jez = vzs[3,:,:].squeeze()*dfac[3] + vzs[1,:,:].squeeze()*dfac[1]
    dne = dns[3,:,:].squeeze()*dfac[3] + dns[1,:,:].squeeze()*dfac[1]

    vex = jex/dne
    vey = jey/dne
    vez = jez/dne

    vex[dne < 0.005] = 0
    vey[dne < 0.005] = 0
    vez[dne < 0.005] = 0

    pxxe = pxx[3,:,:].squeeze()*dfac[3] + pxx[1,:,:].squeeze()*dfac[1]
    pyye = pyy[3,:,:].squeeze()*dfac[3] + pyy[1,:,:].squeeze()*dfac[1]
    pzze = pzz[3,:,:].squeeze()*dfac[3] + pzz[1,:,:].squeeze()*dfac[1]
    pxye = pxy[3,:,:].squeeze()*dfac[3] + pxy[1,:,:].squeeze()*dfac[1]
    pxze = pxz[3,:,:].squeeze()*dfac[3] + pxz[1,:,:].squeeze()*dfac[1]
    pyze = pyz[3,:,:].squeeze()*dfac[3] + pyz[1,:,:].squeeze()*dfac[1]

    pxze = mass[3]*(pxze-vex*jez)*wpewce[0]**2
    pxxe = mass[3]*(pxxe-vex*jex)*wpewce[0]**2
    pxye = mass[3]*(pxye-vex*jey)*wpewce[0]**2
    pyye = mass[3]*(pyye-vey*jey)*wpewce[0]**2
    pyze = mass[3]*(pyze-vey*jez)*wpewce[0]**2
    pzze = mass[3]*(pzze-vez*jez)*wpewce[0]**2

    ### Scaled for 2nd time
    jex = jex*np.sqrt(mass[0])*wpewce[0]
    jey = jey*np.sqrt(mass[0])*wpewce[0]
    jez = jez*np.sqrt(mass[0])*wpewce[0]

    vex = jex/dne
    vey = jey/dne
    vez = jez/dne

    vex[dne < 0.005] = 0
    vey[dne < 0.005] = 0
    vez[dne < 0.005] = 0




    ### Done with hydrogen

    threshold = 1e-5
    pH=(pxxH+pyyH+pzzH)/3
    peH=(pxxeH+pyyeH+pzzeH)/3
    pe=(pxxe+pyye+pzze)/3

    ti = pH/dnH
    te = pe/dne

    tmp=(jeHx**2+jeHy**2+jeHz**2)/2
    ekeH = tmp/dneH
    ekeH[dneH <threshold] = 0

    tmp=(jeOx**2+jeOy**2+jeOz**2)/2
    ekeO = tmp/dneO
    ekeO[dneO < threshold] = 0

    tmp=(jHx**2+jHy**2+jHz**2)/2

    ekH = tmp/dnH
    ekH[dnH < threshold] = 0

    tmp=(jOx**2+jOy**2+jOz**2)/2
    ekO = tmp/dnO
    ekO[dnO < threshold] = 0


    ## circ(A,1,2) = np.roll(A,1,1), circ(A,1,1) = np.roll(A,1,0)

    bx_int = 0.5*(bx + np.roll(bx,1,0))
    ez_int = 0.5*(ez + np.roll(ez,1,0))
    by_int = 0.25*(by + np.roll(by,1,0) + np.roll(by,1,1) + np.roll(by,[1, 1], axis=(0, 1)))


    dx = xe[1]-xe[0]
    dz = ze[1]-ze[0]

    ixm = 10-1
    a = np.zeros((nnz,nnx))
    a[ixm,1] = 0
    for j in range(1,nnz):
        a[j,ixm] = a[j-1,ixm] + dz*bx[j-1,ixm]
    for ix in range(ixm+1,nnx):
        a[:,ix]= a[:,ix-1] - bz[:,ix-1]*dx
    for ix in range(ixm,0,-1):
        a[:,ix]= a[:,ix+1] + bz[:,ix]*dx

    # Proceed with the plotting

    testDict = {}
    testDict.update({"dnH":dnH})#,testDict.update({"dnO":dnO})
    #testDict.update({"dneH":dneH}),testDict.update({"dneO":dneO})

    varDictionary = {}
    varDictionary.update({"dnH":dnH}),varDictionary.update({"dnO":dnO})
    varDictionary.update({"dneH":dneH}),varDictionary.update({"dneO":dneO})

    varDictionary.update({"ex":ex}),varDictionary.update({"bx":bx})
    varDictionary.update({"ey":ey}),varDictionary.update({"by":by})
    varDictionary.update({"ez":ez}),varDictionary.update({"bz":bz})

    varDictionary.update({"jHx":jHx}),varDictionary.update({"jOx":jOx})
    varDictionary.update({"jHy":jHy}),varDictionary.update({"jOy":jOy})
    varDictionary.update({"jHz":jHz}),varDictionary.update({"jOz":jOz})
    varDictionary.update({"jeHx":jeHx}),varDictionary.update({"jeOx":jeOx})
    varDictionary.update({"jeHy":jeHy}),varDictionary.update({"jeOy":jeOy})
    varDictionary.update({"jeHz":jeHz}),varDictionary.update({"jeOz":jeOz})

    varDictionary.update({"vHx":vHx}), varDictionary.update({"vOx":vOx})
    varDictionary.update({"vHy":vHy}), varDictionary.update({"vOy":vOy})
    varDictionary.update({"vHz":vHz}), varDictionary.update({"vOz":vOz})
    varDictionary.update({"veHx":veHx}),varDictionary.update({"veOx":veOx})
    varDictionary.update({"veHy":veHy}),varDictionary.update({"veOy":veOy})
    varDictionary.update({"veHz":veHz}),varDictionary.update({"veOz":veOz})
    varDictionary.update({"vex":vex}),varDictionary.update({"pzze":pzze})
    varDictionary.update({"vey":vey}),varDictionary.update({"bx_int":bx_int})
    varDictionary.update({"vez":vez}),varDictionary.update({"by_int":by_int})
    varDictionary.update({"ez_int":ez_int}), varDictionary.update({"dne":dne})

    varDictionary.update({"pxze":pxze}),varDictionary.update({"pxxe":pxxe})
    varDictionary.update({"pyze":pyze}),varDictionary.update({"pxye":pxye})
    varDictionary.update({"pxxH":pxxH}), varDictionary.update({"pxxO":pxxO})
    varDictionary.update({"pxxeH":pxxeH}), varDictionary.update({"pxxeO":pxxeO})
    varDictionary.update({"pxyeH":pxyeH}), varDictionary.update({"pxyeO":pxyeO})
    varDictionary.update({"pxyH":pxyH}), varDictionary.update({"pxyO":pxyO})
    varDictionary.update({"pzzeH":pzzeH}), varDictionary.update({"pzzeO":pzzeO})
    varDictionary.update({"pzzH":pzzH}), varDictionary.update({"pzzO":pzzO})
    varDictionary.update({"pxzeH":pxzeH}), varDictionary.update({"pxzeO":pxzeO})
    varDictionary.update({"pxzH":pxzH}), varDictionary.update({"pxzO":pxzO})
    varDictionary.update({"pyzeH":pyzeH}), varDictionary.update({"pyzeO":pyzeO})
    varDictionary.update({"pyzH":pyzH}), varDictionary.update({"pyzO":pyzO})
    varDictionary.update({"pyyH":pyyH}), varDictionary.update({"pyyO":pyyO})
    varDictionary.update({"pyyeH":pyyeH}), varDictionary.update({"pyyeO":pyyeO})
    varDictionary.update({"ti":ti}), varDictionary.update({"te":te})

    varDictionary.update({"ekeH":ekeH}),varDictionary.update({"ekeO":ekeO})
    varDictionary.update({"ekH":ekH}),varDictionary.update({"ekO":ekO})
    varDictionary.update({"pO":pO}),varDictionary.update({"peO":peO})
    varDictionary.update({"pH":pH}),varDictionary.update({"peH":peH})



    quiverDictionary = {}
    quiverDictionary.update({"dnH":[vHx,vHz,"Hydrogen"]})
    quiverDictionary.update({"dnO":[vOx,vOz,"Oxygen"]})
    quiverDictionary.update({"dneH":[veHx,veHz,"Hydrogen Electrons"]})
    quiverDictionary.update({"dneO":[veOx,veOz,"Oxygen Electrons"]})

    varDictionaryPanel = {}
    varDictionaryPanel.update({"dnH":dnH}),varDictionaryPanel.update({"dnO":dnO})
    varDictionaryPanel.update({"dneH":dneH}),varDictionaryPanel.update({"dneO":dneO})
    varDictionary.update({"ey":ey}),varDictionary.update({"by":by})

    data = []
    fluxField = a
    time = round(float(file[-9:-4])/wpewce[0]/mass[0],2)

    data.append(varDictionary),data.append(fluxField),data.append(coords), data.append(time), data.append([mass,wpewce])
    data.append(varDictionaryPanel), data.append(quiverDictionary)
    return data
