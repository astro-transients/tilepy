import numpy as np
import healpy as hp

def get_ligo_maps(sky_map, Order, map_names='all'):
    
    un_inds = sky_map['UNIQ']
    
    order = (np.log2(un_inds/4).astype(int)/
             2).astype(int)
    inds = (un_inds - 4*(np.power(4,order))).astype(int)
             
    if Order == 'max':
        Order = np.max(order)

    Nside = int(2**Order)
    Npix = hp.nside2npix(Nside)

    if map_names == 'all':
        keys = ['PROBDENSITY', 'DISTMU', 'DISTSIGMA', 'DISTNORM']
    else:
        keys = map_names
    maps = {}
    
    for k in keys:
        maps[k] = np.zeros(Npix)

    #print np.min(order), np.max(order)

    for ii in range(np.max(order),np.min(order)-1,-1):
    
        nside = 2**ii
        npix = hp.nside2npix(nside)
        bl = (order==ii)
        #print ii
        
        for k in maps.keys():
            
            a = hp.UNSEEN*np.ones(npix)
            a[inds[bl]] = sky_map[k][bl]
            if ii == Order:
                bl_ = (a!=hp.UNSEEN)
                maps[k][bl_] += a[bl_]
                del a
            else:
                a_ = hp.ud_grade(a, nside_out=Nside,\
                    order_in='Nested', order_out='Nested')
                bl_ = (a_!=hp.UNSEEN)
                maps[k][bl_] += a_[bl_]
                del a, a_

    return maps

def uniq2order_ind(uniq):
    order = (np.log2(uniq/4).astype(int)/2).astype(int)
    inds = (uniq - 4*(np.power(4,order))).astype(int)
    return order, inds

def order_inds2uniq(order, inds):
    uniq = 4*(np.power(4,order)).astype(int) + inds
    return uniq

if __name__ == "__main__":

    NewNside = 512

    ra = sys.argv[1]
    dec = sys.argv[2]
    PercentCov = float(sys.argv[3])

    #filename = sys.argv[4]
    #filename = '/Users/hashkar/Desktop/GWfollowup_master_alert/G297595_LALInference.fits'

    skymap_fname='/Users/mseglar/Desktop/AMON_Gamma_O1/LVC/img/2015-12-23T16:13:55.82440.fits'
    sky_tab = Table.read(skymap_fname)
    healpix_skymaps_dict = get_ligo_maps(sky_tab, 'max')
    prob=healpix_skymaps_dict['PROBDENSITY']
    #hp.mollview(healpix_skymaps_dict['PROBDENSITY'],nest=True)
    
    npix = len(prob)
    nside = hp.npix2nside(npix)

    pix_ra1, pix_dec1, area = Get90RegionPixReduced(prob, PercentCov, NewNside)


    coordinates=TransformRADec(ra,dec)
    ra1= coordinates.ra.rad
    dec1 = coordinates.dec.rad


    radecs = co.SkyCoord(pix_ra1, pix_dec1, frame='fk5', unit=(u.deg, u.deg))
    pix_ra11 = radecs.ra.rad
    pix_dec11 = radecs.dec.rad


    Igrb = hp.ang2pix(NewNside, ra1, dec1, lonlat=True)
    Imap = hp.ang2pix(NewNside, pix_ra11, pix_dec11, lonlat=True)


    print(np.isin(Igrb, Imap))

    hp.mollview(prob,title="Inside Region",nest=True)
    hp.visufunc.projplot(coordinates.ra, coordinates.dec, 'r.', lonlat=True)
    hp.graticule()
    plt.show()
