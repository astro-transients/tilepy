#
# Copyright (C) 2016-2024  tilepy developers (Monica Seglar-Arroyo, Halim Ashkar, Fabian Schussler, Mathieu de Bony)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
import os
import time

from astropy.io import fits


##################################################################################################
#                        Read Healpix map from fits file                                         #
##################################################################################################


def GetSkymap(obspar):
    skymap = obspar.skymap

    isURL = False
    # Check if the file is local or not
    if "https" in skymap or "http" in skymap:
        print('it is a url')
        isURL = True
    # Adapt the download of the file from URL to the type of event
    if isURL:
        if obspar.alertType == 'gbmpng':
            fitsMap, filename = GetGBMMap(skymap)
            if fitsMap is None and filename is None:
                print('The localization map is not available, returning.')
                return
            name = skymap.split('/')[-3]
        elif obspar.alertType == 'gbm':
            fitsMap = fits.open(skymap)
            if fitsMap is None:
                print('The localization map is not available, returning.')
                return
            filename = skymap
            name = skymap.split('all_')[1].split('_v00')[0]
        elif obspar.alertType == 'nucascade':
            fitsMap = fits.open(skymap)
            if fitsMap is None:
                print('The localization map is not available, returning.')
                return
            filename = skymap
            name = skymap.split('/')[-1].split('_run')[0]
        else:
            fitsMap, filename = GetGWMap(skymap)
            name = skymap.split('/')[-3]
    else:
        # Open the local file, that is the correct one. There is no further manipulation of the filename needed.
        fitsMap = fits.open(skymap)
        if fitsMap is None:
            print('The localization map is not available, returning.')
            return

        filename = skymap
        if  '/' in filename:
            if  '.' in filename:
                name = filename.split('/')[-1].split('.')[0]
            else:
                name = filename
        else:
            if '.' in filename:
                name = filename.split('.')[0]
            else:
                name = filename
    return fitsMap, filename, name


def GetGBMMap(URL):
    """
    Bottom-level function that takes a url searches for the localisatios maps from the Fermi-GBM database, or waits until it is uplaoded.

    :param URL: the URL of the map
    :type URL: str

    :return: fitsfile, filename
    rtype: fits, str
    """

    filename = URL.split("/")[-1].split(".")[0]
    filename = "./maps/" + filename + ".fit"
    # filename = "glg_healpix_all_bn211130636_2.fit"
    print("The GBM filename is ", filename)
    try:
        fits_map_url_intial = URL
        fits_map_url1 = fits_map_url_intial.split("/")
        fits_map_url2 = fits_map_url_intial.split("_")[-1]
        # fits_map_url1[-1] = ""
        i = 0
        fits_map_url = ""
        for i in range(len(fits_map_url1) - 1):
            fits_map_url += fits_map_url1[i] + "/"
        fits_map_url += "glg_healpix_all" + "_" + \
            fits_map_url2.split(".")[0] + ".fit"

        filename = fits_map_url
    except:
        warn = "Caught exception: "
        print(warn)
        pass

    max_delay = 20
    delay = 0
    d = 0
    while delay == 0:
        delay = 1
        d = d + 1
        try:
            fitsfile = fits.open(filename)
        except:
            print('map is not uploaded yet... Waiting for minute:', d)
            time.sleep(60)
            delay = 0
            if d > max_delay:
                print(
                    f"Waited for {max_delay} minutes... can't wait anymore... I'm leaving")
                fitsfile = None
                filename = None
                break

    return fitsfile, filename


def GetGWMap(URL):
    """
    Bottom-level function that takes a url searches for the localisation maps in multi-order format from the GW database, or waits until it is uplaoded.

    :param URL: the URL of the map
    :type URL: str

    :return: fitsfile, filename
    rtype: fits, str
    """
    filename = URL.split("/")[-1]
    print("The filename is ", filename)
    fits_map_url = URL
    try:
        command = f'curl {fits_map_url} -o {filename}'
        print(command)
        os.system(command)

    except x:
        print('Problem with downloading map from url, it was not multiorder or fits.gz')
        warn = "Caught exception: %s" % x
        print(warn)
        pass

    fitsfile = fits.open(filename)

    return fitsfile, filename
