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


import os
import time
import logging
import traceback
from urllib.parse import urlparse
from urllib.request import urlretrieve

import mhealpy as mh

from astropy.io import fits

##################################################################################################
#                        Read Healpix map from fits file                                         #
##################################################################################################

logger = logging.getLogger(__name__)


class MapReader:
    download_max_nb_try = 20
    time_wait_retry = 60

    def __init__(self, obspar):
        skymap_localisation = obspar.skymap

        # Check if the provided localisation is a url or a path
        parsed_localisation = urlparse(skymap_localisation)
        self.is_remote = parsed_localisation.scheme != ''

        # Download map if necessary
        self.url = None
        self.skymap_filename = skymap_localisation
        if self.is_remote:
            self.url = skymap_localisation
            max_retry = 1
            if 'gbm' in obspar.alertType:
                max_retry = self.download_max_nb_try
            self.skymap_filename = self.downloadMap(download_max_nb_try=max_retry,
                                                    time_wait_retry=self.time_wait_retry)

        # Check the file is available
        if not os.path.isfile(self.skymap_filename):
            if self.is_remote:
                raise Exception('Issue with downloaded file, not existing anymore')
            else:
                raise Exception('Issue with file input, the file doesn\'t exist and is not an url')

        self.fits_map = fits.open(self.skymap_filename)
        self.id_hdu_map = self.getMapHDUId()
        self.name_event = self.getSourceName()
        self.identifyColumns()

    def getMapHDUId(self):
        id_hdu_map = -1
        for i in range(len(self.fits_map)):
            if 'XTENSION' in self.fits_map[i].header and self.fits_map[i].header['XTENSION'] == 'BINTABLE':
                if id_hdu_map != -1:
                    raise Exception('Multiple map detected, please provide a file with only one map')
                else:
                    id_hdu_map = i
        if id_hdu_map == -1:
            raise Exception('No map detected, please provide a file with a map')
        return id_hdu_map

    def identifyColumns(self):
        nb_column = self.fits_map[self.id_hdu_map].header['TFIELDS']
        if nb_column == 1:
            self.id_prob = 0
            self.id_dist_mean = None
            self.id_dist_sigma = None
            self.id_dist_norm = None
        else:
            self.id_dist_mean = None
            self.id_dist_sigma = None
            self.id_dist_norm = None
            for i in range(nb_column):
                if self.fits_map[self.id_hdu_map].header['TTYPE'+str(i)] == 'PROB':
                    self.id_prob = i
                elif self.fits_map[self.id_hdu_map].header['TTYPE'+str(i)] == 'DISTMU':
                    self.id_dist_mean = i
                elif self.fits_map[self.id_hdu_map].header['TTYPE'+str(i)] == 'DISTSIGMA':
                    self.id_dist_sigma = i
                elif self.fits_map[self.id_hdu_map].header['TTYPE'+str(i)] == 'DISTNORM':
                    self.id_dist_norm = i

        if self.id_dist_mean is None or self.id_dist_mean is None or self.id_dist_norm is None:
            self.has3D = False
        else:
            self.has3D = True

    def getSourceName(self):
        """
            Get the source name from the contents of the fits file

            :return: name
        """

        name = 'undefined'
        if 'OBJECT' in self.fits_map[self.id_hdu_map].header:
            name = self.fits_map[self.id_hdu_map].header['OBJECT']
        if 'SENDER' in self.fits_map[self.id_hdu_map].header and self.fits_map[self.id_hdu_map].header['SENDER'] == 'IceCube Collaboration':
            name = str(self.fits_map[self.id_hdu_map].header['RUNID']) + '_' + str(self.fits_map[self.id_hdu_map].header['EVENTID'])

        return name

    def downloadMap(self, download_max_nb_try=1, time_wait_retry=60):
        """
            Download the map from a url

            :param download_max_nb_try: If more than 1, will attempt multiple time to download the file in case of error the first time
            :type download_max_nb_try: int
            :param time_wait_retry: the time to wait between each attempt at downloading the file
            :type time_wait_retry: int

            :return: filename
        """

        if download_max_nb_try < 1:
            download_max_nb_try = 1

        filename = self.url.split("/")[-1]
        logger.info("The filename is {filename}")

        file_exist = os.path.isfile(filename)
        if file_exist:
            logger.warning('The file is already existing, it will be re-downloaded')

        for i in range(download_max_nb_try):
            try:
                urlretrieve(self.url, filename)
                break

            except Exception as e:
                if i == (download_max_nb_try - 1) or str(e) != 'HTTP Error 404: Not Found':
                    logger.error('Issue to download map from url')
                    traceback.print_exc()
                    if not file_exist:
                        raise e
                    else:
                        logger.warning('The already existing file will be used')
                else:
                    logger.info(f'Map not available, waiting for {time_wait_retry} before a new attempt')
                    time.sleep(time_wait_retry)

        return filename

    def getMap(self, mapType):
        if mapType == 'prob':
            mh.HealpixMap.read_map(self.skymap_filename, field=self.id_prob, hdu=self.id_hdu_map, density=True)
        elif mapType == 'distMean' and self.has3D:
            mh.HealpixMap.read_map(self.skymap_filename, field=self.id_dist_mean, hdu=self.id_hdu_map, density=False)
        elif mapType == 'distSigma' and self.has3D:
            mh.HealpixMap.read_map(self.skymap_filename, field=self.id_dist_sigma, hdu=self.id_hdu_map, density=False)
        elif mapType == 'distNorm' and self.has3D:
            mh.HealpixMap.read_map(self.skymap_filename, field=self.id_dist_norm, hdu=self.id_hdu_map, density=False)
        elif not self.has3D and (mapType == 'distMean' or mapType == 'distSigma' or mapType == 'distNorm'):
            raise Exception('No distance information available')
        else:
            raise Exception('Unknown type of map')

    def getDistance(self):
        if self.has3D:
            return self.fits_map[self.id_hdu_map].header['DISTMEAN'], self.fits_map[self.id_hdu_map].header['DISTSTD']
        else:
            raise Exception('No distance information available')
