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


import logging
import os
import time
import traceback
from urllib.parse import urlparse
from urllib.request import urlretrieve

import mhealpy as mh
import astropy.units as u
from astropy.io import fits

##################################################################################################
#                        Read Healpix map from fits file                                         #
##################################################################################################

logger = logging.getLogger(__name__)

__all__ = [
    "MapReader"
]

class MapReader:

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
            max_retry = obspar.downloadMaxRetry
            time_wait_retry = obspar.downloadWaitPeriodRetry
            self.skymap_filename = self.downloadMap(download_max_nb_try=max_retry+1,
                                                    time_wait_retry=time_wait_retry)

        # Check the file is available
        if not os.path.isfile(self.skymap_filename):
            if self.is_remote:
                raise Exception('Issue with downloaded file, not existing anymore')
            else:
                raise Exception('Issue with file input, the file doesn\'t exist and is not an url')

        self.fits_map = fits.open(self.skymap_filename)
        self.id_hdu_map = self.getMapHDUId()
        self.name_event = self.getSourceName()
        if obspar.event_name is not None:
            self.name_event = obspar.event_name
        if obspar.event_name is None: 
            obspar.event_name = self.name_event
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

        # Identify if column numeration start at 0 or 1
        self.offset_column = 1
        if 'TTYPE0' in self.fits_map[self.id_hdu_map].header:
            self.offset_column = 0

        if nb_column == 1:
            self.id_prob = 0
            self.unit_prob = u.dimensionless_unscaled
            self.id_dist_mean = None
            self.unit_dist_mean = None
            self.id_dist_sigma = None
            self.unit_dist_sigma = None
            self.id_dist_norm = None
            self.unit_dist_norm = None
        else:
            self.unit_prob = None
            self.id_dist_mean = None
            self.unit_dist_mean = None
            self.id_dist_sigma = None
            self.unit_dist_sigma = None
            self.id_dist_norm = None
            self.unit_dist_norm = None
            for i in range(self.offset_column, nb_column+self.offset_column):
                columns_name = self.fits_map[self.id_hdu_map].header['TTYPE' + str(i)]
                unit_information = ('TUNIT'+str(i)) in self.fits_map[self.id_hdu_map].header
                if columns_name in ['PROB', 'T', 'PROBABILITY', 'PROBDENSITY']:
                    self.id_prob = i
                    if unit_information:
                        self.unit_prob = u.Unit(self.fits_map[self.id_hdu_map].header['TUNIT' + str(i)])
                    else:
                        columns_name = self.fits_map[self.id_hdu_map].header['TTYPE'+str(i)]
                        if columns_name in ['T', 'PROBABILITY']:
                            self.unit_prob = u.dimensionless_unscaled
                        else:
                            self.unit_prob = u.Unit('sr^-1')
                elif columns_name == 'DISTMU':
                    self.id_dist_mean = i
                    if unit_information:
                        self.unit_dist_mean = u.Unit(self.fits_map[self.id_hdu_map].header['TUNIT' + str(i)])
                    else:
                        self.unit_dist_mean = u.Mpc
                elif columns_name == 'DISTSIGMA':
                    self.id_dist_sigma = i
                    if unit_information:
                        self.unit_dist_sigma = u.Unit(self.fits_map[self.id_hdu_map].header['TUNIT' + str(i)])
                    else:
                        self.unit_dist_sigma = u.Mpc
                elif columns_name == 'DISTNORM':
                    self.id_dist_norm = i
                    if unit_information:
                        self.unit_dist_norm = u.Unit(self.fits_map[self.id_hdu_map].header['TUNIT' + str(i)])
                    else:
                        self.unit_dist_norm = u.Unit('Mpc^-2')

        # Correct unit if needed
        if self.unit_prob.is_equivalent(u.Unit('pix^-1')):
            self.unit_prob = u.dimensionless_unscaled

        # Check if probabilities is probabilities density
        self.prob_density = True
        if self.unit_prob.is_equivalent(u.dimensionless_unscaled):
            self.prob_density = False

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
        # if the event if from LVK and the URL is from GraceDB, get the superevent name
        if "ORIGIN" in self.fits_map[self.id_hdu_map].header.keys():
            if "LIGO/Virgo/KAGRA" in self.fits_map[self.id_hdu_map].header['ORIGIN'] and self.is_remote:
                if "https://gracedb.ligo.org/api/superevents" in self.url:
                    name = self.url.split("/")[5]
        # if the event if from Fermi-GBM, get the GBM name (i.e. replace GRB with bn)
        if "TELESCOP" in self.fits_map["PRIMARY"].header.keys():
            if self.fits_map["PRIMARY"].header["TELESCOP"] == "GLAST":
                name = self.fits_map[self.id_hdu_map].header["OBJECT"].replace("GRB", "bn")

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
        logger.info(f"The filename is {filename}")

        file_exist = os.path.isfile(filename)
        if file_exist:
            logger.warning('The file is already existing, it will be re-downloaded')

        for i in range(download_max_nb_try):
            if 'png' in self.url: #Change Fermi-GBM url(if png) to fit format
                self.url = self.url.replace("png", "fit")
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
            raw_map = mh.HealpixMap.read_map(self.skymap_filename, field=self.id_prob-self.offset_column, hdu=self.id_hdu_map, density=True)
            if not self.prob_density:
                raw_map._data = (raw_map.data * self.unit_prob / raw_map.pixarea()).to_value(u.Unit('sr^-1'))
                raw_map._unit = u.Unit('sr^-1')
            else:
                raw_map._data = (raw_map.data * self.unit_prob).to_value(u.Unit('sr^-1'))
        elif mapType == 'distMean' and self.has3D:
            raw_map = mh.HealpixMap.read_map(self.skymap_filename, field=self.id_dist_mean-self.offset_column, hdu=self.id_hdu_map, density=False)
            raw_map._data = (raw_map.data * self.unit_dist_mean).to_value(u.Unit('Mpc'))
        elif mapType == 'distSigma' and self.has3D:
            raw_map = mh.HealpixMap.read_map(self.skymap_filename, field=self.id_dist_sigma-self.offset_column, hdu=self.id_hdu_map, density=False)
            raw_map._data = (raw_map.data * self.unit_dist_sigma).to_value(u.Unit('Mpc'))
        elif mapType == 'distNorm' and self.has3D:
            raw_map = mh.HealpixMap.read_map(self.skymap_filename, field=self.id_dist_norm-self.offset_column, hdu=self.id_hdu_map, density=False)
            raw_map._data = (raw_map.data * self.unit_dist_norm).to_value(u.Unit('Mpc^-2'))
        elif not self.has3D and (mapType == 'distMean' or mapType == 'distSigma' or mapType == 'distNorm'):
            raise Exception('No distance information available')
        else:
            raise Exception('Unknown type of map')

        return raw_map

    def getDistance(self):
        if self.has3D:
            return self.fits_map[self.id_hdu_map].header['DISTMEAN'], self.fits_map[self.id_hdu_map].header['DISTSTD']
        else:
            raise Exception('No distance information available')
