# -*- coding: UTF-8 -*-
"""
Authors : Julien Pellen

"""
from osgeo import gdal, gdalconst, ogr
import numpy as np
import matplotlib.pyplot as plt
#import rasterio as rst
#from rasterio.merge import merge
#from pyproj import CRS
#import affine
import skimage.filters as sif
import datetime
import csv
import fiona
import re
import sys
import os
import shutil as sh
from math import ceil
from PIL import Image
from pathlib import Path
import pandas as pd
from sklearn.cluster import KMeans
#import rasterio
import json



class Norma_Data:

    def __init__(self, path_img):
        # Check Path
        g = "telenvi"
        if type(path_img) != type(g):
            raise TypeError('Error : The Path must be a string')

        if os.path.exists(path_img) is False:
            raise FileNotFoundError("Error 1 : Path invalid")

        # Loading raster
        ds = gdal.Open(path_img, gdal.GA_Update)

        # Metadata raster
        gt = ds.GetGeoTransform()

        # Raster Projection
        prj = ds.GetProjection()

        # Raster from array
        np_ds = gdal.Open(path_img).ReadAsArray()

        # Raster name
        name = os.path.split(path_img)[1]

        self.objectGDAL = ds
        self.band_Array = np_ds
        self.geodata = gt
        self.shape = np_ds.shape
        self.resolution = gt[1]
        self.name = name
        self.projection = prj
        self.path_file = path_img
        self.band_resample = None  # set by resample()
        self.cm = None  # set by cloud_mask()
        self.array_crop = None  # set by crop_ras() ou crop_ras_shp
        self.median_filter_array = None
        self.gaussian_filter_array = None
        self.fourrier_transform = None
        self.geodata_crop = None
        self.objectGDAL = ds
        self.band_Array = np_ds
        self.geodata = gt
        self.shape = np_ds.shape
        self.resolution = gt[1]
        self.name = name
        self.projection = prj
        self.path_file = path_img
        self.band_resample = None  # set by resample()
        self.cm = None  # set by cloud_mask()
        self.array_crop = None  # set by crop_ras() ou crop_ras_shp
        self.median_filter_array = None
        self.gaussian_filter_array = None
        self.fourrier_transform = None
        self.geodata_crop = None

    def cloud_mask(self):

        ci1 = (self.band_Array[2] + self.band_Array[3] + self.band_Array[4]) / (
                self.band_Array[6] + self.band_Array[12] * 2)

        ci2 = (self.band_Array[2] + self.band_Array[3] + self.band_Array[4] + self.band_Array[6] + self.band_Array[12] +
               self.band_Array[13]) / 6

        classifci1 = np.where(ci1 < 0.1, 5000, 0)

        meanci2 = np.mean(ci2[~np.isnan(ci2)])

        maxci2 = np.max(ci2[~np.isnan(ci2)])

        T2 = meanci2 + (0.1 * (maxci2 - meanci2))

        classifci2 = np.where(ci2 > T2, 5000, 0)

        finalclassif = classifci1 + classifci2

        self.cm = finalclassif

    def resample(self, Xres, Yres, path_out=None, algo=gdalconst.GRA_NearestNeighbour, create_ras=False) -> None:

        if create_ras is True:

            option = gdal.WarpOptions(xRes=Xres, yRes=Yres, resampleAlg=algo)
            ds = self.objectGDAL
            gdal.Warp(path_out + '/' + self.name[:-4] + str("_resample.tif"), ds, options=option)

            self.band_resample = gdal.Open(path_out + '/' + self.name[:-4] + '_resample.tif').ReadAsArray()

        else:

            option = gdal.WarpOptions(xRes=Xres, yRes=Yres, resampleAlg=algo, format="VRT")
            ds = self.objectGDAL
            ds_resample = gdal.Warp(destNameOrDestDS="", srcDSOrSrcDSTab=ds, options=option)
            np_resample = np.array(ds_resample.ReadAsArray())
            self.band_resample = np_resample

    def visualization(self, selectband=False) -> None:
        """
        :param selectband: bool : if your image is a stack you have to select the band to visualize, the parameter must be True. Default = False
        :return: None
        """

        if selectband is False:
            plt.imshow(self.band_Array)
            plt.show()
        else:
            num_band = input('Choose the number (positionnal) band : ')
            plt.imshow(self.band_Array[int(num_band) - 1])
            plt.show()

    def update_stack(self, ds_update, multi_update=False):
        """

        :param multi_update:
        :param ds_update: np.array : A matrix who want add to your raster,
        :return: np.array update : (origin shape (4, 250, 250) -> after method -> shape = (5, 250, 250))
        """

        # Check shape

        if np.shape(self.band_Array[0]) != np.shape(ds_update):
            raise TypeError('Error Shape : Check the shape of your matrix %s and object %s ' % (
                np.shape(self.band_Array[0], np.shape(ds_update))))

        if multi_update is False:

            liste_array_base = []
            for band in range(int(self.shape[0])):
                liste_array_base.append(self.band_Array[band])

            liste_array_base.append(ds_update)
            stack = np.stack(liste_array_base, axis=0)

            # Update attribute
            self.band_Array = stack
            self.shape = np.shape(stack)

        else:
            liste_array_base = []
            for band in range(int(self.shape[0])):
                liste_array_base.append(self.band_Array[band])

            for band in range(int(ds_update.shape[0])):
                liste_array_base.append(ds_update[band])

            stack = np.stack(liste_array_base, axis=0)

            self.band_Array = stack
            self.shape = np.shape(stack)

    def crop_ras(self, zone, create_ras=False):

        """
        else:

            or_x = self.geodata[0]
            print(or_x)
            or_y =self.geodata[3]
            print(or_y)

            width_px = self.geodata[1]
            print(width_px)

            height_px = self.geodata[5]
            print(height_px)

            row1 = int((ymax-or_x)/width_px)

            col1 = int((xmin-or_y)/height_px)

            row2 = int((ymin-or_x)/width_px)

            col2 = int((xmax-or_y)/height_px)

        print(col1, row1, col2-col1+1, row2-row1+1)

        

        if self.shape[0] != 1:
            l_array = []
            for band in range(self.shape[0]):
                ds = self.objectGDAL
                crop = ds.GetRasterBand(band+1).ReadAsArray(col1, row1, col2-col1+1, row2-row1+1)
                l_array.append(crop)

            self.array_crop = l_array
        else:
            ds = self.objectGDAL
            crop = ds.ReadAsArray(abs(col1), abs(row1), abs(col2 - col1 + 1),abs(row2 - row1 + 1)).astype(np.float32)

            self.array_crop = crop
        """

        # Check Dimension raster (CD = Check Dimension)
        # If CD is True we have a raster with multi-band
        CD = None
        if len(self.shape) > 2:
            CD = True
        else:
            CD = False

        if create_ras is True:

            if CD is True:

                option = gdal.TranslateOptions(projWin=zone, bandList=[e + 1 for e in range(self.shape[0])])
                gdal.Translate(destName=self.name[:-4] + '_crop.tif', srcDS=str(self.path_file), options=option)
                ds_crop = gdal.Open(self.name[:-4] + '_crop.tif').ReadAsArray()

                self.array_crop = ds_crop
            else:
                option = gdal.TranslateOptions(projWin=zone, bandList=[1])
                gdal.Translate(destName=self.name[:-4] + '_crop.tif', srcDS=str(self.path_file), options=option)
                ds_crop = gdal.Open(self.name[:-4] + '_crop.tif').ReadAsArray()

                self.array_crop = ds_crop

        else:

            if CD is True:
                print('CD is TRUE')

                option = gdal.TranslateOptions(projWin=zone, bandList=[e + 1 for e in range(self.shape[0])],
                                               format='VRT')
                crop = gdal.Translate(destName="", srcDS=str(self.path_file), options=option)
                crop1 = crop.ReadAsArray()
                ds_crop = np.array(crop1)

                self.array_crop = ds_crop

            else:
                option = gdal.TranslateOptions(projWin=zone, bandList=[1], format='VRT')
                crop = gdal.Translate(destName="", srcDS=str(self.path_file), options=option)
                crop1 = crop.ReadAsArray()
                ds_crop = np.array(crop1)

                self.array_crop = ds_crop

    def crop_ras_shp(self, path_shp, create_ras=False):
        """
        Methode de crop avec en entré un shapefile, Le crop s'effectue sur l'emprise maximale du fichier shape
        La sortie est une étendue rectangulaire qui ne colle pas au fichier shape !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        :param path_shp:
        :param create_ras:
        :return: Les coordonnées à l'origine du découpage xmin et ymax
        """

        with fiona.open(path_shp, "r") as shapefile:
            geoms = [feature['geometry'] for feature in shapefile][0]  # find geometry of the shapefile

        coord = geoms['coordinates'][0]

        print(coord)

        l_x = []
        l_y = []

        for e in coord:
            l_x.append(e[0])
            l_y.append(e[1])

        orga_x = sorted(l_x)
        orga_y = sorted(l_y)

        x_min = orga_x[0]
        x_max = orga_x[-1]
        y_min = orga_y[0]
        y_max = orga_y[-1]

        zone = [x_min, y_max, x_max, y_min]
        print(f'Zone : {zone}')
        CD = None

        if len(self.shape) > 2:
            CD = True
        else:
            CD = False

        if create_ras is True:

            if CD is True:

                option = gdal.TranslateOptions(projWin=zone, bandList=[e + 1 for e in range(self.shape[0])], noData=np.NAN)
                gdal.Translate(destName=self.path_file[:-4] + '_crop.tif', srcDS=str(self.path_file), options=option)
                ds_crop = gdal.Open(self.path_file[:-4] + '_crop.tif').ReadAsArray()

                self.array_crop = ds_crop
                self.geodata_crop = [x_min, self.geodata[1], self.geodata[2], y_max, self.geodata[4], abs(self.geodata[5])]
            else:
                option = gdal.TranslateOptions(projWin=zone, bandList=[1], noData=np.NAN)
                gdal.Translate(destName=self.path_file[:-4] + '_crop.tif', srcDS=str(self.path_file), options=option)
                ds_crop = gdal.Open(self.path_file[:-4] + '_crop.tif').ReadAsArray()

                self.array_crop = ds_crop
                self.geodata_crop = [x_min, self.geodata[1], self.geodata[2], y_max, self.geodata[4], abs(self.geodata[5])]

        else:

            if CD is True:

                option = gdal.TranslateOptions(projWin=zone, bandList=[e + 1 for e in range(self.shape[0])],
                                               format='VRT')
                crop = gdal.Translate(destName="", srcDS=str(self.path_file), options=option)
                crop1 = crop.ReadAsArray()
                ds_crop = np.array(crop1)

                self.array_crop = ds_crop
                self.geodata_crop = [x_min, self.geodata[1], self.geodata[2], y_max, self.geodata[4], abs(self.geodata[5])]


            else:

                option = gdal.TranslateOptions(projWin=zone, bandList=[1], format='VRT')
                crop = gdal.Translate(destName="", srcDS=str(self.path_file), options=option)
                crop1 = crop.ReadAsArray()
                ds_crop = np.array(crop1)

                self.array_crop = ds_crop
                self.geodata_crop = [x_min, self.geodata[1], self.geodata[2], y_max, self.geodata[4], abs(self.geodata[5])]

            return self.array_crop, self.path_file[:-4] + '_crop.tif'

    def update_crs(self, new_proj):
        ds = self.objectGDAL
        ds.SetProjection(new_proj)
        self.objectGDAL = ds
        self.projection = new_proj

    def median_filter_meth(self):

        list_Median = []

        CD = None

        if len(self.shape) > 2:
            CD = True
        else:
            CD = False

        if CD is True:

            for e in range((self.shape[0])):
                filter = sif.median(self.band_Array[e, :, :])
                list_Median.append(filter)
            stack_median = np.stack(list_Median, axis=0)
            self.median_filter_array = stack_median

        else:
            filter = sif.median(self.band_Array)
            self.median_filter_array = filter

    def gaussian_filter_meth(self):

        CD = None

        if len(self.shape) > 2:
            CD = True
        else:
            CD = False
        print(CD)
        Liste_Gaussian = []

        if CD is True:

            for e in range((self.shape[0])):
                filter = sif.gaussian(self.band_Array[e, :, :])
                Liste_Gaussian.append(filter)

            stack_gauss = np.stack(Liste_Gaussian, axis=0)

            self.gaussian_filter_array = stack_gauss

        else:
            filter = sif.gaussian(self.band_Array)


            self.gaussian_filter_array = filter

    def fourier_transform(self, select_band):

        fft_result = np.fft.fft2(self.band_array[select_band, :, :])
        fft_result_abs = np.abs(fft_result)
        filename = self.name[:-4] + '_fft_results.csv'
        filepath = os.path.join(os.path.dirname(self.path_file), filename)
        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Frequency', 'Amplitude'])
            for i in range(len(fft_result_abs)):
                for j in range(len(fft_result_abs[i])):
                    writer.writerow([(i, j), fft_result_abs[i][j]])

        self.fourrier_transform = fft_result_abs

    def img_to_csv(self, name_csv, filter=False):
        img = self
        table = None
        if filter is False:
            if len(img.shape) <= 2:

                band = img.band_Array
                ds_ravel = np.ravel(band)
                # ds_ravel = ds_ravel.reshape(-1, 1)
                dic = {'Bande1': ds_ravel}
                table = pd.DataFrame(data=dic)

            else:

                nb_bande = img.shape[2]
                if nb_bande >= 8:
                    print(img.shape)
                    nb_bande = img.shape[0]

                dic = {}

                for e in range(nb_bande):
                    band = img.band_Array[e, :, :]
                    ds_ravel = np.ravel(band)
                    # ds_ravel = ds_ravel.reshape(-1, 1)
                    dic['Bande{}'.format(e + 1)] = ds_ravel
                print(dic)
                table = pd.DataFrame(data=dic)

            table.to_csv(name_csv)

            return table
        else:
            img.gaussian_filter_meth()
            if len(img.shape) <= 2:

                band = img.gaussian_filter_array
                ds_ravel = np.ravel(band)
                # ds_ravel = ds_ravel.reshape(-1, 1)
                dic = {'Bande1': ds_ravel}
                table = pd.DataFrame(data=dic)

            else:

                nb_bande = img.shape[2]
                if nb_bande >= 8:
                    print(img.shape)
                    nb_bande = img.shape[0]

                dic = {}

                for e in range(nb_bande):
                    band = img.gaussian_filter_array[e, :, :]
                    ds_ravel = np.ravel(band)
                    # ds_ravel = ds_ravel.reshape(-1, 1)
                    dic['Bande{}'.format(e + 1)] = ds_ravel
                print(dic)
                table = pd.DataFrame(data=dic)

            table.to_csv(name_csv)

            return table

    ######## Static Method ################
    @staticmethod
    def clear_df_unnamed(data, name_file):
        for e in data.columns:
            if e[:3] == 'Unn':
                del data[e]
            if e[:2] == 'ID':
                del data[e]
            if e[:4] == 'METH':
                del data[e]
        data.to_csv(name_file)

    @staticmethod
    def stage_rename(pathfile):
        # print('1 PATHFILE /////////',pathfile)

        S2_bi_c_nd = pd.read_csv(pathfile)

        pathfilesplit = pathfile.split('/')[-1]

        # print('PATHFILESPLIT::::::::::::: ', pathfilesplit)

        dict_name_col = {}

        if pathfilesplit.split('_')[-2] or pathfilesplit.split('_')[3] == 'S2':
            for e in S2_bi_c_nd.columns:
                dict_name_col[e] = e.split('_')[0]

        elif pathfilesplit.split('_')[-2] or pathfilesplit.split('_')[-3] == 'PL':
            for e in S2_bi_c_nd.columns:
                dict_name_col[e] = e.split('_')[1]

        da = S2_bi_c_nd.rename(columns=dict_name_col)
        da['ID_ARBRE'] = ['AC25', 'AC08', 'AC07', 'AC06', 'AC20', 'AC23', 'AC26', 'AC29', 'AC04', 'AC14', 'AC27',
                          'AC19', 'AC02',
                          'AC22', 'AC17', 'AC16', 'AC05', 'AC03', 'AC28', 'AC15', 'AC09', 'AC21', 'AC24', 'AC18',
                          'AC10', 'AC01',
                          'AC11', 'AC13', 'AC12', 'FR09', 'FR06', 'FR08', 'FR07', 'FR11', 'FR10', 'FR12', 'FR13',
                          'FR14', 'FR16',
                          'FR15', 'FR19', 'FR18', 'FR17', 'FR21', 'FR20', 'FR24', 'FR29', 'FR25', 'FR02', 'FR01',
                          'FR03', 'FR04',
                          'FR05', 'FR26', 'FR28', 'FR27', 'FR23', 'FR22', 'PL15', 'PL09', 'PL17', 'PL29', 'PL01',
                          'PL25', 'PL26',
                          'PL23', 'PL07', 'PL21', 'PL10', 'PL20', 'PL22', 'PL18', 'PL11', 'PL14', 'PL05', 'PL27',
                          'PL24', 'PL06',
                          'PL08', 'PL19', 'PL13', 'PL12', 'PL16', 'PL28', 'PL02', 'PL03', 'PL04', 'QR11', 'QR15',
                          'QR12', 'QR05',
                          'QR16', 'QR02', 'QR10', 'QR14', 'QR07', 'QR01', 'QR04', 'QR06', 'QR08', 'QR17', 'QR13',
                          'QR09', 'QR03',
                          'QR23', 'QR29', 'QR27', 'QR30', 'QR20', 'QR21', 'QR19', 'QR22', 'QR18', 'QR28', 'QR24',
                          'QR26', 'QR25']

        da.set_index("ID_ARBRE", inplace=True)

        return da

    @staticmethod
    def resample_matrix(mat, Xres, Yres, algo=gdalconst.GRA_NearestNeighbour):
        try:

            option = gdal.WarpOptions(xRes=Xres, yRes=Yres, resampleAlg=algo, format="VRT")
            ds = mat
            ds_resample = gdal.Warp(destNameOrDestDS="A", srcDSOrSrcDSTab=ds, options=option)
            np_resample = np.array(ds_resample.ReadAsArray())
            return np_resample

        except AttributeError:
            print('Attribute Error')
            pass
    """
    @staticmethod
    def test_compa_img(path_img1, path_img2):
        ds_pl = nd(path_img2)
        ds_s = nd(path_img1)

        # Check Projection

        if ds_s.projection != ds_pl.projection:
            raise TypeError('Error projection : The rasters do not have the same projection')

        stack1 = None
        stack2 = None

        # Check dimension raster
        if len(ds_s.shape) > 2:
            stack1 = True
        if len(ds_s.shape) <= 2:
            stack1 = False
        if len(ds_pl.shape) > 2:
            stack2 = True
        if len(ds_pl.shape) <= 2:
            stack2 = False

        # print('Stack1:', stack1, '\n','Stack2:', stack2)

        coor_final_x_s = None
        coor_final_y_s = None

        coor_final_x_pl = None
        coor_final_y_pl = None

        if stack1 is True:
            # print("Stack1 is True")
            coor_final_x_s = ds_s.geodata[0] + ds_s.shape[1] * (ds_s.geodata[1])
            coor_final_y_s = ds_s.geodata[3] + ds_s.shape[2] * (ds_s.geodata[-1])
            # print(coor_final_x_s, coor_final_y_s)
        elif stack1 is False:
            # print("Stack1 is False")
            coor_final_x_s = ds_s.geodata[0] + ds_s.shape[0] * (ds_s.geodata[1])
            coor_final_y_s = ds_s.geodata[3] + ds_s.shape[1] * (ds_s.geodata[-1])
            # print(coor_final_x_s, coor_final_y_s)

        if stack2 is True:
            # print('Stack2 is True')
            coor_final_x_pl = ds_pl.geodata[0] + ds_pl.shape[1] * (ds_pl.geodata[1])
            coor_final_y_pl = ds_pl.geodata[3] + ds_pl.shape[2] * (ds_pl.geodata[-1])
        elif stack2 is False:
            # print("Stack2 is False")
            coor_final_x_pl = ds_pl.geodata[0] + ds_pl.shape[0] * (ds_pl.geodata[1])
            coor_final_y_pl = ds_pl.geodata[3] + ds_pl.shape[1] * (ds_pl.geodata[-1])

        # print('Coordonnées Sentinel-2:', coor_final_x_s, coor_final_y_s)
        # print('Coordonnées Planet:', coor_final_x_pl, coor_final_y_pl)

        # Creation BBOX for Crop raster

        # print("Geodata_Sentinel2 : ", ds_s.geodata[0], ds_s.geodata[3])
        # print("Geodata_Planet : ", ds_pl.geodata[0], ds_pl.geodata[3])

        # Cas où les coordonnées d'origine de l'img2 sont en dehors de l'img1 (orX et orY de img2 sont à l'Ouest-Sud-Ouest de l'img1)
        if ds_s.geodata[0] < ds_pl.geodata[0] and ds_s.geodata[3] < ds_pl.geodata[3]:
            print('SUD-OUEST')
            # Cas où le dernier pixel de l'img2 est dans l'img1
            if coor_final_y_s > coor_final_y_pl:
                print('coor_final is IN')
                bbox = [ds_s.geodata[0], ds_pl.geodata[3], coor_final_x_pl, coor_final_y_pl]
            # Cas où le dernier pixel de l'img2 est en dehors de l'img1
            else:
                print('coor_final is OUT')
                bbox = [ds_s.geodata[0], ds_pl.geodata[3], coor_final_x_pl, coor_final_y_s]

        # Cas où les coordonnées d'origine de l'img2 sont en dehor de l'img1 (orX et orY de l'img2 sont au Nord-Est de l'img1)
        elif ds_s.geodata[0] > ds_pl.geodata[0] and ds_s.geodata[3] > ds_pl.geodata[3]:
            print('NORD-EST')
            # Cas où le dernier pixel est dans l'img1
            if coor_final_y_pl < coor_final_y_s and coor_final_x_pl < coor_final_x_s:
                print('coor_final is IN')
                bbox = [ds_s.geodata[0], ds_pl.geodata[3], coor_final_x_pl, coor_final_y_pl]

            else:
                print('coor_final is OUT')
                bbox = [ds_s.geodata[0], ds_pl.geodata[3], coor_final_x_pl, coor_final_y_s]

        # Cas où les coordonnées d'origine de l'img2 sont en dehor de l'img1 (orX et orY de l'img2 sont au Nord-Ouest de l'img1)
        elif ds_s.geodata[0] > ds_pl.geodata[0] and ds_s.geodata[3] < ds_pl.geodata[3]:

            print('NORD-OUEST')

            if coor_final_y_pl < coor_final_y_s and coor_final_x_pl < coor_final_x_s:
                print('coor_final is IN')
                bbox = [ds_pl.geodata[0], ds_s.geodata[3], coor_final_x_pl, coor_final_y_pl]

            else:
                print('coor_final is OUT')
                bbox = [ds_pl.geodata[0], ds_s.geodata[3], coor_final_x_pl, coor_final_y_s]

        # Cas où les coordonnées d'origine de l'img2 sont en dehor de l'img1 (orX et orY de l'img2 sont au Sud-Est de l'img 1)
        elif ds_s.geodata[0] < ds_pl.geodata[3] and ds_s.geodata[3] > ds_pl.geodata[3]:

            print('SUD-EST')
            if coor_final_y_pl < coor_final_y_s and coor_final_x_pl < coor_final_x_s:
                print('coor_final is IN')
                bbox = [ds_pl.geodata[0], ds_pl.geodata[3], coor_final_x_pl, coor_final_y_pl]

            else:
                print('coor_final is OUT')
                bbox = [ds_pl.geodata[0], ds_pl.geodata[3], coor_final_x_s, coor_final_y_pl]

        ##### Trois coté corrigé mais je sais plus quel coté il manque

        elif ds_s.geodata[0] == ds_pl.geodata[0] and ds_s.geodata[3] == ds_pl.geodata[3]:
            print('Partiulier')

            if ds_s.shape == ds_pl.shape and ds_s.resolution == ds_pl.resolution:
                print('Pas besoin de Crop')

                pass

        ds_s.crop_ras(bbox)
        ds_pl.crop_ras(bbox)
        # cropS2 = np.where(ds_s.array_crop <= 0, np.nan, ds_s.array_crop)
        # crop_pl = np.where(ds_pl.array_crop <= 0, np.nan, ds_pl.array_crop)

        return ds_s.array_crop, ds_pl.array_crop, bbox
    """
    @staticmethod
    def stack_img_s2_resample(path_dos, ktr='FRE', algo=gdalconst.GRA_NearestNeighbour):
        """
        A partir d'un chemin d'accés (répertoire), permet de stack les images présentes dans le répertoire
        Nécessite des fichiers au format .tif.
        IMPORTANT: Cette Fonction marche sur des données S2 provenant du portail THEIA

        :param path_dos: Chemin d'accés vers les images
        :param ktr: Keyword Traitement par exemple pour Sentinel-2 il y a des images FRE ou SRE, ce paramètre permet de sélectionner les images en fonction de ce traitement
        :return: Stack numpy (matrice)
        """
        # Pattern

        p_fre = re.compile(r'(_+[A-Z]{3,4}?_)')
        pattern_band = re.compile(r'(B+\d{1,2})')

        # List Creation
        l_img = [e for e in os.listdir(path_dos) if e[-4:] == '.tif']

        l_clean = [e for e in l_img if
                   e[re.search(p_fre, e).span()[0] + 1:re.search(p_fre, e).span()[0] + 4] == str(ktr)]

        ext_name = []
        l_B10 = []
        l_B20 = []
        for e in l_clean:
            if e[re.search(pattern_band, e).span()[0]:re.search(pattern_band, e).span()[0] + 2] == 'B2':
                l_B10.append(gdal.Open(os.path.join(path_dos, e)).ReadAsArray())
            elif e[re.search(pattern_band, e).span()[0]:re.search(pattern_band, e).span()[0] + 2] == 'B3':
                l_B10.append(gdal.Open(os.path.join(path_dos, e)).ReadAsArray()), ext_name.append(
                    os.path.join(path_dos, e))
            elif e[re.search(pattern_band, e).span()[0]:re.search(pattern_band, e).span()[0] + 2] == 'B4':
                l_B10.append(gdal.Open(os.path.join(path_dos, e)).ReadAsArray())
            elif e[re.search(pattern_band, e).span()[0]:re.search(pattern_band, e).span()[0] + 3] == 'B8.':
                l_B10.append(gdal.Open(os.path.join(path_dos, e)).ReadAsArray())
            elif e[re.search(pattern_band, e).span()[0]:re.search(pattern_band, e).span()[0] + 2] == 'B5':
                l_B20.append(os.path.join(path_dos, e))
            elif e[re.search(pattern_band, e).span()[0]:re.search(pattern_band, e).span()[0] + 2] == 'B6':
                l_B20.append((os.path.join(path_dos, e)))
            elif e[re.search(pattern_band, e).span()[0]:re.search(pattern_band, e).span()[0] + 2] == 'B7':
                l_B20.append((os.path.join(path_dos, e)))
            elif e[re.search(pattern_band, e).span()[0]:re.search(pattern_band, e).span()[0] + 2] == 'B1':
                l_B20.append((os.path.join(path_dos, e)))
            elif e[re.search(pattern_band, e).span()[0]:re.search(pattern_band, e).span()[0] + 3] == 'B8A':
                l_B20.append((os.path.join(path_dos, e)))

        option = gdal.WarpOptions(xRes=10, yRes=10, resampleAlg=algo, format="VRT")
        for e in l_B20:
            ds = gdal.Open(e)
            ds_resample = gdal.Warp(destNameOrDestDS="", srcDSOrSrcDSTab=ds, options=option)
            np_resample = np.array(ds_resample.ReadAsArray())
            l_B10.append(np_resample)

        # Stack numpy
        stack_10 = np.stack(l_B10, axis=0)
        name = str(ext_name[0])
        return [stack_10, name]

    ####### Piste d'amélioration de la fonction stack_img_s2
    """
            l_B10 = [e for e in l_clean if e[re.search(pattern_band, e).span()[0]:re.search(pattern_band, e).span()[0]+2] == 'B2' or 'B3' or 'B4' or 'B8']

            l_B20 = [gdal.Open(os.path.join(path_dos, e)).ReadAsArray() for e in l_clean if e[re.search(pattern_band, e).span()[0]:re.search(pattern_band, e).span()[0]+3] == 'B5.' or 'B6.' or 'B7.' or 'B11' or 'B12' or 'B8A']
            print(l_B10)
            print(len(l_B10))


            # Stack numpy
            stack_10 = np.stack(l_B10, axis=0)
            stack_20 = np.stack(l_B20, axis=0)

            return [stack_10, stack_20]
            """
    """
    if ktr is None:
        l_ds = [gdal.Open(img).ReadAsArray() for img in l_img]
        stack_ds = np.stack(l_ds, axis=0)
        return stack_ds

    else:
        p_fre = re.compile(r'(_+[A-Z]{3,4}?_)')
        l_clean = []
        for img in l_img:

            print(img[re.search(p_fre, img.split('/')[-1]).span()[0]+1:re.search(p_fre, img.split('/')[-1]).span()[0] + 4])
            ds = gdal.Open(img).ReadAsArray()
            l_clean.append(ds)
        #l_clean = [gdal.Open(img).ReadAsArray() for img in l_img if img[re.search(p_fre, img).span()[0] + 1:re.search(p_fre, img).span()[0] + 4] == str(ktr)]
        #print(l_clean)
        stack_ds_clean = np.stack(l_clean, axis=0)
        #return stack_ds_clean
    """

    @staticmethod
    def write_ras_1(ds_np, OUTPUT, Newkey, PathRasRef, dtype=gdalconst.GDT_Int32,
                    setGeotransform=False):  # A adapter !!!!!!!!!!!!!!!!!!!!
        """

        :param ds_np:
        :param OUTPUT:
        :param Newkey:
        :param PathRasRef:
        :param dtype:
        :param setGeotransform:
        :return:
        """

        if setGeotransform is False:

            try:

                ds = gdal.Open(PathRasRef)

                driver = gdal.GetDriverByName('GTiff')

                if len(ds_np.shape) > 2:

                    p = driver.Create(str(OUTPUT + Newkey), ds.RasterXSize, ds.RasterYSize, ds_np.shape[0], dtype)

                    p.SetGeoTransform(ds.GetGeoTransform())

                    p.SetProjection(ds.GetProjection())

                    count = 1

                    for e in range(ds_np.shape[0]):
                        p.GetRasterBand(count).WriteArray(ds_np[(count - 1), :, :])
                        count += 1

                    # p.FlushCache()

                else:

                    p = driver.Create(OUTPUT + Newkey, ds.RasterXSize, ds.RasterYSize, 1, dtype)

                    p.SetGeoTransform(ds.GetGeoTransform())

                    p.SetProjection(ds.GetProjection())

                    p.GetRasterBand(1).WriteArray(ds_np)

                    # p.FlushCache()

            except KeyError:

                print("Problème Key Error")

        else:

            try:

                Geotransform = input('Enter a list of raster features:')

                ds = gdal.Open(PathRasRef)

                driver = gdal.GetDriverByName('GTiff')

                if len(ds_np.shape) > 2:

                    p = driver.Create(OUTPUT + Newkey, ds.RasterXSize, ds.RasterYSize, ds_np.shape[0], dtype)

                    p.SetGeoTransform(Geotransform)

                    p.SetProjection(ds.GetProjection())

                    count = 1

                    for e in range(ds_np.shape[0]):
                        p.GetRasterBand(count).WriteArray(ds_np[(count - 1), :, :])
                        count += 1

                    # p.FlushCache()

                else:

                    p = driver.Create(OUTPUT + Newkey, ds.RasterXSize, ds.RasterYSize, 1, dtype)

                    p.SetGeoTransform(Geotransform)

                    p.SetProjection(ds.GetProjection())

                    p.GetRasterBand(1).WriteArray(ds_np)

                    # p.FlushCache()

            except KeyError:

                print('Key Error')

    @staticmethod
    # Pour créer une image à partir d'une matrice
    def Write_ras(ds_np, path_save, geotranform, projection):

        driver = gdal.GetDriverByName('GTiff')

        if len(ds_np.shape) > 2:

            p = driver.Create(path_save, ds_np.shape[2], ds_np.shape[1], ds_np.shape[0], gdalconst.GDT_Float64)

            p.SetGeoTransform(geotranform)

            p.SetProjection(projection)

            count = 1

            for e in range(ds_np.shape[0]):
                p.GetRasterBand(count).WriteArray(ds_np[(count - 1), :, :])
                count += 1

            # p.FlushCache()

        else:

            p = driver.Create(path_save, ds_np.shape[1], ds_np.shape[0], 1, gdalconst.GDT_Float64)

            p.SetGeoTransform(geotranform)

            p.SetProjection(projection)

            p.GetRasterBand(1).WriteArray(ds_np)

    @staticmethod  # A TESTER !!!!!!!!!
    def get_point_values(path_img, path_shp):

        ds = gdal.Open(path_img)
        gt = ds.GetGeoTransform()
        rb = ds.GetRasterBand(1)

        ds_shp = ogr.Open(path_shp)
        lyr = ds_shp.GetLayer()
        for feat in lyr:
            geom = feat.GetGeometryRef()
            mx, my = geom.GetX(), geom.GetY()
            px = int((mx - gt[0]) / gt[1])
            py = int((my - gt[3]) / gt[5])

            intval = rb.ReadAsArray(px, py, 1, 1)

    @staticmethod
    def look_difference(ds1, ds2, seuil1, seuil2):

        diff = ds1 - ds2
        reclass = np.where(diff > seuil1, 1, diff)
        reclass2 = np.where(reclass == 0, 0, reclass)
        reclass3 = np.where(reclass2 < seuil2, 2, reclass2)

        return reclass3

    # Particular traitement
    """
    @staticmethod
    def mosaic_planet(path_dos, ktr='AnalyticMS_SR_clip.tif'):

        l_planet = [os.path.join(path_dos, e) for e in os.listdir(path_dos) if
                    e[-(len(ktr)):] == str(ktr)]

        motif_date = re.compile(r'(\d{7,8})+_?(\d{5,6})')

        l_try_date = set([e[re.search(motif_date, e).span()[0]:re.search(motif_date, e).span()[0] + 8] +
                          '_' + e[re.search(motif_date, e).span()[0] + 19:re.search(motif_date, e).span()[0] + 23] for
                          e in
                          l_planet])

        dic = {}
        for e in l_try_date:
            dic[e] = list()

        l_try_date = None

        for e in l_planet:

            pos_date_e = re.search(motif_date, e).span()[0]
            date_e = e[pos_date_e:pos_date_e + 8]
            satID_e = e[pos_date_e + 19:pos_date_e + 23]
            fact_e = e[pos_date_e + 9:pos_date_e + 15]

            for i in l_planet:

                pos_date_i = re.search(motif_date, i).span()[0]
                date_i = i[pos_date_i:pos_date_i + 8]
                satID_i = i[pos_date_i + 19:pos_date_i + 23]
                fact_i = i[pos_date_i + 9:pos_date_i + 15]

                if date_e == date_i and fact_i != fact_e and satID_e == satID_i:
                    dic[str(date_e + '_' + satID_e)].append(rst.open(e))
                    dic[str(date_e + '_' + satID_e)].append(rst.open(i))

        for key, values in dic.items():
            dic[key] = list(set(values))

        for key, values in dic.items():
            try:
                name = 'Mosaic_{}'.format(key)

                mosaic, out_trans = merge(values)

                l_coord_Y = [e.meta['transform'][5] for e in values]
                pos_meta_valid_Y = l_coord_Y.index(max(l_coord_Y))

                l_coord_x = [e.meta['transform'][2] for e in values]
                pos_meta_valid_X = l_coord_x.index(min(l_coord_x))

                recp_y = values[pos_meta_valid_Y].meta.copy()
                recp_y_valid = recp_y['transform'][5]

                recp_x = values[pos_meta_valid_X].meta.copy()
                recp_x_valid = recp_x['transform'][2]

                l_coord_x = None
                l_coord_Y = None

                trans_affine = affine.Affine(3.0, 0.0, recp_x_valid, 0.0, -3.0, recp_y_valid)

                out_meta = {"driver": "GTiff",
                            "height": mosaic.shape[1],
                            "dtype": 'uint16',
                            "nodata": 0.0,
                            "width": mosaic.shape[2],
                            "count": 4,
                            "transform": trans_affine,
                            "crs": CRS.from_epsg(32630)}

                with rst.open(name + '.tif', 'w', **out_meta) as dst:
                    dst.write(mosaic)

            except IndexError:
                continue
    """
    @staticmethod
    def histo_array(array, bins=50):

        array_1d = array.ravel()
        plt.hist(array_1d[np.logical_not(np.isnan(array_1d))], bins=bins)
        plt.show()

    # Filter Method
    @staticmethod
    def Median_filter(matrice):
        filter = sif.median(matrice)
        return filter

    @staticmethod
    def Gaussian_filter(matrice):
        filter = sif.gaussian(matrice)
        return filter

    @staticmethod
    def T_S2_P(repS2, repP, ktr='FRE', day=2):
        """

        :param repS2: Repository of Sentinel-2 raster
        :param repP: Repository of PlanetScope raster
        :param ktr: Key traitement Sentinel-2 FRE or SRE
        :param day: daytime difference
        :return: Temporal association dictionary, which contains the paths to the files
        """
        # Check ktr:
        ty = 'STR'
        if type(ktr) != type(ty):
            raise TypeError('Variable ktr must be a string')
            # Check path :
        if type(repP) != type(ty) or type(repS2) != type(ty):
            raise TypeError('Path must be a string')

            # Pattern S2 (Date and TypeCorrection)
        dp_s2 = re.compile(r'(_+\d{7,8}?-)')
        p_fre = re.compile(r'(_+[A-Z]{3,4}?_)')

        # Pattern Planet (Date)
        dp_p = re.compile('(\d{7,8})')

        # Creation list
        ## List of all image in repository Sentinel 2
        l_s2 = [e for e in os.listdir(repS2) if e[-4:] == '.tif']

        ## List of image Sentinel 2 with angle correction (SRE) or not (FRE)
        l_s2_clean = [e for e in l_s2 if e[re.search(p_fre, e).span()[0] + 1:re.search(p_fre, e).span()[0] + 4] == ktr]

        ## List of all image in repository PlanetScope
        l_pl = [e for e in os.listdir(repP) if e[-(len('.tif')):] == '.tif']

        # Dictionnary
        dct_im = {}

        for e in l_s2_clean:
            pos_dat_s2 = re.search(dp_s2, e).span()[0]
            dat_s2 = e[pos_dat_s2 + 1:pos_dat_s2 + 9]
            dc_s2 = datetime.date(int(dat_s2[:4]), int(dat_s2[4:6]), int(dat_s2[6:]))
            print('SENTINEL 2 :', dc_s2)
            for i in l_pl:
                pos_dat_pl = re.search(dp_p, i).span()[0]
                dat_pl = i[pos_dat_pl:pos_dat_pl + 8]
                dc_pl = datetime.date(int(dat_pl[0:4]), int(dat_pl[4:6]), int(dat_pl[6:]))
                print('PLANET : ', dc_pl)
                if dc_pl == dc_s2:
                    dct_im[str(dc_s2) + '_' + str(dc_pl)] = [os.path.join(repS2, e), os.path.join(repP, i)]
                elif datetime.timedelta(day) >= dc_s2 - dc_pl >= datetime.timedelta(-day):
                    dct_im[str(dc_s2) + '_' + str(dc_pl)] = [os.path.join(repS2, e), os.path.join(repP, i)]

        return dct_im

    @staticmethod
    def crop_adapt(nd1, nd2):
        """

        :param path_img1:
        :param path_img2:
        :return:
        """
        #print(path_img2, "\n",path_img1)
        ds_pl = nd1
        ds_s = nd2

        # Check Projection

        if ds_s.projection != ds_pl.projection:
            raise TypeError('Error projection : The rasters do not have the same projection')

        stack1 = None
        stack2 = None

        # Check dimension raster
        if len(ds_s.shape) > 2:
            stack1 = True
        if len(ds_s.shape) <= 2:
            stack1 = False
        if len(ds_pl.shape) > 2:
            stack2 = True
        if len(ds_pl.shape) <= 2:
            stack2 = False

        # print('Stack1:', stack1, '\n','Stack2:', stack2)

        coor_final_x_s = None
        coor_final_y_s = None

        coor_final_x_pl = None
        coor_final_y_pl = None

        if stack1 is True:
            # print("Stack1 is True")
            coor_final_x_s = ds_s.geodata[0] + ds_s.shape[1] * (ds_s.geodata[1])
            coor_final_y_s = ds_s.geodata[3] + ds_s.shape[2] * (ds_s.geodata[-1])
            # print(coor_final_x_s, coor_final_y_s)
        elif stack1 is False:
            # print("Stack1 is False")
            coor_final_x_s = ds_s.geodata[0] + ds_s.shape[0] * (ds_s.geodata[1])
            coor_final_y_s = ds_s.geodata[3] + ds_s.shape[1] * (ds_s.geodata[-1])
            # print(coor_final_x_s, coor_final_y_s)

        if stack2 is True:
            # print('Stack2 is True')
            coor_final_x_pl = ds_pl.geodata[0] + ds_pl.shape[1] * (ds_pl.geodata[1])
            coor_final_y_pl = ds_pl.geodata[3] + ds_pl.shape[2] * (ds_pl.geodata[-1])
        elif stack2 is False:
            # print("Stack2 is False")
            coor_final_x_pl = ds_pl.geodata[0] + ds_pl.shape[0] * (ds_pl.geodata[1])
            coor_final_y_pl = ds_pl.geodata[3] + ds_pl.shape[1] * (ds_pl.geodata[-1])

        # print('Coordonnées Sentinel-2:', coor_final_x_s, coor_final_y_s)
        # print('Coordonnées Planet:', coor_final_x_pl, coor_final_y_pl)

        # Creation BBOX for Crop raster

        # print("Geodata_Sentinel2 : ", ds_s.geodata[0], ds_s.geodata[3])
        # print("Geodata_Planet : ", ds_pl.geodata[0], ds_pl.geodata[3])

        # Cas où les coordonnées d'origine de l'img2 sont en dehors de l'img1 (orX et orY de img2 sont à l'Ouest-Sud-Ouest de l'img1)
        if ds_s.geodata[0] < ds_pl.geodata[0] and ds_s.geodata[3] < ds_pl.geodata[3]:
            print('SUD-OUEST')
            # Cas où le dernier pixel de l'img2 est dans l'img1
            if coor_final_y_s > coor_final_y_pl:
                print('coor_final is IN')
                bbox = [ds_s.geodata[0], ds_pl.geodata[3], coor_final_x_pl, coor_final_y_pl]
            # Cas où le dernier pixel de l'img2 est en dehors de l'img1
            else:
                print('coor_final is OUT')
                bbox = [ds_s.geodata[0], ds_pl.geodata[3], coor_final_x_pl, coor_final_y_s]

        # Cas où les coordonnées d'origine de l'img2 sont en dehor de l'img1 (orX et orY de l'img2 sont au Nord-Est de l'img1)
        elif ds_s.geodata[0] > ds_pl.geodata[0] and ds_s.geodata[3] > ds_pl.geodata[3]:
            print('NORD-EST')
            # Cas où le dernier pixel est dans l'img1
            if coor_final_y_pl < coor_final_y_s and coor_final_x_pl < coor_final_x_s:
                print('coor_final is IN')
                bbox = [ds_s.geodata[0], ds_pl.geodata[3], coor_final_x_pl, coor_final_y_pl]

            else:
                print('coor_final is OUT')
                bbox = [ds_s.geodata[0], ds_pl.geodata[3], coor_final_x_pl, coor_final_y_s]

        # Cas où les coordonnées d'origine de l'img2 sont en dehor de l'img1 (orX et orY de l'img2 sont au Nord-Ouest de l'img1)
        elif ds_s.geodata[0] > ds_pl.geodata[0] and ds_s.geodata[3] < ds_pl.geodata[3]:

            print('NORD-OUEST')

            if coor_final_y_pl < coor_final_y_s and coor_final_x_pl < coor_final_x_s:
                print('coor_final is IN')
                bbox = [ds_pl.geodata[0], ds_s.geodata[3], coor_final_x_pl, coor_final_y_pl]

            else:
                print('coor_final is OUT')
                bbox = [ds_pl.geodata[0], ds_s.geodata[3], coor_final_x_pl, coor_final_y_s]

        # Cas où les coordonnées d'origine de l'img2 sont en dehor de l'img1 (orX et orY de l'img2 sont au Sud-Est de l'img 1)
        elif ds_s.geodata[0] < ds_pl.geodata[3] and ds_s.geodata[3] > ds_pl.geodata[3]:

            print('SUD-EST')
            if coor_final_y_pl < coor_final_y_s and coor_final_x_pl < coor_final_x_s:
                print('coor_final is IN')
                bbox = [ds_pl.geodata[0], ds_pl.geodata[3], coor_final_x_pl, coor_final_y_pl]

            else:
                print('coor_final is OUT')
                bbox = [ds_pl.geodata[0], ds_pl.geodata[3], coor_final_x_s, coor_final_y_pl]

        ##### Trois coté corrigé mais je sais plus quel coté il manque

        elif ds_s.geodata[0] == ds_pl.geodata[0] and ds_s.geodata[3] == ds_pl.geodata[3]:
            print('Partiulier')

            if ds_s.shape == ds_pl.shape and ds_s.resolution == ds_pl.resolution:
                print('Pas besoin de Crop')

                pass



        ds_s.crop_ras(bbox)
        ds_pl.crop_ras(bbox)
        # cropS2 = np.where(ds_s.array_crop <= 0, np.nan, ds_s.array_crop)
        # crop_pl = np.where(ds_pl.array_crop <= 0, np.nan, ds_pl.array_crop)
        return ds_s.array_crop, ds_pl.array_crop, bbox

    @staticmethod
    def seuil_index():
        print("CHEH")

    @staticmethod
    def findfile(name_img, path):
        for dirpath, dirname, filename in os.walk(path):
            if name_img in filename:
                return os.path.join(dirpath, name_img)

    @staticmethod
    def cloud_mask_Sentinel2(ds):

        ci1 = (ds[2, :, :] + ds[3, :, :] + ds[4, :, :]) / (
                ds[6, :, :] + ds[12, :, :] * 2)

        ci2 = (ds[2, :, :] + ds[3, :, :] + ds[4, :, :] + ds[6, :, :] + ds[12, :, :] +
               ds[13, :, :]) / 6

        classifci1 = np.where(ci1 < 0.1, 5000, 0)

        meanci2 = np.mean(ci2[~np.isnan(ci2)])

        maxci2 = np.max(ci2[~np.isnan(ci2)])

        T2 = meanci2 + (0.1 * (maxci2 - meanci2))

        classifci2 = np.where(ci2 > T2, 5000, 0)

        finalclassif = classifci1 + classifci2

        return finalclassif

    @staticmethod
    def seuil_modulable(seuil_inf, seuil_sup, matrice):
        matrice_seuil = np.where(matrice < seuil_inf, seuil_inf, matrice)
        matrice_seuil = np.where(matrice_seuil > seuil_sup, seuil_sup, matrice_seuil)

        return matrice_seuil

    @staticmethod
    def seuil_index_norma(matrice):
        matrice_seuil = np.where(matrice, -1 < -1, matrice)
        matrice_seuil = np.where(matrice_seuil > 1, 1, matrice_seuil)

        return matrice_seuil

    """
    @staticmethod
    def full_index_vegetation_Sentinel2(stack_path, path_enregistrement):
        img = nd(stack_path)
        stack = img.band_Array

        R490 = stack[0, :, :] / 10000  # Bande Bleu
        R490 = seuil_modulable(0, 10000, R490)

        R550 = stack[1, :, :] / 10000  # Bande Verte
        R550 = seuil_index(0, 10000, R550)
        # R550 = np.where(R550 == -10000, np.nan, R550)

        R670 = stack[2, :, :] / 10000  # Bande Rouge
        R670 = seuil_index(0, 10000, R670)
        # R670 = np.where(R670 == -10000, np.nan, R670)

        R700 = stack[3, :, :] / 10000  # Bande RE1
        R700 = seuil_index(0, 10000, R700)
        # R700 = np.where(R700 == -10000, np.nan, R700)

        R740 = stack[4, :, :] / 10000  # Bande RE2
        R740 = seuil_index(0, 10000, R740)
        # R740 = np.where(R740 == -10000, np.nan, R740)

        R780 = stack[5, :, :] / 10000  # Bande RE3
        R780 = seuil_index(0, 10000, R780)
        # R780 = np.where(R780 == -10000, np.nan, R780)

        R800 = stack[6, :, :] / 10000  # Bande NIR
        R800 = seuil_index(0, 10000, R800)
        # R800 = np.where(R800 == -10000, np.nan, R800)

        R860 = stack[7, :, :] / 10000  # Bande RE4
        R860 = seuil_index(0, 10000, R860)
        # R860 = np.where(R860 == -10000, np.nan, R860)

        # Calcul CCI
        CCII = (3 * ((R700 - R670) - (0.2 * (R700 - R550) * (R700 / R670)))) / (
                ((1 + 0.16) * (R800 - R670)) / (R800 + R670 + 0.16))

        # Calcul STVI_norm
        SAT = 0.5 * ((105 * (R700 - R550)) - (145 * (R670 - R550)))
        SRT = 0.5 * ((125 * (R780 - R740)) - (43 * (R860 - R740)))
        STVI_norm = (SAT - SRT) / (SAT + SRT)
        # STVI_norm = seuil_index_norma(STVI_norm)

        # Calcul NDVIRE
        NDVIRE = (R740 - R700) / (R740 + R700)
        # NDVIRE = seuil_index_norma(NDVIRE)

        # Calcul DVI
        DVI = R800 / R670

        # Calcul GLI
        GLI = ((2 * R550) - R670 - R490) / ((2 * R550) + R670 + R490)
        # GLI = seuil_index_norma(GLI)

        # Calcul MCARI2
        MCARI2 = (1.5 * ((2.5 * (R800 - R670) - 1.3 * (R800 - R550))) /
                  (np.sqrt(np.power((2 * R800 + 1), 2) - (6 * R800 - 5 * np.sqrt(R670)) - 0.5)))

        # Calcul IRECI
        IRECI = (R800 - R670) / (R700 / R740)

        # Calcul NDVI
        NDVI = (R800 - R670) / (R800 + R670)
        # NDVI = seuil_index_norma(NDVI)

        # Calcul SeLI
        SelI = (R860 - R700) / (R860 + R700)
        # SelI = seuil_index_norma(SelI)

        # Libération de mémoire

        R490 = None
        R550 = None
        R670 = None
        R700 = None
        R740 = None
        R780 = None
        R800 = None
        R860 = None

        # Matrice des indices organisé en liste pour pouvoir stack
        l_index = [CCII, STVI_norm, NDVIRE, DVI, GLI, MCARI2, IRECI, NDVI, SelI]

        index_stack = np.stack(l_index, axis=0)

        # Ecriture du raster
        return index_stack
    """

class ILearn(Norma_Data):

    def __init__(self, path_img, filter = False):

        img = Norma_Data.__init__(self,path_img)
        table = None
        if filter is False:
            if len(img.shape) <= 2:

                band = img.band_Array
                ds_ravel = np.ravel(band)
                # ds_ravel = ds_ravel.reshape(-1, 1)
                dic = {'Bande1': ds_ravel}
                table = pd.DataFrame(data=dic)


                self.table_img = table
            else:

                nb_bande = img.shape[2]
                if nb_bande >= 8:
                    print(img.shape)
                    nb_bande = img.shape[0]

                dic = {}

                for e in range(nb_bande):
                    band = img.band_Array[e, :, :]
                    ds_ravel = np.ravel(band)
                    # ds_ravel = ds_ravel.reshape(-1, 1)
                    dic['Bande{}'.format(e + 1)] = ds_ravel
                print(dic)
                table = pd.DataFrame(data=dic)



                self.table_img = table
        else:
            img.gaussian_filter_meth()
            if len(img.shape) <= 2:

                band = img.gaussian_filter_array
                ds_ravel = np.ravel(band)
                # ds_ravel = ds_ravel.reshape(-1, 1)
                dic = {'Bande1': ds_ravel}
                table = pd.DataFrame(data=dic)


                self.table_img_filter = table
            else:

                nb_bande = img.shape[2]
                if nb_bande >= 8:
                    print(img.shape)
                    nb_bande = img.shape[0]

                dic = {}

                for e in range(nb_bande):
                    band = img.gaussian_filter_array[e, :, :]
                    ds_ravel = np.ravel(band)
                    # ds_ravel = ds_ravel.reshape(-1, 1)
                    dic['Bande{}'.format(e + 1)] = ds_ravel
                print(dic)
                table = pd.DataFrame(data=dic)


                self.table_img_filter = table
"""
class ImgData:

    def __init__(self, pathdossier, product, zone, lieu):


        xMin, xMax, yMin, yMax = zone

        if product == 'S2':

            try:

                path_band20 = [os.path.join(pathdossier,file) for file in os.listdir(pathdossier)]

                #date = path_band20[1][-48 :-40]

                dic_band = {}

                for e in path_band20 :
                    ds = gdal.Open(e)
                    ds_gt = ds.GetGeoTransform()
                    row1 = int((yMax - ds_gt[3]) / ds_gt[5])
                    col1 = int((xMin - ds_gt[0]) / ds_gt[1])
                    row2 = int((yMin - ds_gt[3]) / ds_gt[5])
                    col2 = int((xMax - ds_gt[0]) / ds_gt[1])
                    nameBand = e[-11 :-8]
                    nArray = np.array(ds.ReadAsArray(col1, row1, col2 - col1 + 1, row2 - row1 + 1).astype(np.float32))
                    dic_band[nameBand] = nArray
                    ds = None



                self.Blue = dic_band["B02"]
                self.Green = dic_band["B03"]
                self.Red = dic_band["B04"]
                self.VRE1 = dic_band["B05"]
                self.VRE2 = dic_band["B06"]
                self.VRE3 = dic_band["B07"]
                self.NIR = dic_band["B08"]
                self.B8A = dic_band["B8A"]
                self.SWIR1 = dic_band["B11"]
                self.SWIR2 = dic_band["B12"]
                self.AOT = dic_band["AOT"]
                self.SCL = dic_band["SCL"]
                self.TCI =dic_band["TCI"]
                self.VIS = dic_band["VIS"]
                self.WVP=dic_band["WVP"]
                self.product = product
                #self.date = date
                self.lieu = lieu


            except NotADirectoryError :
                print("Chemin d'accès vers un répertoire attendue,pour Sentinel2")
                pass

            except AttributeError :
                print(
                    "Vérifier la correspondance entre l'image et les coordonnées géographique du découpage, pour Sentinel2")
                pass

        elif product == "ORTHO":

            try :
                path_ortho = [os.path.join(pathdossier, file) for file in os.listdir(pathdossier)]

                # dic_ortho = {k[-7:-11]:gdal.Open(k).ReadAsArray().astype(np.float16) for k in nameImg}
                #date = input('Rentrer la date d\'acquisition de l\'ortho : ')

                dic_ortho = {}

                for ortho in path_ortho :
                    ds = gdal.Open(ortho)
                    ds_gt = ds.GetGeoTransform()
                    row1 = int((yMax - ds_gt[3]) / ds_gt[5])
                    col1 = int((xMin - ds_gt[0]) / ds_gt[1])
                    row2 = int((yMin - ds_gt[3]) / ds_gt[5])
                    col2 = int((xMax - ds_gt[0]) / ds_gt[1])
                    nArray = np.array(ds.ReadAsArray(col1, row1, col2 - col1 + 1, row2 - row1 + 1).astype(np.float32))
                    dic_ortho['REN'] = nArray
                    ds = None

                self.Blue = dic_ortho['REN'][0]
                self.Green = dic_ortho['REN'][1]
                self.Red = dic_ortho['REN'][2]
                self.Div = dic_ortho['REN'][3]
                self.product = product
                #self.date = date
                self.lieu = lieu

            except NotADirectoryError :
                print("Chemain d'accès vers un répertoire attendue, pour l'orthophotoplan")
                pass

            except AttributeError :
                print(
                    "Vérifier la correspondance entre l'image et les coordonnées géographique du découpage, pour l'orthophotplan")
                pass

            except KeyError :
                print("Vous n'avez pas la bonne option sur les produits satellitaires, pour L'orthophotplan")
                pass

        elif product == 'UA':

            try:
                path_band20 = [os.path.join(pathdossier,file) for file in os.listdir(pathdossier)]

                #date = path_band20[1][-48:-40]

                dic_band = {}
                for e in path_band20:
                    ds = gdal.Open(e)
                    ds_gt = ds.GetGeoTransform()
                    row1 = int((yMax - ds_gt[3]) / ds_gt[5])
                    col1 = int((xMin - ds_gt[0]) / ds_gt[1])
                    row2 = int((yMin - ds_gt[3]) / ds_gt[5])
                    col2 = int((xMax - ds_gt[0]) / ds_gt[1])
                    nameBand = 'UA'
                    nArray = np.array(ds.ReadAsArray(col1, row1, col2 - col1 + 1, row2 - row1 + 1).astype(np.float32))
                    dic_band[nameBand] = nArray
                    ds = None

                # dic = {k[-7:-4]: gdal.Open(k).ReadAsArray().astype(np.float16) for k in path_band20}  Option pour charger toute l'image

                self.classif = dic_band['UA']
                self.product = product
                self.lieu =lieu

            except NotADirectoryError :
                print("Chemin d'accès vers un répertoire attendue,pour Sentinel2")
                pass

            except AttributeError :
                print("Vérifier la correspondance entre l'image et les coordonnées géographique du découpage, pour Sentinel2")
                pass

            except KeyError :
                print("Vous n'avez pas la bonne option sur les produits satellitaires, pour UA")
                pass

        elif product == 'S1':

            try:

                ds = gdal.Open(pathdossier)

                date = pathdossier.split('_')[3][0:8]


                self.Sigma0VH = np.array(ds.GetRasterBand(0))
                self.Sigma0VV = np.array(ds.GetRasterBand(2))
                self.product = product
                self.lieu = lieu
                self.date = date
                self.namefile = pathdossier[:-4]


            except NotADirectoryError :
                print("Chemin d'accès vers un répertoire attendue,pour Sentinel2")
                pass

            except AttributeError :
                print(
                    "Vérifier la correspondance entre l'image et les coordonnées géographique du découpage, pour Sentinel2")
                pass


        elif product == "Hyper":

            try:



                with rasterio.open(pathdossier) as src:

                    out_meta = src.meta

                    l_band = [src.read(e) for e in range(1, len(src.indexes))]

                    stack = np.stack(l_band,axis=0)

                self.ds = src
                self.product = product
                self.stack = stack
                self.shape = np.shape(src)
                self.lieu = lieu
                self.projection = out_meta
                self.resolution = src.transform[0]
                self.namefile = pathdossier[:-4]


            except NotADirectoryError:
                print(f"Chemin d'accès vers un répertoire attendue,\n {pathdossier} semble ne pas convenir")
                pass


    def CropVigSta(self, TailleVignette, Path_Output, Zcouvert):

        Dic = {}
        try:

            ratio_taillevignette = [e for e in range(2000) if e % TailleVignette == 0]

            print(f"Le ratio de la taille vignette est {ratio_taillevignette}")

            if self.product == 'S2':

                ar10 = [self.Blue, self.Green, self.Red, self.NIR,self.VRE1, self.VRE2, self.VRE3, self.B8A, self.SWIR1, self.SWIR2]

                count = 0

                for x in range(ceil(np.shape(self.Red)[0] / (ratio_taillevignette[2] * Zcouvert))):

                    for y in range(ceil(np.shape(self.Red)[1] / (ratio_taillevignette[2] * Zcouvert))):

                        count += 1

                        countband = 2

                        for e in ar10:

                            try :
                                if x == 0:

                                    vign = e[int((x * ratio_taillevignette[2])) : int((x + 1) * ratio_taillevignette[2]), int(((y * ratio_taillevignette[2]) * Zcouvert)) :int(
                                               (y * ratio_taillevignette[2]) * Zcouvert) + (ratio_taillevignette[2])]

                                    if vign.shape[0] == ratio_taillevignette[2] and vign.shape[1] == \
                                            ratio_taillevignette[2] :

                                        new_key = str(self.product) + '_' + str(
                                            vign.shape[0]) + '_' + str('10m') \
                                                  + '_' + 'B' + str(countband) + '_' + (
                                                              4 - len(str(count))) * '0' + str(count) + '_' + str(
                                            self.lieu)

                                        Image.fromarray(vign).save(os.path.join(Path_Output, new_key) + '.tif')

                                        countband += 1

                                    else:
                                        continue
                                else:

                                    vign = e[int((x * ratio_taillevignette[2]) * Zcouvert):int(
                                        ((x * ratio_taillevignette[2]) * Zcouvert) + (ratio_taillevignette[2])),
                                           int(((y * ratio_taillevignette[2]) * Zcouvert)):int(
                                               ((((y)) * ratio_taillevignette[2]) * Zcouvert) + (
                                               ratio_taillevignette[2]))]

                                    if vign.shape[0] == ratio_taillevignette[2] and vign.shape[1] == \
                                            ratio_taillevignette[2] :

                                        new_key = str(self.product) + '_' + str(
                                            vign.shape[0]) + '_' + str('10m') \
                                                  + '_' + 'B' + str(countband) + '_' + (
                                                              4 - len(str(count))) * '0' + str(count) + '_' + str(
                                            self.lieu)

                                        Image.fromarray(vign).save(os.path.join(Path_Output, new_key) + '.tif')

                                        countband += 1
                                    else :
                                        continue

                            except UserWarning :
                                print("Except activé")
                                continue

            elif self.product == 'ORTHO':

                ar = [self.Blue, self.Green, self.Red]

                Count = 0

                for x in range(ceil(np.shape(self.Red)[0] / (ratio_taillevignette[8] * Zcouvert))) :

                    for y in range(ceil(np.shape(self.Red)[1] / (ratio_taillevignette[8] * Zcouvert))) :

                        Count += 1

                        countband = 1

                        for e in ar :

                            if x == 0 :

                                vign = e[int((x * ratio_taillevignette[8])) :int((x + 1) * ratio_taillevignette[8]),
                                       int(((y * ratio_taillevignette[8]) * Zcouvert)) :int(
                                           ((y) * ratio_taillevignette[8]) * Zcouvert) + (ratio_taillevignette[8])]

                                if vign.shape[0] == ratio_taillevignette[8] and vign.shape[1] == ratio_taillevignette[
                                    8] :

                                    new_key = str(self.product) + '_' + str(vign.shape[0]) \
                                              + '_' + str('2m') + '_' + 'B' + str(countband) + '_' + (
                                                          4 - len(str(Count))) * '0' + str(Count) + '_' + str(self.lieu)

                                    Image.fromarray(vign).save(os.path.join(Path_Output, new_key) + '.tif')

                                    countband += 1
                                else :
                                    continue
                            else :

                                vign = e[int((x * ratio_taillevignette[8]) * Zcouvert) :int(
                                    ((x * ratio_taillevignette[8]) * Zcouvert) + (ratio_taillevignette[8])),
                                       int(((y * ratio_taillevignette[8]) * Zcouvert)) :int(
                                           ((((y)) * ratio_taillevignette[8]) * Zcouvert) + (
                                           ratio_taillevignette[8]))]

                                if vign.shape[0] == ratio_taillevignette[8] and vign.shape[1] == ratio_taillevignette[
                                    8] :
                                    new_key = str(self.product) + '_' + str(
                                        vign.shape[0]) + '_' + str('2_m') + '_' + 'B' + str(countband) + '_' + (
                                                          4 - len(str(Count))) * '0' + str(Count) + '_' + str(self.lieu)

                                    Image.fromarray(vign).save(os.path.join(Path_Output, new_key) + '.tif')

                                    countband += 1

            elif self.product == 'UA':

                ar = self.classif

                count = 0

                for x in range(ceil(np.shape(self.classif)[0] / (ratio_taillevignette[8] * Zcouvert))):

                    for y in range(ceil(np.shape(self.classif)[1] / (ratio_taillevignette[8] * Zcouvert))):

                        countband = 1
                        count += 1

                        if x == 0:

                            vign = ar[int((x * ratio_taillevignette[8])): int((x + 1) * ratio_taillevignette[8]), int(((y * ratio_taillevignette[8]) * Zcouvert)):int(
                                       (y * ratio_taillevignette[8]) * Zcouvert) + (ratio_taillevignette[8])]

                            if vign.shape[0] == ratio_taillevignette[8] and vign.shape[1] == ratio_taillevignette[8]:

                                new_key = str(self.product)+ '_' + str(
                                    vign.shape[0]) + '_' + str('2_5m') + '_' + 'B' + str(countband) + '_' + (
                                                  4 - len(str(count))) * '0' + str(count) + '_' + str(self.lieu)

                                Image.fromarray(vign).save(os.path.join(Path_Output, new_key) + '.tif')

                                countband += 1

                            else:
                                pass
                        else:

                            vign = ar[int((x * ratio_taillevignette[8]) * Zcouvert):int(
                                ((x * ratio_taillevignette[8]) * Zcouvert) + (ratio_taillevignette[8])),
                                   int(((y * ratio_taillevignette[8]) * Zcouvert)):int(
                                       ((y * ratio_taillevignette[8]) * Zcouvert) + (ratio_taillevignette[8]))]

                            if vign.shape[0] == ratio_taillevignette[8] and vign.shape[1] == ratio_taillevignette[8]:
                                new_key = str(self.product) + '_' + str(
                                    vign.shape[0]) + '_' + str('10m') + '_' + 'B' + str(countband) + '_' + (
                                                  4 - len(str(count))) * '0' + str(count) + '_' + str(self.lieu)

                                Image.fromarray(vign).save(os.path.join(Path_Output, new_key) + '.tif')

                                countband += 1

            elif self.product == 'S1':

                ar10 = [self.Sigma0VH, self.Sigma0VV]

                count = 0

                for x in range(ceil(np.shape(self.Sigma0VH)[0] / (ratio_taillevignette[2] * Zcouvert))):

                    for y in range(ceil(np.shape(self.Sigma0VH)[1] / (ratio_taillevignette[2] * Zcouvert))):

                        count += 1

                        countband = 0

                        for e in ar10:

                            try:
                                if x == 0:

                                    vign = e[int((x * ratio_taillevignette[2])): int((x + 1) * ratio_taillevignette[2]),
                                           int(((y * ratio_taillevignette[2]) * Zcouvert)):int(
                                               (y * ratio_taillevignette[2]) * Zcouvert) + (ratio_taillevignette[2])]

                                    if vign.shape[0] == ratio_taillevignette[2] and vign.shape[1] == \
                                            ratio_taillevignette[2]:
                                        # Changer Newkey pour modification du nom
                                        new_key = str(self.product) + '_' + str(
                                            vign.shape[0]) + '_' + str('10m') \
                                                  + '_' + 'B' + str(countband) + '_' + (
                                                          4 - len(str(count))) * '0' + str(count) + '_' + str(
                                            self.lieu+str(self.date))

                                        Image.fromarray(vign).save(os.path.join(Path_Output, new_key) + '.tif')

                                        countband += 1

                                    else:
                                        continue
                                else:

                                    vign = e[int((x * ratio_taillevignette[2]) * Zcouvert):int(
                                        ((x * ratio_taillevignette[2]) * Zcouvert) + (ratio_taillevignette[2])),
                                           int(((y * ratio_taillevignette[2]) * Zcouvert)):int(
                                               ((((y)) * ratio_taillevignette[2]) * Zcouvert) + (
                                                   ratio_taillevignette[2]))]

                                    if vign.shape[0] == ratio_taillevignette[2] and vign.shape[1] == \
                                            ratio_taillevignette[2]:

                                        new_key = str(self.product) + '_' + str(
                                            vign.shape[0]) + '_' + str('10m') \
                                                  + '_' + 'B' + str(countband) + '_' + (
                                                          4 - len(str(count))) * '0' + str(count) + '_' + str(
                                            self.lieu + str(self.date))

                                        Image.fromarray(vign).save(os.path.join(Path_Output, new_key) + '.tif')

                                        countband += 1
                                    else:
                                        continue

                            except UserWarning:
                                print("Except activé")
                                continue

            elif self.product == "Hyper":

                count = 0


                for x in range(ceil(np.shape(self.stack[0])[0] / (ratio_taillevignette[1] * Zcouvert))):

                    for y in range(ceil(np.shape(self.stack[0])[1] / (ratio_taillevignette[1] * Zcouvert))):

                        count += 1

                        countband = 1

                        for e in self.stack:

                            try:
                                if x == 0:

                                    vign = e[int((x * ratio_taillevignette[1])): int((x + 1) * ratio_taillevignette[1]),
                                           int(((y * ratio_taillevignette[1]) * Zcouvert)):int(
                                               (y * ratio_taillevignette[1]) * Zcouvert) + (ratio_taillevignette[1])]

                                    #print(vign)

                                    if vign.shape[0] == ratio_taillevignette[1] and vign.shape[1] == \
                                            ratio_taillevignette[1]:

                                        #print(vign)

                                        new_key = str(self.product) + '_' + str(
                                            vign.shape[0]) + '_' +str(int(self.resolution)) \
                                                  + '_' + 'B' + str(countband) + '_' + (
                                                          4 - len(str(count))) * '0' + str(count) + '_' + str(self.lieu)

                                        Dic[new_key] = vign

                                        countband += 1
                                        vign = None

                                    else:
                                        continue
                                else:

                                    vign = e[int((x * ratio_taillevignette[1]) * Zcouvert):int(
                                        ((x * ratio_taillevignette[1]) * Zcouvert) + (ratio_taillevignette[1])),
                                           int(((y * ratio_taillevignette[1]) * Zcouvert)):int(
                                               ((((y)) * ratio_taillevignette[1]) * Zcouvert) + (
                                                   ratio_taillevignette[1]))]

                                    if vign.shape[0] == ratio_taillevignette[1] and vign.shape[1] == \
                                            ratio_taillevignette[1]:

                                        new_key = str(self.product) + '_' + str(
                                            vign.shape[0]) + '_' + str(self.resolution) \
                                                  + '_' + 'B' + str(countband) + '_' + (
                                                          4 - len(str(count))) * '0' + str(count) + '_' + str(
                                            self.lieu)

                                        Dic[new_key] = vign


                                        countband += 1
                                        vign = None

                                    else:
                                        continue

                            except UserWarning:
                                print("Except activé UserWarning !!!")
                                continue



        except NotADirectoryError :
            print("Crop Statique ne marche pas!!!")
            sys.exit(1)
        except AttributeError:

            print("CHEH PAS BON DANS LE DIVERS")
        print(Dic)
        with open('Data_Vignette_Luciole.json','w') as fp:
            json.dump(Dic, fp)
    def tri_vignette(self, path_dossier):

        try:

            name_vign = [fichier for fichier in os.listdir(path_dossier) if fichier[-4:] == ".tif"]

            motif = re.compile(r'_+0{1,5}[0-9]')

            motif_lieu = re.compile(r"_+[A-Z]+[a-z]")

            liste_idscene = []

            for e in name_vign:
                pos = re.search(motif, e).span()[0]
                #print(e[21:24])
                pos_lieu = re.search(motif_lieu, e).span()[0]
                idscene = e[pos + 1:pos_lieu]
                liste_idscene.append(idscene)

            l = list(set(liste_idscene))

            path = [os.path.join(path_dossier, path) for path in l]

            for k in l:
                res = Path(os.path.join(path_dossier, k))
                if res.exists() & res.is_dir():
                    continue
                else:
                    res.mkdir()
                    continue

            regroup_id = {}
            for e in path:
                regroup_id[e] = list()
                for i in name_vign:
                    pos = re.search(motif, i).span()[0]
                    pos_lieu = re.search(motif_lieu, i).span()[0]
                    if i[pos + 1 :pos_lieu] == e[-4:]:
                        regroup_id[e].append(i)

                    # print(regroup_id, '\n', len(regroup_id))

            for i in regroup_id.items():
                for e in range(len(i[1])):
                    pat = i[0]
                    file = i[1][e]
                    sh.move(path_dossier + '/' + file, pat)


        except ValueError:
            print('Pas bon')
            pass
        except NotADirectoryError:
            print("Chemin d'accès vers un répertoire attendue")
            sys.exit(1)
        except FileNotFoundError:
            print('T\'es Mauvais Jack !!!!!')
            pass

    def data_recap(self, path_dossier, seuil):

        file_exist = os.path.exists(path_dossier + '/' + 'bilan.csv')

        if file_exist is False:

            try:
                listpath = []
                for u in [os.listdir(chemin) for chemin in [os.path.join(path_dossier, rep) for rep in os.listdir(path_dossier)]]:
                    for e in u:
                        if e[0:3] == 'IMP':
                            idscene = e[-15:-11]
                            listpath.append(os.path.join(path_dossier + '/' + idscene, e))

                list_file = []
                list_product = []
                list_invproduct = []
                for file in listpath:
                    ds = gdal.Open(file).ReadAsArray().astype(np.float32)
                    img_seuil = np.where(ds > seuil, 255, 0)
                    product = (np.count_nonzero(img_seuil) / (np.shape(img_seuil)[0] * np.shape(img_seuil)[1])) * 100
                    list_file.append(os.path.split(file)[1][-15:-11])
                    list_product.append(product)
                    list_invproduct.append(abs(100 - product))

                dico = {'Identifiant image': list_file, 'Taux d\'imperméabilisation': list_product,
                        'Taux non imperméabilisé': list_invproduct}
                tableur = pd.DataFrame(dico)
                tableur.to_csv(path_dossier + '/' + 'bilan.csv')

            except KeyError:
                print("Problème dans le try !!!")
                pass

        if file_exist is True:
            print('Le fichier bilan.csv existe')
            pass

    def stack_all_scene(self, path_dossier):
        try:

            name_repertoire_scene = [os.path.join(path_dossier, repertoire) for repertoire    in os.listdir(path_dossier)]

            dic = {}

            for e in name_repertoire_scene:

                list_elem = os.listdir(e)

                list_array_file = []

                for i in list_elem:

                    path_file = e+"/"+i

                    ds = gdal.Open(path_file).ReadAsArray()

                    list_array_file.append(ds)

                dic[str(e)] = list_array_file

            name_repertoire_scene = None
            list_elem = None
            list_array_file = None

            dic_stack = {}

            for y, u in dic.items():

                stack = np.stack(u, axis=0)

                print(f"La shape du stack est : {np.shape(stack)}")
                print(f"le Y de dic.items() est : {y}")

                dic_stack[y] = stack

            dic = None

            self.dic_stack = dic_stack

        except ValueError:
            print("ValueErrror t'es Mauvias Jack !!!!")

    def all_path_scene(self, path_vignette):

        repo = [os.path.join(path_vignette,e) for e in os.listdir(path_vignette)]

        dic ={}

        for e in repo:
            list_files = os.listdir(e)
            dic[str(e)] = list_files


        self.path_vign = dic

# =========================== STATIC METHOD ===========================
    @staticmethod
    def class_imper(pathdossier,nclass):
        try:
            data = pd.read_csv(pathdossier+'/'+'bilan.csv')
            num_data = data.to_numpy()
            num_data_clear = num_data[:, 2:]
            classif = KMeans(n_clusters=nclass, random_state=0).fit(num_data_clear)
            predi = pd.DataFrame(classif.predict(num_data_clear))
            data.insert(4, "CLASS",predi)
            data.to_csv(pathdossier +'/'+'bilan.csv')

        except ValueError:
            print("CHEH")
            pass
"""