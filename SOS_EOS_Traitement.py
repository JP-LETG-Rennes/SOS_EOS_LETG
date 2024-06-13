import json
from pathlib import Path
import datetime
import os
import re
import time
from osgeo import gdal, gdalconst
import multiprocessing
from rasterio.crs import CRS
import rasterio as rst
from rasterio import mask as msk
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import fiona
from Prep_Data import Norma_Data as nd
import pandas as pd
import numpy as np
#import seaborn as sb
import decimal
from multiprocessing import Pool

# Fonction


def cropRaster(inShape, inRaster, outRaster=None):
    """Crop/Clip a raster with a given shapefile.

    Args:
        inShape (_str_): pathname of the shapefile (polygon).
        inRaster (_str_): pathname of the raster to crop.
        outRaster (_str_): pathname of the output cropped/clipped raster.

    Returns:
        out_img _array_: array of the cropped/clipped raster.
    """
    # Open shapefile and get geometry feture
    if inShape.endswith(".shp"):  # if it is a shapefile filename
        with fiona.open(inShape, "r") as shapefile:
            geoms = [feature['geometry'] for feature in shapefile]  # find geometry of the shapefil

    elif type(inShape) == dict:
        geoms = inShape['geometry']

    # Crop image
    with rst.open(inRaster) as src:
        out_img, out_transform = msk.mask(src, geoms, invert=False, crop=True)

        if out_img.shape[0] < 10:
            out_img = out_img.reshape(out_img.shape[1], out_img.shape[2])
        out_meta = src.meta.copy()
        img_dtype = out_img.dtype
        out_meta.update({"driver": "GTiff",
                         "dtype": img_dtype,
                         "nodata": 0,
                         "height": out_img.shape[0],
                         "width": out_img.shape[1],
                         "transform": out_transform})

    print('Le Raster est bien Crop')
    return out_img, out_meta


def seuil_index(seuil_inf, seuil_sup, matrice):
    matrice_seuil = np.where(matrice < seuil_inf, seuil_inf, matrice)
    matrice_seuil = np.where(matrice_seuil > seuil_sup, seuil_sup, matrice_seuil)

    return matrice_seuil


def seuil_index_norma(matrice):
    matrice_seuil = np.where(matrice, -1 < -1, matrice)
    matrice_seuil = np.where(matrice_seuil > 1, 1, matrice_seuil)

    return matrice_seuil


def findfile(name_img, path):
    for dirpath, dirname, filename in os.walk(path):
        if name_img in filename:
            return dirpath


def ni4Find(path_dezip):
    list_pre_dos = [e for e in os.listdir(path_dezip) if e[-5:] == ".SAFE"]

    list_dos = [os.path.join(path_dezip, e) for e in os.listdir(path_dezip) if e[-5:] == ".SAFE"]

    list_nf = [str(e.split('_')[-2] + "_" + e.split('_')[2] + "_B02.jp2") for e in list_pre_dos]

    list_res = [findfile(e, i) for e, i in zip(list_nf, list_dos)]

    return list_res


def dbl_sigmoid_function(t, EVI_w, EVI_m, mS, S, mA, A):
    sigma1 = 1. / (1 + np.exp(mA * (t - A)))
    sigma2 = 1. / (1 + np.exp(-mS * (t - S)))
    return EVI_w + (EVI_m - EVI_w) * (sigma1 + sigma2 - 1)


def img_to_csv(filepath, name_csv, filter=False):
    img = nd(filepath)
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
            #print(dic)
            table = pd.DataFrame(data=dic)

        np.savez(name_csv, **dic)

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


def date_Data(path, select_format):
    l = [e for e in os.listdir(path) if e.endswith(select_format)]

    l_tuile_date = []

    for file in l:
        date_str = file.split('_')[1][:8]  # Mise en forme de la date

        tile_date = datetime.datetime(year=int(date_str[0:4]), month=int(date_str[4:6]),
                                      day=int(date_str[6:8]))
        l_tuile_date.append(tile_date)

    return l_tuile_date


def name_extract(path, select_format):
    l = [e for e in os.listdir(path) if e.endswith(select_format)]

    return l


def ravel_img(dos_img, select_format):

    l = [os.path.join(dos_img, e) for e in os.listdir(dos_img) if e.endswith(select_format)]

    l_rav = {}

    for e in l:
        img = nd(e)

        rav = np.ravel(img.band_Array)

        rav = np.where(rav <= 0, np.nan, rav)

        rav = np.where(rav > 1, 1, rav)

        rav = np.where(rav == np.nan, -255, rav)

        name = e.split('/')[-1]

        date_str = name.split('_')[1]

        tile_date = datetime.datetime(year=int(date_str[0:4]), month=int(date_str[4:6]), day=int(date_str[6:8]))

        DOY = tile_date.timetuple().tm_yday

        l_rav["Year"] = int(date_str[0:4])

        l_rav["DOY_" + "0" + str(DOY)] = rav

    return l_rav


def savgol_in_dict(dictio):

    df = pd.DataFrame(dictio)
    dic = dictio.copy()

    # df_no_NaN = df.dropna(axis=0)

    for i in dictio.keys():
        if i[0:3] == 'DOY':
            array_savgol = savgol_filter(df[i], window_length=5, polyorder=1)
            dic['Savgol' + '_' + i] = array_savgol

    df_good = pd.DataFrame(dic)
    df_good['Index'] = df_good.index

    return df_good


def dbl_sigmoid_function(t, EVI_w, EVI_m, mS, S, mA, A):
    sigma1 = 1. / (1 + np.exp(mA * (t - A)))
    sigma2 = 1. / (1 + np.exp(-mS * (t - S)))
    return EVI_w + (EVI_m - EVI_w) * (sigma1 + sigma2 - 1)


def organisation_temporal_dic(temporal_dict):
    l_x = [x for x in temporal_dict.keys() if x[:2] == 'Sa']

    l_x_DOY = [int((x[-3:])) for x in l_x]

    dic = {}

    decimal.getcontext().prec = 100

    for e in (temporal_dict.index):
        l_value = list(temporal_dict.loc[e])

        dic[int(e)] = [l_x_DOY, l_value]

    return dic


def First_deriative(temporal_dict):

    dic = organisation_temporal_dic(temporal_dict)
    print(f"Le Len du dic.items() en entré est {len(dic.items())}")
    t = np.arange(0, 365)

    dic_Final = {}
    np.seterr('ignore')

    try:
        for e in dic.items():
            pix = e[0]
            #print(pix)
            l_x = e[1][0]
            l_y = np.asarray((e[1][1]))
            l_y = l_y * 10000
            min_y = l_y.min()
            max_y = l_y.max()
            popt, pcov = curve_fit(dbl_sigmoid_function, l_x, l_y, p0=[min_y, max_y, 0.029, 2, 0.062, 125], maxfev=200000000)
            dbl_logistic = dbl_sigmoid_function(t, *popt)
            NDVI_w, NDVI_max, mS, S, mA, A = popt
            # print('Parameters : \n- NDVI_w', NDVI_w, '\n- NDVI_max', NDVI_max, '\n- mS : ', mS, '\n- S : ', S, '\n- mA : ', mA, '\n- A :', A)
            # curve_fit_evaluation = pd.DataFrame(columns=['i', 'j', 'Score'])
            try:
                # Calcul a
                produit_a = (-mS * (t - S))
                produit_a = np.asarray(produit_a, dtype=np.float128)
                exp_produit_a = np.exp(produit_a)
                exp_produit_a = mS * exp_produit_a
                diviseur_a = (-mS * (t - S))
                diviseur_a = np.asarray(diviseur_a, dtype=np.float128)
                diviseur_a = np.exp(diviseur_a)
                diviseur_a = 1 + diviseur_a
                diviseur_a = np.square(diviseur_a)
                a = exp_produit_a / diviseur_a
                # print(a)

                # Calcul b
                produit_b = (mA * (t - A))
                produit_b = np.asarray(produit_b, dtype=np.float128)
                exp_produit_b = np.exp(produit_b)
                exp_produit_b = -(mA * exp_produit_b)
                diviseur_b = (mA * (t - A))
                diviseur_b = np.exp(diviseur_b)
                diviseur_b = 1 + diviseur_b
                diviseur_b = np.square(diviseur_b)
                b = exp_produit_b / diviseur_b

            except RuntimeWarning:
                pass

            np.nan_to_num(a, nan=5000)
            np.nan_to_num(b, nan=5000)
            np.nan_to_num(NDVI_max, nan=5000)
            np.nan_to_num(NDVI_w, nan=5000)
            # a = (mS * np.exp(-mS * (t - S))) / np.square(1 + np.exp(-mS * (t - S)))
            # b = -(mA * np.exp(mA * (t - A))) / np.square(1 + np.exp(mA * (t - A)))
            first_df = (NDVI_max - NDVI_w) * (a + b)
            first_df = np.nan_to_num(first_df, nan=-255)

            try:
                EOS_DOY = np.where(first_df == np.nanmin(first_df))[0]

                SOS_DOY = np.where(first_df == np.nanmax(first_df))[0]
            # print(f"Le Start Of Season est : {(SOS_DOY)}")
            # print(f"Le End Of Season est : {(EOS_DOY)}")

            except RuntimeWarning:
                pass
            dic_Final[pix] = [SOS_DOY, EOS_DOY]


    except RuntimeError:
        pass
    except ValueError:
        pass
    except TypeError:
        pass
    return dic_Final


def prepa_Dataframe_multiprocess (dic_temporal, nb_process, name_file_SOS, name_file_EOS, create_raster=True):


    if create_raster is True:
        df = dic_temporal  #pd.read_json(dic_temporal)

    # Prochaine ligne a enlever car juste pour test
    #=============== Enlever ===============

    #df = df.transpose()

    #=======================================

        nb_row = np.shape(df)[0]
        print(nb_row)

        ratio_int = int(nb_row/nb_process)

        l_ratio = [e for e in range(nb_row) if e % ratio_int == 0]

        l_decoup = [e for e in l_ratio]

        l_decoup.append(nb_row)

        l_table = []

        for e in range(len(l_decoup)):

            try:
                #print(f"Couple de découpage : {l_decoup[e]}, {l_decoup[e+1]}")
                table = df[l_decoup[e]:l_decoup[e+1]]
                l_table.append(table)

            except IndexError:
                pass

        pool = Pool(processes=nb_process)

        data_multi = pool.map(First_deriative, l_table)

        agencement_data_multiproces(data_multi, nb_process, name_file_SOS, name_file_EOS)

    else:

        df = dic_temporal  # pd.read_json(dic_temporal)

        # Prochaine ligne a enlever car juste pour test
        # =============== Enlever ===============

        # df = df.transpose()

        # =======================================

        nb_row = np.shape(df)[0]
        print(nb_row)

        ratio_int = int(nb_row / nb_process)

        l_ratio = [e for e in range(nb_row) if e % ratio_int == 0]

        l_decoup = [e for e in l_ratio]

        l_decoup.append(nb_row)

        l_table = []

        for e in range(len(l_decoup)):

            try:
                # print(f"Couple de découpage : {l_decoup[e]}, {l_decoup[e+1]}")
                table = df[l_decoup[e]:l_decoup[e + 1]]
                l_table.append(table)

            except IndexError:
                pass

        pool = Pool(processes=nb_process)

        data_multi = pool.map(First_deriative, l_table)

        return data_multi


def agencement_data_multiproces (dataframePandas, nb_process, name_file_SOS, name_file_EOS):

    ###### Reagencement des données multiprocessing

    #data_parall = pd.read_json(Name_sauvegarde_json)
    #data_parall = dataframePandas.transpose()
    # print(data_parall.columns)

    dic_iter = {}

    for e in range(nb_process):
        dic_iter[e] = f"Iter{e}"

    data_parall = dataframePandas.transpose()
    data_parall = data_parall.rename(columns=dic_iter)

    #print(data_parall)


    data_parall['RES'] = data_parall[data_parall.columns[:]].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1)
    data_parall['Index'] = data_parall.index
    data_res = {"Res": data_parall['RES'], 'Index': data_parall['Index']}


    # Problème l'élément ou la liste qui contient les DOY est mis en str

    # np.savez('Data_Test_parall_BERAMBADI.npz', **data_res)

    # print((data_res['Res']))

    #data = np.load('Data_Test_parall_BERAMBADI.npz', allow_pickle=True)

    data = data_res['Res']  # ATTENTION

    #print(data)

    index = data_res['Index'] # ATTENTION

    pattern_SOS = re.compile(r'(\d{1,3})')

    l_SOS = [int(re.findall(pattern_SOS, line)[0]) for line in data]
    l_EOS = [int(re.findall(pattern_SOS, line)[1]) for line in data]

    dict_good = {"SOS": l_SOS,
           "EOS": l_EOS,
           "Index": index}

    data_fin = pd.DataFrame(dict_good)

    #print(data_fin)

    rav = np.ravel(img_ref.band_Array[:, :])
    dicto = {'Value_img': rav}

    data_tot = pd.DataFrame(dicto)
    data_tot['Index'] = data_tot.index

    #print(data_fin)

    p = data_tot.merge(data_fin, on='Index', how='left')

    SOS = np.asarray(list(p['SOS'])).reshape(img_ref.shape[0], img_ref.shape[1])
    EOS = np.asarray(list(p['EOS'])).reshape(img_ref.shape[0], img_ref.shape[1])

    nd.Write_ras(SOS, name_file_SOS, img_ref.geodata,
                 img_ref.projection)
    nd.Write_ras(EOS, name_file_EOS, img_ref.geodata,
                 img_ref.projection)



if __name__ == '__main__':

    debut = time.time()

    ###### Chemin d'accès

    rep_NDVI_path = "/home/adr2.local/pellen_j/Projet_Laura_Samuel/NDVI/"
    Name_sauvegarde_json = 'Data_Test_parall_BERAMBADI.json'
    img_ref = nd('/home/adr2.local/pellen_j/Projet_Laura_Samuel/NDVI/T43PFP_20210121T052119_Stack_NDVI.tif')
    name_result_SOS = "SOS_Berambadi_2021.tif"
    name_result_EOS = "EOS_Berambadi_2021.tif"

    ###### Processus

    filepath = name_extract(rep_NDVI_path, '.tif')

    du = ravel_img(rep_NDVI_path, select_format='.tif')

    df = savgol_in_dict(du)

    l_keys_dict = df.keys()
    l_select_Sav = [e for e in l_keys_dict if e[:2] == "Sa"]

    df_Savgol = df[l_select_Sav]

    df_Savgol = df_Savgol.dropna(axis=0)
    
    # Partie test sur multiprocessing  
    
    df_final = prepa_Dataframe_multiprocess(dic_temporal=df_Savgol, nb_process=10, name_file_SOS=name_result_SOS, name_file_EOS=name_result_EOS, create_raster=True)




    """
    agencement_data_multiproces(df,10, "SOS_TEST_fonction.tif", "EOS_TEST_fonction.tif")
    df_final = pd.DataFrame(df_final)

    df_final = df_final.transpose()
    df_final = df_final.rename(columns={0: 'SOS', 1: 'EOS'})
    #print(df_final)
    df_final['Index'] = df_final.index
    df_final['Index'] = df_final['Index'].astype(int)
    #print(df_final.dtypes)
    SOS = []
    EOS = []
    Index_li = []

    for e in df_final.values:
        SOS_value = int(e[0][0])
        # print(type(SOS_value),SOS_value)

        EOS_value = int(e[1][0])
        Index = e[2]
        SOS.append(SOS_value)
        EOS.append(EOS_value)
        Index_li.append(Index)

    df_final['SOS'] = SOS
    df_final['EOS'] = EOS

    p = df.merge(df_final, on='Index', how='left')
    print(p)

    SOS_good = p['SOS']
    EOS_good = p['EOS']



    np_SOS = np.asarray(SOS_good).reshape(img_ref.shape[0], img_ref.shape[1])
    np_EOS = np.asarray(EOS_good).reshape(img_ref.shape[0], img_ref.shape[1])

    nd.Write_ras(np_SOS,name_result_SOS, img_ref.geodata,img_ref.projection)
    nd.Write_ras(np_EOS, name_result_EOS, img_ref.geodata, img_ref.projection)

    #plt.imshow(np_EOS, cmap='Greens')
    #plt.show()
    

    ###### Reagencement des données multiprocessing


    data_parall = pd.read_json(Name_sauvegarde_json)
    data_parall = data_parall.transpose()
    #print(data_parall.columns)

    data_parall = data_parall.rename(columns={0: 'Iter0', 1: 'Iter1', 2: 'Iter2', 3: 'Iter3', 4: 'Iter4',
                                              5: 'Iter5', 6: 'Iter6', 7: 'Iter7', 8: 'Iter8', 9: 'Iter9'})

    data_parall['RES'] = data_parall[data_parall.columns[:]].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1)
    data_parall['Index'] = data_parall.index
    data_res = {"Res": data_parall['RES'], 'Index': data_parall['Index']}

    # Problème l'élément ou la liste qui contient les DOY est mis en str



    #np.savez('Data_Test_parall_BERAMBADI.npz', **data_res)

    #print((data_res['Res']))

    data = np.load('Data_Test_parall_BERAMBADI.npz', allow_pickle=True)




    ds = data['Res']
    index = data['Index']

    pattern_SOS = re.compile(r'(\d{1,3})')


    l_SOS = [int(re.findall(pattern_SOS, line)[0]) for line in ds]
    l_EOS = [int(re.findall(pattern_SOS, line)[1]) for line in ds]

    dic = {"SOS":l_SOS,
         "EOS": l_EOS,
         "Index": index}


    data_fin = pd.DataFrame(dic)


    rav = np.ravel(img_ref.band_Array[:,:])
    dicto = {'Value_img': rav}

    data_tot = pd.DataFrame(dicto)
    data_tot['Index'] = data_tot.index

    p = data_tot.merge(data_fin, on='Index', how='left')

    SOS = np.asarray(list(p['SOS'])).reshape(img_ref.shape[0],img_ref.shape[1])
    EOS = np.asarray(list(p['EOS'])).reshape(img_ref.shape[0], img_ref.shape[1])

    nd.Write_ras(SOS, "/home/adr2.local/pellen_j/Projet_Laura_Samuel/NDVI/Test_paral_SOS.tif",img_ref.geodata, img_ref.projection)
    nd.Write_ras(EOS, "/home/adr2.local/pellen_j/Projet_Laura_Samuel/NDVI/Test_paral_EOS.tif", img_ref.geodata, img_ref.projection)




    """


    #p = table_fin.merge(df_final, on='Index', how='left')
    #print(table_fin)








    #loc_start =
    #
    # int(re.findall(pattern_SOS, line)[0])




    #l_SOS = [int(ds[e][2:5]) for e in range(len(ds))]
    #l_EOS = [int(ds[e][9:12]) for e in range(len(ds))]




    #l_EOS = [int(ds[e][])]
    #data_dic = data_parall.to_dict()

    #printds

    #np.savez('Data_Test_parall_BERAMBADI.npz', **data_dic)

    """
    Exemple :
    
    df['ColumnA'] = df[df.columns[1:]].apply(lambda x: ','.join(x.dropna().astype(str)),axis=1)

    data_parall['RES'] = data_parall.iloc[:].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1)
    #data_parall['Index'] = data_parall.index
    #data_parall['Res'] = data_parall[0]+data_parall[1]+data_parall[2]+data_parall[3]+data_parall[4]+data_parall[5]+data_parall[6]+data_parall[7]+data_parall[8]+data_parall[9]
    print(data_parall['Res'])
    """



    fin = time.time()
    print("Temps de Traitement : ", (fin - debut) / 60, "minutes")











