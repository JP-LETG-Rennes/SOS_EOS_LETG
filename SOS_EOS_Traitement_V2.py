import datetime
import os
import re
import time
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from Prep_Data import Norma_Data as nd
from Prep_Data import Norma_Data as nd
import pandas as pd
import numpy as np
import decimal
from multiprocessing import Pool


# Fonction

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

    pattern_date = re.compile(r'(\d{7,8})')

    l_rav = {}

    for e in l:
        img = nd(e)

        rav = np.ravel(img.band_Array)

        rav = np.where(rav <= 0, np.nan, rav)

        rav = np.where(rav > 1, 1, rav)

        rav = np.where(rav == np.nan, -255, rav)

        name = e.split('/')[-1]

        date_str = re.findall(pattern_date, name)[0]

        #date_str = name.split('_')[1]

        tile_date = datetime.datetime(year=int(date_str[0:4]), month=int(date_str[4:6]), day=int(date_str[6:8]))

        DOY = tile_date.timetuple().tm_yday

        l_rav["Year"] = int(date_str[0:4])

        if len(str(DOY)) == 3:
            l_rav["DOY_" + str(DOY)] = rav

        if len(str(DOY)) == 2:
            l_rav["DOY_" + "0" + str(DOY)] = rav

        if len(str(DOY)) == 1:
            l_rav["DOY_" + "00" + str(DOY)] = rav

    return l_rav

def savgol_in_dict(dictio):

    df = pd.DataFrame(dictio)

    dic = dictio.copy()

    # df_no_NaN = df.dropna(axis=0)

    for i in dictio.keys():
        if i[0:3] == 'DOY':
            array_savgol = savgol_filter(df[i], window_length=5, polyorder=1)
            dic['Savgol' + '_' + i] = array_savgol
            #print("C'est le I :", i)

    df_good = pd.DataFrame(dic)
    df_good['Index'] = df_good.index

    return df_good

def organisation_temporal_dic(temporal_dict):

    l_x = [x for x in temporal_dict.keys() if x[:2] == 'Sa']

    l_x_DOY = [int((x[-3:])) for x in l_x] # le x[-3:] génére une erreur lorsque la date est entre le 01 et le 09

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
            popt, pcov = curve_fit(dbl_sigmoid_function, l_x, l_y, p0=[min_y, max_y, 0.029, 2, 0.062, 125],
                                   maxfev=200000000)
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

            #np.nan_to_num(a, nan=5000)
            #np.nan_to_num(b, nan=5000)
            #np.nan_to_num(NDVI_max, nan=5000)
            #np.nan_to_num(NDVI_w, nan=5000)
            # a = (mS * np.exp(-mS * (t - S))) / np.square(1 + np.exp(-mS * (t - S)))
            # b = -(mA * np.exp(mA * (t - A))) / np.square(1 + np.exp(mA * (t - A)))
            first_df = (NDVI_max - NDVI_w) * (a + b)
            #first_df = np.nan_to_num(first_df, nan=-255)

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

def prepa_Dataframe_multiprocess (dic_temporal, nb_process, name_file_SOS, name_file_EOS, img_ref,create_raster=True):


    if create_raster is True:

        df = dic_temporal

        nb_row = np.shape(df)[0]

        ratio_int = int(nb_row/nb_process)

        l_ratio = [e for e in range(nb_row) if e % ratio_int == 0]

        l_decoup = [e for e in l_ratio]

        l_decoup.append(nb_row)

        l_table = []

        for e in range(len(l_decoup)):

            try:
                table = df[l_decoup[e]:l_decoup[e+1]]
                l_table.append(table)

            except IndexError:
                pass

        pool = Pool(processes=nb_process)

        data_multi = pool.map(First_deriative, l_table)

        data_multi = pd.DataFrame(data_multi)

        agencement_data_multiproces(data_multi, nb_process, name_file_SOS, name_file_EOS, img_ref)

    else:

        df = dic_temporal  # pd.read_json(dic_temporal)

        nb_row = np.shape(df)[0]

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

def agencement_data_multiproces (dataframePandas, nb_process, name_file_SOS, name_file_EOS, img_ref):

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

    img_ref = nd(img_ref)

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

def lecture_config(configuration_path):

    # Récupération des paramètres contenus dans le fichier conf.json
    with open(configuration_path, "r") as fileConf:
        conf = json.load(fileConf)

    rep_NDVI_path = conf["rep_NDVI_path"]
    img_ref = conf["img_ref"]
    name_result_SOS = conf["name_result_SOS"]
    name_result_EOS = conf["name_result_EOS"]
    nb_process = conf["nb_process"]

    return rep_NDVI_path, img_ref, name_result_SOS, name_result_EOS, nb_process


if __name__ == '__main__':

    debut = time.time()

    ############## Chemin d'accès ##############
    rep_NDVI_path, img_ref, name_result_SOS, name_result_EOS, nb_process = lecture_config(config_SOS_EOS)

    ############## Processus ##############

    filepath = name_extract(rep_NDVI_path, '.tif')
 
    du = ravel_img(rep_NDVI_path, select_format='.tif')

    df = savgol_in_dict(du)

    l_keys_dict = df.keys()

    l_select_Sav = [e for e in l_keys_dict if e[:2] == "Sa"]

    df_Savgol = df[l_select_Sav]

    df_Savgol = df_Savgol.dropna(axis=0)
    
    # Partie sur multiprocessing
    
    df_final = prepa_Dataframe_multiprocess(dic_temporal=df_Savgol, nb_process=nb_process, name_file_SOS=name_result_SOS, name_file_EOS=name_result_EOS, img_ref=img_ref, create_raster=True)

    fin = time.time()

    print("Temps de Traitement : ", (fin - debut) / 60, "minutes")












