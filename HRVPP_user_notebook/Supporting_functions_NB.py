#general packages
import subprocess
import requests
import time
import os
import logging
import concurrent.futures
from tqdm.auto import tqdm
from datetime import datetime, date


#raster/s3 packages
import boto3
from rasterio.session import AWSSession
from rasterio.windows import Window
from rasterio.features import geometry_window, geometry_mask
import rasterio

#Geometry packages
import shapely
from shapely.geometry import Point, shape
from pyproj import Transformer, CRS
import shapefile
from affine import Affine


#matrix packages
import numpy as np
import pandas as pd
import geopandas as gpd

#plotting packages
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from ipyleaflet import *
from pandas.plotting import register_matplotlib_converters



###########################
# CLASSES
##########################
          
class data_extraction:
    def __init__(self, header_file, credentials, config, catalogue, version, tile, xoff, yoff, N_pixel_window, shape, years, fill_value=0):
        self.header_file = header_file
        self.credentials = credentials
        self.config = config
        self.catalogue = catalogue
        self.version = version
        self.tile = tile
        self.xoff = xoff
        self.yoff = yoff
        self.years = years
        self.N_pixel_window = N_pixel_window
        if N_pixel_window ==0 :
            self.shape = shape
            
        
        self.rio_gdal_options = {
                    'GDAL_HTTP_HEADER_FILE' : self.header_file,
                    'GDAL_DISABLE_READDIR_ON_OPEN': 'EMPTY_DIR',
                    'GDAL_HTTP_MERGE_CONSECUTIVE_RANGES': 'YES',
                    'CPL_VSIL_CURL_ALLOWED_EXTENSIONS': 'tif',
                    'VSI_CACHE': 'FALSE',
                    'VSI_CACHE_SIZE': 0,
                    'GDAL_CACHEMAX': 0}

        self.loader = HTTPParallelLoader(rio_gdal_options=self.rio_gdal_options,
                     max_workers=15, fill_value=fill_value)
        
        
        
        
    def get_token(self):
        command = [f'curl -s -d "client_id={self.config.oidc_client_id}"']
        command.append('-d "client_secret=%s"'%"uC2pdaxozGkI8iJh7TqfJKzoKfIa")
        command.append('-d "username=%s"'%self.credentials[0])
        command.append('-d "password=%s"'%self.credentials[1])
        command.append('-d "grant_type=password"')
        command.append(f'"{self.config.oidc_token_endpoint}"')
        out = eval(subprocess.check_output(" ".join(command), shell=True))
        if 'error' in out:
           logging.error("Following error occured when getting the token: {}".format(out))
        return out['access_token']
    
    
    def savetoken(self):
        with open(self.header_file,'w') as header:
            header.write(f"Authorization: Bearer {self.get_token()}") 

    

class VI_class(data_extraction):
    def __init__(self, header_file, credentials ,config, catalogue, version, tile, xoff, yoff, N_pixel_window, shape, option_PPI, option_NDVI, years, fill_value=-32768):
        self.option_PPI = option_PPI
        self.option_NDVI = option_NDVI
        self.dict_VI_TS = {}
        self.dict_no_data_VI = {}
        self.dict_no_data_QF2 = {}
        self.dict_VI_df = {}

        super().__init__(header_file,credentials,config,catalogue, version,  tile, xoff, yoff, N_pixel_window, shape, years, fill_value=fill_value)
        
        
    ## FUNCTION TO BUILD_UP VI FILELIST FROM WEKEO:
    def build_vi_list(self, VI_name='PPI'):
        vi_list = list(self.catalogue.get_products(f"copernicus_r_utm-wgs84_10_m_hrvpp-vi_p_2017-now_{self.version}", 
                                       tileId=self.tile,
                                       productType=VI_name,
                                       start=date(int(self.years[0]), 1, 1),
                                       end=date(int(self.years[-1]), 12, 31)))

        if len(vi_list) == 0:
            log.error("No data found for the years given, please check your years to be extracted.")
        df_vi = pd.DataFrame(vi_list, columns=['DL'])
        df_vi['path']= df_vi.applymap(lambda x: x.data[0].href)
        df_temp = df_vi["path"].str.split("_", n=3, expand=True)
        df_temp.columns = ['version','vi', 'date', 'name']
        df_vi['year'] = df_temp['date'].str[0:4]
        df_vi['month'] = df_temp['date'].str[4:6]
        df_vi['day'] = df_temp['date'].str[6:8]
        df_vi['time'] =  pd.to_datetime(df_vi[['year', 'month', 'day']])
        df_temp = None

        return df_vi
    
        ## FUNCTION TO EXTRACT THE DATA FROM WEKEO:
    def extract_vi_data(self, df_vi, VI_name):
        #Change for easier typeing
        N_pixel_window = self.N_pixel_window
        xoff = self.xoff
        yoff = self.yoff

        self.savetoken()  
        print('EXTRACTING raw {} timeseries (this can take a while ... several minutes)'.format(VI_name))

        start = time.time()

        product_locations = df_vi.path
        qf_locations = [file.replace('_{}'.format(VI_name), '_QFLAG2') for file in product_locations]
        date_times = df_vi.time
        with  rasterio.Env(**self.rio_gdal_options) as env:
                        file= product_locations[1]
                        with rasterio.open(file, 'r') as ds:
                            no_data_VI = ds.nodata
                            scale = ds.scales
                            if self.N_pixel_window == 0:
                                geom, transform = get_reproj_geom_from_shape_file(self.shape, self.tile, ds)
                                window = geometry_window(ds,[geom])
                                polygon_mask = geometry_mask(geometries=[geom],
                                                       out_shape=(window.height, 
                                                                  window.width,
                                                                 ),
                                                       transform=transform * Affine.translation(window.col_off, window.row_off),
                                                       all_touched=False,
                                                       invert=False)
                                N_pixel_window = (window.width + window.height) / 2
                            else: 
                                N_pixel_window = self.N_pixel_window 
                                window = Window(int(xoff / 10 - N_pixel_window / 2), int(yoff / 10 - N_pixel_window / 2), N_pixel_window, N_pixel_window)
                                
                    
                            
                            
        data = np.zeros((len(product_locations),window.height,window.width))
        qf_array = np.zeros((len(product_locations),window.height,window.width))
        max_length = 300

        for chunk in range(0,len(product_locations)//max_length + 1):
            intermed = time.time()
            self.savetoken() 
            data[chunk*max_length: min((chunk+1)*max_length,len(product_locations)),:,:] = np.asarray(self.loader.load_arrays(product_locations[chunk*max_length: min((chunk+1)*max_length,len(product_locations))],window)).astype(np.float).reshape(np.shape(data[chunk*max_length:min((chunk+1)*max_length,len(product_locations)),:,:]))
            self.savetoken() 
            qf_array[chunk*max_length: min((chunk+1)*max_length,len(product_locations)),:,:] = np.asarray(self.loader.load_arrays(qf_locations[chunk*max_length: min((chunk+1)*max_length,len(product_locations))], window)).reshape(np.shape(data[chunk*max_length: min((chunk+1)*max_length,len(product_locations)),:,:]))
            print('EXTRACTION VI CHUNK {}/{} FINISHED in {} mins\n'.format(chunk+1,len(product_locations)//max_length+1,(time.time()-intermed)/60))
        
        data[data!=no_data_VI]=data[data!=no_data_VI]*scale
        no_data_QF2 = 65535

        if N_pixel_window > 1:
            loc_cloudiness = np.where(qf_array != 1, 1,0)
            if self.N_pixel_window == 0:
                data = np.ma.masked_where( loc_cloudiness | (data==no_data_VI) | np.tile(polygon_mask, (np.shape(data)[0],1,1)),data)
                pixel_area = np.sum(polygon_mask)
            else:
                data = np.ma.masked_where(loc_cloudiness | (data==no_data_VI),data)
                pixel_area = (N_pixel_window^2)
            
            data2 = np.ma.mean(data, axis= (1,2))
            data2 = np.ma.filled(data2,no_data_VI)
            cloudiness = (np.sum(loc_cloudiness, axis=(1,2)) / pixel_area) > (0.25)          
            data2[cloudiness] = no_data_VI
            clouds = data2.copy()
            clouds[cloudiness] = no_data_QF2
            clouds[~cloudiness] = 1
            VI_TS = [[date_times[i], data2[i], clouds[i]] for i in range(0, len(date_times))]


        else:
            VI_TS = [[date_times[i], data[i], qf_array[i]] for i in range(0, len(date_times))]

        
        
        return VI_TS, no_data_VI, no_data_QF2


    
        ## WRAPPER FUNCTION TO EXTRACT THE VI DATA AND STORE IT IN A DICTIONARY FOR EACH VI OF INTEREST:
    def load_and_store_vi_data(self):
        N_pixel_window = self.N_pixel_window
        xoff = self.xoff
        yoff = self.yoff

        print('** Build VI file list **')
        if self.option_NDVI and self.option_PPI:
            VIs_extract = ['PPI', 'NDVI']
        elif self.option_NDVI and not self.option_PPI:
            VIs_extract = ['NDVI']
        else:
            VIs_extract = ['PPI']

        for VI in VIs_extract:       
            df_vi = self.build_vi_list(VI_name=VI)
            VI_TS, no_data_VI, no_data_QF2 = self.extract_vi_data(df_vi, VI)
            self.dict_VI_TS.update({VI: VI_TS})
            self.dict_no_data_VI.update({VI: no_data_VI})
            self.dict_no_data_QF2.update({VI: no_data_QF2})

    ## FUNCTION TO WRITE OUT THE EXTRACTED VI(s) INTO A DATAFRAME:

    def VI_to_df(self):
        for key in list(self.dict_VI_TS.keys()):
            VI_TS = self.dict_VI_TS.get(key)
            no_data_VI = self.dict_no_data_VI.get(key)
            no_data_QF2 = self.dict_no_data_QF2.get(key)
            df_VI = pd.DataFrame(VI_TS, dtype=float)
            df_VI.columns = ['date', '{}'.format(key), 'QF2']
            # determine weights used for curve fitting
            df_VI['weight'] = 0
            df_VI['QF2'] = df_VI['QF2'].astype(np.uint16)  # V091 correct weight calculation
            df_VI.loc[(df_VI['QF2'] & 1 == 1), 'weight'] = 2
            df_VI.loc[(((df_VI['QF2'] & 1 == 1) & (df_VI['QF2'] & 1024 == 1024)) | (
                        (df_VI['QF2'] & 1 == 1) & (df_VI['QF2'] & 2048 == 2048))), 'weight'] = 1
            df_VI.loc[(df_VI['QF2'] == no_data_QF2), 'weight'] = 0
            df_VI[df_VI['{}'.format(key)] == no_data_VI] = np.NaN
            df_VI.set_index('date', inplace=True)
            df_VI = df_VI.dropna().sort_index()
            self.dict_VI_df.update({key: df_VI})
            
    def plot_VI(self):
        for VI in list(self.dict_VI_df.keys()):
                df_VI = self.dict_VI_df.get(VI)
                if VI == 'PPI':
                    colors = {0: 'gainsboro', 1: 'darkslategray', 2: 'black'}
                else:
                    colors = {0: 'palegreen', 1: 'limegreen', 2: 'darkgreen'}

                df_VI['color'] = df_VI['weight'].apply(lambda x: colors[x])
                for g, b in df_VI.groupby(by='color'):
                    if g == colors[2]:
                        fill = 'full'
                        w = '-100%'
                    elif g == colors[1]:
                        fill = 'bottom'
                        w = '-50%'
                    else:
                        fill = 'none'
                        w = '-0%'
                    if VI == 'PPI':
                        ax1.plot(b.index, b['{}'.format(VI)], color=g, linestyle='', fillstyle=fill,
                                 label='{}'.format('raw_{}'.format(VI) + w), marker='o')
                    elif len(list(dict_VI_df.keys())) > 1 and VI == 'NDVI':
                        ax2.plot(b.index, b['{}'.format(VI)], color=g, linestyle='', fillstyle=fill,
                                 label='{}'.format('raw_{}'.format(VI) + w), marker='o')
                    else:
                        ax1.plot(b.index, b['{}'.format(VI)], color=g, linestyle='', fillstyle=fill,
                                 label='{}'.format('raw_{}'.format(VI) + w), marker='o')
        return


class phenology_class(data_extraction):
    def __init__(self,header_file, credentials, config, catalogue, version, tile, xoff, yoff, N_pixel_window, shape, years):        
        self.phen = []
        self.df_VPP =  pd.DataFrame()
        self.ST_TS = []  
        self.df_ST = pd.DataFrame()

        super().__init__(header_file,credentials,config,catalogue, version, tile, xoff, yoff, N_pixel_window, shape, years, fill_value=0)
        
    def extract_phenology_data(self):   
        #for easier typeing
        years = self.years
        version= self.version
        xoff = self.xoff
        yoff = self.yoff
        
        DD = ['01', '11', '21']
        MM = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
        VPP_params = ['SOSD', 'MAXD', 'EOSD']
        seasons = ['s1', 's2']
        

        print(f'EXTRACTING ST DATA FOR YEARS: {years[0]} - {years[-1]}')
        self.savetoken()
        start = time.time()
        products = list(self.catalogue.get_products(f"copernicus_r_utm-wgs84_10_m_hrvpp-st_p_2017-now_{version}", 
                               start=date(int(years[0]), 1, 1), 
                               end=date(int(years[-1]),12, 31), 
                               tileId=self.tile,                                                                                 
                               productType="PPI"))
        product_locations = [x.data[0].href for x in products]
        date_times = [x.beginningDateTime for x in products]



        with  rasterio.Env(**self.rio_gdal_options) as env:
                file= product_locations[0]
                with rasterio.open(file, 'r') as ds:
                    self.no_data_ST = ds.nodata
                    scale = ds.scales
                    if self.N_pixel_window == 0:
                        geom, transform = get_reproj_geom_from_shape_file(self.shape, self.tile, ds)
                        window = geometry_window(ds,[geom])
                        polygon_mask = geometry_mask(geometries=[geom],
                                               out_shape=(window.height, 
                                                          window.width,
                                                         ),
                                               transform=transform * Affine.translation(window.col_off, window.row_off),
                                               all_touched=False,
                                               invert=False)
                       
                        N_pixel_window = (window.width + window.height) / 2
                    else: 
                        N_pixel_window = self.N_pixel_window 
                        window = Window(int(xoff / 10 - N_pixel_window / 2), int(yoff / 10 - N_pixel_window / 2), N_pixel_window, N_pixel_window)
                    

        data = np.asarray(self.loader.load_arrays(product_locations , window)).astype(np.float)
        data[data!=self.no_data_ST] = data[data!=self.no_data_ST]*scale

        if N_pixel_window > 1:
            if self.N_pixel_window == 0:
                data = np.ma.masked_where( (data==self.no_data_ST) | np.tile(polygon_mask, (np.shape(data)[0],1,1)),data)
            else:
                data = np.ma.masked_where((data==self.no_data_ST),data)
            data2 = np.ma.mean(data, axis= (1,2))
            data = np.ma.filled(data2,self.no_data_ST)

        self.ST_TS.extend([[date_times[i], data[i]] for i in range(0, len(date_times))])


        print(f'EXTRACTING VPP DATA FOR YEARS: {years[0]} - {years[-1]}')        
        products = []            
        for VPP_param in VPP_params:         
            pre_products = list(self.catalogue.get_products(f"copernicus_r_utm-wgs84_10_m_hrvpp-vpp_p_2017-now_{self.version}", 
                   start=date(int(years[0]), 1, 1), 
                   end=date(int(years[-1]),12, 31), 
                   tileId=self.tile, 
                   productType=VPP_param))
            products.extend(pre_products)


        product_locations = [x.data[0].href for x in products]
        product_year = [x.beginningDateTime.year for x in products]
        product_season = [f"{x.properties['productInformation']['productType']}_{x.properties['productInformation']['productGroupId']}" for x in products ]
        data = np.asarray(self.loader.load_arrays(product_locations, window))



        with  rasterio.Env(**self.rio_gdal_options) as env:
            file= product_locations[0]
            with rasterio.open(file, 'r') as ds:
                self.no_data_VPP = ds.nodata

        if N_pixel_window > 1:
            if self.N_pixel_window == 0:
                data = np.ma.masked_where((data==self.no_data_VPP) | np.tile(polygon_mask, (np.shape(data)[0],1,1)),data)
            else:
                 data = np.ma.masked_where((data==self.no_data_VPP),data)
      
            VPP_data = (data%1000) + (data//1000)*365
            stats_site = np.floor(np.ma.mean(VPP_data, axis= (1,2))-1)                
            data2 = ((stats_site//365) *1000  +  stats_site%365 + 1)

            data = np.ma.filled(data2,self.no_data_VPP)

        self.phen.extend([[product_year[i], data[i],product_season[i]] for i in range(0, len(product_year))])      

        end = time.time()
        print('... in {} mins'.format((end-start)/60))
        print('EXTRACTION FINISHED')

    def VPP_to_df(self): 
        #for easier typeing
        years = self.years

        # The date axis on which all the values will be plotted
        idx_years = pd.date_range('{}-01-01'.format(years[0]), '{}-12-31'.format(years[-1]))

        # VPP extracts to dataframe
        df_VPP =  pd.DataFrame(self.phen, dtype=float)
        df_VPP.columns = ['year','transition_date', 'phenophase']
        self.df_VPP = df_VPP
        df_ST = pd.DataFrame(self.ST_TS, dtype = float)
        df_ST.columns = ['date','ST']

        # reindex dataframe
        df_ST.set_index('date',inplace  = True)
        self.df_ST = df_ST.reindex(idx_years)
    
        ## FUNCTION TO PRINT THE DATES OF THE PHENOLOGY ESTIMATES
    def print_VPP_params(self):
        #for easier typeing
        years = self.years
        df_VPP = self.df_VPP
        no_data_VPP = self.no_data_VPP
        
        if len(np.unique(df_VPP.transition_date))==1:
            logging.warning("no seasons have been found")
            return

        for year in df_VPP.year.unique():

            SOS_s1 = df_VPP[((df_VPP.year == year) & (df_VPP.phenophase == 'SOSD_s1'))]['transition_date'].values
            SOS_s2 = df_VPP[((df_VPP.year == year) & (df_VPP.phenophase == 'SOSD_s2'))]['transition_date'].values
            MAX_s1 = df_VPP[((df_VPP.year == year) & (df_VPP.phenophase == 'MAXD_s1'))]['transition_date'].values
            MAX_s2 = df_VPP[((df_VPP.year == year) & (df_VPP.phenophase == 'MAXD_s2'))]['transition_date'].values
            EOS_s1 = df_VPP[((df_VPP.year == year) & (df_VPP.phenophase == 'EOSD_s1'))]['transition_date'].values
            EOS_s2 = df_VPP[((df_VPP.year == year) & (df_VPP.phenophase == 'EOSD_s2'))]['transition_date'].values

            if SOS_s1 != int(no_data_VPP):
                print('\n\n sos DOY (first season) for {} at : '.format(str(int(year))), end='')
                print(datetime.strptime(str(int(SOS_s1[-1])), '%y%j'), end='')
            if MAX_s1 != int(no_data_VPP):
                print('\n max DOY (first season) for {} at : '.format(str(int(year))), end='')
                print(datetime.strptime(str(int(MAX_s1[-1])), '%y%j'), end='')
            if EOS_s1 != int(no_data_VPP):
                print('\n eos DOY (first season) for {} at: '.format(str(int(year))), end='')
                print(datetime.strptime(str(int(EOS_s1[-1])), '%y%j'), end='')

            if SOS_s2 != int(no_data_VPP):
                print('\n\n sos DOY (second season) for {} at : '.format(str(int(year))), end='')
                print(datetime.strptime(str(int(SOS_s2[-1])), '%y%j'), end='')
            if MAX_s2 != int(no_data_VPP):
                print('\n max DOY (second season) for {} at: '.format(str(int(year))), end='')
                print(datetime.strptime(str(int(MAX_s2[-1])), '%y%j'), end='')
            if EOS_s2 != int(no_data_VPP):
                print('\n eos DOY (second season) for {} at : '.format(str(int(year))), end='')
                print(datetime.strptime(str(int(EOS_s2[-1])), '%y%j'), end='')
                

                
        
###################################################################
#PARALLELIZATION CLASSES
###################################################################
class ParallelLoader:

    def __init__(self,
                 max_workers=20,
                 progressbar=False,
                 rio_gdal_options=None,
                 fill_value=0,
                 random_order=False):

        self._max_workers = max_workers
        self._progressbar = progressbar

        if rio_gdal_options is None:
            rio_gdal_options = {'GDAL_CACHEMAX': 40}

        self._rio_gdal_options = rio_gdal_options
        self._fill_value = fill_value
        self._random_order = random_order

    def _load_array_windows(self, fname, window):

        with rasterio.Env(**self._rio_gdal_options):
            with rasterio.open(fname) as src:

                arr = src.read(window=window, boundless=True,
                               fill_value=self._fill_value)

        arr = np.squeeze(arr)
        return arr



    def load_arrays(self, filenames, windows):

        def f_handle(filename):
            return self._load_array_windows(filename, windows)

        if self._random_order:
            ids = list(range(len(filenames)))
            ids_random = ids.copy()
            random.shuffle(ids_random)
            ids_map = [ids_random.index(i) for i in ids]
            filenames = [filenames[i] for i in ids_random]

        arrs_list = self._run_parallel(f_handle,
                                       filenames,
                                       self._max_workers,
                                       threads=True,
                                       progressbar=self._progressbar)

        if self._random_order:
            arrs_list = [arrs_list[i] for i in ids_map]

        return arrs_list

    @staticmethod
    def _run_parallel(f, my_iter, max_workers, threads=True, progressbar=True):

        if threads:
            Pool = concurrent.futures.ThreadPoolExecutor
        else:
            Pool = concurrent.futures.ProcessPoolExecutor

        with Pool(max_workers=max_workers) as executor:
            print(executor)
            if progressbar:
                results = list(executor.map(f, my_iter))
            else:
                results = list(executor.map(f, my_iter))

        return results



class HTTPParallelLoader(ParallelLoader):
    def __init__(self,
                 rio_gdal_options=None,
                 max_workers=40,
                 fill_value=0):


        if rio_gdal_options is None:
            rio_gdal_options = {
                'GDAL_DISABLE_READDIR_ON_OPEN': 'EMPTY_DIR',
                'GDAL_HTTP_MERGE_CONSECUTIVE_RANGES': 'YES',
                'CPL_VSIL_CURL_ALLOWED_EXTENSIONS': 'tif',
                'VSI_CACHE': 'FALSE',
                'VSI_CACHE_SIZE': 0,
                'GDAL_CACHEMAX': 0}

        super().__init__(rio_gdal_options=rio_gdal_options,
                         max_workers=max_workers,
                         fill_value=fill_value)

                
######################################
#    FUNCTIONS (MOSTLY PLOTTING)     #
######################################
path_shapes = '/home/jovyan'



def latlon_to_S2(lat, lon):
    proj_wgs84 = CRS('EPSG:4326')
    point_ll = Point(lon, lat)
    shp_tiles = os.path.join(path_shapes,'S2_tiles_EEA39.shp')
    if not os.path.exists(shp_tiles):
        logging.error('Shapefile S2_tiles_EEA39.shp needs to be located at {}'.format(path_shapes))
        raise
    with shapefile.Reader(shp_tiles) as oVector:
        all_shapes = oVector.shapes()
        all_records = oVector.records()
        for idx, polygon in enumerate(all_shapes):
            if point_ll.within(shapely.geometry.shape(polygon)):
                tile = all_records[idx][0]
                # print(all_records[idx])
                proj_utm = CRS('EPSG:326' + str(tile[0:2]))
                transformer = Transformer.from_crs(proj_wgs84, proj_utm)
                utm_x, utm_y = transformer.transform(lat, lon)
                # print(str(utm_x) + ' , ' + str(utm_y))
                x_offset = utm_x - all_records[idx][2]
                y_offset = all_records[idx][3] - utm_y
                print("Point found in tile {} at col {} line {}".format(tile, int(x_offset / 10), int(y_offset / 10)))
                if (x_offset < 0) or (y_offset < 0):
                    continue
                return (tile, x_offset, y_offset)
    logging.warning('No tile found, please check your coordinates')
    return (None, -1, -1)

def get_reproj_geom_from_shape_file(shape, tile, dest_raster):
    gdf = gpd.read_file(os.path.join(path_shapes,shape))
    geom = gdf['geometry'][0]
    proj_wgs84 = CRS('EPSG:4326')
    proj_utm = CRS('EPSG:326' + str(tile[0:2]))
    transformer = Transformer.from_crs(proj_wgs84, proj_utm , always_xy=True)
    reproj_geom = shapely.ops.transform(transformer.transform,geom)
    
    return reproj_geom, dest_raster.transform
   
    


## FUNCTION TO PLOT AN INTERACTIVE MAP IN THE NOTEBOOK
def create_interactive_map(center=[50.38611400625332, 12.362440398248282]):
    map_layer = basemaps.Esri.WorldImagery  # basemaps.OpenStreetMap.Mapnik, basemaps.OpenTopoMap   
    shp_tiles = os.path.join(path_shapes,r'S2_tiles_EEA39.shp')  #full EEA39 coverage
    if not os.path.exists(shp_tiles):
        logging.error('Shapefile S2_tiles_EEA39.shp needs to be located at {}'.format(path_shapes))
    df_tiles = gpd.read_file(shp_tiles)
    geo_data_tiles = GeoData(geo_dataframe=df_tiles,
                                  style={'color': 'white', 'fillColor': 'white', 'opacity': 1, 'weight': 1.9,
                                         'dashArray': '2', 'fillOpacity': 0.1},
                                  hover_style={'fillColor': 'white', 'fillOpacity': 0.03},
                                  name='Tiles')

    # Set center point of the map and the zoom level
    center = center
    zoom = 4

    # Set up the map
    m = Map(center=center, zoom=zoom)
    WorldImagery_layer = basemap_to_tiles(basemaps.Esri.WorldImagery)
    m.add_layer(WorldImagery_layer)
   
    m.add_layer(geo_data_tiles)
    m.add_control(LayersControl())
    return m


    
### PLOT MAXV BOX AROUND POINT OF INTEREST
def get_centroid(shp_file):
    gdf = gpd.read_file(os.path.join(path_shapes,shp_file))
    centroid = gdf['geometry'][0].centroid
    
    return latlon_to_S2(centroid.coords[0][1], centroid.coords[0][0])

def plot_box(header_file, credentials, config, catalogue, version, tile, xoff, yoff, parameter = 'MAXV', year=2018, box=200):
    data_ex = data_extraction(header_file, credentials, config, catalogue, version, tile, xoff, yoff, N_pixel_window = box, shape=None, years=[str(year)], fill_value=-32768)
    
    data_ex.savetoken()
    
   
    product = list(data_ex.catalogue.get_products(f"copernicus_r_utm-wgs84_10_m_hrvpp-vpp_p_2017-now_{data_ex.version}", 
                           start=date(year, 1, 1), 
                           end=date(year, 1, 1), 
                           tileId=data_ex.tile,
                           productGroupId = 's1',                                                     
                           productType=parameter))[0]

    file = product.data[0].href
    with  rasterio.Env(**data_ex.rio_gdal_options) as env:
        with rasterio.open(file, 'r') as ds:
            try:
                array = ds.read(1, window=Window(int(xoff / 10 - box / 2), int(yoff / 10 - box / 2), box, box))
            except:
                loggin.error(f"Error in reading {file}")
                
    fig, ax = plt.subplots(1)
    ax.imshow(array)
    plt.axvline(x=int(box / 2), color='red')
    plt.axhline(y=int(box / 2), color='red')
    plt.title("VPP-MAXV for {} at {}".format(data_ex.tile, year))
    plt.show(fig)


    ### PLOT THE TIMESAT PHENOLOGY PARAMERTERS ON THE GRAPH
def plot_VPP_ST(phenology_data, VI_data, pt, outdir = None, save_plot = False):
    register_matplotlib_converters()
    #change for easier typing
    df_ST = phenology_data.df_ST
    dict_VI_df = VI_data.dict_VI_df
    df_VPP = phenology_data.df_VPP
    no_data_VPP = phenology_data.no_data_VPP
    
    if len(np.unique(df_VPP.transition_date))==1:
        logging.warning("no seasons have been found")
        return
    
    # Interpolate the ST to a daily-scale (for plotting. Interpolation is done based on the given lenght interval ('daily')  
    df_ST_smoothed = df_ST.asfreq('D').interpolate(method='time')  
    
    
    plt.rcParams['figure.figsize'] = [20, 12]
    fig, ax1 = plt.subplots(figsize=(15, 10))
    if (VI_data.option_NDVI) :
        ax2 = ax1.twinx()

    ax1.plot(df_ST.index, df_ST['ST'], color='blue', label='{}'.format('ST_10_daily'), marker='o')
    ax1.plot(df_ST_smoothed.index, df_ST_smoothed['ST'], color='blue', label='{}'.format('ST_daily_interpolated'))
    if VI_data.option_PPI or VI_data.option_NDVI :
        for VI in list(dict_VI_df.keys()):
            if (VI == 'PPI') and not VI_data.option_PPI:
                ax1.set_ylabel(f'ST value')
                logging.info('PPI data is available but option_PPI = False: PPI will not be plotted')
                continue
            if (VI == 'NDVI') and not VI_data.option_NDVI:
                logging.info('NDVI data is available but option_NDVI = False: NDVI will not be plotted')
                continue
            df_VI = dict_VI_df.get(VI)
            if VI == 'PPI':
                logging.info('Plotting PPI')
                colors = {0: 'gainsboro', 1: 'darkslategray', 2: 'black'}
            else:
                logging.info('Plotting NDVI')
                colors = {0: 'palegreen', 1: 'limegreen', 2: 'darkgreen'}

            df_VI['color'] = df_VI['weight'].apply(lambda x: colors[x])
            for g, b in df_VI.groupby(by='color'):
                if g == colors[2]:
                    fill = 'full'
                    w = '-100%'
                elif g == colors[1]:
                    fill = 'bottom'
                    w = '-50%'
                else:
                    fill = 'none'
                    w = '-0%'
                if VI == 'PPI':
                    ax1.plot(b.index, b['{}'.format(VI)], color=g, linestyle='', fillstyle=fill,
                             label='{}'.format('raw_{}'.format(VI) + w), marker='o')
                    ax1.set_ylabel('ST/PPI value')
                elif len(list(dict_VI_df.keys())) > 1:
                    ax2.plot(b.index, b['{}'.format(VI)], color=g, linestyle='', fillstyle=fill,
                             label='{}'.format('raw_{}'.format(VI) + w), marker='o')
                    ax2.set_ylabel(f'{VI} value')
                    
                else:
                    ax2.plot(b.index, b['{}'.format(VI)], color=g, linestyle='', fillstyle=fill,
                             label='{}'.format('raw_{}'.format(VI) + w), marker='o')
                    ax1.set_ylabel(f'ST value')
                    ax2.set_ylabel(f'{VI} value')

    for year in df_VPP.year.unique():
        SOS_s1 = df_VPP[((df_VPP.year == year) & (df_VPP.phenophase == 'SOSD_s1'))]['transition_date'].values
        SOS_s2 = df_VPP[((df_VPP.year == year) & (df_VPP.phenophase == 'SOSD_s2'))]['transition_date'].values
        MAX_s1 = df_VPP[((df_VPP.year == year) & (df_VPP.phenophase == 'MAXD_s1'))]['transition_date'].values
        MAX_s2 = df_VPP[((df_VPP.year == year) & (df_VPP.phenophase == 'MAXD_s2'))]['transition_date'].values
        EOS_s1 = df_VPP[((df_VPP.year == year) & (df_VPP.phenophase == 'EOSD_s1'))]['transition_date'].values
        EOS_s2 = df_VPP[((df_VPP.year == year) & (df_VPP.phenophase == 'EOSD_s2'))]['transition_date'].values

        if SOS_s1 != int(no_data_VPP):
            SOS_s1_to_datetime = datetime.strptime(str(int(SOS_s1[0])), '%y%j')
            ax1.axvline(x=SOS_s1_to_datetime, color='yellow')

        if SOS_s2 != int(no_data_VPP):
            SOS_s2_to_datetime = datetime.strptime(str(int(SOS_s2[0])), '%y%j')
            ax1.axvline(x=SOS_s2_to_datetime, color='black')

        if MAX_s1 != int(no_data_VPP):
            MAX_s1_to_datetime = datetime.strptime(str(int(MAX_s1[0])), '%y%j')
            ax1.axvline(x=MAX_s1_to_datetime, color='green')

        if MAX_s2 != int(no_data_VPP):
            MAX_s2_to_datetime = datetime.strptime(str(int(MAX_s2[0])), '%y%j')
            ax1.axvline(x=MAX_s2_to_datetime, color='black')

        if EOS_s1 != int(no_data_VPP):
            EOS_s1_to_datetime = datetime.strptime(str(int(EOS_s1[0])), '%y%j')
            ax1.axvline(x=EOS_s1_to_datetime, color='brown')

        if EOS_s2 != int(no_data_VPP):
            EOS_s2_to_datetime = datetime.strptime(str(int(EOS_s2[0])), '%y%j')
            ax1.axvline(x=EOS_s2_to_datetime, color='black')

    # where some data has already been plotted to ax
    handles, labels = ax1.get_legend_handles_labels()
    if VI_data.option_NDVI :
        handles_ax2, labels_ax2 = ax2.get_legend_handles_labels()
        handles = handles + handles_ax2

    # manually define seasonal
    handles.append(mpatches.Patch(color='yellow', linestyle='-', linewidth=1, fill=False, label='SOS_s1'))
    handles.append(mpatches.Patch(color='green', linestyle='-', linewidth=1, fill=False, label='MAX_s1'))
    handles.append(mpatches.Patch(color='brown', linestyle='-', linewidth=1, fill=False, label='EOS_s1'))
    handles.append(mpatches.Patch(color='black', linestyle='-', linewidth=1, fill=False, label='dates_s2'))
    # plot the legend
    ax1.legend(handles=handles, loc='upper right')
    ax1.set_xlabel('Date')
    if phenology_data.N_pixel_window ==0:
        plt.title(
            'HRVPP series {} for shapefile {}, tile:{}'.format(phenology_data.version, phenology_data.shape, phenology_data.tile,
                                                                                                   ))
        plt.tight_layout()

        if save_plot:
            plt.savefig(os.path.join(outdir,'HRVPP_{}_shapefile_{}'.format(phenology_data.version, phenology_data.shape)+ '.png'))
        

    else:
        plt.title(
            'HRVPP series {} for lat:{}, lon:{}, boxsize:{}\n tile:{}, row:{}, col:{}'.format(phenology_data.version, pt[0],
                                                                                                    pt[1], pt[2], phenology_data.tile,
                                                                                                    int(phenology_data.yoff/10), int(phenology_data.xoff/10)))
        plt.tight_layout()

        if save_plot:
            plt.savefig(os.path.join(outdir,'HRVPP_{}_lat_{}_lon_{}_tile_{}_row_{}_col_{}'.format(phenology_data.version, str(pt[0]),
                                                                                                    str(pt[1]), phenology_data.tile,
                                                                                                    str(int(phenology_data.yoff/10)), str(int(phenology_data.xoff/10)))+ '.png'))
        



      