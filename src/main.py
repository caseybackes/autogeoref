'''
References: 
    https://www.matecdev.com/posts/landsat-sentinel-aws-s3-python.html
    https://stackoverflow.com/questions/24346872/python-equivalent-of-a-given-wget-command
    http://kapadia.github.io/usgs/reference/api.html

Notes:
    1. ca-certificates were installed with 'brew install ca-certificates' to resolve the error:
        "rasterio.errors.RasterioIOError: CURL error: error setting certificate file: /etc/ssl/certs/ca-certificates.crt"

'''
from json import load
# from lib2to3.pgen2 import driver
import satsearch
from pystac_client import Client
import wget
import os
import boto3
import rasterio as rio
from pyproj import Transformer
from rasterio.features import bounds
from rasterio.plot import show
import matplotlib.pyplot as plt
import numpy as np

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('roi',help='short name for roi in roi dir', type = str, default="")
parser.add_argument('-l','--limit',help='limit landsat download count', type = int, default=10)
args = parser.parse_args()

# Set the SSL certificat location
os.environ['CURL_CA_BUNDLE'] = '/etc/ssl/certs/cacert.pem' # see note 1 above

# Define the FOV of the footprint of your recent collect
def fov_extent(center_fov_coordinate_pair):
    '''Takes a tuple of coordinate pairs of center of fov on the ground. 
    Returns widest possible extent of image footprint in geojson format.'''
    # we'll do this when we have some fake gps data, or after we build a 
    # way to simulate fake gps data. 
    pass


def getSubset(aws_session,geotiff_file, bbox):
    with rio.Env(aws_session):
        with rio.open(geotiff_file) as geo_fp:
            # Calculate pixels with PyProj 
            # breakpoint() # look at the .crs component and get the shape and look at a few values
            Transf = Transformer.from_crs("epsg:4326", geo_fp.crs) 
            lat_north, lon_west = Transf.transform(bbox[3], bbox[0])
            lat_south, lon_east = Transf.transform(bbox[1], bbox[2]) 
            x_top, y_top = geo_fp.index( lat_north, lon_west )
            x_bottom, y_bottom = geo_fp.index( lat_south, lon_east )
            # Define window in RasterIO
            window = rio.windows.Window.from_slices( ( x_top, x_bottom ), ( y_top, y_bottom ) )                
            # Actual HTTP range request
            subset = geo_fp.read(1, window=window)
            meta = geo_fp.meta
    return subset,meta


def plotNDVI(nir,nir_meta,red,red_meta,filename):
    if (np.min(nir)<0 or np.min(red) < 0):
        nir = nir+1
        red = red +1
    
    ndvi = (nir-red)/(nir+red)
    ndvi[ndvi>1] = 1
    #### If we want to visualize the image, do so here
    # plt.imshow(ndvi)
    # plt.savefig(filename)
    # plt.close()
    return None

def plotRGB(red,blue,green):
    ''' show the RGB true color composite image from the R,G,B bands'''
    pass

def fetch_landsat(fov_center_coord_pair, fov_extent):
    '''get the latest landsat scene at the coordinate pair center location. build RGB composite. 
    inputs: 
        fov_center_coord_pair(tuple)
        - center of your satellite's field of view(fov) as given by platform ephemeris. 
    '''




    # Initiate pystac client
    LandsatSTAC = Client.open("https://landsatlook.usgs.gov/stac-server", headers=[])





    # Locate geojson file for roi
    '''   
    Needs to be automatically generated based on platform ephemeris, but for now 
    is saved from having used geojson.io'''
    try:
        file_path = list(filter(lambda x: (".geojson" in x and args.roi in x), os.listdir('../rois/')))[0]
    except IndexError as e: 
        print(f"\n\n{'*'*30}   ERROR   {'*'*30}\nNo such results for the given roi of: '{args.roi}'. Exiting.\n\n")
        exit()

    file_content = load(open(os.path.join('../rois',file_path)))
    geometry = file_content["features"][0]["geometry"]
    bbox = bounds(geometry)
    timeRange = '2021-06-01/2022-06-01'
    LandsatSearch       = LandsatSTAC.search ( 
        intersects      = geometry,               # from geojson defined ROI
        datetime        = timeRange,              # time range of interest
        query           = ['eo:cloud_cover95'],   # min acceptable cloud cover 
        collections     = ["landsat-c2l2-sr"] )   # the product to search for (c212-sr = collection 2 level 1 )
    # Convert each returned item to dictionary
    Landsat_items = [i.to_dict() for i in LandsatSearch.get_items()]





    # Initiate AWS Session to use for Landsat file download
    print("Creating AWS Session...\n")
    aws_session = rio.session.AWSSession(boto3.Session(), requester_pays=True)






    # For each item in the collection get the S3 link for R,G,B,NIR,SW1,SW2,and Date
    for i,item in enumerate(Landsat_items[0:args.limit]):
        blue_s3     = item['assets']['blue']['alternate']['s3']['href']
        green_s3    = item['assets']['green']['alternate']['s3']['href']
        red_s3      = item['assets']['red']['alternate']['s3']['href']
        nir_s3      = item['assets']['nir08']['alternate']['s3']['href']
        swir16_s3   = item['assets']['swir16']['alternate']['s3']['href']
        swir22_s3   = item['assets']['swir22']['alternate']['s3']['href']
        date        = item['properties']['datetime'][0:10]

        print("Landsat item number " + str(i) + "/" + str(len(Landsat_items)) + " " + date)
        
        red,red_meta        = getSubset(aws_session, red_s3, bbox)
        nir,nir_meta        = getSubset(aws_session, nir_s3, bbox)
        blue,blue_meta      = getSubset(aws_session, blue_s3, bbox)
        green,green_meta    = getSubset(aws_session, green_s3, bbox)
        swir16,swir16_meta  = getSubset(aws_session, swir16_s3, bbox)
        swir22,swir22_meta  = getSubset(aws_session, swir22_s3, bbox)
        # breakpoint() # looking for crs of the bands, just one but any one will do. 
        rgb_landsat_image = rio.open(
            f"../data/landsat_{args.roi}_{date}_Scene{i}.tif","w",
            driver="GTiff", # save as geotiff
            height=red_meta['height'], # pixels tall
            width=red_meta['width'], # pixels wide
            count = 3, # number of bands for each "RGB" image 
            nodata=red_meta['nodata'],
            dtype = red_meta['dtype'],
            crs = red_meta['crs'],
            transform=red_meta['transform']
        )
        # Using rasterio to write the downloaded bands to tiffs in ../data
        rgb_landsat_image.write(red,1)
        rgb_landsat_image.write(green,2)
        rgb_landsat_image.write(blue,3)
        # rgb_landsat_image.write(nir,4)
        # rgb_landsat_image.write(swir16,5)
        # rgb_landsat_image.write(swir22,6)

        rgb_landsat_image.close()


    return Landsat_items
items = fetch_landsat(None, None)


'''


# Take a list of numbers. 
my_list = [12, 65, 54, 39, 102, 339, 221, 50, 70, ]
  
# use anonymous function to filter and comparing 
# if divisible or not
result = list(filter(lambda x: (x % 13 == 0), my_list)) 

'''