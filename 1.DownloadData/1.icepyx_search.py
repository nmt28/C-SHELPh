from icepyx import is2class as ipd
#from icepyx import icesat2data as ipd
import os
import shutil

'''
short_name = the dataset of interest, known as its "short name". See https://nsidc.org/data/icesat-2/data-sets for a list of the available datasets.
spatial extent = a region of interest to search within. This can be entered as a bounding box, polygon vertex coordinate pairs, or a polygon geospatial file (currently shp, kml, and gpkg are supported).
bounding box: Given in decimal degrees for the lower left longitude, lower left latitude, upper right longitude, and upper right latitude
polygon vertices: Given as longitude, latitude coordinate pairs of decimal degrees with the last entry a repeat of the first.
polygon file: A string containing the full file path and name.
date_range = the date range for which you would like to search for results. Must be formatted as a set of 'YYYY-MM-DD' strings.
'''
#bounding box
short_name = 'ATL08'
# Bermuda
#spatial_extent = [-65, 32, -64.5, 32.5]
# Puerto Rico
spatial_extent = [-81.513, 25.078, -80.885, 25.928]
date_range = ['2018-01-01','2021-01-21']
# default version is most recent if left blank
version='003'
start_time='00:00:00'
end_time='23:59:59'

'''PASS IN A KML INSTEAD'''
#polygon geospatial file
#short_name = 'ATL06'
#spatial_extent = './supporting_files/data-access_PineIsland/glims_polygons.kml'
#date_range = ['2019-02-22','2019-02-28']

region_a = ipd.Icesat2Data(short_name, spatial_extent, date_range, start_time, end_time, version)


print('\n print avail granules = ', region_a.avail_granules())
print('\n')

Cont = input("Continue...?")
#Number = input("How many...?")
if Cont == "yes":
    # EarthData credentials
    earthdata_uid = 'nathanmthomas'
    email = 'nathan.m.thomas@nasa.gov'
    session=region_a.earthdata_login(earthdata_uid, email)


    page_size = 10 # This is the number of granules we will request per order.
    page_num = 3 # Determine the number of pages based on page size and the number of granules available. If no page_num is specified, this calculation is done automatically to set page_num, which then provides the number of individual orders we will request given the number of granules.
    #request_mode = 'async'
    #agent = 'NO'
    #include_meta = 'Y'


    region_a.build_reqconfig_params('download')#page_size=int(Number))

    region_a.order_granules(session)

    print(region_a.avail_granules())


    #Must be empty
    path = '/Users/nmthoma1/Documents/Research/ICESat2/MangroveHeight/ATL08/'
    region_a.download_granules(session, path)

    #Clean up Outputs folder by removing individual granule folders
    '''
    for root, dirs, files in os.walk(path, topdown=False):
        for file in files:
            try:
                shutil.move(os.path.join(root, file), path)
            except OSError:
                pass

    for root, dirs, files in os.walk(path):
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    '''
else:
    print("Closed")



