"""
# script for downloading IS2 data with ICEPYX
# built from example at: https://icepyx.readthedocs.io/en/latest/getting_started/example_link.html
"""

import icepyx as ipx

def download_atl03(short_name, spatial_extent, date_range):
    """ Function to search and download IS2"""
    region_a = ipx.Query(short_name, spatial_extent, date_range)

    print('\n print avail granules = ', region_a.avail_granules(), '\n')

    cont = input("Continue...? (yes/no): ")
    if cont == "yes":
        # EarthData credentials
        earthdata_uid = 'bhylee'
        email = 'brianlee52@ucsb.edu'
        session=region_a.earthdata_login(earthdata_uid, email)
        page_size = 10 # This is the number of granules we will request per order.
        # Determine the number of pages based on page size and the number of granules available.
        #If no page_num is specified, this calculation is done automatically to set page_num,
        #which then provides the number of individual orders we will request given the number
        #of granules.
        page_num = 3
        #request_mode = 'async'
        #agent = 'NO'
        #include_meta = 'Y'
        #Must be empty
        path = '/Users/brian.h.lee/Desktop/icesat/data/raw'
        region_a.download_granules(session, path, page_size, page_num)
    elif cont == "no":
        print("answer 'yes' to continue")
    else:
        print('Answer "yes" or "no"')

def main():
    """ run functions """
    #short_name = the dataset of interest, known as its "short name".
    #See https://nsidc.org/data/icesat-2/data-sets for a list of the available datasets.
    #spatial extent = a region of interest to search within. This can be entered as a bounding box,
    #polygon vertex coordinate pairs, or a polygon geospatial file (currently shp, kml,
    #and gpkg are supported).
    #bounding box: Given in decimal degrees for the lower left longitude,
    #lower left latitude, upper right longitude, and upper right latitude
    #polygon vertices: Given as longitude, latitude coordinate pairs
    #of decimal degrees with the last entry a repeat of the first.
    #polygon file: A string containing the full file path and name.
    #date_range = the date range for which you would like to search for results.
    #Must be formatted as a set of 'YYYY-MM-DD' strings.
    short_name = 'ATL03'
    spatial_extent = [-82.14258, 21.60278, -81.47461, 22.18276]
    #spatial_extent = './supporting_files/data-access_PineIsland/glims_polygons.kml'
    date_range = ['2021-02-05','2021-02-07']

    download_atl03(short_name, spatial_extent, date_range)

if __name__=='__main__':
    main()
