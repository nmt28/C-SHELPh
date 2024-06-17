#!/usr/bin/env python
# coding: utf-8

'''
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'''


#####
# This group of functions processes ICESAT2 data and creates a bathymetric model.
# To do this, it follows a number of steps in the form of functions, including:
# 1. Reading data (ReadATL03())
# 2. Orthometrically correcting the dataset (OrthometricCorrection())
# 3. pulling down the data segment ID (getAtl03SegID())
# 4. Bin the data along latitudinal and height gradients (bin_data())
# 5. Calculate sea height (get_sea_height())
# 6. Get water temperature (get_water_temp())
# 7. Correct bathymetry surface for refraction (RefractionCorrection())
# 8. Calculate bathymetry height (get_bath_height())
# 9. produce figures (produce_figures())
#####

import cshelph.run_cshelph
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, default=None, required=True, help="Specify the input ICESAT H5 file")
    parser.add_argument("-l", "--laser", type=str, default=None, required=True, help="Specify the ICESAT-2 laser number (1, 2 or 3)")
    parser.add_argument("-th", "--thresh", type=int, default=None, help="Specify the threshold percentage")
    parser.add_argument("-tl", "--threshlist", nargs='+', default=None, help="Specify a list of thresholds to trial e.g., -tl 20 25 30h")
    parser.add_argument("-o", "--output", type=str, required = False, help="Specify the output location")
    parser.add_argument("-lr", "--lat_res", type=float, default = 10, help="Specify the latitudinal resoltuion (normally 10)")
    parser.add_argument("-hr", "--h_res", type=float, default = 0.5, help="Specify the height resolution (normally 0.5)")
    parser.add_argument("-wt", "--water_temp", type=float, default = None, required = False, help="Specify the water temperature in degrees C")
    parser.add_argument("-slat", "--start_lat", type=float, required = False, help="Specify the start latitude")
    parser.add_argument("-elat", "--end_lat", type=float, required = False, help="Specify the stop latitude")
    parser.add_argument("-minb", "--min_buffer", type=float, default = -40, required = False, help="Specify the stop latitude")
    parser.add_argument("-maxb", "--max_buffer", type=float, default=5, required = False, help="Specify the stop latitude")
    parser.add_argument("-sb", "--surface_buffer", type=float, default=-0.5, required = False, help="Specify the point at which sea surface points are excluded")
    
    args = parser.parse_args()

    if (args.thresh is None) and (args.threshlist is None):
        raise Exception('Must provide either a threshold percentage or a list of thresholds to test.')
    
    
    
    print('''THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.''')

    cshelph.run_cshelph.run_cshelp(args.input, args.laser, args.thresh, args.threshlist, args.output, args.lat_res, args.h_res, args.water_temp, args.start_lat, args.end_lat, args.min_buffer, args.max_buffer, args.surface_buffer)
    
if __name__ == '__main__':
    main()
