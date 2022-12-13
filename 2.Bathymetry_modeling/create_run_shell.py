'''
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'''

import os
import subprocess
import glob

thresh_array = [30,35,40,45,50,55,60,65,70,75]

laser_num = [1, 2, 3]

dir = 'dir to h5 files'

files = glob.glob(dir + '/*.h5')

txt_file = 'run_bathy.sh'
txt = open(txt_file, 'w+')


for file in files:
    for l in laser_num:
        for t in thresh_array:
            txt.write('python run_bathy.py -i ' + str(file) + ' -l ' + str(l) + ' -th ' + str(t) + '\n')
