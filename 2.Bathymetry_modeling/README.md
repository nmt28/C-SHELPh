'''
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'''

Here are the primary scripts for constructing C-SHELPh commands:

Bathy\_utils.py: Contains all of the functions that C-SHELPh needs

run\_bathy.py: Calls the functions from Bathy\_utils.py in order to extract bathymetry photons

create\_run\_shell.py: A helper script to iterate over a range of threshold values and IS2 lasers. Creates run\_bathy.sh.

Note: The refraction correction currently uses a simplificatio of the incidence angle as outlined by Parrish et al, 2019. For greater precision, the curvature of the Earth needs to be accounted for along teh ICESat-2 track. This will be integrated into C-SHELPh as time permits.
