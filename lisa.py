from os import environ; environ['USE_PYGEOS'] = '0'
from geopandas import read_file
from matplotlib.pyplot import savefig
from esda.moran import Moran, Moran_Local
from splot.esda import plot_moran, lisa_cluster
from libpysal.weights import DistanceBand, min_threshold_distance

# p value threshold for LISA outputs
significance = 0.05

# convert LISA output to useful labels
quadList = ["NA", "HH", "LH", "LL", "HL"]

def getQuadrants(qs, sigs, acceptable_sig):
    """
    * Return list of quadrant codes depending upon specified significance level
    """
    # return quad code rather than number
    out = []
    for q in range(len(qs)):
        # overrride non-significant values as N/A
        if sigs[q] < acceptable_sig:
            out.append(quadList[qs[q]])
        else:
            out.append(quadList[0])
    return out

# open file
tanks = read_file("./out/tank_wetness.shp")

# calculate min distance
centroids = tanks.geometry.centroid
min_threshold = min_threshold_distance(list(zip(centroids.x, centroids.y)))
print(f'min threshold distance: {min_threshold}')

# calculate and row standardise weights matrix
W = DistanceBand.from_dataframe(tanks, min_threshold)
W.transform = 'r'
print(f"max neighbours: {W.max_neighbors}, min neighbours: {W.min_neighbors}")

# calculate and report global Moran's I (Null Hypothesis: CSR)
mi = Moran(
    tanks['evans_coef'], 
    W, 
    permutations=9999
    )
print(f"\nGlobal Moran's I Results (evans_coef)")
print("I:\t\t\t", mi.I)					   # value of Moran's I
print("Expected I:\t\t", mi.EI)			   # expected Moran's I
print("Simulated p:\t\t", mi.p_sim, "\n")  # simulated p

# plot output
plot_moran(mi, zstandard=True, figsize=(10,4))
savefig(f"./moran.png")


# only do LISA if significant
if mi.p_sim < 0.005:

    # calculate local I (Rate Adjusted for Population)
    lisa = Moran_Local(
        tanks['evans_coef'], 
        W, 
        transformation='R', 
        permutations=9999)

    # update GeoDataFrame
    tanks['Morans_I'] = lisa.Is                                        # value of Moran's I
    tanks['sig'] = lisa.p_sim                                          # simulated p
    tanks['quadrant'] = getQuadrants(lisa.q, lisa.p_sim, significance) # quadrant (HH, HL, LH, LL)

    # plot local moran figure
    lisa_cluster(lisa, tanks, significance)
    savefig(f"./lisa.png", bbox_inches='tight')
    
    # esport file
    tanks.to_file('tanks_lisa.shp')
print("done")