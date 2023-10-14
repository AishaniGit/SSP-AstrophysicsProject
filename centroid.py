import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
def findCentroid(fits_file, target_x, target_y, radius=3, sky_radius=5):
# YOUR CODE HERE
    data = fits.getdata(fits_file)
    plt.imshow(data, origin = "lower", cmap = "gray_r", vmin = 2500, vmax = 3000)
    plt.show()
    x = target_x
    y = target_y
    r = sky_radius
    slice = data[target_y-sky_radius:target_y+sky_radius, target_x-sky_radius:target_x+sky_radius]
    plt.imshow(slice)
    plt.show()
    ySum = 0
    xSum = 0
    skyweights = 0
    weights = 0
    for o in range(target_y-radius-sky_radius, target_y+radius+sky_radius+1):
        skyweights += sum(data[o,target_x-radius-sky_radius:target_x+radius+sky_radius+1])
    
    for x in range(target_y-radius, target_y+radius+1):
        weights += sum(data[x,target_x-radius:target_x+radius+1])
    avg = (skyweights - weights) / (((radius*2 + 1) + sky_radius*2)**2 - (radius*2 + 1)**2) #fix hardcode
    print(avg)

    for r in range(target_y-radius, target_y+radius+1):
        ySum += sum(data[r, target_x-radius:target_x+radius+1])  * r 
    y_centroid = ySum / weights
    for c in range(target_x-radius,target_x+radius+1):
        xSum += sum(data[target_y-radius:target_y+radius+1, c]) * c 
    x_centroid = xSum / weights
    return (x_centroid, y_centroid)

    
    #return x_centroid, y_centroid 
print(findCentroid("D:\\sampleimage.fits", 351, 154, 3, 5))
#centroid_x, centroid_y = findCentroid("sampleimage.fits", 351,
#154, 3, 5)

#if abs(centroid_x - 350.9958) < 1e-3 and abs(centroid_y - 153.9955) < 1e-3:
 #   print("centroid calculation CORRECT")
#else:
 #   print("centroid calculation INCORRECT, expected (350.9958, 153.9955), got ({}, {))".format()
        
    

