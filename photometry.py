import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
def photometry(image,x, y, ap_radius, in_an, out_an, read_noise, dark_current):
    data = fits.getdata(image)
    plt.imshow(data, origin = "lower", cmap = "gray_r", vmin = 500, vmax = 3000)
    plt.show()
    centroid_values = findCentroid(image, x, y, radius = ap_radius)
    X = centroid_values[0]
    Y = centroid_values[1]
    ADU = 0
    nap = 0
    nan = 0
    skyADU = 0
    for i in range(int(Y) - ap_radius - 1, int(Y) + ap_radius + 2):
        for j in range(int(X) - ap_radius - 1, int(X) + ap_radius + 2):
            if ((i - Y)**2 + (j - X)**2)**(1/2) < ap_radius:
                ADU += data[i,j]
                nap += 1
    skyADUlist = []
    for i in range(int(Y) - out_an - 1, int(Y) + out_an + 2):
        for j in range(int(X) - out_an - 1, int(X) + out_an+ 2):
            if ((i - Y)**2 + (j - X)**2)**(1/2) < out_an and  ((i - Y)**2 + (j - X)**2)**(1/2) > in_an:
                nan += 1
                skyADUlist.append(data[i,j])
    skyArray = np.array(skyADUlist)
    Skye = np.median(skyArray)
    Se = ADU - Skye*nap

    De = dark_current
    R2 = read_noise**2
    SNR = Se / ((Se + nap*(1 + (nap / nan))*(Skye + De + R2))**(1/2))

    uncertainty = 1.0875 / SNR

    minst = -2.5*np.log10(Se)
    print("Centroid:", X, Y)
    print("Signal", Se)
    print("SNR", SNR)
    print("minst:", minst)
    print("Error:", uncertainty)
    


def findCentroid(fits_file, target_x, target_y, radius=3, sky_radius=5):
    data = fits.getdata(fits_file)
    plt.imshow(data, origin = "lower", cmap = "gray_r", vmin = 500, vmax = 3000)
    plt.show()
    x = target_x
    y = target_y
    r = sky_radius
    ySum = 0
    xSum = 0
    skyweights = 0
    weights = 0
    for o in range(target_y-radius-sky_radius, target_y+radius+sky_radius+1):
        skyweights += sum(data[o,target_x-radius-sky_radius:target_x+radius+sky_radius+1])
    
    for x in range(target_y-radius, target_y+radius+1):
        weights += sum(data[x,target_x-radius:target_x+radius+1])
    avg = (skyweights - weights) / (((radius*2 + 1) + sky_radius*2)**2 - (radius*2 + 1)**2)

    for r in range(target_y-radius, target_y+radius+1):
        ySum += sum(data[r, target_x-radius:target_x+radius+1])  * r 
    y_centroid = (ySum - avg*len(range(target_y-radius, target_y+radius+1))) / weights
    for c in range(target_x-radius,target_x+radius+1):
        xSum += sum(data[target_y-radius:target_y+radius+1, c]) * c 
    x_centroid = (xSum - avg*len(range(target_x-radius,target_x+radius+1)))/ weights
    return (x_centroid, y_centroid)

photometry("D:\\aptest.fit", 490, 293, 5, 8, 13, 11, 10)
