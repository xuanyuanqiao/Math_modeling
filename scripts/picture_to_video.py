import cv2 as cv
import os

fps = 5.0
#The frequence you want, please use float.
picture_location = './figure/cis/'
#Your picture folder location, note that there is no '/' at the end of path.
#The out put video will bo also saved in this folder.

filename = picture_location + '/video' + '.avi'
filelist = os.listdir(picture_location)
picture_list = []
for file in filelist:
    if ".png" in file:
        picture_list.append(picture_location + '/' + file)
picture_list.sort()
fourcc = cv.VideoWriter_fourcc(*'XVID')
print(picture_list)

start = 0
for picture in picture_list:
    img = cv.imread(picture)
    if start == 0:
        size = img.shape
        out = cv.VideoWriter(filename, fourcc, fps, (int(size[1]), int(size[0])))
        start = 1
    out.write(img)
out.release()
print("Done!")