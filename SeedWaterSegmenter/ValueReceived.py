# Value Received is an enhancer for a matplotlib image axis that displays the value as well as the x,y
import matplotlib.pyplot as plt

def GetReportPixel(arr):
    def report_pixel(x,y):
        s=arr.shape
        x=int(x+0.5)
        y=int(y+0.5)
        if 0<x<s[1] and 0<y<s[0]:
            return "value=" + str(arr[y,x]) + "  x=" + str(x) + " y=" + str(y)
        else:
            return "x=" + str(x) + " y=" + str(y)
    return report_pixel
def ValueReceived(axis,arr):
    axis.format_coord = GetReportPixel(arr)

def imshow_valuereceived(arr,*args,**kwds):
    plt.imshow(arr,*args,**kwds)
    plt.gca().format_coord = GetReportPixel(arr)

imshow_vr = imshow_valuereceived # synonym

# Most usage would be:
# from ValueReceived import imshow_vr
# (use imshow_vr just like imshow)
