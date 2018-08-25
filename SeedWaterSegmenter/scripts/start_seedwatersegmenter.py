from SeedWaterSegmenter.SeedWaterSegmenter import *

InitializeMPL()
app = SegmenterApp(0)

if sys.platform=='darwin':
    # Workaround for Mac bug
    # Create and immediately destroy a file dialog to prevent a crash
    dlg = wx.FileDialog(None,"")
    wx.FutureCall(1, dlg.Destroy)
    dlg.ShowModal()

f = (wx.FileSelector('Select a tiff or gif file') if len(sys.argv)<=1 else
     ' '.join(sys.argv[1:]) if CAT_SPACED_ARGS else
     sys.argv[1])

if not os.path.exists(f):
    print('The specified path does not exist.')
    print(f)
    exit()

d = '' if CAT_SPACED_ARGS or len(sys.argv)<=2 else sys.argv[2]

fig1,ax1,fig2,ax2 = InitMPL()
app.frame.OnInit(filename=f,saveDir=d,fig1=fig1,ax1=ax1,fig2=fig2,ax2=ax2)
app.MainLoop()
