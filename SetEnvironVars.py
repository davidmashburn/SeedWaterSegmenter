import os

d = 'C:/WHEREVER YOU PUT THE FOLDER...'

os.environ['MYPYTHON'] = d
os.environ['MYPYREX'] = os.path.join(os.environ['MYPYTHON'],'Pyx')
