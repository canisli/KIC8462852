import os

# input the list outputted by ap.py baddies = ['./data/elp0m411-kb55-20190827-0154-e91.fits.fz',
# './data/elp0m411-kb55-20190920-0154-e91.fits.fz', './data/tfn0m410-kb23-20190822-0138-e91.fits.fz',
# './data/tfn0m410-kb23-20190831-0066-e91.fits.fz', './data/tfn0m414-kb81-20190916-0134-e91.fits.fz',
# './data/tfn0m410-kb23-20190822-0141-e91.fits.fz', './data/elp0m411-kb55-20190924-0085-e91.fits.fz',
# './data/ogg0m406-kb27-20190815-0145-e91.fits.fz', './data/elp0m411-kb55-20190920-0144-e91.fits.fz',
# './data/elp0m411-kb55-20190823-0103-e91.fits.fz', './data/elp0m411-kb55-20190823-0154-e91.fits.fz',
# './data/tfn0m410-kb23-20190830-0070-e91.fits.fz', './data/elp0m411-kb55-20190920-0152-e91.fits.fz',
# './data/elp0m411-kb55-20190830-0249-e91.fits.fz', './data/elp0m411-kb55-20190823-0157-e91.fits.fz',
# './data/tfn0m410-kb23-20190823-0068-e91.fits.fz', './data/tfn0m414-kb81-20190916-0132-e91.fits.fz',
# './data/tfn0m410-kb23-20190916-0038-e91.fits.fz', './data/ogg0m406-kb27-20190815-0146-e91.fits.fz',
# './data/elp0m411-kb55-20190924-0086-e91.fits.fz', './data/ogg0m406-kb27-20190815-0143-e91.fits.fz',
# './data/tfn0m410-kb23-20190831-0065-e91.fits.fz', './data/elp0m411-kb55-20190827-0157-e91.fits.fz',
# './data/elp0m411-kb55-20190823-0155-e91.fits.fz', './data/ogg0m406-kb27-20190815-0144-e91.fits.fz',
# './data/elp0m411-kb55-20190823-0102-e91.fits.fz', './data/ogg0m406-kb27-20190921-0159-e91.fits.fz',
# './data/tfn0m410-kb23-20190822-0140-e91.fits.fz', './data/tfn0m414-kb81-20190916-0135-e91.fits.fz',
# './data/elp0m411-kb55-20190924-0084-e91.fits.fz', './data/elp0m411-kb55-20190920-0155-e91.fits.fz',
# './data/elp0m411-kb55-20190827-0155-e91.fits.fz', './data/elp0m411-kb55-20190830-0251-e91.fits.fz',
# './data/elp0m411-kb55-20190827-0156-e91.fits.fz', './data/elp0m411-kb55-20190920-0143-e91.fits.fz',
# './data/elp0m411-kb55-20190823-0104-e91.fits.fz', './data/tfn0m410-kb23-20190831-0064-e91.fits.fz',
# './data/elp0m411-kb55-20190924-0087-e91.fits.fz', './data/elp0m411-kb55-20190823-0101-e91.fits.fz',
# './data/tfn0m414-kb81-20190916-0133-e91.fits.fz', './data/elp0m411-kb55-20190823-0156-e91.fits.fz',
# './data/elp0m411-kb55-20190920-0153-e91.fits.fz', './data/elp0m411-kb55-20190830-0248-e91.fits.fz']

bad_files = ['./2018/ogg0m406-kb27-20180308-0144-e91.fits.fz', './2018/tfn0m414-kb99-20180319-0165-e91.fits.fz',
           './2018/ogg0m406-kb27-20180308-0143-e91.fits.fz']

for f in bad_files:
    command = 'mv ' + " " + f + " ./data/bad_files"
    os.system(command)
