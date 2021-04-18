import gams
from gams import GamsWorkspace
from os import getcwd
from os.path import getsize

fileName = 'NTimeReg.csv'
# workspace de gams
ws = GamsWorkspace(getcwd())

with open(fileName, 'w') as reg:
    reg.write(f'past,fo,time,size\n')
    for past in range(1, 24+1):
        dbFile = f'NS_past{past}_delay0.gdx'
        print(f'Processing {dbFile}')
        db = ws.add_database_from_gdx(dbFile)
        fo = db['foval'].first_record().value
        time = db['time_elapsed'].first_record().value
        size = getsize(dbFile)
        reg.write(f'{past},{fo},{time},{size}\n')
input('...')
