import os
def createDir(path):
    try: 
        os.mkdir(path)
    except OSError as error:
        pass
