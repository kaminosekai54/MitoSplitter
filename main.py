# Import of needed package
from functions import *



def main():
    setup()
    writeLog("./logs/", "Starting treatement",firstTime=True)
    run()


if __name__ == '__main__':
    main()