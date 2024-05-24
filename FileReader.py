
import numpy as np


def CoffMaker(GGM):
    file = open(GGM, 'r')
    line = file.readline()
    while True:
        line = file.readline()
        if 'end_of_head' in line:
            break

    while True:
        line = file.readline()
        if not line:
            break
        if 'd' in line:
            line = line.replace('d', 'e')
        numbers = [float(num) for num in line[3:].split()]


    file.close()


    num = int(numbers[0]) +1
    anm = np.zeros((num,num))
    bnm = np.zeros((num,num))

    GGM = 'EGM2008.gfc'
    file = open(GGM, 'r')
    line = file.readline()
    while True:
        line = file.readline()
        if 'end_of_head' in line:
            break

    while True:
        line = file.readline()
        if not line:
            break
        if 'd' in line:
            line = line.replace('d', 'e')
        numbers = [float(num) for num in line[3:].split()]

        anm[int(numbers[0]),int(numbers[1])] = numbers[2]
        bnm[int(numbers[0]),int(numbers[1])] = numbers[3]
        

    file.close()
    return anm , bnm , num