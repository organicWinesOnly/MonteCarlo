import numpy as np
import csv
import matplotlib.pyplot as mp

timess = np.arange(4000, 6999)
value = []
with open('corr2.4.csv') as f:
    reader = csv.reader(f)
    for line in reader:
        value.extend(line)

val_prime = list(map(lambda x: float(x), value))
mp.plot(val_prime)
mp.xlim(0, 10)
mp.show()
