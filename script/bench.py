#!/usr/bin/env python3

import os, time

TIMES = 5
SIZE = 2000

ponyflags =  " --ponyminthreads 2"

sum = 0
maxt = 0
mint = 0
i = 0
while i < TIMES:
    a = time.time_ns()
    r = os.system("./ponymath" + ponyflags)
    b = time.time_ns()
    elapsed = round((b - a) / 10 ** 6)

    if r == 0:
        print("#" + str(i+1) + " lap: " + str(elapsed) + " ms elapsed")
        sum += elapsed
        if maxt < elapsed:
            maxt = elapsed
        if (mint > elapsed) or (i == 1):
            mint = elapsed
    else:
        print("killed")
    
    i += 1

av = round(sum / TIMES)
print("min " + str(mint) + " ms\tav " + str(av) + " ms\tmax " + str(maxt) + " ms")

GALL = 2 * (SIZE) ** 3 / (1024) ** 3
gmin = round(1000 * GALL / maxt, 2)
gav = round(1000 * GALL / av, 2)
gmax = round(1000 * GALL / mint, 2)
print("min " + str(gmin) + " GFLOPS\tav " + str(gav) +
      " GFLOPS\tmax " + str(gmax) + " GFLOPS")
