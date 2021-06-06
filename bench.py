#!/usr/bin/env python3

import os, time

TIMES = 5
SIZE = 1000

ponyflags =  " --ponyminthreads 2 --ponysuspendthreshold 1000"
ponyflags += " --ponynoyield --ponynoblock --ponypin"

sum = 0
maxt = 0
mint = 0
i = 0
while i < TIMES:
    a = time.clock_gettime_ns(time.CLOCK_PROCESS_CPUTIME_ID)
    r = os.system("./cmatrix" + ponyflags)
    b = time.clock_gettime_ns(time.CLOCK_PROCESS_CPUTIME_ID)
    elapsed = round((b - a) / 1000)

    if r == 0:
        print("#" + str(i+1) + " lap: " + str(elapsed) + " ms elapsed")
        i += 1
        sum += elapsed
        if maxt < elapsed:
            maxt = elapsed
        if (mint > elapsed) or (i == 1):
            mint = elapsed
    else:
        print("killed")

av = round(sum / TIMES)
print("min " + str(mint) + " ms\tav " + str(av) + " ms\tmax " + str(maxt) + " ms")

GALL = 2 * (SIZE) ** 3 / (1024) ** 3
gmin = round(1000 * GALL / maxt, 2)
gav = round(1000 * GALL / av, 2)
gmax = round(1000 * GALL / mint, 2)
print("min " + str(gmin) + " GFLOPS\tav " + str(gav) +
      " GFLOPS\tmax " + str(gmax) + " GFLOPS")
