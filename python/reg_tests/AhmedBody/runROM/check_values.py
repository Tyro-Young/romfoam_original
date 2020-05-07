#!/usr/bin/env python
"""
Check if the values predicted by ROM matches the reference value
"""

import sys

f=open('objFuncs.dat','r')
lines=f.readlines()
f.close()


line=lines[0]
cols=line.split()
CD=float(cols[1])

Ref=float(sys.argv[1])

Tol=float(sys.argv[2])

Err=abs(CD-Ref)/abs(Ref)

if Err>Tol:
    print("Failed!")
else:
    print("Success!")

