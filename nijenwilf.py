# Nijenwilf.py, written 21 May 2021 by Vince Vatter

# The code uses the Nijenhuis--Wilf recurrence to generate what we call
# "weakly indecomposable matchings", which are called "irreducible chord diagrams"
# in the OEIS and are counted by sequence A000699.

# It then filters these to find what we call the "strongly indecomposable matchings".
# These were not counted in the OEIS at the time this code was written.

from itertools import product
from tqdm import tqdm

def nijenwilf(mo, mi, k):
  newmo = (mo[0]+k,) + tuple(map(lambda x: x+len(mi) if x>mo[0] else x, mo[1:]))
  newmi = tuple(map(lambda x: x+mo[0]+1 if x>=k else x+mo[0], mi))
  return(newmo[:mo[0]] + newmi[:k] + (0,) + newmi[k:] + newmo[mo[0]+1:])

def weakindecomps(n):
  if n == 1:
    yield tuple([1,0])
  WI = set()
  for m in range(1,n):
    for (mo, mi, k) in product(weakindecomps(m), weakindecomps(n-m), range(1,2*(n-m))):
      yield(nijenwilf(mo, mi, k))

def isstrongindecomp(mat):
  # Check for intervals of size 2 first
  for i in range(len(mat)-1):
    if abs(mat[i] - mat[i+1]) == 1:
      return(False)
  for ln in range(4, int(len(mat)/2)):
    # ln will be the length of the interval
    for i in range(len(mat)-ln+1):
      Ends = [mat[j] for j in range(i,i+ln)]
      if max(Ends) - min(Ends) + 1 == ln:
        return(False)
  return(True)

def strongindecomps(n):
  print("Counting strongly indecomposable matchings with", n, "edges.")
  for mat in tqdm(weakindecomps(n), total=numweak(n), smoothing=0):
    if isstrongindecomp(mat):
      yield(mat)

def numweak(n):
  if n == 1:
    return(1)
  return((n-1) * sum(numweak(k)*numweak(n-k) for k in range(1,n)))

def numstrong(n):
  return(sum(1 for _ in strongindecomps(n)))