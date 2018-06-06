import random
from functools import reduce
n = 2000
mod = 9001
T = 2000
degg = random.randint(1,n-1)
degh = random.randint(1,n-1)
coeffs = reduce(lambda x,y: str(x)+" " +str(y), list(map(lambda x: random.randint(0,mod),list(range(0,degg)))))
print(T)
print(f"{degg} {mod} {coeffs}")
coeffs = reduce(lambda x,y: str(x)+" " +str(y), list(map(lambda x: random.randint(0,mod),list(range(0,degh)))))
print(f"{degh} {mod} {coeffs}")
coeffs = reduce(lambda x,y: str(x)+" " +str(y), list(map(lambda x: random.randint(0,mod),list(range(0,n)))))
print(f"{n} {mod} {coeffs}")
for i in range(1,T):
    degg = random.randint(1,n-1)
    coeffs = reduce(lambda x,y: str(x)+" " +str(y), list(map(lambda x: random.randint(0,mod),list(range(0,degg)))))
    print(f"{degg} {mod} {coeffs}")
   
