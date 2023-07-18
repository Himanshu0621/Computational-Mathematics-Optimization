#!/usr/bin/env python
# coding: utf-8

# ![Newton's%20Method.png](attachment:Newton's%20Method.png)

# In[21]:


import numpy as np
import sympy as smp
import math


# In[22]:


x = 0
f = (np.exp(-x) - 10**(-9))


# (i) |f (x
# k
# )| < εa, εa = 10<sup>-3</sup>

# In[23]:


for i in range(100):
    fn = (np.exp(-x) - 10**(-9))
    x_new = x - ((np.exp(-x) - 10**(-9))/(-np.exp(-x)))
    if abs(x_new - x) < 0.001:
        break
    x = x_new
print(fn, i)
        


# 2. |f (x
# k
# )| < εr
# |f (x
# 0
# )| + εa, εa = 10<sup>-3</sup>, εr = 10<sup>-10</sup>

# In[24]:


for i in range(1000):
    fn = (np.exp(-x) - 10**(-9))
    x_new = x - ((np.exp(-x) - 10**(-9))/(-np.exp(-x)))
    if x_new < 10**(-10) * abs(x) + (10**(-3)):
        break
    x = x_new
print(fn, i)
        


# 3. x
# k+1 − x
# k
# | < εa

# In[25]:


for i in range(100):
    fn = (np.exp(-x) - 10**(-9))
    x_new = x - ((np.exp(-x) - 10**(-9))/(-np.exp(-x)))
    if abs(x_new - x) < 10**(-3):
        break
    x = x_new
print(fn, i)
        


# 4. |x
# k+1 − x
# k
# | < εr
# |f (x
# 0
# )| + εa

# In[26]:


for i in range(100):
    fn = (np.exp(-x) - 10**(-9))
    x_new = x - ((np.exp(-x) - 10**(-9))/(-np.exp(-x)))
    if abs(x_new - x) < 10**(-10) * abs(x) + (10**(-3)):
        break
    x = x_new
print(fn, i)
        


# 2.>

# In[27]:


x0, x1, x2 = smp.symbols('x0, x1, x2')


# In[28]:


def f(x):
    return (x1 - 1)**4 + 2*((x1 - 1)**2)*(x2 - 1)**2 + (x2 + 1)**4 - 2*((x2 - 1)**2) - (2*x2 + 1)**2 + 1


# In[29]:


dfdx = smp.diff(f(x), x1)
dfdx


# In[30]:


dfdb = smp.diff(f(x), x2)
dfdb


# In[31]:


Ja = smp.Matrix([dfdx, dfdb])
Ja


# In[32]:


Hess = smp.Matrix([[Ja.diff(x1) , Ja.diff(x2)]])
Hess


# In[33]:


x0 = (1.21, -1.15)
for i in range(10000):
    x=x0
    fun = Hess.subs(x1, x0[0]).subs(x2, x0[1]).inv() * Ja.subs(x1, x0[0]).subs(x2, x0[1])
    x0 = [e1 - e2 for (e1, e2) in zip(x0, fun)]
    
    delx = [e1 - e2 for (e1, e2) in zip(x0, x)]
    absDelx = [abs(elem) for elem in delx]
    if max(absDelx) < 10**(-4):
        minimafound = True
        x = i
        break
        
if  minimafound:
    print(x0)


# er = 0.6

# In[35]:


x0 = (1.21, -1.15)
er = 0.6
for i in range(10000):
    x=x0
    fun = Hess.subs(x1, x0[0]).subs(x2, x0[1]).inv() * Ja.subs(x1, x0[0]).subs(x2, x0[1])
    x0 = [e1 - e2 for (e1, e2) in zip(x0, fun)]
    
    delx = [e1 - e2 for (e1, e2) in zip(x0, x)]
    absDelx = [abs(elem) for elem in delx]
    
    vdel = Ja.subs(x1, x0[0]).subs(x2, x0[1])
    nor_Ja = vdel[0]**2 + vdel[1]**2
    nor_Ja = math.sqrt(nor_Ja)
    
    
    if max(absDelx) < er * nor_Ja + 10**(-4):
        minimafound = True
        x = i
        break
        
if  minimafound:
    print(x0)   


# er = 0.32

# In[36]:


x0 = (1.21, -1.15)
er = 0.32
for i in range(10000):
    x=x0
    fun = Hess.subs(x1, x0[0]).subs(x2, x0[1]).inv() * Ja.subs(x1, x0[0]).subs(x2, x0[1])
    x0 = [e1 - e2 for (e1, e2) in zip(x0, fun)]
    
    delx = [e1 - e2 for (e1, e2) in zip(x0, x)]
    absDelx = [abs(elem) for elem in delx]
    
    vdel = Ja.subs(x1, x0[0]).subs(x2, x0[1])
    nor_Ja = vdel[0]**2 + vdel[1]**2
    nor_Ja = math.sqrt(nor_Ja)
    
    
    if max(absDelx) < er * nor_Ja + 10**(-4):
        minimafound = True
        x = i
        break
        
if  minimafound:
    print(x0)   


# er = 0.0069

# In[38]:


x0 = (1.21, -1.15)
er = 0.0069
for i in range(10000):
    x=x0
    fun = Hess.subs(x1, x0[0]).subs(x2, x0[1]).inv() * Ja.subs(x1, x0[0]).subs(x2, x0[1])
    x0 = [e1 - e2 for (e1, e2) in zip(x0, fun)]
    
    delx = [e1 - e2 for (e1, e2) in zip(x0, x)]
    absDelx = [abs(elem) for elem in delx]
    
    vdel = Ja.subs(x1, x0[0]).subs(x2, x0[1])
    nor_Ja = vdel[0]**2 + vdel[1]**2
    nor_Ja = math.sqrt(nor_Ja)
    
    
    if max(absDelx) < er * nor_Ja + 10**(-4):
        minimafound = True
        x = i
        break
        
if  minimafound:
    print(x0)   


# er = 0.0001

# In[39]:


x0 = (1.21, -1.15)
er = 0.0001
for i in range(10000):
    x=x0
    fun = Hess.subs(x1, x0[0]).subs(x2, x0[1]).inv() * Ja.subs(x1, x0[0]).subs(x2, x0[1])
    x0 = [e1 - e2 for (e1, e2) in zip(x0, fun)]
    
    delx = [e1 - e2 for (e1, e2) in zip(x0, x)]
    absDelx = [abs(elem) for elem in delx]
    
    vdel = Ja.subs(x1, x0[0]).subs(x2, x0[1])
    nor_Ja = vdel[0]**2 + vdel[1]**2
    nor_Ja = math.sqrt(nor_Ja)
    
    
    if max(absDelx) < er * nor_Ja + 10**(-4):
        minimafound = True
        x = i
        break
        
if  minimafound:
    print(x0)   


# In[ ]:




