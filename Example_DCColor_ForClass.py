#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Test domain coloring example
#
# Using DColor by Anthony Hernandez, 2021 (accessed 12/29/2022)
#  See https://github.com/hernanat/dcolor
#
# I (AKB) made some modifications, see ## GITHUB_LINK_HERE


# In[ ]:


import numpy as np
import sys
sys.path.append('dcolor-master')
import dcolor
dc = dcolor.DColor(xmin=-2,xmax=2,ymin=-2,ymax=2)


# In[ ]:


dc.plot(lambda z : z,title='z')
dc.plot(lambda z : np.sin(z),title='sin(z)')


# In[ ]:


# Some options: include a grid
dc.plot(lambda z: z*np.sin(z*np.pi/2),grid=True)


# In[ ]:


# Different color schemes
dc.plot(lambda z: z*np.sin(z*np.pi/2),grid=True,cscheme='p')

"""
'h':  A. Hernandez's scheme (makes everything look like trippy glazed donuts)
'p': plain phase plot
'm': phase + modulus 
'c': phase+conformal grid 
'd': Standard domain coloring 
'e': enhanced domain coloring 
 """


# In[ ]:


# Choose a different grid size
dc2 = dcolor.DColor(xmin=-4,xmax=4,ymin=-2,ymax=2)

def sampleF(z):
    y = 1j*np.exp(np.pi*z)
    return y

dc2.plot(lambda z: 1j*sampleF(z))


# In[ ]:


# Visualize how 


# In[ ]:





# In[ ]:


### Below: leave this for later. 


# In[ ]:


## Some functions implemented for PPGUI

def PolyTrunc(z,k):
    # PolyTrunc: visualize what happens at the edge of the ROC
    #
    # 1/(1-z), truncated at order k
  
    #If necessary...
    k = np.round(k);

    w = 1;

    for k1 in np.arange(k):
        w = w + (z**(k1+1))
    
    return w
    
dc2.plot(lambda z: PolyTrunc(z,20))


# In[ ]:


def EssSingTrunc( z,k ):
    # EssSingTrunc: visualize what happens near an essential singularity
    #
    #   exp(1/z), truncated at order k
    #  

    # If necessary...
    k = np.round(k);

    w = 1;

    for k1 in np.arange(k):
        w = w + 1/z**(k1+1)/np.math.factorial((k1+1))
    return w

dc3 = dcolor.DColor(xmin=-0.5,xmax=0.5,ymin=-0.5,ymax=0.5)
dc3.plot(lambda z: EssSingTrunc(z,20))


# In[ ]:


dc3.plot(lambda z: EssSingTrunc(z,20), cscheme ='d')


# In[ ]:




