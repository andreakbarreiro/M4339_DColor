#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.colors import hsv_to_rgb

class DColor:
    def __init__(self, samples=1000, xmin=-8, xmax=8, ymin=-8, ymax=8):
        mpl.rcParams['toolbar'] = 'None'
        self._samples = samples
        #axes
        self._xmin = xmin
        self._xmax = xmax
        self._ymin = ymin
        self._ymax = ymax
        self.makeDomain()

    def makeDomain(self):
        """Create the domains for Real (x) and Imaginary (y) values respectively"""
        x = np.linspace(self._xmin, self._xmax, self._samples)
        y = np.linspace(self._ymin, self._ymax, self._samples)
        self.xx, self.yy=np.meshgrid(x,y)

    def makeColorModel(self, zz, cscheme='h'):
        """Create the HSV color model for the function domain that will be plotted
        
            cscheme options:
                'h':  A. Hernandez's scheme (makes everything look like trippy glazed donuts)
                'p':  phase only
        """
        
        # DO NOT normalize, this will not work if your results do not include the 
        #   positive real axis where 0->1      
        #H = self.normalize(np.angle(zz) % (2. * np.pi)) # Hue determined by arg(z)
        
        # Instead, divide by 2pi
        H = (np.angle(zz) % (2. * np.pi))/(2. * np.pi) # Hue determined by arg(z)
        r = np.log2(1. + np.abs(zz))                   # Modulus on log scale
                                                       # i.e. map [0, Inf]-> [0, Inf], log scale
        if (cscheme == 'h'):
            S = (1. + np.abs(np.sin(2. * np.pi * r))) / 2.  # Each time r increases by an integer, 
            V = (1. + np.abs(np.cos(2. * np.pi * r))) / 2.  # cycle through light/dark cycle
        elif cscheme == 'p':                                                      
            S = np.ones(r.shape) 
            V = np.ones(r.shape)
        #temp = self.sawfct(r,1,0.7,1)    
        #V = temp     
        
        return H,S,V

    def makeColorModel_RGB(self, zz, cscheme='h'):
        """Create the RGB color model for the function domain that will be plotted
        
            cscheme options: (*** means NOT IMPLEMENTED)
                'h':  A. Hernandez's scheme (makes everything look like trippy glazed donuts)
                'p': plain phase plot
                'm': phase + modulus 
                'c': phase+conformal grid 
                'd': Standard domain coloring 
                'e': enhanced domain coloring 
        """
        if (cscheme in {'h','p'}):
            H,S,V = self.makeColorModel(zz,cscheme)
            rgb   = hsv_to_rgb(np.dstack((H,S,V)))
        elif cscheme in {'m','c'}:
            # Get baseline hue
            H,S,V = self.makeColorModel(zz,'p')
            rgb   = hsv_to_rgb(np.dstack((H,S,V)))          
           
            phaseres = 20
            black   = self.sawfct(np.log(np.abs(zz)),2*np.pi/phaseres,0.75,1)
            if (cscheme == 'c'):
                blackp = self.sawfct(H,1/phaseres,0.75,1)
                black = black*blackp
            
            R = np.multiply(black,rgb[:,:,0])
            G = np.multiply(black,rgb[:,:,1])
            B = np.multiply(black,rgb[:,:,2])
            
            rgb = np.dstack((R,G,B)) 
            rgb = self.BrightenRGB(rgb,.1)
            
        elif cscheme in {'d','e'}:
            # Get baseline hue
            H,S,V = self.makeColorModel(zz,'p')
            rgb   = hsv_to_rgb(np.dstack((H,S,V)))          
           
            logf   = np.log(np.abs(zz));
            phaseres = 20
            if (cscheme == 'e'):
                black   = self.sawfct(logf,2*np.pi/phaseres,0.75,1)
                R = np.multiply(black,rgb[:,:,0])
                G = np.multiply(black,rgb[:,:,1])
                B = np.multiply(black,rgb[:,:,2])
                rgb = np.dstack((R,G,B))
                
            # Now a modulus-dependendent brightening factor         
            bright = (logf>=0)*((1.-(1./(1+logf)))**2)-(logf<0)*(1.-(1./(1-logf)))**2;
            # Both *, / should act ELEMENT-WISE
    
            rgb = self.BrightenRGB(rgb,bright)
                
        # List in GUI_PhasePlot
        # color_schemes = {'p', 'm', 'c', 'd', 'e', 'u', 'v', 'a', 'b', 'x', 'y'};

        # Matching strings??
        #{'plain phase plot','phase + modulus',...
        #  'phase + conformal grid','standard domain coloring', ...
        #  'enhanced domain coloring', 'polar chessboard ', ...
        #  'cartesian chessboard', 'alternating b&w phase', ...
        #  'alternating b&w modulus', ...
        #  'b&w stripes (real part)', 'b&w stripes (imag part)'},
            
        return rgb    
    
    """ From PPGUI:   
        various color schemes for use with phase plots
 
          Usage: rgb = colscheme(f,cs,t,pres)
 
          to get schemes bases on standard hsv coloring, specify string cs as follows: 
 
          a - alternating black and white phase
          b - alternating black and white modulus
      c - phase plot with conformal polar grid
      d - standard domain coloring
      e - enhanced domain coloring
      f - like 'y' but white and blue, especially for stream lines
      i - like 'y' with spacing depending on size of Im f 
      j - colored phase plot with specific phase jumps
      l - like 'y' but spacing in integer fractions of 2 pi
      m - colored phase plot with module jumps
      n - like 'c' - with brighter color for background 
      p - proper phase plot
      q - phase plot colored in steps 
      r - conformal cartesian grid
      s - conformal polar grid
      t - polar chessboard - light gray
      u - polar chessboard
      v - cartesian chessboard
      w - cartesian chessboard - light gray
      x - black and white stripes corresponding to real part
      y - black and white stripes corresponding to imaginary part
 
   """     
    
    def BrightenRGB(self,rgb,bright):
        
        # If bright is a number
        if np.isscalar(bright):
            bright = bright*np.ones(rgb[:,:,0].shape)
        
        R0 = rgb[:,:,0]
        R = (bright>=0)*((1-bright)*R0 + bright*np.ones(R0.shape)) + (bright<0)*((1+bright)*R0)
      
        G0 = rgb[:,:,1]
        G = (bright>=0)*((1-bright)*G0 + bright*np.ones(G0.shape)) + (bright<0)*((1+bright)*G0)
        
        B0 = rgb[:,:,2]
        B = (bright>=0)*((1-bright)*B0 + bright*np.ones(B0.shape)) + (bright<0)*((1+bright)*B0)
        
        rgb = np.dstack((R,G,B)) 
        return rgb
        """                
        function RGB = BrightenRGB(RGB,bright)
% modification of color scheme
% bright - scalar value or field  
% between 0 and 1 for brightening
% between -1 and 0 for darkening 

if size(size(bright))==1
    bright = bright * ones(size(RGB(:,:,1)));
end
    
RGB(:,:,1) = (bright>=0).* ...
  ((1-bright).*RGB(:,:,1) + bright.*ones(size(RGB(:,:,1)))) ...
  + (bright<0).*((1+bright).*RGB(:,:,1));
 % 

RGB(:,:,2) = (bright>=0).* ...
    ((1-bright).*RGB(:,:,2) + bright.*ones(size(RGB(:,:,2)))) ...
   + (bright<0).*((1+bright).*RGB(:,:,2));
 %  

RGB(:,:,3) = (bright>=0).* ...
    ((1-bright).*RGB(:,:,3) + bright.*ones(size(RGB(:,:,3)))) ...
    + (bright<0).*((1+bright).*RGB(:,:,3));

end
"""
    
    def normalize(self, arr):
        """Used for normalizing data in array based on min/max values"""
        arrMin = np.min(arr)
        arrMax = np.max(arr)
        arr = arr - arrMin
        return arr / (arrMax - arrMin)
    
    def plot(self, f, xdim=8, ydim=8, plt_dpi=100,title='',grid=False,cscheme='h'):
        """Plot a complex-valued function
            Arguments:
            f -- a (preferably) lambda-function defining a complex-valued function
            Keyword Arguments:
            xdim -- x dimensions
            ydim -- y dimensions
            plt_dpi -- density of pixels per inch
            
            NEW ARGS:
            title -- duh
            grid -- include grid lines?
            cscheme -- color scheme. See makeColorModel_RGB for definitions
        """
        zz=f(self.z(self.xx,self.yy))
        
        # H,S,V = self.makeColorModel(zz,cscheme)
        # rgb = hsv_to_rgb(np.dstack((H,S,V)))
        
        # We now have a wrapper function to return rgb.
        # Most coloring schemes are most naturally defined by R,G,B; 
        rgb = self.makeColorModel_RGB(zz,cscheme)
        
        fig = plt.figure(figsize=(xdim, ydim), dpi=plt_dpi)
        ax = fig.gca()
        val = str('Re(z) : xmin=')
        val = val + str(self._xmin) + ", xmax=" + str(self._xmax)
        ax.set_xlabel(val)
        val = str('Im(z) : ymin=')
        val = val + str(self._ymin) + ", ymax=" + str(self._ymax)
        ax.set_ylabel(val)
        
        # Get actual dimensions on axes
        # Reverse roles of ymin, ymax
        extent = self._xmin, self._xmax, self._ymax, self._ymin
       
        ax.imshow(rgb, extent=extent)
        
        ax.invert_yaxis() # make CCW orientation positive
        ax.get_xaxis().set_visible(True)
        ax.get_yaxis().set_visible(True)
        ax.set_title(title)
        if grid:
            ax.grid(color='k')
        plt.show()
        
        

    def z(self, x, y):
        """return complex number x+iy
            If inputs are arrays, then it returns an array with corresponding x_j+iy_j values
        """
        return x+1j*y
    
    
    # Added by AKB
    # I used this to test how argument was being normalized
    def evalArg(self, f, xdim=8, ydim=8):
        zz=f(self.z(self.xx,self.yy))
        H = self.normalize(np.angle(zz) % (2. * np.pi)) # Hue determined by arg(z)
        return zz, H
    
    
    ## auxiliary function: sawfct

    #    function y= sawfct(x,dx,a,b)
    #
    #%y = sawfct(x,dx,a,b)
    #% saw tooth function on R with period dx onto [a,b]
    #
    #x = x/dx-floor(x/dx);
    #y = a + (b-a)*x;
    #
    #end

    def sawfct(self,x,dx,a,b):
    #% saw tooth function on R with period dx onto [a,b]

        x = x/dx-np.floor(x/dx)
        y = a + (b-a)*x
        return y

