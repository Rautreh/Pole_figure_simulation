# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 20:08:57 2023

@author: rautr
"""

                #color = ['#EE3377', '#004488', '#997700', '#009988',  '#CC3311', '#66CCEE']
                #color = ['r', 'r','#004488','r', '#997700', '#997700', '#009988', '#009988', '#CC3311', '#CC3311', '#66CCEE', '#EE3377']
                
                #color = ['r', 'r','#004488','r', '#997700', '#997700', '#CC3311', 'g', '#66CCEE', '#009988', '#66CCEE', '#EE3377']

from crystal import Crystal

#We first define a crystal named c1, composed of GaN in the cubic system 
#zincblende), of the class Crystal. 
#with Miller indices (001)[110]
c1 = Crystal(mat='GaN', system = 'hex', z=[0,0,1], x = [1,0,0])

#We add the 113 pole
#c1.rotate_plot([0,0,1], 90)
x=0
x=54.7
c1.add_pole(ref=[1,0,3], twin_axis='y', twin_angle=54.7)

vr = c1.rotate_vector(c1.z_vector, c1.y_vector, -x)
ivr = c1.vector_to_indices(vr)

