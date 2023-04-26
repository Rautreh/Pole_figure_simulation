from crystal import Crystal


#We first define a crystal named c1, composed of GaN in the cubic system (zincblende), of the class Crystal. 
#with Miller indices (001)[110]
c1 = Crystal(mat='GaN', system = 'cub', z=[0,0,1], x = [1,1,0])
c1.add_pole([1,1,3], None, None)
c1.save_fig('113 GaN')


c2 = Crystal(mat='GaN', system = 'cub', z=[0,0,1], x = [1,1,0])

#We add 
c2.add_pole(ref=[1,1,3], twin_axis=[1,-1,0], twin_angle=109.5)
c2.add_pole(ref=[1,1,3], twin_axis=[1,-1,0], twin_angle=-109.5)
c2.add_pole(ref=[1,1,3], twin_axis=[1,1,0], twin_angle=109.5)
c2.add_pole(ref=[1,1,3], twin_axis=[1,1,0], twin_angle=-109.5)
c2.save_fig('Twinned GaN')



