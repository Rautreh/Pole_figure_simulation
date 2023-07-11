from crystal import Crystal

#We first define a crystal named c1, composed of GaN in the cubic system 
#zincblende), of the class Crystal. 
#with Miller indices (001)[110]
c1 = Crystal(mat='GaN', system = 'hex', z=[0,0,1], x = [1,1,0])

#We add the 113 pole
c1.add_pole([1,0,0])
#and save it under the name "113 GaN"
#c1.save_fig('113 GaN')

#We define another crystal named c2, composed of GaN in the cubic system 
#(zincblende), of the class Crystal. 
#with Miller indices (001)[110]
c2 = Crystal(mat='GaN', system = 'cub', z=[0,0,1], x = [1,1,0])

#We add the four twins about the <110> axes by +- 109.5 degrees.
c2.add_pole(ref=[1,1,3], twin_axis=[1,-1,0], twin_angle=109.5)
c2.add_pole(ref=[1,1,3], twin_axis=[1,-1,0], twin_angle=-109.5)
c2.add_pole(ref=[1,1,3], twin_axis=[1,1,0], twin_angle=109.5)
c2.add_pole(ref=[1,1,3], twin_axis=[1,1,0], twin_angle=-109.5)

#and save the file under the name "Twinned GaN".
#c2.save_fig('Twinned GaN')

c3= Crystal(mat='Ga2O3', system = 'mon', z=[4,0,0], x = [1,1,0])

c3.add_pole([1, 1, 0], None, None)
c3.add_pole([1, 1, 0], [4,0,0], 90)
c3.add_pole([1, 1, 0], [4,0,0], 180)
c3.add_pole([1, 1, 0], [4,0,0], 270)

c3.add_pole([4, 0,-1], None, None)
c3.add_pole([4, 0,-1], [4,0,0], 90)
c3.add_pole([4, 0,-1], [4,0,0], 180)
c3.add_pole([4, 0,-1], [4,0,0], 270)


c4= Crystal(mat='Ga2O3', system = 'mon', z=[-1,1,2], x = [4,0,0])
c4.rotate_plot([-1,1,2], 45)
c4.add_pole([1, 1, 0], None, None)
c4.add_pole([1, 1, 0], [-1,1,2], 90)
c4.add_pole([1, 1, 0], [-1,1,2], 180)
c4.add_pole([1, 1, 0], [-1,1,2], 270)

c4.add_pole([4, 0,-1], None, None)
c4.add_pole([4, 0,-1], [-1,1,2], 90)
c4.add_pole([4, 0,-1], [-1,1,2], 180)
c4.add_pole([4, 0,-1], [-1,1,2], 270)

c4.add_pole([4, 0, 0], None, None)
c4.add_pole([4, 0, 0], [-1,1,2], 90)
c4.add_pole([4, 0, 0], [-1,1,2], 180)
c4.add_pole([4, 0, 0], [-1,1,2], 270)



