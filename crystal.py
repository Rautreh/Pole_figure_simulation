import os
from collections import Counter
from itertools import product, permutations
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import lines
import pandas as pd
from database import database

class Crystal:
    def __init__(self, mat='GaN', system='cub', z=[0,0,1], x = [1,0,0], approximation='x'):
        self.mat = mat
        self.system = system
        self.z_indices = np.array([z])
        self.x_indices = np.array([x])
        self.approximation = approximation
        self.rotation = False
        self.rotation_axis = None
        self.rotation_angle = None
        size = 35
        self.markers = {'1': {'color': 'r', 'marker': 'o', 'size': size}, 
                        '2':{'color': 'b', 'marker': 's', 'size': size},
                        '3':{'color':'m', 'marker': '^', 'size': size}, 
                        '4':{'color':'g', 'marker': 'v', 'size': size},
                        '5': {'color': 'r', 'marker': 'o', 'size': size}, 
                        '6':{'color': 'b', 'marker': 's', 'size': size}, 
                        '7':{'color':'m', 'marker': '^', 'size': size}, 
                        '8':{'color':'g', 'marker': 'v', 'size': size},
                        '9':{'color':'m', 'marker': '^', 'size': size}, 
                        '10':{'color':'g', 'marker': 'v', 'size': size},
                        '11': {'color': 'r', 'marker': 'o', 'size': size}, 
                        '12':{'color': 'b', 'marker': 's', 'size': size}, 
                        '13':{'color':'m', 'marker': '^', 'size': size}, 
                        '14':{'color':'g', 'marker': 'v', 'size': size},
                        '15': {'color': 'r', 'marker': 'o', 'size': size}, 
                        '16':{'color': 'b', 'marker': 's', 'size': size}, 
                        '17':{'color':'m', 'marker': '^', 'size': size}, 
                        '18':{'color':'g', 'marker': 'v', 'size': size},
                        '19':{'color':'m', 'marker': '^', 'size': size}, 
                        '20':{'color':'g', 'marker': 'v', 'size': size}}
        self.dfs = {}
        self.get_database()
        self.lattice_parameters()
        self.orient()
        self.fig, self.ax = self.set_plot()
        
    def get_database(self):
        self.data_base_latice_paramaeters = database
        
        self.real_spacing = np.array([
            [self.data_base_latice_paramaeters[self.mat][self.system]['a'],
             self.data_base_latice_paramaeters[self.mat][self.system]['b'],
             self.data_base_latice_paramaeters[self.mat][self.system]['c']]])
        
        self.real_ang = np.array(
            [self.data_base_latice_paramaeters[self.mat][self.system]['alpha'],
             self.data_base_latice_paramaeters[self.mat][self.system]['beta'],
             self.data_base_latice_paramaeters[self.mat][self.system]['gamma']])
       
        self.real_ang_rad = np.deg2rad(self.real_ang)
        
    def lattice_parameters(self):
        
        cos = np.cos(self.real_ang_rad)
        sin = np.sin(self.real_ang_rad)
        
        if ((cos[1]*cos[2]-cos[0])/(sin[1]*sin[2]) > 1 or 
        (cos[2]*cos[0]-cos[1])/(sin[2]*sin[0]) > 1 or 
        (cos[0]*cos[1]-cos[2])/(sin[0]*sin[1]) > 1):
            raise ValueError('Please make sure that:\n0 < (cos(alpha_j)*cos' 
                             '(alpha_k)-cos(alpha_i))/(sin(alpha_j)*sin(alpha_k))'
                             ' < pi\nfor i =/= j =/= k, with alpha and beta the'
                             ' real and reciprocal angles')
        
        self.recip_ang_rad = np.arccos([
                            (cos[1]*cos[2]-cos[0])/(sin[1]*sin[2]),
                            (cos[2]*cos[0]-cos[1])/(sin[2]*sin[0]),
                            (cos[0]*cos[1]-cos[2])/(sin[0]*sin[1])])
        
        self.recip_ang = np.rad2deg(self.recip_ang_rad)
        #in form such that matrix multiplication is correct
        print(self.recip_ang_rad)
        self.real_vectors = self.real_spacing * np.array([
                            [1, cos[2],   cos[1]],
                            [0, sin[2],   -sin[1]*np.cos(self.recip_ang_rad[0])],
                            [0, 0,        sin[1]*np.sin(self.recip_ang_rad[0])]])
        
        self.real_vectors_inv = np.linalg.inv(self.real_vectors)
        
        #transposed so that each element of real_vectors is a real vector a, b and c
        real_vectors = np.transpose(self.real_vectors)
        Volume = np.linalg.det(self.real_vectors)
        #transposed back to matrix multiply with miller indices
        
        self.recip_vectors = np.transpose(np.array([
                            np.cross(real_vectors[1], real_vectors[2]),
                            np.cross(real_vectors[2], real_vectors[0]),
                            np.cross(real_vectors[0], real_vectors[1])])/Volume)
        
        self.b_inverse = np.linalg.inv(self.recip_vectors)
        print(self.recip_vectors)
        
    def indices_to_vectors(self, indices):
        return np.sum(indices * self.recip_vectors, axis=1)
    
    def vector_to_indices(self, vector):
        return np.sum(vector * self.b_inverse, axis=1)
    
    def interplanar_distance(self, indices):
        H = np.sum(indices * self.recip_vectors, axis=1)
        if np.linalg.norm(H) == 0:
            return 0
        # Interplanar distance
        d_hkl = 1 / np.linalg.norm(H)
        return d_hkl
        
    def orient(self):
        self.z_vector = self.indices_to_vectors(self.z_indices)
        self.x_vector = self.indices_to_vectors(self.x_indices)
        
        self.y_vector = np.cross(self.z_vector, self.x_vector)
        
        angle_zx = self.angle_between_vectors(self.z_vector, self.x_vector)

        if  angle_zx == 0 or angle_zx == 180:
            raise ValueError('z and x directions are the same.')
        
        if self.approximation != 'z':
            self.x_vector = self.rotate_vector(self.x_vector, 
                                               self.y_vector, 90 - angle_zx )
        else:
            self.z_vector = self.rotate_vector(self.z_vector, 
                                               self.y_vector, angle_zx - 90 )
        
    def angle_between_vectors(self, v1, v2):
        a = np.dot(v1, v2)/(np.linalg.norm(v1) * np.linalg.norm(v2))
        a_angle = np.rad2deg(np.arccos(np.round(a, 5)))
        return a_angle

    def rotate_vector(self, v, u, ang):
        #v = v/np.linalg.norm(v)
        u = u/np.linalg.norm(u)
        angle = np.deg2rad(ang)
        if angle == 0:
            return v

        u = np.sin(angle/2)*u
        
        u0 = np.cos(angle/2)
        
        u_conj = -u
        
        uv = u0*v+np.cross(u,v)
        uv0 = -np.dot(u,v)
        
        vr = uv0*u_conj + u0*uv + np.cross(uv, u_conj)
        
        return vr
    
    def planes_in_family(self, ref, twin=False, twin_axis=None, twin_angle=None):
        '''
        Parameters
        ----------
        ref : List
            DESCRIPTION. Reflection hkl
        phase : string, optional
            DESCRIPTION. 'c' for cubic, 'h' for hexagonal. The default is 'c'.

        Returns
        -------
        planes list of tuples: 
            DESCRIPTION. List of all planes in the family.

        '''
        planes_in_family = []
        plane_in_df = False
        df = pd.DataFrame()
        alphas = []
        betas = []
        if len(self.dfs) != 0:
            for key in self.dfs:
                if (self.dfs[key]['Twin_axis'] == twin_axis and 
                    self.dfs[key]['Twin_angle'] == twin_angle):
            
                    if self.dfs[key]['ref'] == ref:
                        return
                    elif tuple(ref) in self.dfs[key]['df']['Planes'].tolist():
                        return
            for key in self.dfs:    
                if self.dfs[key]['ref'] == ref:
                    planes_in_family = self.dfs[key]['df']['Planes']
                    plane_in_df = True
                    break
                
                elif tuple(ref) in self.dfs[key]['df']['Planes'].tolist():
                    planes_in_family = self.dfs[key]['df']['Planes']
                    plane_in_df = True
                    break
                
        if self.system == 'hex' and not plane_in_df:
            refh = np.array([ref[0], ref[1], -ref[0] - ref[1], ref[2]])
            planesh = []
            perm = set([*Counter(permutations(refh[0:3]))] + 
                       [*Counter(permutations(-refh[0:3]))])
            for i in perm:
                lst = [*i]
                lst.append(ref[-1])
                planesh.append(lst)
                if ref[-1] != 0:
                    lst = [*i]
                    lst.append(-ref[-1])
                    planesh.append(lst)
            planes_in_family = []

            for plane in planesh:
                planeh = plane[0:2] + plane[3:]
                planes_in_family.append(planeh)  # remove i from hkil
                


        elif not plane_in_df:
            d_hkl = self.interplanar_distance(ref)
            # Takes indices of plane
            indices = [*Counter(ref)]  
            # Adds its negative indices
            indices.extend(-np.array(indices))  
            # Obtains all possible combinations
            plane_combination = [*Counter([*product(indices, repeat=3)])] 
            for plane in plane_combination:
                d_hkl_pl = self.interplanar_distance(plane)
                # If interplanar distance is the same, belongs to the same family
                if abs((d_hkl-d_hkl_pl)/d_hkl) < 0.001: 
                    planes_in_family.append(plane)
        
        for plane in planes_in_family:
            alpha, beta = self.calc_alpha_beta(
                        plane, twin, twin_axis, 
                        twin_angle, self.rotation, 
                        self.rotation_axis, self.rotation_angle)
            
            alphas.append(alpha)
            betas.append(beta)
        df['Planes'] = planes_in_family
        df['Alpha'] = alphas
        df['Beta'] = betas
        
        df_nr = str(len(self.dfs)+1)
        self.dfs[df_nr] = {'ref':'', 'Twin': twin, 'Twin_axis':None, 'Twin_angle':None}
        
        self.dfs[df_nr]['ref'] = ref
        self.dfs[df_nr]['Twin_axis'] = twin_axis
        self.dfs[df_nr]['Twin_angle'] = twin_angle
        self.dfs[df_nr]['df'] = df               
        
    def clear_plot(self):
        self.dfs = {}
        self.PF_plot()

    def calc_alpha_beta(self, plane, twin, twin_axis, twin_angle, rotation, 
                        rotation_axis, rotation_angle):
        
        plane_vector = self.indices_to_vectors(plane)
        if twin:
            twin_vector = self.indices_to_vectors(twin_axis)
            plane_vector = self.rotate_vector(plane_vector, twin_vector, 
                                              twin_angle)
            
        if rotation:
            rotation_vector = self.indices_to_vectors(rotation_axis)
            plane_vector = self.rotate_vector(plane_vector, rotation_vector, 
                                              rotation_angle)
            
        alpha = self.angle_between_vectors(self.z_vector, plane_vector)
        proj_v = self.proj(plane_vector, self.z_vector)
        # To see on which side of the Y plane is.
        v1dot = np.dot(proj_v, self.y_vector) 
        beta = self.angle_between_vectors(proj_v, self.x_vector)
        
        if v1dot < 0:
            beta = 360 - beta
        
        if alpha > 90: #Not on the visible side, appears on the opposite side.
            beta += 180 
            alpha = 180 - alpha
            
            if beta > 360:
                beta -= 360
                
        if alpha == 0:
            beta = 0
        
        alpha = np.round(alpha, 2)
        beta = np.round(beta, 2)
            
        return alpha, beta
    
    def proj(self, vector, plane_vector):
        v = np.array(vector)
        plane = np.array(plane_vector)
        
        proj = (v - np.dot(v, plane)/np.dot(plane, plane) * plane)
        return proj
    
    def add_pole(self, ref, twin_axis, twin_angle):
        if twin_axis == None:
            twin = False
        else:
            twin = True
        self.planes_in_family(ref, twin=twin, twin_axis=twin_axis, twin_angle=twin_angle)
        self.PF_plot()

    def PF_plot(self, ref=None, sim=True, stereographic=False, scale='log'):
        '''
        Parameters
        ----------
        alpha : TYPE
            DESCRIPTION.
        beta : TYPE
            DESCRIPTION.
        sim : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        fig : TYPE
            DESCRIPTION.
        ax : TYPE
            DESCRIPTION.

        '''
        
        if ref != None:
            self.planes_in_family(ref)
        
        xlim = -1.4
        xlim2 = xlim - 0.24
        ylim = xlim + 0.05
        number_pf_plots = 1
        for key in self.dfs:
            df = self.dfs[key]['df']
            color = self.markers[str(number_pf_plots)]['color']
            marker = self.markers[str(number_pf_plots)]['marker']
            size = self.markers[str(number_pf_plots)]['size']
            number_pf_plots += 1
            b_dirmuth = np.deg2rad(df.Beta.to_numpy())
            dip = df.Alpha.to_numpy()
            # Rigaku data is given reversed
            reverse = not sim
            if reverse:
                x = (90 - dip)/90 * np.cos(b_dirmuth)
                y = (90 - dip)/90 * np.sin(b_dirmuth)
            else:
                x = (dip)/90 * np.cos(b_dirmuth)
                y = (dip)/90 * np.sin(b_dirmuth)
            if sim:
                a = self.ax.scatter(x, y, c=color,marker=marker, s=size)
            
        self.add_artists(self.ax, xlim, ylim, xlim2)

    def rotate_plot(self, rotation_axis, rotation_angle):
        self.rotation = True
        self.rotation_axis = rotation_axis
        self.rotation_angle = rotation_angle
        for key in self.dfs:
            alphas = []
            betas = []
            for ref in self.dfs[key]['df']['Planes']:
                alpha, beta = self.calc_alpha_beta(ref, self.dfs[key]['Twin'], self.dfs[key]['Twin_axis'], self.dfs[key]['Twin_angle'], True,  rotation_axis, rotation_angle)
                alphas.append(alpha)
                betas.append(beta)
            self.dfs[key]['df']['Alpha'] = alphas
            self.dfs[key]['df']['Beta'] = betas
        self.PF_plot()
        
    def set_plot(self):

        fig, ax = plt.subplots(dpi=600)

        ax.tick_params(
            axis='both',
            which='both',
            bottom=False,
            top=False,
            left=False,
            labelbottom=False,
            labelleft=False)
        return fig, ax
    
    def add_artists(self, ax, xlim, ylim, xlim2, stereographic=False):
        circ = plt.Circle((0, 0), 1.0, facecolor='none', edgecolor='black')
        self.ax.add_patch(circ)

        self.ax.axis('equal')  # equal aspect ratio
        self.ax.axis('off')  # remove the box
        self.ax.add_artist(lines.Line2D([xlim + 0.05, xlim + 0.05], [-1, 1], color='black', linewidth=1))
        self.ax.set_xlim(-1.5, 1.3)
        self.ax.set_ylim(-1.3, 1.3)
        ax.add_artist(lines.Line2D([xlim, ylim], [-1, -1], color='black', linewidth=1))
        ax.add_artist(lines.Line2D([xlim, ylim], [1, 1], color='black', linewidth=1))
        ax.add_artist(lines.Line2D([xlim + 0.03, ylim], [-0.33, -0.33], color='black', linewidth=1))
        ax.add_artist(lines.Line2D([xlim + 0.03, ylim], [0.33, 0.33], color='black', linewidth=1))
        ax.add_artist(lines.Line2D([xlim + 0.03, ylim], [-0.66, -0.66], color='black', linewidth=1))
        ax.add_artist(lines.Line2D([xlim + 0.03, ylim], [0.66, 0.66], color='black', linewidth=1))
        ax.add_artist(lines.Line2D([xlim, ylim], [0, 0], color='black', linewidth=1))

        ax.text(xlim2, -1.03, '90$^\circ$')
        ax.text(xlim2, 0.97, '90$^\circ$')
        ax.text(xlim2 + 0.05, -0.03, '0$^\circ$')

        ax.text(xlim2, -0.36, '30$^\circ$')
        ax.text(xlim2, 0.30, '30$^\circ$')
        ax.text(xlim2, -0.69, '60$^\circ$')
        ax.text(xlim2, 0.63, '60$^\circ$')

        if stereographic:
            r = np.tan((np.pi / 4.0) - (np.deg2rad(60) / 2)) * np.sin(90)
        else:
            r = 30/90

        for x in range(0, 360, 45):
            xcoords = (r*np.sin(np.deg2rad(x)), np.sin(np.deg2rad(x)))
            ycoords = (r*np.cos(np.deg2rad(x)), np.cos(np.deg2rad(x)))
            if x == 0:
                txtcoords = (1.13*np.cos(np.deg2rad(x)), 1.13*np.sin(np.deg2rad(x)))
            else:
                txtcoords = (1.16*np.cos(np.deg2rad(x)), 1.16*np.sin(np.deg2rad(x)))
            ax.text(txtcoords[0], txtcoords[1], str(x) + '$^\circ$', ha='center', va='center')

            ax.add_artist(lines.Line2D(xcoords, ycoords, color='black', linewidth=0.25))

        for x in range(0, 90, 30):
            if stereographic:
                dip = np.deg2rad(90 - x)
                ycir = np.tan((np.pi / 4.0) - (dip / 2)) * np.cos(0)
                circ = plt.Circle((0, 0), ycir, facecolor='none', edgecolor='black', linewidth=0.25, linestyle='dashed')
            else:
                circ = plt.Circle((0, 0), x/90, facecolor='none', edgecolor='black', linewidth=0.25, linestyle='dashed')
            ax.add_patch(circ)
        return
    
    def save_fig(self, name):
        saves_dir = os.getcwd() + '/SavedFigures' 
        if not os.path.exists(saves_dir):
            os.makedirs(saves_dir)
        self.fig.savefig(saves_dir + '/' + name)
        print('Saved figure: ' + saves_dir + '/' + name + '.png')
    