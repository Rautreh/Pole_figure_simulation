import os
from collections import Counter
from itertools import product, permutations

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import lines
import pandas as pd

from database import database

class Crystal:
    def __init__(self, mat='GaN', system='cub', z=[0,0,1], x = [1,0,0], aproximation='x'):
        self.mat = mat
        self.system = system
        self.z_indices = np.array([z])
        self.x_indices = np.array([x])
        self.aproximation = aproximation
        self.rotation = False
        self.rotation_axis = None
        self.rotation_angle = None
        size = 35
        self.markers = {'1': {'color': 'r', 'marker': 'o', 'size': size}, '2':{'color': 'b', 'marker': 's', 'size': size}, 
                        '3':{'color':'m', 'marker': '^', 'size': size}, '4':{'color':'g', 'marker': 'v', 'size': size},
                        '5': {'color': 'r', 'marker': 'o', 'size': size}, '6':{'color': 'b', 'marker': 's', 'size': size}, 
                        '7':{'color':'m', 'marker': '^', 'size': size}, '8':{'color':'g', 'marker': 'v', 'size': size},
                        '9':{'color':'m', 'marker': '^', 'size': size}, '10':{'color':'g', 'marker': 'v', 'size': size},
                        '11': {'color': 'r', 'marker': 'o', 'size': size}, '12':{'color': 'b', 'marker': 's', 'size': size}, 
                        '13':{'color':'m', 'marker': '^', 'size': size}, '14':{'color':'g', 'marker': 'v', 'size': size},
                        '15': {'color': 'r', 'marker': 'o', 'size': size}, '16':{'color': 'b', 'marker': 's', 'size': size}, 
                        '17':{'color':'m', 'marker': '^', 'size': size}, '18':{'color':'g', 'marker': 'v', 'size': size},
                        '19':{'color':'m', 'marker': '^', 'size': size}, '20':{'color':'g', 'marker': 'v', 'size': size}}
        self.dfs = {}
        self.get_database()
        self.lattice_parameters()
        self.orient()
        self.fig, self.ax = self.set_plot()
        
        
    def get_database(self):
        self.data_base_latice_paramaeters = database
        
        
    def lattice_parameters(self):
        self.a = np.array([[self.data_base_latice_paramaeters[self.mat][self.system]['a'],
                  self.data_base_latice_paramaeters[self.mat][self.system]['b'],
                  self.data_base_latice_paramaeters[self.mat][self.system]['c']]])
        
        self.alpha = np.array([self.data_base_latice_paramaeters[self.mat][self.system]['alpha'],
                      self.data_base_latice_paramaeters[self.mat][self.system]['beta'],
                      self.data_base_latice_paramaeters[self.mat][self.system]['gamma']])
       
        self.alpha_rad = np.deg2rad(self.alpha)
        cos = np.cos(self.alpha_rad)
        sin = np.sin(self.alpha_rad)
        
        if (cos[1]*cos[2]-cos[0])/(sin[1]*sin[2]) > 1 or (cos[2]*cos[0]-cos[1])/(sin[2]*sin[0]) > 1 or (cos[0]*cos[1]-cos[2])/(sin[0]*sin[1]) > 1:
            raise ValueError('Please make sure that:0 < (cos(alphaⱼ)*cos(alphaₖ)-cos(alphaᵢ))/(sin(alphaⱼ)*sin(alphaₖ)) < pi for i =/= j =/= k, with alpha and beta the real and reciprocal angles')
        
        self.beta_rad = np.arccos([(cos[1]*cos[2]-cos[0])/(sin[1]*sin[2]),
                                   (cos[2]*cos[0]-cos[1])/(sin[2]*sin[0]),
                                   (cos[0]*cos[1]-cos[2])/(sin[0]*sin[1])])
        
        self.beta = np.rad2deg(self.beta_rad)
        #in form such that matrix multiplication is correct
        print(self.beta_rad)
        self.a_vectors = self.a * np.array([[1, cos[2],   cos[1]],
                                            [0, sin[2],   -sin[1]*np.cos(self.beta_rad[0])],
                                            [0, 0,        sin[1]*np.sin(self.beta_rad[0])]])
        
        self.a_vectors_inv = np.linalg.inv(self.a_vectors)
        
        #transposed so that each element of a_vectors is a real vector a, b and c
        a_vectors = np.transpose(self.a_vectors)
        Volume = np.linalg.det(self.a_vectors)
        #transposed back to matrix multiply with miller indices
        self.b_vectors = np.transpose(np.array([np.cross(a_vectors[1], a_vectors[2]),
                                   np.cross(a_vectors[2], a_vectors[0]),
                                   np.cross(a_vectors[0], a_vectors[1])])/Volume)
        
        self.b_inverse = np.linalg.inv(self.b_vectors)
        print(self.b_vectors)
        
    def indices_to_vectors(self, indices):
        return np.sum(indices * self.b_vectors, axis=1)
    
    def vector_to_indices(self, vector):
        return np.sum(vector * self.b_inverse, axis=1)
    
    def interplanar_distance(self, indices):
        H = np.sum(indices * self.b_vectors, axis=1)
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
        
        if self.aproximation != 'z':
            self.x_vector = self.rotate_vector(self.x_vector, self.y_vector, 90 - angle_zx )
        else:
            self.z_vector = self.rotate_vector(self.z_vector, self.y_vector, angle_zx - 90 )
        
    def angle_between_vectors(self, v1, v2):
        a = np.dot(v1, v2)/(np.linalg.norm(v1) * np.linalg.norm(v2))
        a_angle = np.rad2deg(np.arccos(np.round(a, 5)))
        return a_angle

    def rotate_vector(self, v1, v2, ang):
        v1 = v1/np.linalg.norm(v1)
        v2 = v2/np.linalg.norm(v2)
        angle = np.deg2rad(ang)
        if angle == 0:
            return v1

        v2 = np.sin(angle/2)*v2

        q = [np.cos(angle/2), v2[0], v2[1], v2[2]]

        q_1 = [np.cos(angle/2), -v2[0], -v2[1], -v2[2]]
        v = [0, v1[0], v1[1], v1[2]]
        Qr = self.quaternion_multiply(q, v)

        Qr = self.quaternion_multiply(Qr, q_1)

        return np.array(Qr[1:])
    
    def quaternion_multiply(self, Q0, Q1):
        """
        Multiplies two quaternions.

        Input
        :param Q0: A 4 element array containing the first quaternion (q01,q11,q21,q31) 
        :param Q1: A 4 element array containing the second quaternion (q02,q12,q22,q32) 

        Output
        :return: A 4 element array containing the final quaternion (q03,q13,q23,q33) 

        """
        # Extract the values from Q0
        w0 = Q0[0]
        x0 = Q0[1]
        y0 = Q0[2]
        z0 = Q0[3]

        # Extract the values from Q1
        w1 = Q1[0]
        x1 = Q1[1]
        y1 = Q1[2]
        z1 = Q1[3]

        # Computer the product of the two quaternions, term by term
        Q0Q1_w = w0 * w1 - x0 * x1 - y0 * y1 - z0 * z1
        Q0Q1_x = w0 * x1 + x0 * w1 + y0 * z1 - z0 * y1
        Q0Q1_y = w0 * y1 - x0 * z1 + y0 * w1 + z0 * x1
        Q0Q1_z = w0 * z1 + x0 * y1 - y0 * x1 + z0 * w1

        # Create a 4 element array containing the final quaternion
        final_quaternion = np.array([Q0Q1_w, Q0Q1_x, Q0Q1_y, Q0Q1_z])

        # Return a 4 element array containing the final quaternion (q02,q12,q22,q32)

        return final_quaternion
    def calculate_PF_data(ref):
        pass
    
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
                if self.dfs[key]['Twin_axis'] == twin_axis and self.dfs[key]['Twin_angle'] == twin_angle:
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
            perm = set([*Counter(permutations(refh[0:3]))] + [*Counter(permutations(-refh[0:3]))])
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
            indices = [*Counter(ref)]  # Takes indices of plane
            indices.extend(-np.array(indices))  # Adds its negative indices
            plane_combination = [*Counter([*product(indices, repeat=3)])]  # Obtains all possible combinations
            for plane in plane_combination:
                d_hkl_pl = self.interplanar_distance(plane)
                if abs((d_hkl-d_hkl_pl)/d_hkl) < 0.001: # If interplanar distance is the same, belongs to the same family
                    planes_in_family.append(plane)
        
        for plane in planes_in_family:
            alpha, beta = self.calc_alpha_beta(plane, twin, twin_axis, twin_angle, self.rotation, self.rotation_axis, self.rotation_angle)
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

    def calc_alpha_beta(self, plane, twin, twin_axis, twin_angle, rotation, rotation_axis, rotation_angle):
        plane_vector = self.indices_to_vectors(plane)
        if twin:
            twin_vector = self.indices_to_vectors(twin_axis)
            plane_vector = self.rotate_vector(plane_vector, twin_vector, twin_angle)
            
        if rotation:
            rotation_vector = self.indices_to_vectors(rotation_axis)
            plane_vector = self.rotate_vector(plane_vector, rotation_vector, rotation_angle)
            
        alpha = self.angle_between_vectors(self.z_vector, plane_vector)
        proj_v = self.proj(plane_vector, self.z_vector)
        
        v1dot = np.dot(proj_v, self.y_vector) # To see on which side of the Y plane is.
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
        plane_vector = np.array(plane_vector)
        proj_v_on_plane = v - np.dot(v, plane_vector)/np.dot(plane_vector, plane_vector) * plane_vector
        return proj_v_on_plane
    
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
            
            # under construction:
            # else:
            #     if scale == 'log':
            #         a = self.ax.tricontourf(x, y, df.Intensity, cmap=cmaps.precip3_16lev,
            #                            c_dir=colors.Logc_dir(vmin=df.Intensity.min(), vmax=df.Intensity.max()))
            #     else:
            #         a = self.ax.tricontourf(x, y, df.Intensity, cmap=cmaps.precip3_16lev)
    
            
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
    