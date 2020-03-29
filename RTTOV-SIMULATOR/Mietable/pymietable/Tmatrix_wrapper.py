# -*- coding: utf-8 -*-

'''
@Description: Tmatrix_wrapper.py: A wrapper to read, processing and plot the cascading
directory database produced by II-Tmatrix
@Author: Hejun Xie
@Date: 2020-03-24 21:35:41
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-29 11:14:18
'''

# global import
import os
import numpy as np
import sys
import copy

# local import
import pymietable.coatratio
from pymietable.utils import readtable, float_index, get_dim

# some global conatnsts
C         = 299792458.  # [m/s]

class OptNode():
    'class to handle a folder containing optical properties produced by II-Tmatrix'

    def __init__(self, folder, num_sca_angles, pmtype = 0, isctype = 0, random_orientation = True,
                spherecoat=False, passive=False, melting=True, **kwargs):
        
        self.folder = folder
        self.num_sca_angles = num_sca_angles
        self.random_orientation = random_orientation
        self.passive = passive
        self.spherecoat = spherecoat
        self.axis = kwargs

        self.integrity = True

        # get alias for axes
        if 'frequency' in self.axis.keys() and self.axis['frequency'] is not None:
            self.axis['lambda'] =  C / self.axis['frequency'] * 1e-9 * 1e+3 # [mm]
        elif 'lambda' in self.axis.keys() and self.axis['lambda'] is not None:
            self.axis['frequency'] = C / self.axis['lambda'] * 1e-9 * 1e+3 # [GHz] 

        if 'xsize' in self.axis.keys() and self.axis['xsize'] is not None:
            self.axis['Dmax'] = self.axis['xsize'] * (self.axis['lambda']  / (2 * np.pi)) # [mm]
        elif 'Dmax' in self.axis.keys() and self.axis['Dmax'] is not None:
            self.axis['xsize'] = self.axis['Dmax'] * ((2 * np.pi) / self.axis['lambda']) # [-]

        if melting and 'coat_ratio' in self.axis.keys() and self.axis['coat_ratio'] is not None:
            self.axis['wc'] = coatratio.get_water_fraction(self.axis['coat_ratio'],
            spherecoat=self.spherecoat, AR=self.axis['aspect_ratio'])    

        # print(self.axis)

        self.file_phase_matrix = os.path.join(folder, 'phase_matrix.dat')
        self.file_isca = os.path.join(folder, 'isca.dat')
        self.file_par = os.path.join(folder, 'par.dat')
        
        self.sca_angles = np.linspace(0, 180, self.num_sca_angles)
        
        self.phase_matrix = None

        if not self.passive:
            self.mueller_matrix_b = None
            self.mueller_matrix_f = None
            self.jones_matrix_f = None 

        # phase matrix
        try:
            self.read_phase_matrix(pmtype)
            self.read_isca(isctype)
            self.read_par()
        except ValueError:
            self.integrity = False
            return

        if not self.passive:
            self.set_mueller()
            self.set_jones()

    def read_phase_matrix(self, pmtype):
        
        with open(self.file_phase_matrix, 'r') as fin:
            if pmtype == 0:
                NUM_PHASE_ELE = 10
            elif pmtype == 1:
                NUM_PHASE_ELE = 6
            
            N = (NUM_PHASE_ELE + 1) * self.num_sca_angles
            phase_matrix = readtable(fin, NUM_PHASE_ELE + 1, N, datatype='float64')
            phase_matrix = phase_matrix.reshape((self.num_sca_angles, NUM_PHASE_ELE + 1))
        
        self.phase_matrix = np.zeros((4, 4, self.num_sca_angles), dtype='float64')

        # put the data in right position in phase matrix
        self.phase_matrix[0, 0, :] = phase_matrix[:, 1]
        self.phase_matrix[0, 1, :] = phase_matrix[:, 2]
        self.phase_matrix[1, 1, :] = phase_matrix[:, 3]
        self.phase_matrix[2, 2, :] = phase_matrix[:, 4]
        self.phase_matrix[3, 3, :] = phase_matrix[:, 5]
        self.phase_matrix[2, 3, :] =  - phase_matrix[:, 6]
        if pmtype == 0:
            self.phase_matrix[0, 2, :] =  phase_matrix[:, 7]
            self.phase_matrix[0, 3, :] =  phase_matrix[:, 8]
            self.phase_matrix[1, 2, :] =  phase_matrix[:, 9]
            self.phase_matrix[1, 3, :] =  phase_matrix[:, 10]

        if self.random_orientation == True:
            self.phase_matrix[1, 0, :] = self.phase_matrix[0, 1, :]
            self.phase_matrix[3, 2, :] = - self.phase_matrix[2, 3, :]

    def read_isca(self, isctype):
        with open(self.file_isca, 'r') as fin:
            if isctype == 0:
                self.Cext = float((fin.readline().split())[-1])
                self.Csca = float((fin.readline().split())[-1])

                self.Csca = self.Csca * (self.axis['lambda'] / (2 * np.pi)) ** 2 # [mm^2]
                self.Cext = self.Cext * (self.axis['lambda'] / (2 * np.pi)) ** 2 # [mm^2]

            elif isctype == 1:
                self.Qext = float((fin.readline().split())[-1])
                self.Qsca = float((fin.readline().split())[-1])
                self.area = float((fin.readline().split())[-1])
                self.g = float((fin.readline().split())[-1])

                self.Cext = self.Qext * self.area
                self.Csca = self.Qsca * self.area
    
    def read_par(self):
        with open(self.file_par, 'r') as fin:
            fin.readline()
            self.Ssec = float(fin.readline().split()[1])
            
        self.Ssec = self.Ssec * (self.axis['lambda'] / (2 * np.pi)) ** 2
        self.Veq = np.power(self.Ssec / np.pi, 1.5) * (4*np.pi/3) 

    def set_mueller(self):
        self.mueller_matrix_b = (self.phase_matrix * self.Csca / (4 * np.pi)) [:, :, -1]
        self.mueller_matrix_f = (self.phase_matrix * self.Csca / (4 * np.pi)) [:, :, 0]
    
    def set_jones(self):
        self.jones_matrix_f = np.zeros((2, 2), dtype='complex')
        mu_f = self.mueller_matrix_f
        Shh_f_i = (2 * np.pi / self.axis['lambda']) * self.Cext / (4 * np.pi)
        Shh_f_sq = 0.5 * (mu_f[0, 0] + mu_f[1, 1] + mu_f[0, 1] + mu_f[1, 0])
        Shh_f_r = np.sqrt(Shh_f_sq - Shh_f_i ** 2)

        self.jones_matrix_f[0, 0] = Shh_f_r + Shh_f_i * 1j
        self.jones_matrix_f[1, 1] = Shh_f_r + Shh_f_i * 1j

    def get_sigma(self):
        mu_b = self.mueller_matrix_b 
        sigmah = 2 * np.pi * (mu_b[0, 0] + mu_b[1, 1] + mu_b[0, 1] + mu_b[1, 0])
        sigmav = 2 * np.pi * (mu_b[0, 0] + mu_b[1, 1] - mu_b[0, 1] - mu_b[1, 0])

        return sigmah, sigmav
    
    def get_att_factor(self):
        return self.jones_matrix_f[0, 0].imag / self.axis['Dmax']
    
    def get_phsft_factor(self):
        return self.jones_matrix_f[0, 0].real / self.axis['Dmax']
    
    def get_axis_sector(self):
        Ar = self.axis['aspect_ratio']
        axis_sector = 8 * Ar  / (4 + Ar ** 2)
        return axis_sector * (self.axis['lambda'] / (2 * np.pi)) ** 2
    
    def get_eqv_sector(self):
        if self.spherecoat:
            return self.axis['Dmax'] ** 2 * np.pi
        else:
            Ar = self.axis['aspect_ratio']
            x = 2 / np.sqrt(4 + Ar ** 2)
            y = x * Ar
            top = (3 * np.sqrt(3) / 2) * x ** 2
            side = x * y    
            eqv_sector = (2 * top + 6 * side) / 4

            return eqv_sector * (self.axis['lambda'] / (2 * np.pi)) ** 2 

    def get_sigma_factor(self):
        sigmah, sigmav = self.get_sigma()
        eqv_sector  = self.get_eqv_sector()
        return sigmah / eqv_sector, sigmav / eqv_sector

class OptDB():
    'class to handle a root folder containing a database of optical properties produced by II-Tmatrix'

    def __init__(self, root_folder, casca_ls, num_sca_angles, pmtype = 0, isctype = 0,random_orientation = True,
                spherecoat=False, dmnt_axis_ls=None, passive=False, melting=True, **kwargs):
        '''
        Input:
            root_folder: The root folder of the OptDB
            casca_ls: The cascading directory list (Outer --> Inner)
            pm_type: The phase matrix type generated by II-Tmatrix
            isc_type: The isca type generated by II-Tmatrix
            random_orientation: (bool)
            spherecoat: (bool) True if the water coat is sphere; False if the water coat is isomorphic with the core
            dmnt_axis: A list assigning the axis that should be substituted as domninator axis,
                       ex: [[casca_axis1, dmnt_axis1, dmnt_dim1], [casca_axis2, dmnt_axis2, dmnt_dim2], ...]
            kwargs: mutual dictionary for the axis of instances of OptNode in Nodes
        '''

        self.root_folder = root_folder
        self.root_segments = len(root_folder.split('/'))
        self.num_sca_angles = num_sca_angles
        self.sca_angles = np.linspace(0, 180, self.num_sca_angles)
        self.Node_dic_mutual = kwargs
        self.pmtype = pmtype
        self.isctype = isctype
        self.random_orientation = random_orientation
        self.melting = melting
        self.passive = passive
        self.spherecoat = spherecoat
        
        self.casca_ls = casca_ls
        self.casca_dim = {}
        for casca in casca_ls:
            self.casca_dim[casca] = None
        
        # make dominator axis
        self.dmnt_ls = copy.deepcopy(casca_ls)
        self.dmnt_dim = {}
        if dmnt_axis_ls is not None:
            for dmnt_axis in dmnt_axis_ls:
                index = self.dmnt_ls.index(dmnt_axis[0])
                self.dmnt_ls[index] = dmnt_axis[1]
                self.dmnt_dim[dmnt_axis[1]] = dmnt_axis[2]

        for dmnt in self.dmnt_ls:
            if dmnt not in self.dmnt_dim.keys():
                self.dmnt_dim[dmnt] = None
        
        self.Nodes = list()

        self.num_segments = len(casca_ls) + self.root_segments

        self.read_DB()

        # optical property arrays
        self.op_ls = None
        self.op_dim = None
        if not self.passive:
            self.sigma_factor = None
            self.att_factor   = None
            self.phsft_factor = None
        else:
            self.Cext = None
            self.Csca = None
            self.g = None
        self.phase_matrix = None
        self.Veq = None

        # # some aux dimension variables, may be deprecated later
        self.Dmax = None
        if self.melting:
            self.wc = None
    
    def read_DB(self):
        '''
        Read the cascading directory

        Input: 
            casca_ls: The cascading directory list (Outer --> Inner)
            casca_dim: A dictionary describe the dimension of cascading directory
            Node_dic_mutual: The base of the dictionary for OptNode
            root_folder: The root folder of the OptDB
            root_segments: Number of sements of root_folder

            random_orientation: bool
            pmtype: The phase matrix type generated by II-matrix
        
        Ouput:
            Nodes: A list of OptNode instances
        '''

        self.casca_dim[self.casca_ls[0]] = get_dim(self.root_folder)
        if self.dmnt_dim[self.dmnt_ls[0]] is None:
            self.dmnt_dim[self.dmnt_ls[0]] = get_dim(self.root_folder)

        for root, dirs, files in os.walk(self.root_folder):
            for dir in dirs:
                folder = os.path.join(root, dir)

                # check if it is an terminal dir
                if len(folder.split('/')) != self.num_segments:
                    # get the dimension of cascading folder
                    if self.casca_dim[self.casca_ls[len(folder.split('/')) - self.root_segments]] is None:
                        self.casca_dim[self.casca_ls[len(folder.split('/')) - self.root_segments]] \
                            = get_dim(folder)
                    
                    if self.dmnt_dim[self.dmnt_ls[len(folder.split('/')) - self.root_segments]] is None:
                        self.dmnt_dim[self.dmnt_ls[len(folder.split('/')) - self.root_segments]] \
                            = get_dim(folder)

                    continue

                # generate Node_dic
                Node_dic = copy.deepcopy(self.Node_dic_mutual)
                folder_segments = folder.split('/')
                for i, casca in enumerate(self.casca_ls):
                    Node_dic[casca] = float(folder_segments[i + self.root_segments].split('_')[-1])

                # some sanity check for Node_dic: Invalid coat ratio
                if self.melting:
                    if Node_dic['coat_ratio'] is not None and not ( 0 <= Node_dic['coat_ratio'] <= 1):
                        continue

                Node = OptNode(folder, self.num_sca_angles, pmtype=self.pmtype, isctype=self.isctype, random_orientation=self.random_orientation,
                                spherecoat=self.spherecoat, passive=self.passive, melting=self.melting, **Node_dic)

                if Node.integrity:
                    self.Nodes.append(Node)

    def set_optical_property_array(self, op_ls):
        '''
        get optical properties from the database instance
        
        Input:
            op_ls: The cascading dimensions list to generate the optical property arrays
        
        Output:
            The optical property arrays in the instance of OptDB
            ex: sigma_factor, att_factor, phsft_factor, phase_matrix
        '''

        self.op_ls = op_ls
        self.op_dim = {}
        
        # some aux dimension variables, may be deprecated later
        if self.melting:
            if 'wc' not in self.dmnt_ls:
                self.wc = np.zeros(len(self.casca_dim['coat_ratio'])) if 'coat_ratio' in self.casca_dim else None
            else:
                self.wc = np.array(self.dmnt_dim['wc'])
        
        if 'xsize' in self.casca_dim:
            self.Dmax = np.zeros(len(self.casca_dim['xsize']))
        elif 'Dmax' in self.casca_dim:
            self.Dmax = np.asarray(self.casca_dim['Dmax'], dtype='float32')

        # make slices
        op_shape_ls = list()
        for op in self.op_ls:
            self.op_dim[op] = self.dmnt_dim[op]
            op_shape_ls.append(len(self.dmnt_dim[op]))
        
        op_shape_tuple = tuple(op_shape_ls)
        op_shape_ls.extend([4, 4, self.num_sca_angles])
        pm_shape_tuple = tuple(op_shape_ls)

        if not self.passive:
            self.sigma_factor = np.zeros(op_shape_tuple, dtype='float64') * np.nan
            self.att_factor = np.zeros(op_shape_tuple, dtype='float64') * np.nan
            self.phsft_factor = np.zeros(op_shape_tuple, dtype='float64') * np.nan
        else:
            self.Cext = np.zeros(op_shape_tuple, dtype='float64') * np.nan
            self.Csca = np.zeros(op_shape_tuple, dtype='float64') * np.nan
            self.g = np.zeros(op_shape_tuple, dtype='float64') * np.nan
        self.phase_matrix = np.zeros(pm_shape_tuple, dtype='float64') * np.nan
        self.Veq = np.zeros(op_shape_tuple, dtype='float64') * np.nan

        for Node in self.Nodes:
            
            slice_ls = list()
            for op in op_ls:
                index = float_index(self.op_dim[op], Node.axis[op])
                slice_ls.append(slice(index, index+1))
            slice_tuple = tuple(slice_ls)

            if not self.passive:
                self.sigma_factor[slice_tuple], _ = Node.get_sigma_factor()
                self.att_factor[slice_tuple] = Node.get_att_factor()
                self.phsft_factor[slice_tuple] = Node.get_phsft_factor()
            else:
                self.Cext[slice_tuple] = Node.Cext
                self.Csca[slice_tuple] = Node.Csca
                self.g[slice_tuple] = Node.g
            
            self.phase_matrix[slice_tuple] = Node.phase_matrix
            self.Veq[slice_tuple] = Node.Veq

            # some aux dimension variables, may be deprecated later
            if self.melting and 'coat_ratio' in self.casca_dim and 'wc' not in self.dmnt_ls:
                iCR = self.casca_dim['coat_ratio'].index(Node.axis['coat_ratio'])
                self.wc[iCR] = Node.axis['wc']
            
            if 'xsize' in self.casca_dim:
                isize = self.casca_dim['xsize'].index(Node.axis['xsize'])
                self.Dmax[isize] = Node.axis['Dmax']

# some unit test code
if __name__  == '__main__':
    a = [1.1, 1.2, 1.3]
    b = 1.20005
    print(float_index(a, b))
