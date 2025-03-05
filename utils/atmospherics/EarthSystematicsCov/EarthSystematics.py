from numpy import *
from pandas import DataFrame
import matplotlib.pyplot as plt
import sys

# Constants for unit conversion 
km_cm = 1e5
g_kg = 1e-3
cm2_m2 = 1e-4

class EarthSystematics:
    def __init__(self, filename):
        """
        A class to find correlations between Earth Model Systematics and help manage covariance objects for MaCh3.
        """
        
        self.model_filename = None
        self.yaml_filename= None
        self.NLayers = None
        self.ParNames = []
        self.radii = None
        self.layer_widths = None
        self.rhos = None
        self.density_weights = None
        self.radii_sigma = None
        self.density_weights_sigma = None
        self.set_of_radii_pulls = []
        self.set_of_weight_pulls = []
        self.N_accepted_pulls = 0
        self.df = None
        self.corr = None
        self.cov = None
        self.EarthM_Tolerance = 1e-3
        self.EarthI_Tolerance = 1e-3
        self.isPoly = False
        self.PolyCoeffs_a = None
        self.PolyCoeffs_b = None
        self.PolyCoeffs_c = None
        self.EarthMs = []
        self.EarthIs = []

        # Earth's parameters
        
        # Ref: H. Moritz, Geodetic Reference System 1980, Journal of Geodesy, vol. 74, pp. 128–133, 2000
        # Earth Radius
        self.R = 6371 #km

        # Ref: A. M. Dziewonski and D. L. Anderson, Preliminary Reference Earth Model, Phys. Earth Planet. Interiors, vol. 25, pp. 297–356, 1981
        # Earth Mass 
        self.M = 5.9722*1e24 #kg
        # Earth Moment of Inertia (I/MR^2)
        self.I = 0.3308

        print("\n####   Initializing EarthSystematics   ####\n")
        print("Retrieving information from the file: ", filename)

        self.model_filename = filename

        File = loadtxt(filename)
        self.radii = File[:, 0]
        self.layer_widths = self.radii[1:] - self.radii[:-1]
                
        # Check which kind of table we have
        if(File.shape[-1]==3): 
            self.isPoly=False
            self.rhos = File[:, 1]
            if self.radii[0]==0:
                self.radii = self.radii[1:]
                self.rhos = self.rhos[1:]
            self.NLayers = len(self.radii)
            print('Found %d layers in the format [radius, rho].\n'%(self.NLayers))
        elif(File.shape[-1]==5): 
            self.isPoly=True
            self.PolyCoeffs_a = File[:, 1]
            self.PolyCoeffs_b = File[:, 2]
            self.PolyCoeffs_c = File[:, 3]
            if self.radii[0]==0:
                self.radii = self.radii[1:]
                self.PolyCoeffs_a = self.PolyCoeffs_a[1:]
                self.PolyCoeffs_b = self.PolyCoeffs_b[1:]
                self.PolyCoeffs_c = self.PolyCoeffs_c[1:]
            self.NLayers = len(self.radii)
            print('Found %d layers in the format [radius, a, b, c].\n'%(self.NLayers))
        else: raise ValueError("Unrecognized Earth model table format: %s."%filename)

        self.density_weights = ones(self.NLayers)

        # Check if provided model is sane
        M = self.calculate_M()
        I = self.calculate_I()/M/(self.radii[-1]*km_cm)**2
        if abs(M-self.M)/self.M > self.EarthM_Tolerance or abs(I - self.I) > self.EarthI_Tolerance:
            print("M calculated from %s: %.3e.\t\t Expected: %.3e"%(filename, M, self.M))
            print("I/MR^2 calculated from %s: %.3f.\t Expected: %.3f\n"%(filename, I, self.I))
            raise ValueError("Please make sure the input model is compatible with the Earth's mass and moment of inertia.")
        
        for iPar in range(self.NLayers):
            self.ParNames.append("r_%d"%(iPar+1))

        for iPar in range(self.NLayers):
            self.ParNames.append("w_%d"%(iPar+1))
        
        print("Loaded model with %d parameters. Please set the uncertainties."%(self.NLayers *2))


    def set_N_accepted_pulls(self, N):
        """
        Sets the number of accepted pulls (not the actual number of pulls!)

        Args:
            N (int): the number of desired accepted pulls
        """
        self.N_accepted_pulls = int(N)

    def set_sigma_r_array(self, sigma_arr):
        """
        Sets the radii uncertainties (in km) from an array.

        Args:
            sigma_arr (array(float)): array with radii sigmas
        """
        if len(sigma_arr)!=self.NLayers:
            raise ValueError("Incompatible size in set_sigma_r_array(). Expected: %d. Got: %d"%(self.NLayers, len(sigma_arr)))
        
        self.radii_sigma = sigma_arr

    def set_sigma_w_array(self, sigma_arr):
        """
        Sets the denisity weights uncertainties from an array.

        Args:
            sigma_arr (array(float)): array with weights sigmas
        """
        if len(sigma_arr)!=self.NLayers:
            raise ValueError("Incompatible size in set_weights_sigma(). Expected: %d. Got: %d"%(self.NLayers, len(sigma_arr)))
        self.density_weights_sigma = sigma_arr

    def set_sigma(self, ParName, sigma):
        """
        Sets a specific uncertainty. Radii must be given in km.

        Args:
            ParName (string): name of the parameter whose uncertainty is to be set.
            sigma (float): uncertainty.
        """
        i=-1
        for j, p in enumerate(self.ParNames): 
            if p == ParName: 
                i = j
                break
        if i==-1:
            raise ValueError("Provided parameter name not found. Please make sure it is in format 'r_i' for radii and 'w_i' for weights, with i={1, ..., %d}"%self.NLayers)
        
        if 'r' in ParName: self.radii_sigma[i] = sigma
        elif 'w' in ParName: self.density_weights_sigma[i-self.NLayers] = sigma

    def set_sigmas_uniform(self, r, w):
        """
        Sets a uniform uncertainty for the radii and the weights using a fraction of the layer widths and the densities.

        Args:
            r (float): weights uncertainty (0<r<1).
            w (float): weights uncertainty (0<w<1).
        """
        # Set density weights sigmas
        self.density_weights_sigma = w*ones(self.NLayers)
        
        # Set radii sigmas
        self.radii_sigma = r*self.layer_widths
        self.radii_sigma[-1] = 10

    def print_systematics(self):
        """
        Prints out the systematics.
        """
        
        print("\n------------------------------------")
        for iPar in range(self.NLayers):
            print("%s = %.1f +- %.1f"%(self.ParNames[iPar], self.radii[iPar], self.radii_sigma[iPar]))
        print("------------------------------------")
        for iPar in range(self.NLayers):
            print("%s = %.2f +- %.2f"%(self.ParNames[iPar+self.NLayers], self.density_weights[iPar], self.density_weights_sigma[iPar]))
        print("------------------------------------\n")

    def set_M_tolerance(self, t):
        """
        Sets the allowed tolerance for the error in the calculated Earth's mass.

        Args:
            t (float): relative error allowed in the Earth's mass
        """
        self.EarthM_Tolerance = t

    def set_I_tolerance(self, t):
        """
        Sets the allowed tolerance for the error in the calculated Earth's moment of inertia.

        Args:
            t (float): error allowed in the Earth's I/MR^2
        """
        self.EarthM_Tolerance = t

    def make_pulls(self, ProgressBar=True):
        """
        Generates new models from the loaded model and uncertainties.

        Args:
            ProgressBar (bool): turn on/off the progress bar.
        """
        
        print("Using mass tolerance of %.1e"%self.EarthM_Tolerance)
        print("Using moment of inertia tolerance of %.1e\n"%self.EarthI_Tolerance)

        set_of_radii_pulls = []
        set_of_weight_pulls = []

        print("Generating new models...")

        if self.N_accepted_pulls <=0: 
            raise ValueError("Please set the desired number of accepted pulls: set_N_accepted_pulls(int)")
        
        if ProgressBar==True:
            counter = 0
            div = 50
            while(len(set_of_radii_pulls)<self.N_accepted_pulls):
                tmp_radii = self.pull_radii()
                tmp_weights = self.pull_weights()
                if self.test_M_I(tmp_radii, tmp_weights)==True:
                    set_of_radii_pulls.append(array(tmp_radii))
                    set_of_weight_pulls.append(array(tmp_weights))
                    sys.stdout.write('\r')
                    sys.stdout.write("[%-50s] %.2f%%" % ('='*int((counter+1)/(self.N_accepted_pulls)*div), ((counter+1)*100/self.N_accepted_pulls)))
                    sys.stdout.flush()
                    counter+=1
            print('\n')
        else: 
            while(len(set_of_radii_pulls)<self.N_accepted_pulls):
                tmp_radii = self.pull_radii()
                tmp_weights = self.pull_weights()
                if self.test_M_I(tmp_radii, tmp_weights)==True:
                    set_of_radii_pulls.append(array(tmp_radii))
                    set_of_weight_pulls.append(array(tmp_weights))

        set_of_radii_pulls = array(set_of_radii_pulls)
        set_of_weight_pulls = array(set_of_weight_pulls)
        
        self.set_of_radii_pulls = set_of_radii_pulls
        self.set_of_weight_pulls = set_of_weight_pulls

    def save_pulls(self, savename_r='', savename_w=''):
        """
        Saves the generated accepted radii and weights pulls to a npy file.

        Args:
            savename_r (string): the name of the file for the radii pulls
            savename_w (string): the name of the file for the weight pulls
        """

        if savename_r=='': 
            savename_r='SetPulls_Radii_%d'%(self.N_accepted_pulls)
            print('savename not provided for radii pulls. Saving to %s.npy'%savename_r)
        if savename_w=='': 
            savename_w='SetPulls_Weights_%d'%(self.N_accepted_pulls)
            print('savename not provided for weights pulls. Saving to %s.npy'%savename_w)

        save('SavedPulls/%s'%savename_r, self.set_of_radii_pulls, allow_pickle=True)
        save('SavedPulls/%s'%savename_w, self.set_of_weight_pulls, allow_pickle=True)

    def load_pulls(self, filename_r='', filename_w='', N=-1):
        """
        Loads a previously generated accepted radii pulls from a npy file.

        Args:
            filename_r (string): the name of the file for the radii pulls
            filename_w (string): the name of the file for the weight pulls
            N (int): the number of pulls to be used. Defaults to use everything
        """

        if filename_r=='': 
            filename_r='SavedPulls/SetPulls_Radii_%d.npy'%(self.N_accepted_pulls)
            print('filename not provided for radii pulls. Loading %s.npy'%filename_r)
        if filename_w=='': 
            filename_w='SavedPulls/SetPulls_Weights_%d.npy'%(self.N_accepted_pulls)
            print('filename not provided for weights pulls. Loading %s.npy'%filename_w)
        
        r_pulls = load(filename_r)
        w_pulls = load(filename_w)

        if N==-1:
            self.N_accepted_pulls += len(w_pulls[:, 0])
            for i in range(self.N_accepted_pulls):
                self.set_of_radii_pulls.append(r_pulls[i, :])
                self.set_of_weight_pulls.append(w_pulls[i, :])
            
        else:
            N=int(N)
            self.N_accepted_pulls += N
            for i in range(self.N_accepted_pulls):
                self.set_of_radii_pulls.append(r_pulls[i, :])
                self.set_of_weight_pulls.append(w_pulls[i, :])
                
        self.set_of_radii_pulls = array(self.set_of_radii_pulls)
        self.set_of_weight_pulls = array(self.set_of_weight_pulls)
        self.NLayers = len(self.set_of_weight_pulls[0])

    def find_correlations(self):
        """
        Generates dataframe and calculates covariance and correlation matrices.
        """
        data_dict = {}

        for iPar in range(self.NLayers):
            namestring = 'r_%d'%(iPar+1)
            data_dict[namestring] = self.set_of_radii_pulls[:,iPar].tolist()

        for iPar in range(self.NLayers):
            namestring = 'w_%d'%(iPar+1)
            data_dict[namestring] = self.set_of_weight_pulls[:,iPar].tolist()

        df = DataFrame(data_dict)

        self.df = df
        self.corr = df.corr(method='pearson')
        self.cov = df.cov()
        print("\nCorrelations found!")
        if(self.is_posdef(self.corr.to_numpy())==False): print('\n***********************\nWarning: the correlation matrix is not positive definite!\n***********************\n')
        else: print("Correlation matrix was found to be positive definite.")

    def display_correlation_matrix(self):
        """
        Display the correlation matrix in a readable format.
        """
        print("\nCorrelation Matrix:")
        print(self.corr)
        print("\n")

    def get_correlation_matrix(self):
        """
        Returns the correlation matrix as a numpy array

        Returns:
            array(float) correlation matrix
        """
        return self.corr.to_numpy

    def display_covariance_matrix(self):
        """
        Display the covariance matrix in a readable format.
        """
        print("\nCovariance Matrix:")
        print(self.cov)
        print("\n")

    def add_osccov_to_yaml(self, yamlfile_source, yamlfile_dest=''):
        """
        Adds the generated systematics and correlations to a yaml file in the expected order.

        Args:
            yamlfile_source: yaml config file model from which the old objects will be loaded.
            yamlfile_dest: yaml config in which in which the old and new objects will be saved. If not provided, the new yaml will be saved to yamlfile_source + TAG referring to the model used in the systematics.
        """

        if yamlfile_dest=='':
            tag = self.model_filename[self.model_filename.rfind('/')+1:self.model_filename.rfind('.')]
            output_filename = yamlfile_source[:yamlfile_source.find('.yaml')]+'_%s.yaml'%tag
        else: output_filename=yamlfile_dest

        with open(yamlfile_source, 'r') as infile:
            file_content = infile.readlines()

        file_content.extend(["\n ############# Earth Model Systematics #############"])
        for iPar in range(2*self.NLayers):
            new_lines = self.make_covobj(iPar)
            print("Adding Covariance Object: ")
            for l in new_lines: print(l[:-1])
            file_content.extend(new_lines)

        with open(output_filename, 'w') as outfile:
            outfile.writelines(file_content)

        print(f"Successfully added the objects to the new file. File saved as: {output_filename}")

    def plot_M_I_dist(self, savename=''):
        """
        Plots the distributions of masses and moments of inertial from the generated pulls.

        Args:
            savename (string): name of the file to be saved
        """

        if self.EarthMs==[]:
            for iPull in range(len(self.set_of_weight_pulls)):
                weights_arr = self.set_of_weight_pulls[iPull]
                radii_arr = self.set_of_radii_pulls[iPull]
                if radii_arr[0]!=0: radii_arr = concatenate((array([0.]), radii_arr), axis=None)

                if self.isPoly==True:
                    radii3 = (radii_arr*km_cm)**3
                    radii4 = radii3*(radii_arr*km_cm)
                    radii5 = radii4*(radii_arr*km_cm)
                    radii6 = radii5*(radii_arr*km_cm)
                    radii7 = radii6*(radii_arr*km_cm)
                    M = 4*pi*(weights_arr*self.PolyCoeffs_a*1e-10*(radii5[1:]-radii5[:-1])/5 + weights_arr*self.PolyCoeffs_b*1e-5*(radii4[1:]-radii4[:-1])/4 + weights_arr*self.PolyCoeffs_c*(radii3[1:]-radii3[:-1])/3).sum()*g_kg
                    I = 8/3*pi*(weights_arr*self.PolyCoeffs_a*1e-10*(radii7[1:]-radii7[:-1])/7 + weights_arr*self.PolyCoeffs_b*1e-5*(radii6[1:]-radii6[:-1])/6 + weights_arr*self.PolyCoeffs_c*(radii5[1:]-radii5[:-1])/5).sum()*g_kg*g_kg*g_kg*cm2_m2*km_cm*km_cm/(M*radii_arr[-1]*radii_arr[-1]*km_cm*km_cm)
                    self.EarthMs.append(M)
                    self.EarthIs.append(I)
                else:
                    radii3 = (radii_arr*km_cm)**3
                    radii5 = radii3*radii_arr*radii_arr*km_cm*km_cm
                    M = (4/3*pi*(radii3[1:] - radii3[:-1])*self.rhos*weights_arr*g_kg).sum()
                    I = (4/5*pi*(radii5[1:] - radii5[:-1])*self.rhos*weights_arr*g_kg).sum()/(M*radii_arr[-1]*radii_arr[-1]*km_cm*km_cm)
                    self.EarthMs.append(M)
                    self.EarthIs.append(I)

        _, ax = plt.subplots(1, 2, figsize=(15, 6))

        h0 = ax[0].hist(self.EarthMs, bins=100)
        ax[0].set_title('Mass')
        ax[0].set_xlabel(r'$M$ [kg]')
        ax[0].set_ylabel('Counts')
        ax[0].vlines(self.M, 0, 1.1*h0[0].max(), color='red')

        h1 = ax[1].hist(self.EarthIs, bins=100)
        ax[1].set_title('Moment of Inertia')
        ax[1].set_xlabel(r'$I/MR^2$')
        ax[1].set_ylabel('Counts')
        ax[1].vlines(self.I, 0, 1.1*h1[0].max(), color='red')

        plt.suptitle('Accepted Models')
        plt.tight_layout()
        if savename=='': savename='Accepted_MI.png'
        plt.savefig(savename)
        print('Saved %s'%savename)

    def plot_2D_distribution(self, var1, var2, nbins=50, savename='', ax=None, cmap='Greys'):
        """
        Plots the 2D distribution of a pair of variables.

        Args:
            var1 (string): variable in the x-axis
            var2 (string): variable in the y-axis
            nbins (int): number of bins in each dimension
            ax (matplotlib.pyplot.axes): axis used to plot (optional)
            cmap: color map
            savename: save file name
        """
        if ax==None: 
            fig, ax = plt.subplots(figsize=(6.5, 5))
            StandAloneFig = True

        i_var1 = int(var1[var1.find('_')+1:])-1
        i_var2 = int(var2[var2.find('_')+1:])-1

        if 'r' in var1: 
            x = self.set_of_radii_pulls[:,i_var1]
            x0 = self.radii[i_var1]
            sigmax = self.radii_sigma[i_var1]
        elif 'w' in var1: 
            x = self.set_of_weight_pulls[:,i_var1]
            x0 = self.density_weights[i_var1]
            sigmax = self.density_weights_sigma[i_var1]

        if 'r' in var2: 
            y = self.set_of_radii_pulls[:,i_var2]
            y0 = self.radii[i_var2]
            sigmay = self.radii_sigma[i_var2]
        elif 'w' in var2: 
            y = self.set_of_weight_pulls[:,i_var2]
            y0 = self.density_weights[i_var2]
            sigmay = self.density_weights_sigma[i_var2]    

        histo = ax.hist2d(x, y,  bins=[linspace(x0-4*sigmax, x0+4*sigmax, nbins), linspace(y0-4*sigmay, y0+4*sigmay, nbins)], cmap=cmap)
        if StandAloneFig==True: fig.colorbar(histo[3], ax=ax)
        ax.plot(x0, y0, 'r*')
        ax.set_xlim(x0-4*sigmax, x0+4*sigmax)
        ax.set_ylim(y0-4*sigmay, y0+4*sigmay)
        ax.set_xlabel(r'$%s$'%var1)
        ax.set_ylabel(r'$%s$'%var2)

        if StandAloneFig==True:
            plt.tight_layout()
            if savename=='': savename='2D_Distributions_%s_%s.png'%(var1, var2)
            plt.savefig(savename)
            print('Saved %s'%savename)

    def plot_2D_weights(self, savename='', display_avg_rho=False):
        """
        Plots the 2D distributions of weights from the generated pulls.

        Args:
            savename (string): name of the file to be saved
            display_avg_rho (bool): weather or not to plot the weights as actual weights or in terms of the average layer densities
        """
        if display_avg_rho==True:
            avg_rhos = self.calculate_avg_rho()
            _, ax = plt.subplots(self.NLayers, self.NLayers, figsize=(11, 10))

            for i in range(self.NLayers-1):
                ax[i+1, 0].set_ylabel(r'$\langle \rho_%d \rangle$  [g/cm$^{-3}$]'%(i+1))
                for j in range(self.NLayers-1-i):
                    ax[i, j+1+i].set_axis_off()
            
            #Plot diagonal elements
            for i in range(self.NLayers):
                ax[i,i].hist(avg_rhos[i]*self.set_of_weight_pulls[:,i], 50, color='red')
                ax[i,i].set_ylabel('Counts')
                x0 = avg_rhos[i]*self.density_weights[i]
                sigmax = avg_rhos[i]*self.density_weights_sigma[i]
                ax[i,i].set_xlim(x0-4*sigmax, x0+4*sigmax)
                ax[self.NLayers-1, i].set_xlabel(r'$\langle \rho_%d \rangle$ [g/cm$^{-3}$]'%(i+1))

            # Plot off-diagonal elements
            for j in range(self.NLayers-1):
                for i in range(self.NLayers-1-j):
                    x0 = avg_rhos[j]*self.density_weights[j]
                    sigmax = avg_rhos[j]*self.density_weights_sigma[j]
                    y0 = avg_rhos[1+i+j]*self.density_weights[1+i+j]
                    sigmay = avg_rhos[1+i+j]*self.density_weights_sigma[1+i+j]
                    ax[1+i+j,j].hist2d(avg_rhos[j]*self.set_of_weight_pulls[:,j], avg_rhos[1+i+j]*self.set_of_weight_pulls[:,1+i+j],  bins=[linspace(x0-4*sigmax, x0+4*sigmax, 50), linspace(y0-4*sigmay, y0+4*sigmay, 50)], cmap='Reds')
                    ax[1+i+j,j].set_xlim(x0-4*sigmax, x0+4*sigmax)
                    ax[1+i+j,j].set_ylim(y0-4*sigmay, y0+4*sigmay)
                    ax[1+i+j,j].plot(avg_rhos[j]*self.density_weights[j], avg_rhos[1+i+j]*self.density_weights[1+i+j], 'g*')

            plt.tight_layout()
            if savename=='': savename='2D_Distributions_Densities.png'
            plt.savefig(savename)
            print('Saved %s'%savename)

        else:
            _, ax = plt.subplots(self.NLayers, self.NLayers, figsize=(11, 10))

            for i in range(self.NLayers-1):
                ax[i+1, 0].set_ylabel(r'$%s$'%(self.ParNames[self.NLayers+i+1]))
                for j in range(self.NLayers-1-i):
                    ax[i, j+1+i].set_axis_off()
            
            #Plot diagonal elements
            for i in range(self.NLayers):
                ax[i,i].hist(self.set_of_weight_pulls[:,i], 50, color='red')
                ax[i,i].set_ylabel('Counts')
                x0 = self.density_weights[i]
                sigmax = self.density_weights_sigma[i]
                ax[i,i].set_xlim(x0-4*sigmax, x0+4*sigmax)
                ax[self.NLayers-1, i].set_xlabel(r'$%s$'%(self.ParNames[self.NLayers+i]))

            # Plot off-diagonal elements
            for j in range(self.NLayers-1):
                for i in range(self.NLayers-1-j):
                    x0 = self.density_weights[j]
                    sigmax = self.density_weights_sigma[j]
                    y0 = self.density_weights[1+i+j]
                    sigmay = self.density_weights_sigma[1+i+j]
                    ax[1+i+j,j].hist2d(self.set_of_weight_pulls[:,j], self.set_of_weight_pulls[:,1+i+j],  bins=[linspace(x0-4*sigmax, x0+4*sigmax, 50), linspace(y0-4*sigmay, y0+4*sigmay, 50)], cmap='Reds')
                    ax[1+i+j,j].set_xlim(x0-4*sigmax, x0+4*sigmax)
                    ax[1+i+j,j].set_ylim(y0-4*sigmay, y0+4*sigmay)
                    ax[1+i+j,j].plot(self.density_weights[j], self.density_weights[1+i+j], 'g*')

            plt.tight_layout()
            if savename=='': savename='2D_Distributions_Weights.png'
            plt.savefig(savename)

    def plot_2D_radii(self, savename=''):
        """
        Plots the 2D distributions of radii from the generated pulls.

        Args:
            savename (string): name of the file to be saved
        """
        _, ax = plt.subplots(self.NLayers, self.NLayers, figsize=(11, 10))

        for i in range(self.NLayers-1):
            ax[i+1, 0].set_ylabel(r'$%s$ [km]'%(self.ParNames[i+1]))
            for j in range(self.NLayers-1-i):
                ax[i, j+1+i].set_axis_off()
        
        #Plot diagonal elements
        for i in range(self.NLayers):
            ax[i,i].hist(self.set_of_radii_pulls[:,i], 50, color='blue')
            ax[i,i].set_ylabel('Counts')
            x0 = self.radii[i]
            sigmax = self.radii_sigma[i]
            ax[i,i].set_xlim(x0-4*sigmax, x0+4*sigmax)
            ax[self.NLayers-1, i].set_xlabel(r'$%s$ [km]'%(self.ParNames[i]))

        # Plot off-diagonal elements
        for j in range(self.NLayers-1):
            for i in range(self.NLayers-1-j):
                x0 = self.radii[j]
                sigmax = self.radii_sigma[j]
                y0 = self.radii[1+i+j]
                sigmay = self.radii_sigma[1+i+j]
                ax[1+i+j,j].hist2d(self.set_of_radii_pulls[:,j], self.set_of_radii_pulls[:,1+i+j],  bins=[linspace(x0-4*sigmax, x0+4*sigmax, 50), linspace(y0-4*sigmay, y0+4*sigmay, 50)], cmap='Blues')
                ax[1+i+j,j].set_xlim(x0-4*sigmax, x0+4*sigmax)
                ax[1+i+j,j].set_ylim(y0-4*sigmay, y0+4*sigmay)
                ax[1+i+j,j].plot(self.radii[j], self.radii[1+i+j], 'r*')

        plt.tight_layout()
        if savename=='': savename='2D_Distributions_Radii.png'
        plt.savefig(savename)
        print('Saved %s'%savename)

    def plot_2D_all(self, savename='', display_avg_rho=False):
        """
        Plots the 2D distributions of radii x weights from the generated pulls.

        Args:
            savename (string): name of the file to be saved
            display_avg_rho (bool): weather or not to plot the weights as actual weights or in terms of the average layer densities
        """
        fig, ax = plt.subplots(self.NLayers, self.NLayers, figsize=(11, 10))

        if display_avg_rho==True:
            avg_rhos = self.calculate_avg_rho()
            for i in range(self.NLayers):
                ax[i, 0].set_ylabel(r'$\langle \rho_%d \rangle$  [g/cm$^{-3}$]'%(i+1))
                ax[self.NLayers-1,i].set_xlabel(r'$%s$ [km]'%(self.ParNames[i]))
                for j in range(self.NLayers):
                    x0 = self.radii[j]
                    sigmax = self.radii_sigma[j]
                    y0 = avg_rhos[i]*self.density_weights[i]
                    sigmay = avg_rhos[i]*self.density_weights_sigma[i]
                    ax[i,j].hist2d(self.set_of_radii_pulls[:,j], avg_rhos[i]*self.set_of_weight_pulls[:,i],  bins=[linspace(x0-4*sigmax, x0+4*sigmax, 50), linspace(y0-4*sigmay, y0+4*sigmay, 50)], cmap='Purples')
                    ax[i,j].set_xlim(x0-4*sigmax, x0+4*sigmax)
                    ax[i,j].set_ylim(y0-4*sigmay, y0+4*sigmay)
                    ax[i,j].plot(self.radii[j], avg_rhos[i]*self.density_weights[i], 'y*')
                    if savename=='': savename='2D_Distributions_All_rho.png'
        else:
            for i in range(self.NLayers):
                ax[i, 0].set_ylabel(r'$%s$'%(self.ParNames[self.NLayers+i]))
                ax[self.NLayers-1,i].set_xlabel(r'$%s$ [km]'%(self.ParNames[i]))
                for j in range(self.NLayers):
                    x0 = self.radii[j]
                    sigmax = self.radii_sigma[j]
                    y0 = self.density_weights[i]
                    sigmay = self.density_weights_sigma[i]
                    ax[i,j].hist2d(self.set_of_radii_pulls[:,j], self.set_of_weight_pulls[:,i],  bins=[linspace(x0-4*sigmax, x0+4*sigmax, 50), linspace(y0-4*sigmay, y0+4*sigmay, 50)], cmap='Purples')
                    ax[i,j].set_xlim(x0-4*sigmax, x0+4*sigmax)
                    ax[i,j].set_ylim(y0-4*sigmay, y0+4*sigmay)
                    ax[i,j].plot(self.radii[j], self.density_weights[i], 'y*')
                    if savename=='': savename='2D_Distributions_All.png'

        plt.tight_layout()
        plt.savefig(savename)
        print('Saved %s'%savename)

    def plot_correlation_matrix(self, savename=''):
        """
        Plots the correlation matrix.

        Args:
            savename (string): name of the file to be saved
        """
        fig, ax = plt.subplots(figsize=(6.5, 5))
        h = ax.imshow(self.corr.to_numpy(), vmin=-1, vmax=1, cmap='PiYG')
        fig.colorbar(h, ax=ax)
        ax.set_xticks(linspace(0, 2*self.NLayers-1, 2*self.NLayers), [r'$%s$'%par for par in self.ParNames])
        ax.set_yticks(linspace(0, 2*self.NLayers-1, 2*self.NLayers), [r'$%s$'%par for par in self.ParNames])
        ax.set_title('Correlation Matrix')
        plt.tight_layout()
        if savename=='': savename='CorrelationMatrix.png'
        plt.savefig(savename)
        print('Saved %s'%savename)

    def plot_covariance_matrix(self, savename=''):
        """
        Plots the correlation matrix.

        Args:
            savename (string): name of the file to be saved
        """
        fig, ax = plt.subplots(figsize=(6.5, 5))
        M = self.cov.to_numpy()
        vabs = max(M.max(), -M.min())
        h = ax.imshow(M, vmin=-vabs, vmax=vabs, cmap='PuOr')
        fig.colorbar(h, ax=ax)
        ax.set_xticks(linspace(0, 2*self.NLayers-1, 2*self.NLayers), [r'$%s$'%par for par in self.ParNames])
        ax.set_yticks(linspace(0, 2*self.NLayers-1, 2*self.NLayers), [r'$%s$'%par for par in self.ParNames])
        ax.set_title('Covariance Matrix')

        plt.tight_layout()
        if savename=='': savename='CovarianceMatrix.png'
        plt.savefig(savename)
        print('Saved %s'%savename)

    #############################
    # "Private" methods
        
    def calculate_avg_rho(self):
        #Calculates the average densities per layer.
        #Returns:
        #    array(float) with average density layers
 

        if self.isPoly==False: return self.rhos
        else:
            avg_rhos = []
            rad_lim = [0.] + self.radii.tolist()
            for i in range(self.NLayers):
                interval= self.layer_widths[i]
                integral_up = self.PolyCoeffs_a[i]*(rad_lim[i+1]**3)/3 + self.PolyCoeffs_b[i]*(rad_lim[i+1]**2)/2 + self.PolyCoeffs_c[i]*(rad_lim[i+1])
                integral_lo = self.PolyCoeffs_a[i]*(rad_lim[i]**3)/3 + self.PolyCoeffs_b[i]*(rad_lim[i]**2)/2 + self.PolyCoeffs_c[i]*(rad_lim[i])
                avg_rhos.append((integral_up - integral_lo)/interval)
            return avg_rhos

    def calculate_M(self):
        #Calculates the mass using the provided density model.
        
        if self.isPoly == True:
            if self.radii[0]!=0: radii_arr = concatenate((array([0.]), self.radii), axis=None)
            else: radii_arr = self.radii
            radii3 = (radii_arr*km_cm)**3
            radii4 = radii3*(radii_arr*km_cm)
            radii5 = radii4*(radii_arr*km_cm)
            return 4*pi*(self.PolyCoeffs_a*1e-10*(radii5[1:]-radii5[:-1])/5 + self.PolyCoeffs_b*1e-5*(radii4[1:]-radii4[:-1])/4 + self.PolyCoeffs_c*(radii3[1:]-radii3[:-1])/3).sum()*g_kg

        else:
            if self.radii[0]!=0: radii_arr = concatenate((array([0.]), self.radii), axis=None)
            else: radii_arr = self.radii
            radii3 = (radii_arr*km_cm)**3
            return (4/3*pi*(radii3[1:] - radii3[:-1])*self.rhos*g_kg).sum()

    def calculate_I(self):
        #Calculates the moment of inertia using the provided density model.
        
        if self.isPoly == True:
            if self.radii[0]!=0: radii_arr = concatenate((array([0.]), self.radii), axis=None)
            else: radii_arr = self.radii
            radii5 = (radii_arr*km_cm)**5
            radii6 = radii5*(radii_arr*km_cm)
            radii7 = radii6*(radii_arr*km_cm)
            return 8/3*pi*(self.PolyCoeffs_a*1e-10*(radii7[1:]-radii7[:-1])/7 + self.PolyCoeffs_b*1e-5*(radii6[1:]-radii6[:-1])/6 + self.PolyCoeffs_c*(radii5[1:]-radii5[:-1])/5).sum()*g_kg*g_kg*g_kg*cm2_m2*km_cm*km_cm
        else:
            if self.radii[0]!=0: radii_arr = concatenate((array([0.]), self.radii), axis=None)
            else: radii_arr = self.radii
            radii5 = (radii_arr*km_cm)**5
            return (4/5*pi*(radii5[1:] - radii5[:-1])*self.rhos*g_kg).sum()

    def test_M_I(self, radii_arr, weights_arr):
        # Calculates the mass and moment of inertia using the new set of parameters from a pull.

        # Returns:
        #     (bool) True if they pass the thresholds; False if they don't.
        

        if radii_arr[0]!=0: radii_arr = concatenate((array([0.]), radii_arr), axis=None)

        if self.isPoly==True:
            radii3 = (radii_arr*km_cm)**3
            radii4 = radii3*(radii_arr*km_cm)
            radii5 = radii4*(radii_arr*km_cm)
            radii6 = radii5*(radii_arr*km_cm)
            radii7 = radii6*(radii_arr*km_cm)

            M = 4*pi*(weights_arr*self.PolyCoeffs_a*1e-10*(radii5[1:]-radii5[:-1])/5 + weights_arr*self.PolyCoeffs_b*1e-5*(radii4[1:]-radii4[:-1])/4 + weights_arr*self.PolyCoeffs_c*(radii3[1:]-radii3[:-1])/3).sum()*g_kg
            I = 8/3*pi*(weights_arr*self.PolyCoeffs_a*1e-10*(radii7[1:]-radii7[:-1])/7 + weights_arr*self.PolyCoeffs_b*1e-5*(radii6[1:]-radii6[:-1])/6 + weights_arr*self.PolyCoeffs_c*(radii5[1:]-radii5[:-1])/5).sum()*g_kg*g_kg*g_kg*cm2_m2*km_cm*km_cm/(M*radii_arr[-1]*radii_arr[-1]*km_cm*km_cm)
            
        else:
            radii3 = (radii_arr*km_cm)**3
            radii5 = radii3*radii_arr*radii_arr*km_cm*km_cm

            M = (4/3*pi*(radii3[1:] - radii3[:-1])*self.rhos*weights_arr*g_kg).sum()
            I = (4/5*pi*(radii5[1:] - radii5[:-1])*self.rhos*weights_arr*g_kg).sum()/(M*radii_arr[-1]*radii_arr[-1]*km_cm*km_cm)

        if abs(M-self.M)/self.M < self.EarthM_Tolerance and abs(I - self.I) < self.EarthI_Tolerance: 
            self.EarthMs.append(M)
            self.EarthIs.append(I)
            return True
        else: return False

    def pull_radii(self):
        # Pulls a new set of radii from the provided center values and standard deviations.
      
        new_pull = []
        for iPull in range(self.NLayers):
            new_pull.append(random.normal(self.radii[iPull], self.radii_sigma[iPull]))
        return new_pull
    
    def pull_weights(self):
        #Pulls a new set of weights from the provided center values and standard deviations.
 
        new_pull = []
        for iPull in range(self.NLayers):
            new_pull.append(random.normal(self.density_weights[iPull], self.density_weights_sigma[iPull]))
        return new_pull
   
    def is_posdef(self, matrix):
        # Checks if a matrix is positive definite.

        # Args:
        #     matrix: Matrix to be evaluated
        # Returns: (bool)
  
        return all(linalg.eigvals(matrix) > 0)

    def make_covobj(self, iPar):
        # Make a covariance object to be added to a yaml file.

        # Args:
        #     iPar: the index of the parameter from 0 to (2*NLayers - 1)

        if iPar<self.NLayers:
            error = self.radii_sigma[iPar]
            value = self.radii[iPar]
        else: 
            error = self.density_weights_sigma[iPar-self.NLayers]
            value = 1.

        corr_matrix = self.corr.to_numpy()
        new_lines = ["\n- Systematic:\n", "    Names:\n", "      FancyName: %s\n"%(self.ParNames[iPar]), 
                     "      ParameterName: %s\n"%(self.ParNames[iPar]), "    DetID: [\"FD\",\"ATM\"]\n", "    Error: %.3f\n"%(error), "    FlatPrior: false\n", 
                     "    ParameterBounds: [0, 999]\n", "    ParameterGroup: Osc\n", "    ParameterValues:\n", "      Generated: %.2f\n"%value, 
                     "      PreFitValue: %.2f\n"%value, "    StepScale:\n", "      MCMC: 1.0\n", "    Type: Osc\n", "    Correlations:\n"]
        for j, entry in enumerate(corr_matrix[iPar]): new_lines.append("    - %s: %.5f\n"%(self.ParNames[j], entry))
        new_lines.append("\n")
        
        return new_lines
