import numpy as np
import math
import glob
from skimage import data, filters, color, morphology, measure
from skimage.segmentation import flood, flood_fill

def combine_3D_data_for_timestep(folder,step,gNx,gNy,gNz):


    files=glob.glob(folder+'/time_step-'+str(step)+'/*.raw')
    phi1_global = np.zeros((gNx,gNy,gNz))
    u_face_global = np.zeros((gNx,gNy,gNz))
    v_face_global = np.zeros((gNx,gNy,gNz))
    w_face_global = np.zeros((gNx,gNy,gNz))

    l=0
    
    if not files:
        # Either raise an error or return some defaults.
        raise ValueError(f"No .raw files found in {folder}/time_step-{step}")
    
    for file in files:

        data = np.fromfile(file,dtype=float)
        np.size(data)
        
        #removing the last 18 elements and splitting the variables
        n_elem = 18
        n_vars = 15
        u_face,v_face,w_face,c_pos,c_neg,phi1,phi2,rho,p,mu,k,esp,rel_perm,K,psi = np.split(data[:len(data)-n_elem],n_vars)
        info = data[len(data)-n_elem:]
    
        # Setting information
        sim_time = info[0];
        Ndim = info[1];
        Nx=int(info[2]);
        Ny=int(info[3]);
        Nz=int(info[4]);
        xcoord=int(info[5]);
        ycoord=int(info[6]);
        zcoord=int(info[7]);
        pad=int(info[8]);
        dx=info[9];
        dy=info[10];
        dz=info[11];
        gNx=int(info[12]);
        gNy=int(info[13]);
        gNz=int(info[14]);
        Nxnp=Nx-2*pad;
        Nynp=Ny-2*pad;
        Nznp=Nz-2*pad;
        
        #local coordinates
        xi = xcoord*Nxnp;
        yi = ycoord*Nynp;
        zi = zcoord*Nznp;
        xf = xcoord*Nxnp+Nxnp;
        yf = ycoord*Nynp+Nynp;
        zf = zcoord*Nznp+Nznp;
        
        phi1_global[xi:xf,yi:yf,zi:zf] = phi1.reshape((Nx,Ny,Nz))[pad:Nxnp+pad,pad:Nynp+pad,pad:Nznp+pad]
        u_face_global[xi:xf,yi:yf,zi:zf] = u_face.reshape((Nx,Ny,Nz))[pad:Nxnp+pad,pad:Nynp+pad,pad:Nznp+pad]
        v_face_global[xi:xf,yi:yf,zi:zf] = v_face.reshape((Nx,Ny,Nz))[pad:Nxnp+pad,pad:Nynp+pad,pad:Nznp+pad]
        w_face_global[xi:xf,yi:yf,zi:zf] = w_face.reshape((Nx,Ny,Nz))[pad:Nxnp+pad,pad:Nynp+pad,pad:Nznp+pad]
    
        #local array
        del(u_face)
        del(v_face)
        del(w_face)
        del(phi1)
        del(phi2)
        del(rho)
        del(c_pos)
        del(c_neg)
        del(p)
        del(mu)
        del(K)
        del(psi)
        del(k)
        del(esp)
        del(rel_perm)
        
        l=l+1
        
    return phi1_global,dx,sim_time,Nx

def axes_from_eigvals(eigvals, tol=1e-12):
    s1, s2, s3 = sorted(eigvals, reverse=True)
    wA =  s1 + s2 - s3
    wB =  s1 - s2 + s3
    wC = -s1 + s2 + s3

    #round tiny negatives to zero, reject true negatives
    wA, wB, wC = [0 if abs(w) < tol else w for w in (wA, wB, wC)]
    if min(wA, wB, wC) < 0:
        return None

    return np.sqrt(5/2*wA), np.sqrt(5/2*wB), np.sqrt(5/2*wC)

#Proposed Approach: Spherical Approximation of Surface Area with Hyperbolic Tangent Function interface
#Returns: Sum of the volume of the drops in the flow field
#phi_c: Cutoff Phi value to use for droplet masking
#flow: the 3D array of the flow field

def calculateDropletVolumes(phi_c,flow):
    #initialize Variables
    volume = 0
    
    #calculate Integral of the Hyperbolic Tangent Function (edge thickness)
    c = phi_c + 10**-100
    c = -1*c/(c-1)
    t = dx*np.log(np.cosh(0.5*np.log(c))) + 0.5*dx*np.log(c) + dx*np.log(2)
    
    #mask and label the regions
    mask = flow > phi_c
    labels, numDrops = measure.label(mask, return_num=True)
    p = 1.6075
    #calculate average volume of droplets
    V_avg = sum(sum(sum(flow)))/numDrops*dx**3
    
    #measure properties
    props = measure.regionprops(labels, flow, cache=True)  
    
    volume = np.zeros(labels.max())
    for index in range(0, labels.max()):
        #flood fill volume
        vol = props[index].area*dx**3
        axes = axes_from_eigvals(props[index].inertia_tensor_eigvals)
        a, b, c = axes          # A is the semi-major *diameter*, a = A/2
        # Thomsen/Michon surface-area approximation
        mean_term = ((a * b) ** p + (a * c) ** p + (b * c) ** p) / 3.0
        ellipsoid_SA = 4.0 * np.pi * mean_term**(1.0 / p) * dx**2
        
        volume[index] = vol*props[index].intensity_mean + t * ellipsoid_SA
    
    return volume

#Global grid size
Nx=256
Ny=256
Nz=256

folder = './200_87_256/forced-two-phase'
first = 120000
last = 176000
outputiter = 2000
eps = 0.51

steps = np.arange(first + outputiter, last, outputiter)

pdf_series = []

for i in steps:
    #Put Timestep into the Parameters - Update with path to data
    phi1,dx,sim_time,_ = combine_3D_data_for_timestep(folder,i,Nx,Ny,Nz)
    dx = eps*dx
    volume = calculateDropletVolumes(0.5,phi1)
    pdf_series.append(volume)

np.savez(folder + "/dropletVolumes.npz", *pdf_series)

