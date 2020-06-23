#!/usr/bin/env python
# coding: utf-8

# In[3]:


#get_ipython().run_line_magic('load_ext', 'autotime')
import psi4
import yaml
import argparse
import numpy as np

# In[4]:


from joblib import Parallel, delayed


# In[ ]:





# In[8]:


def get_psi4_dftenergy(ts):
        #psi4.core.Molecule.set_multiplicity(psi4.geometry(ts), 2)
        #psi4.core.Molecule.set_molecular_charge(psi4.geometry(ts), 0)
        psi4.set_memory('8500 MB')
        psi4.set_options({'reference': 'uhf'})
        #psi4.geometry(ts)
        try:
           psi4.geometry(ts)
           return psi4.energy('wb97x-d/def2-svp')
        except:
           return np.nan


# In[6]:


def generate_grid(ymlfile=None):
    print('ymal file:', ymlfile)
    with open(ymlfile, 'r') as outfile:
        test_dict = yaml.load(outfile, Loader=yaml.FullLoader)
    
    energies=Parallel(n_jobs=-1)(delayed(get_psi4_dftenergy)(test_dict[k]) for k, v in test_dict.items())
    #energies=range(len(test_dict.keys()))
    energy_dic = dict()
    for kv, En in zip(test_dict.items(), energies):
        energy_dic[kv[0]] = En
    with open('en_'+ymlfile, 'w') as outfile:
        yaml.dump(energy_dic, outfile, default_flow_style=False)


# In[12]:


if __name__ == '__main__':
    """generate 2d scan energies"""
    parser = argparse.ArgumentParser()
    parser.add_argument('ymlfile', type=str, help='yaml file name do nor pass full path')
    args = parser.parse_args()
    generate_grid(ymlfile=args.ymlfile)


# In[ ]:




