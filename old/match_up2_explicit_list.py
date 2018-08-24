#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from astropy.cosmology import WMAP7 as cosmo
import sys
import time
from scipy.interpolate import interp1d 
import dtk
pd.set_option('display.expand_frame_repr', False)
h=0.702
total_start = time.time()
param = dtk.Param(sys.argv[1])

lc_loc = param.get_string("lc_loc")
gltcs_loc = param.get_string("gltcs_loc")
shear_loc = param.get_string("shear_loc")
gio_loc   = param.get_string("gio_loc")
use_shear = param.get_bool("use_shear")
gltcs_file_num = param.get_int("gltcs_file_num")
gltcs_file_step = param.get_int("gltcs_file_step")
gltcs_steps     = param.get_int_list("gltcs_steps")
gltcs_internal  = param.get_int_list("gltcs_internal")
gltcs_str_z     = param.get_string_list("gltcs_str_z")
gltcs_file_list = param.get_string_list("gltcs_file_list")
nan_str_z       = param.get_string_list("nan_str_z")
use_mr_cut      = False;
mr_cut          = param.get_float("mr_cut")
tmp_to_disk     = param.get_bool("tmp_to_disk")
tmp_file_loc    = param.get_string("tmp_file_loc")
output          = param.get_string("output")
box             = param.get_bool("box")
stepz = dtk.StepZ(200,0,500)
zs = np.linspace(0,1.5,1000)
z_to_dl = interp1d(zs,cosmo.luminosity_distance(zs))
step2nan_str_z = {}
for i in range(0,len(nan_str_z)):
    step2nan_str_z[gltcs_steps[i]]=nan_str_z[i]
def load_bins(dic, loc,step,var_name,dtype,as_name=None,z_str=None):
    file_name = loc.replace("${step}",str(step)).replace("${var_name}",var_name)
    if(z_str != None):
        file_name = file_name.replace("${nan_str_z}",z_str)
    data = np.fromfile(file_name,dtype=dtype)
    if(as_name == None):
        dic[var_name]=data
    else:
        dic[as_name]=data

def add_step(dic,step):
    size = dic[dic.keys()[0]].size
    data = np.ones(size,dtype="int")*step;
    dic["step"]=data
    print "\tadded steps. data size:",data.size

def pandas_to_hdf5(df,output):
    hfile = h5py.File(output,'w')
    columns = list(df.columns.values)
    print("pandas to hdf5")
    for col in columns:
        #print("col",col)
        hfile[col] = df[col][:]
    hfile.flush()
    hfile.close()

def pandas_from_hdf5(output):
    hfile = h5py.File(output,'r')
    keys = hfile.keys()
    dic = {}
    for key in keys:
        dic[key]=hfile[key].value
    return pd.DataFrame(dic)

def load_shear(file_loc,step):
    dic = {}
    file_loc = file_loc.replace("${step}",str(step))
    #    data = np.loadtxt(file_loc,skiprows=1,dtype='i8,f4,f4,f4,f4,f4,f4,f4')
    nodeIndex,tx1,tx2,zs,s1,s2,k0,m0 = np.loadtxt(file_loc,skiprows=1,unpack=True,dtype={'names': ('nodeIndex', 'tx1', 'tx2','zs','s1','s2','k0','m0'), 'formats': ('i8',       'f4',  'f4','f4','f4','f4','f4','f4')})
    dic['nodeIndex']=nodeIndex
    dic['shear1']=s1
    dic['shear2']=s2
    dic['k0']    =k0
    dic['m0']    =m0
    dic['step']  =np.ones(nodeIndex.size,dtype='i4')*step
    print "loaded: ", nodeIndex.size
    return dic

def load_gio(file_loc,step,box=False):
    start = time.time()
    print "\treading in gio..."
    file_step_loc = file_loc.replace("${step}",str(step))
    dic = {}
    dic['nodeIndex']   = dtk.gio_read(file_step_loc,'nodeIndex')
    dic['hostIndex']   = dtk.gio_read(file_step_loc,'hostIndex')
    dic['hostHaloTag'] = dtk.gio_read(file_step_loc,"host_halo_tag")
    dic['hostHaloMass']= dtk.gio_read(file_step_loc,'host_halo_mass')/h
    dic['infallHaloTag']=dtk.gio_read(file_step_loc,'infall_halo_tag')
    dic['infallHaloMass']=dtk.gio_read(file_step_loc,'infall_halo_mass')/h
    dic['step']=np.ones(dic['nodeIndex'].size,dtype='i4')*step
    dic['placementType']=dtk.gio_read(file_step_loc,'positionType')
    if box:
        dic['x'] = dtk.gio_read(gio_loc.replace('${step}',str(step)),'x')
        dic['y'] = dtk.gio_read(gio_loc.replace('${step}',str(step)),'y')
        dic['z'] = dtk.gio_read(gio_loc.replace('${step}',str(step)),'z')
        dic['vx'] = dtk.gio_read(gio_loc.replace('${step}',str(step)),'vx')
        dic['vy'] = dtk.gio_read(gio_loc.replace('${step}',str(step)),'vy')
        dic['vz'] = dtk.gio_read(gio_loc.replace('${step}',str(step)),'vz')
        dic['redshift'] = np.ones(pos_dict['nodeIndex'].size)*stepz.get_z(step)

    print "\t\tdone. time:",time.time()-start
    return dic

def cat_hdf5(file_loc,nums,output):
    keys = h5py.File(file_loc.replace("${num}","0"),'r').keys()
    hfile= h5py.File(output,'w')
    i_key_max = len(keys)
    for i_key, key in enumerate(keys):
        data = []
        print "%.2f  %d/%d"%(float(i_key)/float(i_key_max),i_key,i_key_max)
        for i in range(0,nums):
            htmp = h5py.File(file_loc.replace("${num}",str(i)),'r')
            data.append(htmp[key][:])
        data_tot = np.concatenate(data)
        hfile[key]=data_tot
        hfile.flush()

dfs = []
total_start = time.time()
if not box:
    print "Loading LC"
    for step in gltcs_steps:
        start = time.time()
        print "working on step: ", step
        lc_dic = {}
        load_bins(lc_dic,lc_loc,step,"id",dtype='i8',as_name='nodeIndex')
        load_bins(lc_dic,lc_loc,step,"phi",dtype='f4')
        load_bins(lc_dic,lc_loc,step,"theta",dtype='f4')
        load_bins(lc_dic,lc_loc,step,"redshift",dtype='f4')
        load_bins(lc_dic,lc_loc,step,"x",dtype='f4')
        load_bins(lc_dic,lc_loc,step,"y",dtype='f4')
        load_bins(lc_dic,lc_loc,step,"z",dtype='f4')
        load_bins(lc_dic,lc_loc,step,"vx",dtype='f4')
        load_bins(lc_dic,lc_loc,step,"vy",dtype='f4')
        load_bins(lc_dic,lc_loc,step,"vz",dtype='f4')
        if(use_shear):
            load_bins(lc_dic,shear_loc,step,'sr1',dtype='f4',as_name='shear1',
                      z_str=step2nan_str_z[step])
            load_bins(lc_dic,shear_loc,step,'sr2',dtype='f4',as_name='shear2',
                      z_str=step2nan_str_z[step])
            load_bins(lc_dic,shear_loc,step,'kr0',dtype='f4',as_name='k0',
                      z_str=step2nan_str_z[step])
            load_bins(lc_dic,shear_loc,step,'xr1',dtype='f4',as_name='theta_lensed',
                      z_str=step2nan_str_z[step])
            load_bins(lc_dic,shear_loc,step,'xr2',dtype='f4',as_name='phi_lensed',
                      z_str=step2nan_str_z[step])
            load_bins(lc_dic,shear_loc,step,"mra",dtype='f4',as_name='m0',
                      z_str=step2nan_str_z[step])
        add_step(lc_dic,step)
        df1 = pd.DataFrame(lc_dic)
        gio_df = pd.DataFrame(load_gio(gio_loc,step))
        print "\t",df1.shape,gio_df.shape, "merged into",
        df = pd.merge(df1,gio_df,on=['nodeIndex','step'])
        print df.shape
        dfs.append(df)
else: # a snapshot box simluation
    for step in gltcs_steps:
        print "working on step: ", step
        # pos_dict = {}
        # pos_dict['x'] = dtk.gio_read(gio_loc.replace('${step}',str(step)),'x')
        # pos_dict['y'] = dtk.gio_read(gio_loc.replace('${step}',str(step)),'y')
        # pos_dict['z'] = dtk.gio_read(gio_loc.replace('${step}',str(step)),'z')
        # pos_dict['vx'] = dtk.gio_read(gio_loc.replace('${step}',str(step)),'vx')
        # pos_dict['vy'] = dtk.gio_read(gio_loc.replace('${step}',str(step)),'vy')
        # pos_dict['vz'] = dtk.gio_read(gio_loc.replace('${step}',str(step)),'vz')
        # pos_dict['nodeIndex'] = dtk.gio_read(gio_loc.replace('${step}',str(step)),'nodeIndex')
        # pos_dict['step'] = np.ones(pos_dict['nodeIndex'].size,dtype=int)*step
        # pos_dict['redshift'] = np.ones(pos_dict['nodeIndex'].size)*stepz.get_z(step)
        # df1 = pd.DataFrame(pos_dict)
        gio_df = pd.DataFrame(load_gio(gio_loc,step,box=True))
        # df = pd.merge(df1,gio_df,on=['nodeIndex','step'])
        # dfs.append(df)
        dfs.append(gio_df)

print "\t done. ",time.time()-start
total_end = time.time()
print "done loading LC", total_end - total_start
pos_df = pd.concat(dfs,ignore_index=True)
print "done making pd dataframe",time.time() -total_end
print "\t", pos_df.keys()
print pos_df.head()
print pos_df.shape

gal_cat = dtk.Catalog(gltcs_loc)
gal_cat.add_steps(gltcs_steps,in_file_steps=gltcs_internal)
gal_cat.add_var_step_replace("${step_var_name}",gltcs_steps,gltcs_str_z)

hfile = h5py.File(gltcs_file_list[0],'r')
keys_tmp = hfile["/Outputs/Output32/nodeData/"].keys()
keys = []

#don't copy these table that have these words
    avoids = ["AngularMomentum","hotHalo","basic","blackHoleCount","blackHoleJetPower","blackHoleRadiativeEfficiency","blackHoleSpin","darkMatterProfileScale","indicesBranchTip","position","satellite","siblingIndex","merge","merging","Gas"]
if box: # the snapshot box has less data columns than the lightcone because
    avoids += ["SED","Dust","dust","age","black","Continuum","Line","Metal","total"]
print "box: ", box
print avoids
for key in keys_tmp:
    if any(x in  key for x in avoids):
        print "throw",key
    elif((":contam_nitrogenII6584" in key) and ("balmerAlpha6563" not in key)):
        print "throw", key
    else:
        key_infile = key.replace(":z0.0000",":z${step_var_name}")
        key = key.replace(":z0.0000","")
        gal_cat.add_var("/Outputs/Output${step}/nodeData/"+key_infile,as_name=key)
        print "keep",key
        keys.append(key)


final_create = True
results = []
tmp_file_num =0
print gltcs_file_num
current_file_list = []
for i in range(0,len(gltcs_file_list)):
    start = time.time()
    gal_cat.add_subfile(i)
    current_file_list.append(gltcs_file_list[i])
    if(i%gltcs_file_step!=gltcs_file_step-1 and not (i==len(gltcs_file_list)-1)):
        #print "time: ",(time.time()-start)
        continue
    print "\n\n\n=============================================="
    print "working on files: ",i
    print "=============================================="
    print current_file_list
    gal_cat.set_explicit_files(current_file_list)
    current_file_list = []
    start_read = time.time()
    gal_cat.read_hdf5_no_step_file()
    gal_cat.clear_subfiles()
    dfs = []
    for step in gltcs_steps:
        gal_dic = {}
        for key in gal_cat[step].keys():
            gal_dic[key] = gal_cat[step][key]
        gal_dic["step"]=np.ones(gal_dic[gal_dic.keys()[0]].size,dtype='i4')*step
        dfs.append(pd.DataFrame(gal_dic))

    gal_df = pd.concat(dfs,ignore_index=True)
    gal_df = pd.merge(pos_df,gal_df,on=["step","nodeIndex"])

  
    if(tmp_to_disk):
        print "writing to subfile to disk"
        pandas_to_hdf5(gal_df,tmp_file_loc.replace("${num}",str(tmp_file_num)))
        tmp_file_num += 1
    else:
        results.append(gal_df)

    end = time.time()
    print "time: ",(end-start)

print "merging everything together..."

if(tmp_to_disk):
    print "here we go to cat"
    cat_hdf5(tmp_file_loc,tmp_file_num,output)
else:
    total_result = pd.concat(results,ignore_index=True)
    print lc_df.shape, "vs", total_result.shape 
    print "diff", lc_df.shape[0]-total_result.shape[0]
    pandas_to_hdf5(total_result,output)

total_end = time.time()
print "done. Time:",(total_end - total_start)

import subprocess 
p1 = subprocess.Popen("./add_metadata.py "+sys.argv[1], shell=True)
p1.wait()
p2 = subprocess.Popen("./rename.py "+sys.argv[1], shell=True)
p2.wait()
