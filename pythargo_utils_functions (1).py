path_main = '/home/581/ao2582/SONIF/data/ARGO/'
#path_sat = '/home/aoetjens/PhD/Data/chla_oc/monthly/' # satellite products
#path = '/home/aoetjens/PhD/ARGO/Profiles_SO/'
path_g = '/../../../g/data/jk72/ao2582/'
path = path_g + 'Profiles/'
path_sonif ='/home/581/ao2582/SONIF/data/ARGO/' # path to floats
path_argo =  '/../../../g/data/jk72/ao2582/Profiles/'

import sys 
import os
sys.path.append(os.path.abspath(path_main))
# Import from other scipts
from pythargo_utils import *

print('interp_argo -- imported')

print('interp_argo -- imported')
def interp_argo(argo_df):
    argo_df = argo_df[argo_df['mask_doxy'] == 1]
    print('mask_doxy')
    argo_df = argo_df[argo_df['mask_nit'] == 1]
    argo_df = argo_df[argo_df['mask_bbp'] == 1]
    print('mask_bbp')

    # list of wmo's that satisfy criterion, this can e more elegant....
    files = []
    for wmo in argo_df['wmo'].drop_duplicates().values:
        files.append(glob.glob(path+ str(wmo)+'_Sprof.nc')[0])
    print(path)
        
    print('Floats with DOXY and NITRATE and BBP #: ' + str(len(files)))

    fnames = files
    # this is the pressure coords you want to interpolate your data on.
    pres_itp = np.arange(0,1501,1) 
    
    # create a dataset for each float
    var_list = ['TEMP_ADJUSTED','PSAL_ADJUSTED','CHLA_ADJUSTED', 'NITRATE_ADJUSTED', 'DOXY_ADJUSTED']
    print(var_list)
    ds_dict = {}
    for ifn,fn in enumerate(fnames[:]):
        
        # could have this as looping over dedicated wmo's instead
        if len(ds_dict) == 0:
            old_i = 0
        else:
            old_i = ds_dict[list(ds_dict.keys())[-1]].n_prof.values[-1]+1
        
        ds = xr.open_dataset(fn)
        wmo = int(ds['PLATFORM_NUMBER'].data[0])
        print(str(wmo))
        PRES_ADJUSTED = ds['PRES_ADJUSTED'].transpose('N_LEVELS','N_PROF').values # put the levels to the first axis, and cyc number to the second
        n_cyc = PRES_ADJUSTED.shape[1] # this is the number of cycles, is it also N_PROFS? Should be right? Yes
        print(fn, 'start loop over n_cyc')
        new_cycle = True
        for i in range(n_cyc):
            k = i + old_i
    
            ds_new = xr.Dataset({})
            pres_old = PRES_ADJUSTED[:,i]
            # check if cycle is within geographic boundary
            nprof = ds['N_PROF'].values[i] #ds['CYCLE_NUMBER'].values[i]
            cycle_number = ds['CYCLE_NUMBER'].values[i]
    
            if nprof in argo_df[argo_df['wmo']==str(wmo)]['cyc'].values:
                
                
                # should this be N_PROF instead? Need to double check!
                lat = ds['LATITUDE'].values[i]
                lon = ds['LONGITUDE'].values[i]
                time = ds['JULD'].values[i]
                
                ### BBP computation
                
                # 1. select float and profile
                dbbp = argo(wmo,path).sel(N_PROF = nprof)
                # 2. get rid of nan's and check if all nan
                var_qc = dbbp['BBP700_ADJUSTED' + '_QC'].values
                mask_qc = (var_qc==b'1') | (var_qc==b'2') | (var_qc==b'5') | (var_qc==b'8')
                mask = (~np.isnan(dbbp.BBP700_ADJUSTED)) & (dbbp.BBP700_ADJUSTED < 0.04) & mask_qc
                
                if (mask).sum() > 7:
                   # bbp700 = dbbp.BBP700_ADJUSTED.where(mask,np.nan).rolling(N_LEVELS = 7,center = True, min_periods=1).min(skipna=True).rolling(N_LEVELS = 7,center = True, min_periods=1).max(skipna=True)
                   # bbp700 = dbbp.BBP700_ADJUSTED[mask].rolling(N_LEVELS = 7,center = True, min_periods=1).min(skipna=True).rolling(N_LEVELS = 7,center = True, min_periods=1).max(skipna=True)
                    bbp700 = dbbp.BBP700_ADJUSTED[mask].rolling(N_LEVELS = 10,center = True, min_periods=1).median(skipna=True)
                    # 3. interpolate data onto uniform pressure grid
                    f = interpolate.interp1d(dbbp.PRES_ADJUSTED[mask].data,bbp700.data,bounds_error=False,fill_value='nan')
                    bbs_itp = f(pres_itp)
                    
                    spikes = dbbp.BBP700_ADJUSTED[mask].data - bbp700.data
                    f = interpolate.interp1d(dbbp.PRES_ADJUSTED[mask].data,spikes,bounds_error=False,fill_value='nan')
                    bbl_itp = f(pres_itp)
                    
                else:
                    bbs_itp = np.nan * np.ones_like(pres_itp)
                    bbl_itp = np.nan * np.ones_like(pres_itp)
                    
                bbs_itp = bbs_itp.reshape(len(bbs_itp),1)
                bbl_itp = bbl_itp.reshape(len(bbl_itp),1)
                
                
                ds_new['bbp700_bs'] = xr.DataArray(data = bbs_itp.data,
                                               dims = ['pressure','n_prof'],
                                               coords = {'pressure': (['pressure'],pres_itp),
                                                         'n_prof':(['n_prof'],[k]),
                                                         'cycle_number':(['n_prof'],[cycle_number]),
                                                         'lon':(['n_prof'],[lon]),
                                                         'lat':(['n_prof'],[lat]),
                                                         'time':(['n_prof'],[time]),
                                                         'wmo':(['n_prof'],[wmo])})
                    
                ds_new['bbp700_bl'] = xr.DataArray(data = bbl_itp.data,
                                               dims = ['pressure','n_prof'],
                                               coords = {'pressure': (['pressure'],pres_itp),
                                                         'n_prof':(['n_prof'],[k]),
                                                         'cycle_number':(['n_prof'],[cycle_number]),
                                                         'lon':(['n_prof'],[lon]),
                                                         'lat':(['n_prof'],[lat]),
                                                         'time':(['n_prof'],[time]),
                                                         'wmo':(['n_prof'],[wmo])})
                
                ### CHL computation
                
                # 2. get rid of nan's and check if all nan
               # mask = (~np.isnan(dbbp.CHLA_ADJUSTED)) and (dbbp.CHLA_ADJUSTED > 0)
               # 2. get rid of nan's and values < 0
                var_qc = dbbp['CHLA_ADJUSTED' + '_QC'].values
                mask_qc = (var_qc==b'1') | (var_qc==b'2') | (var_qc==b'5') | (var_qc==b'8') 
                mask = (dbbp.CHLA_ADJUSTED > -0.001) & mask_qc
                
                if (mask).sum() > 7:
                    chla = dbbp.CHLA_ADJUSTED[mask].rolling(N_LEVELS = 7,center = True, min_periods=1).median(skipna=True)
                    # 3. interpolate data onto uniform pressure grid
                    f = interpolate.interp1d(dbbp.PRES_ADJUSTED[mask].data,chla.data,bounds_error=False,fill_value='nan')
                    chla_itp = f(pres_itp)
                    
                   
                else:
                    chla_itp = np.nan * np.ones_like(pres_itp)
                
                    
                chla_itp = chla_itp.reshape(len(chla_itp),1)
               
                
                ds_new['chla'] = xr.DataArray(data = chla_itp.data,
                                               dims = ['pressure','n_prof'],
                                               coords = {'pressure': (['pressure'],pres_itp),
                                                         'n_prof':(['n_prof'],[k]),
                                                         'cycle_number':(['n_prof'],[cycle_number]),
                                                         'lon':(['n_prof'],[lon]),
                                                         'lat':(['n_prof'],[lat]),
                                                         'time':(['n_prof'],[time]),
                                                         'wmo':(['n_prof'],[wmo])})
                
                
                ### Interpolation for other variables
                
                for var in var_list:
                    if var in ds.keys():
                        var_old = ds[var].transpose('N_LEVELS','N_PROF').values[:,i]
                        
                    
                        var_qc = ds[var + '_QC'].transpose('N_LEVELS','N_PROF').values[:,i]
                        mask_qc = (var_qc==b'1') | (var_qc==b'2') | (var_qc==b'5') | (var_qc==b'8')
                        mask_valid = np.isfinite(var_old) & np.isfinite(pres_old) & mask_qc 
                            
                        if mask_valid.sum() > 1:
                            variable_valid = var_old[mask_valid]
                            pres_old_valid = pres_old[mask_valid]
                                
                            f = interpolate.interp1d(pres_old_valid,variable_valid,bounds_error=False,fill_value='nan')
                            var_itp = f(pres_itp)
                        else:
                               # print('nans!')
                            var_itp = np.nan * np.ones_like(pres_itp)
                        var_itp = var_itp.reshape(len(var_itp),1)
                    else:
                        print('No Nitrate')
                        var_itp = np.nan * np.ones_like(pres_itp)
                        var_itp = var_itp.reshape(len(var_itp),1)
               
                    ds_new[var] = xr.DataArray(data = var_itp,
                                               dims = ['pressure','n_prof'],
                                               coords = {'pressure': (['pressure'],pres_itp),
                                                         'n_prof':(['n_prof'],[k]),
                                                         'cycle_number':(['n_prof'],[cycle_number]),
                                                         'lon':(['n_prof'],[lon]),
                                                         'lat':(['n_prof'],[lat]),
                                                         'time':(['n_prof'],[time]),
                                                         'wmo':(['n_prof'],[wmo])})
                  
                    
                    # if i == 0:
                    if new_cycle:
                        DS = ds_new.copy()
                        new_cycle=False
                    else:
                        DS = xr.merge([DS,ds_new])
        ds_dict[wmo] = DS
        
    return ds_dict


print('bgcalc_argo -- imported')
def bgcalc_argo(ds_dict,fosor):
    ds_copy = {}

    # loop over all floats
    for wmo in list(ds_dict.keys()):
        print(wmo)
        
        # creat copy of ds_dict and only keep float entries where bbp and chl are not all nan
        ds_new_1 = ds_dict[wmo].copy().dropna(dim='n_prof',how='all', subset =['bbp700_bs'])
        ds_new = ds_new_1.dropna(dim='n_prof',how='all', subset =['chla'])
        # initialize arrays
        dyn_height = np.array([])
        t_zone = np.array([])
        mixed_layer_depth = np.array([])
        z_euph = np.array([])
        z_euph_argo = np.array([])
        z_prod = np.array([])
        poc_koest = np.array([])
        poc_koest_pi = np.array([])
        poc_koest_xi = np.array([])
        iNO3_prod = np.array([])
        iNO3_max = np.array([])
        iNO3_min = np.array([])
        iNO3_50 = np.array([])
        eke = np.array([])
        koest23_xi = np.array([])
        ice = np.array([])
        fill_npp = np.array([])
        par = np.array([])
        inpp = np.array([])
        cflux = np.array([])
        chl_sat = np.array([])
        zone_label = np.array([])

            
        SP = ds_new.PSAL_ADJUSTED
        p = ds_new.pressure
        pressure_reshaped = ds_new['pressure'].values[:, np.newaxis]
        lat = ds_new.lat
        lat_reshaped = ds_new['lat'].values[np.newaxis,:]
        lon = ds_new.lon
        t = ds_new.TEMP_ADJUSTED
            
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        CT = gsw.CT_from_t(SA, t, p)

        spice = gsw.spiciness1(SA, CT).data
       # Add spice as a new data variable to the xarray dataset
        ds_new['spice'] = (('pressure', 'n_prof'), spice)    
        
        nsquared = gsw.Nsquared(SA, CT, pressure_reshaped, lat=lat_reshaped)[0]
        zeros_row = np.zeros((1, len(ds_new.n_prof)))
        # Append the zeros row to N_squared along the first axis
        N_squared = np.concatenate((nsquared, zeros_row), axis=0)
        ds_new['n2'] = (('pressure', 'n_prof'), N_squared)
        
        AOU =  gsw.O2sol_SP_pt(SP, t).data - ds_new.DOXY_ADJUSTED.data
        ds_new['aou'] = (('pressure', 'n_prof'), AOU)


        density = gsw.density.sigma0(SA, CT)
        ds_new['rho'] = (('pressure', 'n_prof'), density.data)
        delta_density = density - density.isel(pressure=10)
        ds_new['del_rho'] = (('pressure', 'n_prof'), delta_density.data)
            
        bbp470 = ds_new['bbp700_bs'].data * (470/700)**(-1)
        Cphyto = 12128 * bbp470 + 0.59
        ds_new['cphyto'] = (('pressure', 'n_prof'), Cphyto)
        
        #spikes = ds['BBP700_ADJUSTED'].data - ds['BBP700_ADJUSTED'].where(ds['BBP700_ADJUSTED']<0.004).rolling(pressure=42,center=True,min_periods=1).max(skipna=True).data
        # POC based on Johnson
        ds_new['poc_bs'] =  (('pressure', 'n_prof'), 31200* ds_new['bbp700_bs'].data + 3.04 )        
        ds_new['poc_bl'] =  (('pressure', 'n_prof'), 31200* ds_new['bbp700_bl'].data + 3.04)
        
        
        ds_new['month'] = ('n_prof', ds_new.time.dt.month.data)
        ds_new['sonif'] = ('n_prof', [fosor]*len(ds_new.n_prof))

       # initialize arrays for computations ona prof by prof basis
     
        # loop over all profiles for each wmo
        for nprof in ds_new.n_prof:
            ds = ds_new.sel(n_prof=nprof)
            
            SP = ds.PSAL_ADJUSTED
            p = ds.pressure
            lat = ds.lat
            lon = ds.lon
            t = ds.TEMP_ADJUSTED
            bbp =  ds.bbp700_bs
            ## 25.10.2024 aoetjens: set values within uncertainty to zero (below 200m)
            chla = ds.chla.where((ds.chla>0.04) | (np.isnan(ds.chla)) | (ds.pressure < 200),0)
            argo_time = ds.time
                   
          # Compute mixed payer depth according to pressure criterion
            mixed_layer_depth = np.append(mixed_layer_depth, np.min(ds.pressure.where(ds.del_rho > 0.03), axis = 0))
          # compute euphotic zone 
            chl_surf = ds.CHLA_ADJUSTED.sel(pressure = np.arange(0,15,1)).mean(skipna=True).data/2 # scaled by two to match satellite
            log10zeu = 1.524 - 0.436 * np.log10(chl_surf) - 0.0145 * np.log10(chl_surf)**2 + 0.0186 * np.log10(chl_surf)**3
            z_euph_argo = np.append(z_euph_argo, 10**(log10zeu))
          # compute productive zone
            z_prod = np.append(z_prod, np.nanmax([np.min(ds.pressure.where(ds.del_rho > 0.03), axis = 0),10**(log10zeu),0]))
    
          # Compute dynamic height
          # dyn m = m/9.81
            SA = gsw.SA_from_SP(SP, p, lon, lat)
            CT = gsw.CT_from_t(SA, t, p)
            
           # dyn_height =np.append(dyn_height, gsw.geo_strf_dyn_height(SA[500:1501], CT[500:1501], p[500:1501], p_ref=1500, interp_method='pchip')[0]/9.81)
            dyn_height =np.append(dyn_height, gsw.geo_strf_dyn_height(SA[50:1001], CT[50:1001], p[50:1001], p_ref=1000, interp_method='pchip')[0]/9.81)
          #  print(dyn_height)
        
            if dyn_height[-1] > 0.87:
                label = 'saz'
            elif dyn_height[-1] > 0.62 and dyn_height[-1] < 0.87:
                label = 'pfz'
            elif dyn_height[-1] > 0.45 and dyn_height[-1] < 0.62:
                label = 'az'
            elif dyn_height[-1] < 0.45:
                label = 'sz'
            else:
                label = 'nz' # nan, no dynamic height value
                
            zone_label = np.append(zone_label, label)
       

        # Append profs to dataset
        ds_new['dyn_height'] = ('n_prof', dyn_height)
        ds_new['zone'] = ('n_prof', zone_label)
        ds_new['mld'] = ('n_prof', mixed_layer_depth)
        ds_new['zeuph_argo'] = ('n_prof', z_euph_argo)
        ds_new['z_prod'] = ('n_prof', z_prod)
        
        ds_copy[wmo] = ds_new

    return ds_copy



def create_argo_df(chl_region, path, wmos = []):
    import glob
    print(path)
    files = []
   
    # if there is a list of floats we read in those
    if len(wmos) > 0:
        for wmo in wmos:
            print(wmo)
            files.append(glob.glob(path + str(wmo) +'_Sprof.nc')[0])
    # if no floats provided we search he whole box region for available floats        
    else:
        for file in glob.glob(path + '*_Sprof.nc'):
            files.append(file)
            
        
        
    
    print('Going through ARGO profiles ....')
    wmo = np.array([])
    cyc = np.array([])
    mask = np.array([])
    lon = np.array([])
    lat = np.array([])
    chl = np.array([])
    mon = np.array([])
    yer = np.array([])
    mask_doxy = np.array([])
    mask_par = np.array([])
    mask_nit =np.array([])
    mask_bbp =np.array([])
    res_nitrate = np.array([])
    res_chl = np.array([])
    res_temp =np.array([])
    
    
    for i in range(len(files)):
        dq = xr.open_dataset(files[i])
        
        platform_number = int(dq['PLATFORM_NUMBER'].data[0])
    
        # check fo available parameter and flag them with 0/1
        if 'DOXY' in dq.data_vars:
            doxy = 1
        else:
            doxy = 0
        if 'DOWNWELLING_PAR' in dq.data_vars:
            par = 1
        else:
            par = 0
        if 'NITRATE_ADJUSTED' in dq.data_vars:
            nit = 1
        else:
            nit = 0 
        if 'BBP700_ADJUSTED'  in dq.data_vars:
            bbp = 1
        else:
            bbp = 0 
           
        # read in variables for ? TODO What do we need this for actually?   
        timestamp = pd.to_datetime(dq['JULD'].data)
        chla = dq['CHLA_ADJUSTED'].data
        pres = dq['PRES_ADJUSTED'].data
        # read in coordinates to check for in-bloom/out-bloom condition
        cycn =  dq.N_PROF.data #dq['CYCLE_NUMBER'].data 
        lonargo = dq['LONGITUDE'].data
        latargo = dq['LATITUDE'].data
        # close the curent float xarray 
      #  dq.close()
      #  del(dq)
        print(platform_number)
       # loop over all cycles and check if coordinates are inside box
        for dip in cycn:
            
           # create input for box_region  
           pll = geometry.Point(lonargo[cycn == dip], latargo[cycn == dip])

           if chl_region.contains(pll)[0]:    
               # check resolution in upper 30 dbar
                reschl = len(dq.CHLA_ADJUSTED.where(dq.PRES_ADJUSTED < 30).sel(N_PROF = dip)[~np.isnan(dq.CHLA_ADJUSTED.where(dq.PRES_ADJUSTED < 30).sel(N_PROF = dip))])
                if nit == 1:
                    resnitrate = len(dq.NITRATE_ADJUSTED.where(dq.PRES_ADJUSTED < 30).sel(N_PROF = dip)[~np.isnan(dq.NITRATE_ADJUSTED.where(dq.PRES_ADJUSTED < 30).sel(N_PROF = dip))])
                else:
                    resnitrate = 0
                restemp = len(dq.TEMP_ADJUSTED.where(dq.PRES_ADJUSTED < 30).sel(N_PROF = dip)[~np.isnan(dq.TEMP_ADJUSTED.where(dq.PRES_ADJUSTED < 30).sel(N_PROF = dip))])
                # store 
                if (reschl > 0): 
                
                    year = timestamp[cycn == dip].year[0]
                    month = f'{timestamp[cycn == dip].month[0]:02}'    
                    wmo = np.append(wmo,str(platform_number))
                    
                    cyc = np.append(cyc, int(round(dip)))
                    lat = np.append(lat, latargo[cycn == dip])
                    lon = np.append(lon, lonargo[cycn == dip])
                    mon = np.append(mon, int(month))
                    yer = np.append(yer, int(year))
                    
                    mask_doxy = np.append(mask_doxy, doxy)
                    mask_nit = np.append(mask_nit, nit)
                    mask_par = np.append(mask_par, par)
                    mask_bbp = np.append(mask_bbp, bbp)
              
    
                    res_temp = np.append(res_temp, restemp)
                    res_chl = np.append(res_chl, reschl)
                    res_nitrate = np.append(res_nitrate, resnitrate)
    
    argo_df = pd.DataFrame({
        'wmo': wmo,
        'cyc': np.array(cyc, dtype=int),
        'lat' : lat,
        'lon' : lon,
        'monargo' : mon,
        'yeargo' : yer,
        'mask_doxy' :mask_doxy,
        'mask_par' :mask_par,
        'mask_nit' :mask_nit,
        'mask_bbp' :mask_bbp,
        'res_temp' :res_temp,
        'res_nitrate' :res_nitrate,
        'res_chl' :res_chl,
        })
    
    print('DONE')
    print('Profiles #: ' + str(len(argo_df['wmo'])))
    print('Floats #: ' + str(len(argo_df['wmo'].drop_duplicates())))
    
    return argo_df
   

print('build_argo -- imported')

def build_argo(fosor,argo_df, to_netcdf = False):
    # SIC : selection, interpolating, calculating
    
    '''
    calls the functions:
        pythargo_region
        interp_argo
        bgcalc_argo
        
    Returns:
        xarray Dataset  with interpolated BGC parameter for all profiles in selected region
        dict of all selected floats with parameter availability
    '''
    print(fosor)
    # import region, argo_df#, chla_oc #, chl_region
    # build one pythargo_region file for all the regions
    print('argo profiles -- selected')
    argo_dict = interp_argo(argo_df)
    print('argo profiles -- interpolated')
    argo_DS = bgcalc_argo(argo_dict,fosor)
    print('bgc parameter -- computed')
    ds_floats = xr.concat([*argo_DS.values()],dim='n_prof')
    if to_netcdf == True:
        ds_floats.to_netcdf('SONIF_'+str(fosor)+ '.nc', encoding={'time':{'dtype':'object'},'sonif':{'dtype':'str'},'zone':{'dtype':'str'}})
        print('data stored as ' + 'SONIF_'+str(fosor)+ '.nc')
    return ds_floats, argo_dict


def ssh_zones(path_g, ds_cbpm):
    # Path to AVISO MADT files
    file_list = sorted(glob.glob(path_g + 'AVISO/*_madt_*_m*.nc'))
    
    # Function to load and clean each dataset
    def preprocess(ds):
        # Drop 'nv' dimension if it exists
        ds = ds.drop_dims('nv', errors='ignore')
        return ds
    
    # Open all files as a single xarray.Dataset along the time dimension
    ds_all = xr.open_mfdataset(file_list, 
                               concat_dim='time', 
                               combine='nested', 
                               preprocess=preprocess)
    
    ds_all = ds_all.rename({'time': 'month'})
    ds_all = ds_all.assign_coords(month=np.arange(1, len(ds_all.month) + 1))
    # Confirm structure
    print(ds_all)
    
    # drop data above 40South
    ds_ssh = ds_all.where(ds_all.latitude < -40, drop=True)
    
    # shuffle longitude to -180 to 180
    ds_ssh = ds_ssh.assign_coords(
        longitude=(((ds_ssh.longitude + 180) % 360) - 180)
    ).sortby('longitude')
    
    
    # create mask from ssh contour regions
    levs = [-1.1, -0.65, -0.1, 1]  # ssh
    az_mask = (ds_ssh.adt >= -1.1) & (ds_ssh.adt < -0.65)& (ds_ssh.latitude > -70)
    pfz_mask = (ds_ssh.adt >= -0.65) & (ds_ssh.adt < -0.1)& (ds_ssh.latitude > -70)
    saz_mask = (ds_ssh.adt >= -0.1) & (ds_ssh.latitude < -40) & (ds_ssh.latitude > -70)
    
    # match grid to carbon data
    ssh_interp = ds_ssh.adt.interp(
        latitude=ds_cbpm.lat,
        longitude=ds_cbpm.lon,
        method="nearest"
    )
    
    
    az_mask = (ssh_interp >= -1.1) & (ssh_interp < -0.65)& (ssh_interp.latitude > -70)
    pfz_mask = (ssh_interp >= -0.65) & (ssh_interp < -0.1)& (ssh_interp.latitude > -70)
    saz_mask = (ssh_interp >= -0.1) & (ssh_interp.latitude < -40) & (ssh_interp.latitude > -70)
    return ds_ssh, az_mask, pfz_mask, saz_mask



def load_carbon(path_g):
    filepath = glob.glob(path_g + 'AQUA_MODIS/Climatology 1.mat')[0]
    
    with h5py.File(filepath, 'r') as f:
        # Load data
        lat = f['Lat'][:]
        lon = f['Lon'][:]
        time = f['time_steps'][:].squeeze()
        carbon = f['Carbon'][:]
    
        # Convert from HDF5 format (which stores Fortran-style) to NumPy
        lat = np.squeeze(lat)
        lon = np.squeeze(lon)
        carbon = np.transpose(carbon, (2, 1, 0))  # Now: (time, lat, lon)
    
        # Extract unique 1D lat/lon
        lat_1d = np.unique(lat)
        lon_1d = np.unique(lon)
    
    # Build xarray Dataset assuming carbon shape is (time, lat, lon)
    ds_cbpm = xr.Dataset(
        {
            'carbon': (['month', 'lat', 'lon'], carbon)
        },
        coords={
            'month': time,
            'lat': lat_1d,
            'lon': lon_1d
        }
    )
    
    
    ds_cbpm = ds_cbpm.where(ds_cbpm.lat < -40, drop=True)
    
    # integrate over area, do it the correct way
    
    R = 6.371e6  # Earth radius in meters [m]
    carbon = ds_cbpm
    
    # Compute spacing in radians
    dlat = np.deg2rad(np.abs(carbon.lat[1] - carbon.lat[0]))
    dlon = np.deg2rad(np.abs(carbon.lon[1] - carbon.lon[0]))
    
    # Compute area grid (2D, lat x lon)
    lat_radians = np.deg2rad(carbon.lat)
    area_1d = R**2 * dlat * dlon * np.cos(lat_radians)
    
    # Expand to 2D (lat x lon)
    area_2d = xr.DataArray(
        np.outer(area_1d, np.ones_like(carbon.lon)),
        coords={"lat": carbon.lat, "lon": carbon.lon},
        dims=["lat", "lon"]
    )
    
    
    ds_cbpm['tot_carbon'] = ds_cbpm.carbon * area_2d    #[mg/m3 * m2 = mg/m] for each pixel
    
    print('created ds_cbpm.tot_carbon [mg/m3 * m2 = mg/m] ')
    
    return ds_cbpm




##% BGC specific functions

def Koest23_modelB_700p(bbp, chla, alignment=1, bias='power', alpha=0.125):
    '''
    from Koestner et al. (2024), Github
    converted to Python
    '''
    if len(bbp[~np.isnan(bbp)]) == 0 or len(chla[~np.isnan(chla)]) == 0:
        raise ValueError("bbp and chla arrays must not be empty.")

    
    # Set any negative values to 0
    bbp = np.maximum(bbp, 0)
    chla = np.maximum(chla, 0)

    bbp = bbp * alignment  # alignment adjustment for use with different bbp sensor

    comp = chla / bbp
    comp[comp == 0] = np.min(comp[comp > 0])  # fix all 0 chla values to minimum comp value
    comp[comp > 2000] = 2000  # fix large comp values

    # New model B coefficients from Koestner et al. 2023
    k_modelB = [52.8187501942431, 0.135288289126603, 0.884851394851513, 0.226810214797258]

    # Calculate POC from multivariable model
    pocB = (k_modelB[0] * (bbp ** k_modelB[1]) * 
            (comp ** k_modelB[2]) * 
            (comp ** (k_modelB[3] * np.log10(bbp))))

    # Continuous bias correction for modeled POC < bias_lim
    if bias.startswith('pow'):
        bias_modelB_coef = [1.46918996207386, -0.734453171035830]
        bias_lim = 36.8
        pocB[pocB < bias_lim] = (pocB[pocB < bias_lim] ** bias_modelB_coef[0]) * 10 ** bias_modelB_coef[1]
    elif bias.startswith('lin'):
        bias_modelB_coef = [1.42263442605374, -14.7264359822160]
        bias_lim = 34.8
        pocB[pocB < bias_lim] = (pocB[pocB < bias_lim] * bias_modelB_coef[0]) + bias_modelB_coef[1]

    # Prediction interval stuff
    S = np.array([
        [0.0248767643354040, 0.00967335656242312, -0.0105363562351115, -0.00404078515530673],
        [0.00967335656242312, 0.00405822189121109, -0.00402605221279587, -0.00166729112989270],
        [-0.0105363562351115, -0.00402605221279587, 0.00483637706215940, 0.00180593541611042],
        [-0.00404078515530673, -0.00166729112989270, 0.00180593541611042, 0.000727608368784931]
    ])  # covariance matrix of coefficients

    mse = 0.0311478519483073  # mean square error from regression model

    df = 407 - 4 - 2  # take away 6 for model coefficients (4 + 2 for bias correction)
    t = stats.t.ppf(1 - alpha, df)

    # Calculate prediction bounds log10(y +/- e), where e is prediction bound defined
    # as e = t * sqrt(MSE + xSx'); x is the row vector of the design
    # matrix for new predictors (i.e., [1 log10(bbp) log10(comp) log10(comp)*log10(bbp)])
    x = np.column_stack((np.ones(len(bbp)), np.log10(bbp), np.log10(comp), np.log10(comp) * np.log10(bbp)))

    e = np.zeros(len(bbp))
    for i in range(len(bbp)):
        e[i] = t * np.sqrt(mse + x[i, :] @ S @ x[i, :])

    # Rename final data
    POC = pocB
    POC[POC < 0] = 0  # fix negative values to be 0
    pi = np.column_stack((np.log10(pocB) - e, np.log10(pocB) + e))
    PI = 10 ** pi
  #  print(comp)
    return POC, PI, comp
