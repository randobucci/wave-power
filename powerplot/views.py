from django.shortcuts import get_object_or_404, render 
from django.http import HttpResponse
from django.http import Http404

# Imports required for most plots
import os
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
import math
import numpy as np
from datetime import datetime
import calendar
#- Get all the cdippy modules
#bin_dir = '../cdippy'
from cdippy.cdipnc import *
from cdippy.utils import *
from cdippy.stndata import *


# Default settings
default_x_inch = 9 # Width of image in inches
default_y_inch = 5 # Height of image in inches

def index(request):
    return HttpResponse('Plots are: wave_power')

# Polar Spectrum using netcdf model directional spectrum data
def wave_power(request):
    # Set plot parameters: stn, datearg, typearg
    stn = request.GET.get('stn','100')
    id = stn
    stn = stn+'p1'                               #- 100p1
    datearg = request.GET.get('date','201612')
    typearg = request.GET.get('type','frequency')

    #- Initialize plot types
    plot_types = ["frequency","percent","value"]

    #- Wave power density (wave energy flux): watts per meter of wave crest width
    #- Units are kW/m
    def calc_power(Hs,Te):
        #- Define constants for power calculation
        rho = 1025          #- Density of seawater kg/m^3
        g = 9.8067    #- Acceleration due to gravity m/sec^2
        Hs = np.array(Hs)
        Te = np.array(Te)
        p0 = ((rho*g**2)/(64*math.pi)) * (Te * Hs**2)
        return p0/1000    #- return kW/m

    def calc_energy(Hs):
        g = 9.8067    #- Acceleration due to gravity m/sec^2
        rho = 1025          #- Density of seawater kg/m^3
        e0 = (1/16)*rho*g*(Hs**2)
        return e0 

    #- Check datearg for range 201612-201701
    datestr = datearg.split('-')
    if(len(datestr) == 1):          #- Deal with start/end dates
        startstr = datearg
        endstr = datearg
    else:
        startstr = datestr[0]   
        endstr = datestr[1]   

    #- Get end month/day
    start_yr = startstr[0:4]
    end_yr = endstr[0:4]
    if(len(startstr) == 4):                   #- 2016
        start_mo = '01'
        end_mo = '12'
        start_day,end_day = calendar.monthrange(int(end_yr),int(end_mo))
        startstr = start_yr+'-'+start_mo
        start = startstr+'-01 00:00:00'
        endstr = end_yr+'-'+end_mo
        end = endstr+'-'+str(end_day)+' 23:59:59'
    elif(len(startstr) == 6):                   #- 201612
        start_mo = startstr[4:6]
        end_mo = endstr[4:6]
        start_day,end_day = calendar.monthrange(int(end_yr),int(end_mo))
        startstr = start_yr+'-'+start_mo
        start = startstr+'-01 00:00:00'
        endstr = end_yr+'-'+end_mo
        end = endstr+'-'+str(end_day)+' 23:59:59'
    elif(len(startstr) == 8):                 #- 20161201
        start_mo = startstr[4:6]
        end_mo = endstr[4:6]
        startstr = start_yr+'-'+start_mo+'-'+startstr[6:8]
        start = startstr+' 00:00:00'
        endstr = end_yr+'-'+end_mo+'-'+endstr[6:8]
        end = endstr+' 23:59:59'
    else:
        start = '2016-12-01 00:00:00'
        end = '2016-12-31 23:59:59'

    s_dt = datetime.strptime(start, '%Y-%m-%d %H:%M:%S')    #- 20161201 00:00:00
    e_dt = datetime.strptime(end, '%Y-%m-%d %H:%M:%S')      #- 20161231 23:59:59
    in_date = datearg 

    #- Get data from netcdf files using cdippy
    params = 'waveHs,waveTp,waveDp,waveEnergyDensity'
    params = params.split(',')

    #- Create Station Data Object
    stn_data = StnData(stn)
    meta = stn_data.get_stn_meta() 
    meta_name = meta['metaStationName']         #- 'TORREY PINES OUTER, CA BUOY - 100p1
    stn_name,stn_number = meta_name.split(',')  #- ['TORREY PINES OUTER, CA BUOY', '100p1']

    #- Get station data w/ dict_keys(['waveTp', 'waveTime', 'waveHs', 'waveDp'])
    data = stn_data.get_series(start, end, params)
    if not data:
        raise Http404('Error: no data stn '+stn)

    # Prepare data to show gaps where there is no data
    index_name = 'waveTime'
    if len(data) == 0:
        raise Http404('Error: Zero length data records')

    time_count = len(data['waveTime'])
    a = Archive(stn)
    a.set_request_info(start,end,vrs=['waveFrequency','waveBandwidth'])
    d = a.get_request()
    bw = d['waveBandwidth']
    freqs = d['waveFrequency']
    a0 = data['waveEnergyDensity']
    time_count,freq_count = a0.shape

    #- From f90 wc5_calculate_moments(wset,m0,m1,m4,n1)
    m0 = np.zeros(time_count)
    m1 = np.zeros(time_count)
    m2 = np.zeros(time_count)
    m4 = np.zeros(time_count)
    n1 = np.zeros(time_count)
    hm0 = np.zeros(time_count)
    for i in range(time_count):
        for j in range(freq_count):
            bener = a0[i,j] * bw[j]
            m0[i] = m0[i] + bener
            m1[i] = m1[i] + bener * freqs[j]
            m2[i] = m2[i] + bener * freqs[j]**2.0
            m4[i] = m4[i] + bener * freqs[j]**4.0
            n1[i] = n1[i] + bener * freqs[j]**(-1.0)
            hm0[i] = 4*(math.sqrt(m0[i]))

    #- Convert to wavetime to datetime objects
    wT = [timestamp_to_datetime(x) for x in data['waveTime']]
    #- Convert to datetime to string for printing
    print_time = [datetime_to_format(x) for x in wT]
    start_str = wT[0].strftime("%m/%d/%y")
    end_str = wT[time_count-1].strftime("%m/%d/%y")
    num_str = 'Number of samples: '+str(time_count)

    te_array = []                     #- wave te array
    for i in range(time_count):
        waveTe = n1[i] / m0[i]
        te_array.append(waveTe)

    hs_array = data['waveHs']
    power_array = calc_power(hs_array,te_array)
    #energy_array = calc_energy(hs_array)
    total_power = np.sum(power_array)

    num_recs = len(hs_array)

    #- Set up bins for Te and Hs
    min_hs = 0.0
    max_hs = math.ceil(max(hs_array))
    if (max_hs < 5):
        max_hs = 5.0
    num_bins = (max_hs - min_hs) * 2 + 1
    bins_hs = np.linspace(min_hs,max_hs,num_bins)
    bins_te = np.linspace(3,21,19)
    yv,xv = np.meshgrid(bins_hs,bins_te)
    Nx = len(bins_te)
    Ny = len(bins_hs)
    count_array = np.zeros((Nx,Ny))
    percent_array = np.zeros((Nx,Ny))
    cellpower_array = np.zeros((Nx,Ny))
    cellpercent_array = np.zeros((Nx,Ny))
 
    #- Loop through bins and compute counts in each
    #- bin, not including endpoints
    bin_i = 0
    for j in range(len(bins_hs)):
        for i in range(len(bins_te)):
            if( (i == Nx-1) & (j == Ny-1)):
                num_vals, = np.where( (te_array > bins_te[i]) & 
                                      (hs_array > bins_hs[j]) )
                cell_power = power_array[num_vals]
            elif( (i == Nx-1) & (j < Ny-1) ):
                num_vals, = np.where( (te_array > bins_te[i]) & 
                                      (hs_array > bins_hs[j])&(hs_array <= bins_hs[j+1])  )
                cell_power = power_array[num_vals]
            elif( (i < Nx-1) & (j == Ny-1) ):
                num_vals, = np.where( (te_array > bins_te[i])&(te_array <= bins_te[i+1]) & 
                                      (hs_array > bins_hs[j]) )
                cell_power = power_array[num_vals]
            elif( (i < Nx-1) & (j < Ny-1) ):
                num_vals, = np.where( (te_array > bins_te[i])&(te_array <= bins_te[i+1]) & 
                                      (hs_array > bins_hs[j])&(hs_array <= bins_hs[j+1])  )
                cell_power = power_array[num_vals]
            else:
                num_vals = []
                cell_power = []
            count_array[i,j] = len(num_vals)
            percent_array[i,j] = (100*len(num_vals)) / (num_recs)
            cell_power = np.sum(cell_power)
            cellpower_array[i,j] = np.sum(cell_power)
            cellpercent_array[i,j] = (100*cell_power) / total_power
            bin_i = bin_i + 1

    # Set the font dictionaries (for plot title and axis titles)
    title_font = {'size':'16', 'color':'black', 'weight':'normal',
                  'verticalalignment':'bottom'} # Bottom vertical alignment for more space
    #axis_font = {'fontname':'Arial', 'size':'9'}
    axis_font = {'size':'9'}
    #cell_font = {'fontname':'Arial', 'size':'8'}
    cell_font = {'size':'8'}

    plot_types = ["frequency","percent","value"]
    #- Loop through plot types and create 3 figures:
    #-   (1) # of occurences in each bin
    #-   (2) Total Power in each cell
    #-   (3) Percent of total power in each cell
    for plot_type in (plot_types):
        out_file = './'+stn+'_'+in_date+'_'+plot_type+'.png'
        #title_str = 'Wave Power Matrix Station '+id+' '+in_date
        title_str = 'CDIP '+id+' - '+stn_name+' BUOY'
        if (plot_type == "frequency"):
            plot_array = percent_array
            scale_label = '% Occurence'
            title_str = title_str+'\nWave Power Matrix: Frequency of Occurence'
            str_format = '%.2f'
            str_format2 = '%d'
        elif (plot_type == "percent"):
            plot_array = cellpercent_array
            scale_label = '% of Total Power'
            title_str = title_str+'\nWave Power Matrix: % of Total Power '+str(int(total_power))+' kW/m'
            str_format = '%.2f'
            str_format2 = '%d'
        elif (plot_type == "value"):
            plot_array = cellpower_array
            scale_label = 'Wave Power (kW/m)'
            title_str = title_str+'\nWave Power Matrix:  Wave Power kW/m'
            str_format = '%.e'
            str_format2 = '%.e'


        vmax = int(np.max(plot_array ))
        #- Following will plot green-yellow heatmap w/ Te on top and Hs on left side
        cmap = LinearSegmentedColormap.from_list('mycmap', [(0 / vmax, 'white'),
                                                        (0.00001, 'green'),
                                                        (vmax*0.5 / vmax, 'yellow'),
                                                        (vmax / vmax, 'red')])
        cmaplist = [cmap(i) for i in range(cmap.N)]
        #- Force first colors to be white
        #cmaplist[0] = (1,1,1,1.0)
        cmap = cmap.from_list('Custom_cmap',cmaplist,cmap.N*8)
        #- define the bins and normalize
        bounds = np.linspace(0,vmax,(vmax-0)+1)
        norm = mpl.colors.BoundaryNorm(bounds,cmap.N)

        #- Create the figure
        width,height=10,5
        fig, ax = plt.subplots(1, 1, figsize=(width, height), dpi=100)
        nudge = 0.05

        #- Colorize the matrix
        heatmap = ax.pcolor(xv,yv,plot_array,cmap=cmap)

        #- Label each cell w/ value
        data = plot_array
        for j in range(count_array.shape[1]-1):
            for i in range(count_array.shape[0]-1):
                x = xv[i,j]
                y = yv[i,j]
                if (data[i,j] > 0):
                    plt.text(x + 0.45, y + 0.25, str_format % data[i, j],
                             horizontalalignment='center',
                             verticalalignment='center', **cell_font )

        #- colorbar legend
        cbar = plt.colorbar(heatmap,format=str_format2)
        cbar.ax.set_ylabel(scale_label)

        #- Label the date range and # of samples in upper right
        right = 1.255
        top = 1.15
        ur_str = num_str+'\n'+start_str+' to '+end_str
        label_font = {'size':'11'}
        ax.text(right, top, ur_str,
           horizontalalignment='right',
           verticalalignment='top', color='blue',
           transform=ax.transAxes,**axis_font)

        #- Format the plot axes
        ax.grid(True, which='minor',axis='both',linestyle='-',color='0.45')
        ax.set_xticks(bins_te,minor=True)
        ax.set_yticks(bins_hs,minor=True)
        ax.set_xlim(bins_te[0], bins_te[-1])
        ax.set_ylim(bins_hs[0], bins_hs[-1])
        plt.xticks(bins_te)
        plt.xlabel('Energy Period Te (s)')
        plt.ylabel('Wave Hs (m)')
        plt.title(title_str,**title_font)
        plt.tight_layout()
        response=HttpResponse(content_type='image/png')
        fig.savefig(response,format='png')
        return response

    #return HttpResponse(stn_name+stn_number)

